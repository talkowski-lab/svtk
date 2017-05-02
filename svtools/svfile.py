"""
svfile.py

Wrap the pysam API to permit clustering of standardized SV VCF records.

Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
Distributed under terms of the MIT license.
"""

import numpy as np
from .utils import recip, make_bnd_alt
from .genomeslink import GSNode


class SVFile(object):
    def __init__(self, vcf):
        """
        Wrapper for standardized VCF files.

        Parameters
        ----------
        vcf : pysam.VariantFile
        """
        self.reader = vcf
        self.filename = vcf.filename.decode('utf-8')
        self.samples = list(self.reader.header.samples)

        # Confirm all standard INFO fields are present
        required_info = 'SVTYPE CHR2 END STRANDS SVLEN SOURCES'.split()
        for info in required_info:
            if info not in self.reader.header.info.keys():
                msg = "Required INFO field {0} not found in file {1}"
                msg = msg.format(info, self.filename)
                raise KeyError(msg)

        # Unfortunately no way to index into "source" metadata record
        # via pysam API, must manually check all header records
        self.sources = None
        for hrec in self.reader.header.records:
            if hrec.key == 'source':
                self.sources = hrec.value.split(',')

        if self.sources is None:
            msg = "Source not specified in header of {0}"
            msg = msg.format(self.filename)
            raise KeyError(msg)

    def fetch(self, chrom, start=None, end=None):
        """
        Fetch only calls from specified region of SVFile.

        Requires Tabix index. Updates SVFile in place.

        Parameters
        ----------
        chrom : str
        start : int, optional
        end : int, optional
            Required if start specified
        """
        if start is not None and end is None:
            msg = 'Start {}:{} specified but no end coordinate provided'
            msg = msg.format(chrom, start)
            raise ValueError(msg)

        # First check if VCF is empty
        try:
            pos = self.reader.tell()
            next(self.reader)
            self.reader.seek(pos)
        except StopIteration:
            self.reader = iter(())
            return

        # Then check if index is present
        try:
            self.reader = self.reader.fetch(chrom, start, end)
        except ValueError:
            msg = 'No index found for {0}'.format(self.filename)
            raise FileNotFoundError(msg)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        record = next(self.reader)
        return SVRecord(record)


class SVRecord(GSNode):
    """
    Clusterable VCF record.
    """

    def __init__(self, record):
        """
        record : pysam.VariantRecord
            Must specify 'CHR2' and 'END' in INFO
        """

        self.record = record
        self.sources = record.info['SOURCES']

        chrA = record.chrom
        posA = record.pos
        chrB = record.info['CHR2']
        posB = record.info['END']
        name = record.id

        super().__init__(chrA, posA, chrB, posB, name)

    def clusters_with(self, other, dist, frac=0.0, match_strands=False):
        svtype_match = self.svtype == other.svtype

        if match_strands:
            strand_match = self.info['STRANDS'] == other.info['STRANDS']
        else:
            strand_match = True

        if self.is_tloc or other.is_tloc:
            overlap = True
        else:
            overlap = recip(self.posA, self.posB, other.posA, other.posB, frac)

        return(super().clusters_with(other, dist) and
               svtype_match and overlap and strand_match)

    @property
    def svtype(self):
        """
        Returns
        -------
        svtype : str
            One of {DEL, DUP, INV, BND}
        """
        return self.record.info['SVTYPE']

    @property
    def is_tloc(self):
        return self.chrA != self.chrB

    def __hash__(self):
        return id(self)


class SVRecordCluster:
    def __init__(self, records):
        self.records = records

    def sources(self):
        """
        Return list of source algorithms in clustered records.

        By default, searches for both a SOURCE and a SOURCES INFO field

        Returns
        -------
        sources : list of str
        """
        call_sources = set()
        for record in self.records:
            sources = record.record.info.get('SOURCES')
            if sources:
                call_sources = call_sources.union(sources)

        return sorted(call_sources)

    def merge_record_data(self, new_record):
        """
        Aggregate metadata (coordinates, alts, and INFO) of clustered records.

        * Original record IDs are preserved in a new INFO field

        Parameters
        ----------
        new_record : pysam.VariantRecord
            Blank record to fill with aggregated data

        Returns
        -------
        new_record : pysam.VariantRecord
            VCF record populated with aggregate cluster data
        """

        # Secondary records have duplicate POS/END/INFO
        # TODO: move secondary filtering elsewhere

        if len(self.records) == 0:
            return None

        base_record = self.records[0]

        new_record.chrom = base_record.chrA
        new_record.ref = base_record.record.ref
        new_record.info['SVTYPE'] = base_record.svtype
        new_record.info['CHR2'] = base_record.chrB

        # TODO: Check if strands should be merged
        new_record.info['STRANDS'] = base_record.record.info['STRANDS']

        # Merge coordinates
        POS, END, CIPOS, CIEND = self.merge_pos()
        new_record.pos = POS
        new_record.info['END'] = END
        new_record.info['CIPOS'] = CIPOS
        new_record.info['CIEND'] = CIEND

        # Assign alts, updating translocation alt based on merged coordinates
        if new_record.info['SVTYPE'] == 'BND':
            strands = new_record.info['STRANDS']
            alt = make_bnd_alt(base_record.chrB, END, strands)
            new_record.alts = (alt, )
        else:
            new_record.alts = base_record.record.alts

        # SVLEN for intra-chromosomal is -1
        if base_record.is_tloc:
            new_record.info['SVLEN'] = -1
        else:
            new_record.info['SVLEN'] = END - POS

        # QUAL, FILTER currently unused
        new_record.filter.add('PASS')

        # Report cluster RMSSTD
        new_record.info['RMSSTD'] = self.rmsstd

        # List of aggregate sources
        new_record.info['SOURCES'] = self.sources()

        return new_record

    def merge_record_formats(self, new_record, sourcelist, single_source=False,
                             call_sources=False):
        """
        Aggregate sample genotype data across records.

        1) Set GT to 0/1 for samples called in any record, 0/0 otherwise.
        2) For each provided source, set a corresponding FORMAT field to 1
           for samples called by the source. If the FORMAT field is available
           in the record being merged, use per-sample data. If the FORMAT field
           is not available, use the record's SOURCE or SOURCES to determine
           support for all samples called in that record.

        Parameters
        ----------
        new_record : pysam.VariantRecord
            Record to populate with FORMAT data
        sourcelist : list of str
            List of all sources to add a FORMAT field for
        call_sources : bool, optional
            If True, each record is already annotated with FORMAT data
            indicating the supporting source algorithms for each sample.
            If False, use each record's SOURCE key to derive a supporting
            algorithm for all samples called in the record.

        Returns
        -------
        new_record : pysam.VariantRecord
            Populated record
        """

        # Seed with null values
        for sample in new_record.samples:
            new_record.samples[sample]['GT'] = (0, 0)

            for source in sourcelist:
                new_record.samples[sample][source] = 0

        # Update with called samples
        # TODO: optionally permit ./. instead of rejecting
        # I think that was an issue with one caller, maybe handle in preproc
        null_GTs = [(0, 0), (None, None), (0, ), (None, )]
        for record in self.records:
            for sample in record.record.samples:
                gt = record.record.samples[sample]['GT']

                # Skip samples without a call
                if gt in null_GTs:
                    continue

                # Otherwise call the sample in the new record
                new_record.samples[sample]['GT'] = (0, 1)

                # If call sources are already provided, add them
                if call_sources:
                    for source in sourcelist:
                        # Get the sample's current source call and the
                        # record's source call for that sample
                        call = new_record.samples[sample].get(source)
                        src_call = record.record.samples[sample].get(source)

                        # If sample was already called by the source, or if
                        # it was called by the source in the current record,
                        # call it in the new record
                        call = call or src_call
                        new_record.samples[sample][source] = call
                # Otherwise add the record's SOURCE
                else:
                    for source in record.record.info['SOURCES']:
                        new_record.samples[sample][source] = 1

        return new_record

    def merge_pos(self):
        """
        Compute aggregate POS/END of clustered SVRecords.

        Defaults to computing median of POS and END over constituent
        records. CIPOS/CIEND are calculated as (MIN - MEDIAN, MAX - MEDIAN)
        over POS and END.

        Returns
        -------
        POS : int
        END : int
        CIPOS : list of int
        CIEND : list of int
        """

        # Position bounds
        MIN_POS = min(rec.posA for rec in self.records)
        MAX_POS = max(rec.posA for rec in self.records)
        MIN_END = min(rec.posB for rec in self.records)
        MAX_END = max(rec.posB for rec in self.records)

        POS = int(np.median([rec.posA for rec in self.records]))
        END = int(np.median([rec.posB for rec in self.records]))
        CIPOS = [MIN_POS - POS, MAX_POS - POS]
        CIEND = [MIN_END - END, MAX_END - END]

        return POS, END, CIPOS, CIEND

    @property
    def rmsstd(self):
        """
        Root-mean-square standard deviation of cluster coordinates
        """

        if hasattr(self, '_rmsstd'):
            return self._rmmstd

        starts = np.array([record.posA for record in self.records])
        ends = np.array([record.posB for record in self.records])

        def _meanSS(X):
            mu = np.mean(X)
            return np.sum((X - mu) ** 2) / len(X)

        SS = _meanSS(starts) + _meanSS(ends)
        self._rmmstd = np.sqrt(SS)

        return self._rmmstd
