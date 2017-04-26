"""
svfile.py

Wrap the pysam API to permit clustering of standardized SV VCF records.

Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
Distributed under terms of the MIT license.
"""

import numpy as np
from pysam import VariantFile
from .utils import recip, make_bnd_alt
from .genomeslink import GSNode


class UnsupportedFiletypeError(Exception):
    """Unsupported or unexpected filetype"""


class UnboundedFetchError(Exception):
    """Only start position of tabix fetch region specified, no end"""


class TabixNotFoundError(Exception):
    """No Tabix index found"""


class PyVCFParsingError(Exception):
    """Unexpected datatype from PyVCF"""


class MissingInfoError(Exception):
    """Missing required INFO field"""


class MissingSourceError(Exception):
    """Source not specified in VCF header"""


class SVFile(object):
    def __init__(self, filename):
        """
        Wrapper for standardized VCF files.
        """
        self.filename = filename

        self.reader = VariantFile(filename)

        # Confirm all standard INFO fields are present
        required_info = 'SVTYPE CHR2 END STRANDS SVLEN SOURCE'.split()
        for info in required_info:
            if info not in self.reader.header.info.keys():
                msg = "Required INFO field {0} not found in file {1}"
                msg = msg.format(info, filename)
                raise MissingInfoError(msg)

        # Unfortunately no way to index into "source" metadata record
        # via pysam API, must manually check all header records
        self.source = None
        for hrec in self.reader.header.records:
            if hrec.key == 'source':
                self.source = hrec.value
        if self.source is None:
            msg = "Source not specified in header of {0}".format(filename)
            raise MissingSourceError(msg)

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

        Raises
        ------
        UnboundedFetchError
        TabixNotFoundError
        """
        if start is not None and end is None:
            msg = 'Start {}:{} specified but no end coordinate provided'
            msg = msg.format(chrom, start)
            raise UnboundedFetchError(msg)

        try:
            self.reader = self.reader.fetch(chrom, start, end)
        except ValueError:
            raise TabixNotFoundError(self.filename)

    @property
    def samples(self):
        """
        Returns
        ------
        samples : list of str
        """
        if hasattr(self.reader.header, 'samples'):
            return list(self.reader.header.samples)
        else:
            return []

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
        self.source = record.info['SOURCE']

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

    @staticmethod
    def merge_pos(records):
        """
        Compute aggregate POS/END of clustered SVRecords.

        Defaults to computing median of POS and END over constituent
        records. CIPOS/CIEND are calculated as (MIN - MEDIAN, MAX - MEDIAN)
        over POS and END.

        Parameters
        ----------
        records : list of SVRecord

        Returns
        -------
        POS : int
        END : int
        CIPOS : list of int
        CIEND : list of int
        """

        # Position bounds
        MIN_POS = min(rec.posA for rec in records)
        MAX_POS = max(rec.posA for rec in records)
        MIN_END = min(rec.posB for rec in records)
        MAX_END = max(rec.posB for rec in records)

        POS = int(np.median([rec.posA for rec in records]))
        END = int(np.median([rec.posB for rec in records]))
        CIPOS = [MIN_POS - POS, MAX_POS - POS]
        CIEND = [MIN_END - END, MAX_END - END]

        return POS, END, CIPOS, CIEND

    @classmethod
    def merge(cls, new_record, records, sources):
        """
        Aggregate clustered records.

        * Original record IDs are preserved in a new INFO field

        Parameters
        ----------
        new_record : pysam.VariantRecord
            Blank record to fill with aggregated data
        records : list of SVRecord
            Clustered records
        sources : list of str
            List of all sources
        """

        # Secondary records have duplicate POS/END/INFO
        # TODO: move secondary filtering elsewhere

        if len(records) == 0:
            return None

        base_record = records[0]

        new_record.chrom = base_record.chrA
        new_record.ref = base_record.record.ref
        new_record.info['SVTYPE'] = base_record.svtype
        new_record.info['CHR2'] = base_record.chrB

        # TODO: Check if strands should be merged
        new_record.info['STRANDS'] = base_record.record.info['STRANDS']

        # Merge coordinates
        POS, END, CIPOS, CIEND = cls.merge_pos(records)
        new_record.pos = POS
        new_record.info['END'] = END
        new_record.info['CIPOS'] = CIPOS
        new_record.info['CIEND'] = CIEND

        call_sources = sorted(set([r.record.info['SOURCE'] for r in records]))
        new_record.info['SOURCES'] = call_sources

        # Assign alts, updating translocation alt based on merged coordinates
        if new_record.info['SVTYPE'] == 'BND':
            strands = new_record.info['STRANDS']
            alt = make_bnd_alt(base_record.chrB, END, strands)
            new_record.alts = (alt, )
        else:
            new_record.alts = base_record.alts

        # SVLEN for intra-chromosomal is -1
        if base_record.is_tloc:
            new_record.info['SVLEN'] = -1
        else:
            new_record.info['SVLEN'] = END - POS

        # QUAL, FILTER currently unused
        new_record.filter.add('PASS')

        # Seed with null values
        for sample in new_record.samples:
            new_record.samples[sample]['GT'] = (0, 0)
            for source in sources:
                new_record.samples[sample][source] = 0

        # Update with called samples
        # TODO: optionally permit ./. instead of rejecting
        # I think that was an issue with one caller, maybe handle in preproc
        null_GTs = [(0, 0), (None, None), (0, ), (None, )]
        for record in records:
            for sample in record.samples:
                gt = record.samples[sample]['GT']
                if gt not in null_GTs:
                    new_record.samples[sample]['GT'] = (0, 1)

                    source = record.info['SOURCE']
                    new_record.samples[sample][source] = 1

        return SVRecord(record)

    def __hash__(self):
        return id(self)
