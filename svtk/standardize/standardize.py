# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
"""
standardize.py

Standardize a VCF of SV calls.

Each record corresponds to a single SV breakpoint and will have the following
INFO fields, with specified constraints:
  SVTYPE:  SV type [DEL,DUP,INV,BND]
  CHR2:    Secondary chromosome [Must be lexicographically greater than CHROM]
  END:     SV end position (or position on CHR2 in translocations)
  STRANDS: Breakpoint strandedness [++,+-,-+,--]
  SVLEN:   SV length (-1 if translocation)
"""


from svtk.utils import make_bnd_alt, NULL_GT


def any_called(record):
    null_GTs = [(0, 0), (None, None), (0, ), (None, )]

    def _is_called(sample):
        return record.samples[sample]['GT'] not in null_GTs

    return any([_is_called(sample) for sample in record.samples])


class VCFStandardizer:
    subclasses = {}

    def __init__(self, raw_vcf, std_vcf, prefix=None, min_size=50,
                 include_reference_sites=False, call_null_sites=False):
        """
        Standardize a VCF.

        Parameters
        ----------
        raw_vcf : pysam.VariantFile
            Input VCF.
        std_vcf : pysam.VariantFile
            Standardized VCF.
        prefix : str, optional
            ID prefix to assign to standardized records
        min_size : int, optional
            Minimum SV size to report
        include_reference_sites : bool, optional
            Keep sites where no sample has a variant call (i.e. all 0/0 or ./.)
        call_null_sites : bool, optional
            Convert null genotypes (./.) to variant calls (0/1)
        """
        self.raw_vcf = raw_vcf
        self.std_vcf = std_vcf
        self.prefix = prefix
        self.min_size = min_size
        self.include_reference_sites = include_reference_sites
        self.call_null_sites = call_null_sites

    @classmethod
    def register(cls, source):
        def decorator(subclass):
            cls.subclasses[source] = subclass
            return subclass
        return decorator

    @classmethod
    def create(cls, source, *args):
        if source not in cls.subclasses:
            msg = 'No standardizer defined for {0}'.format(source)
            raise ValueError(msg)
        return cls.subclasses[source](*args)

    def filter_vcf(self):
        """
        Filter records.

        Records are excluded if:
        1) They exist on contigs excluded from new file
        2) They are tagged as SECONDARY
        """
        for record in self.raw_vcf:
            if record.chrom not in self.std_vcf.header.contigs:
                continue
            # Filter on chr2 if a breakend
            if '[' in record.alts[0] or ']' in record.alts[0]:
                chr2, end = parse_bnd_pos(record.alts[0])
                if chr2 not in self.std_vcf.header.contigs:
                    continue

            # Skip SECONDARY events to avoid double counting
            if 'SECONDARY' in record.info.keys():
                continue

            yield record

    def standardize_vcf(self):
        """
        Standardize every record in a VCF.

        Any filtering of records should be implemented in this method.

        Yields
        ------
        std_rec : pysam.VariantRecord
            Standardized records
        """
        idx = 1
        for record in self.filter_vcf():
            std_rec = self.std_vcf.new_record()
            std_rec = self.standardize_record(std_rec, record)

            # Apply size filter (but keep breakends (SVLEN=-1))
            if 0 < std_rec.info['SVLEN'] < self.min_size:
                continue

            # Exclude sites with no called samples unless requested otherwise
            if not any_called(std_rec) and not self.include_reference_sites:
                continue

            # Filter unstranded breakpoints
            if std_rec.info['STRANDS'] not in '++ +- -+ --'.split():
                continue

            # Assign new variant IDs
            if self.prefix is not None:
                std_rec.id = '{0}_{1}'.format(self.prefix, idx)
                idx += 1

            yield std_rec

    def standardize_record(self, std_rec, raw_rec):
        """
        Create a standardized copy of a VCF record.

        Parameters
        ----------
        record : pysam.VariantRecord

        Returns
        -------
        std_rec : pysam.VariantRecord
        """

        # Construct a new record and copy basic VCF fields
        std_rec = self.std_vcf.new_record()
        std_rec.chrom = raw_rec.chrom
        std_rec.pos = raw_rec.pos
        std_rec.id = raw_rec.id
        std_rec.ref = raw_rec.ref
        std_rec.alts = raw_rec.alts

        # Strip filters
        std_rec.filter.add('PASS')

        # Standardize the required INFO fields
        std_rec = self.standardize_info(std_rec, raw_rec)
        std_rec = self.standardize_alts(std_rec, raw_rec)
        std_rec = self.standardize_format(std_rec, raw_rec)

        return std_rec

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize VCF record INFO.

        When implementing this function, assume the basic data fields (CHROM,
        POS, ID, REF, ALTS, and FILTER) have been copied directly from the
        original record.

        The default implementation is a placeholder for testing and should be
        overridden with the logic for each algorithm's formatting.

        Parameters
        ----------
        std_rec : pysam.VariantRecord
            Standardized record with populated basic data.
        raw_rec : pysam.VariantRecord
            Raw record to be standardized.
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']
        std_rec.info['CHR2'] = raw_rec.chrom
        std_rec.stop = raw_rec.pos + 1
        std_rec.info['SVLEN'] = 0
        std_rec.info['ALGORITHMS'] = ['source']

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Copy desired FORMAT fields to new record.

        By default, copy GT and tag source FORMAT appropriately.
        """

        source = std_rec.info['ALGORITHMS'][0]

        # Add per-sample genotypes (ignoring other FORMAT fields)
        for sample in raw_rec.samples:
            gt = raw_rec.samples[sample]['GT']
            if self.call_null_sites:
                if gt == (None, None):
                    gt = (0, 1)
                if gt == (None,):
                    gt = (1,)
            std_rec.samples[sample]['GT'] = gt

            if gt not in NULL_GT:
                std_rec.samples[sample][source] = 1
            else:
                std_rec.samples[sample][source] = 0

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize ALT field.

        Default behavior is to standardize BND alt to VCF spec and leave
        other SVTYPE alts untouched.
        """

        # Standardize tloc ALT after SVTYPE and CHR2/END are standardized
        if std_rec.info['SVTYPE'] == 'BND':
            alt = make_bnd_alt(std_rec.info['CHR2'], std_rec.stop,
                               std_rec.info['STRANDS'])
            stop = std_rec.stop
            std_rec.alts = (alt, )
            std_rec.stop = stop

        return std_rec


def parse_bnd_pos(alt):
    """
    Parses standard VCF BND ALT (e.g. N]1:1000]) into chrom, pos

    Parameters
    ----------
    alt : str
        VCF-formatted BND ALT

    Returns
    -------
    chrom : str
    pos : int
    """
    alt = alt.strip('ATCGN')
    # Strip brackets separately, otherwise GL contigs will be altered
    alt = alt.strip('[]')

    # HLA contigs include colons, so be careful when parsing
    data = alt.split(':')
    chr2 = ':'.join(data[:-1])
    end = int(data[-1])

    return chr2, end


def parse_bnd_strands(alt):
    """
    Parses standard VCF BND ALT (e.g. N]1:1000]) for strandedness

    Note about parsing strands from BND ALT:
    t[p[ piece extending to the right of p is joined after t (+-)
    t]p] reverse comp piece extending left of p is joined after t (++)
    ]p]t piece extending to the left of p is joined before t (-+)
    [p[t reverse comp piece extending right of p is joined before t (--)

    Parameters
    ----------
    alt : str
        VCF-formatted BND ALT

    Returns
    -------
    strands : str
        ++,+-,-+,--
    """
    if alt.endswith('['):
        return '+-'
    elif alt.endswith(']'):
        return '++'
    elif alt.startswith(']'):
        return '-+'
    elif alt.startswith('['):
        return '--'
