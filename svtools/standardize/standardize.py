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

Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
Distributed under terms of the MIT license.
"""


class VCFStandardizer:
    subclasses = {}

    def __init__(self, raw_vcf, std_vcf):
        """
        Standardize a VCF.

        Parameters
        ----------
        raw_vcf : pysam.VariantFile
            Input VCF.
        std_vcf : pysam.VariantFile
            Standardized VCF.
        """
        self.raw_vcf = raw_vcf
        self.std_vcf = std_vcf

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

    def standardize(self):
        """
        Yields
        ------
        std_rec : pysam.VariantRecord
            Standardized records
        """
        for record in self.filter_vcf():
            yield self.standardize_record(record)

    def standardize_record(self, record):
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
        std_rec = self.copy_basics(std_rec, record)

        # Standardize the required INFO fields
        std_rec = self.standardize_info(std_rec, record)

        # Standardize tloc ALT after SVTYPE and CHR2/END are standardized
        if std_rec.info['SVTYPE'] == 'BND':
            alt = make_bnd_alt(std_rec.info['CHR2'], std_rec.info['END'],
                               std_rec.info['STRANDS'])
            std_rec.alts = (alt, )

        # Add per-sample genotypes (ignoring other FORMAT fields)
        for sample in record.samples:
            std_rec.samples[sample]['GT'] = record.samples[sample]['GT']

        return std_rec

    def copy_basics(self, std_rec, raw_rec):
        """
        Copy basic record data (CHROM, POS, ID, REF, ALT).

        Reset FILTER to PASS.
        """
        std_rec.chrom = raw_rec.chrom
        std_rec.pos = raw_rec.pos
        std_rec.id = raw_rec.id
        std_rec.ref = raw_rec.ref
        std_rec.alts = raw_rec.alts

        # Strip filters
        std_rec.filter.add('PASS')

    def filter_vcf(self):
        """
        Filter records in raw VCF.

        By default, no filtering is applied.
        """

        for record in self.raw_vcf:
            yield record

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
        std_rec.info['END'] = raw_rec.pos + 1
        std_rec.info['SVLEN'] = 0
        std_rec.info['SOURCE'] = 'source'

        return std_rec


def standardize_vcf(raw_vcf, std_vcf, std_fn=None, filter_fn=None):
    """
    Iterator over the construction of new standardized records.

    Parameters
    ----------
    raw_vcf : pysam.VariantFile
        Input VCF.
    std_vcf : pysam.VariantFile
        Output VCF. Required to construct new VariantRecords.
    std_fn : (pysam.VariantRecord, pysam.VariantRecord) -> pysam.VariantRecord
        Standardization function for converting each record
    filter_fn : (pysam.VariantFile -> iter of pysam.VariantRecord)
        Optional filtering of raw VCF

    Yields
    ------
    std_rec : pysam.VariantRecord
        Standardized VCF record.
    """
    if filter_fn is not None:
        raw_vcf = filter_fn(raw_vcf)

    for raw_rec in raw_vcf:
        std_rec = std_vcf.new_record()
        std_rec = standardize_record(raw_rec, std_rec)
        yield std_rec


# TODO: pass standardization function instead of source name and allow
# user specification of function with file at command line
def standardize_record(raw_rec, std_rec, source='delly'):
    """
    Copies basic record data and standardizes INFO/FORMAT fields.

    Parameters
    ----------
    raw_rec : pysam.VariantRecord
        VCF record to standardize.
    std_rec : pysam.VariantRecord
        Empty VariantRecord constructed from new VariantFile.

    Returns
    -------
    std_rec : pysam.VariantRecord
        New VariantRecord with standardized data filled in.
    """

    # Copy basic record data

    # Standardize INFO fields, and update basic data as necessary
    #  if source == 'delly':
    #      std_rec = standardize_wham(raw_rec, std_rec)
    #      std_rec = standardize_delly(raw_rec, std_rec)

    return std_rec


def make_bnd_alt(chrom, pos, strands):
    """
    Make ALT for BND record in accordance with VCF specification.
    """

    p = '{0}:{1}'.format(chrom, pos)

    if strands == '++':
        fmt = 'N]{0}]'
    elif strands == '+-':
        fmt = 'N[{0}['
    elif strands == '-+':
        fmt = ']{0}]N'
    elif strands == '--':
        fmt = '[{0}[N'

    return fmt.format(p)
