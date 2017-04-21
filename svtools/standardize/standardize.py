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

    def standardize_vcf(self):
        """
        Standardize every record in a VCF.

        Any filtering of records should be implemented in this method.

        Yields
        ------
        std_rec : pysam.VariantRecord
            Standardized records
        """
        for record in self.raw_vcf:
            std_rec = self.std_vcf.new_record()
            yield self.standardize_record(std_rec, record)

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

        # Standardize tloc ALT after SVTYPE and CHR2/END are standardized
        if std_rec.info['SVTYPE'] == 'BND':
            alt = make_bnd_alt(std_rec.info['CHR2'], std_rec.info['END'],
                               std_rec.info['STRANDS'])
            std_rec.alts = (alt, )

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
        std_rec.info['END'] = raw_rec.pos + 1
        std_rec.info['SVLEN'] = 0
        std_rec.info['SOURCE'] = 'source'

        return std_rec

    def standardize_format(self, std_rec, raw_rec):
        """
        Copy desired FORMAT fields to new record.

        By default, only GT is copied.
        """

        # Add per-sample genotypes (ignoring other FORMAT fields)
        for sample in raw_rec.samples:
            std_rec.samples[sample]['GT'] = raw_rec.samples[sample]['GT']

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
    chr2, end = alt.split(':')
    end = int(end)
    return chr2, end


def parse_bnd_strands(alt):
    """
    Parses standard VCF BND ALT (e.g. N]1:1000]) for strandedness

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
