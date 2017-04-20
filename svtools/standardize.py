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

from svtools.utils import is_smaller_chrom


def standardize_vcf(raw_vcf, std_vcf):
    """
    Iterator over the construction of new standardized records.

    Arguments
    ---------
    raw_vcf : pysam.VariantFile
        Input VCF.
    std_vcf : pysam.VariantFile
        Output VCF. Required to construct new VariantRecords.

    Yields
    ------
    std_rec : pysam.VariantRecord
        Standardized VCF record.
    """
    for raw_rec in raw_vcf:
        std_rec = std_vcf.new_record()
        std_rec = standardize_record(raw_rec, std_rec)
        yield std_rec


def standardize_record(raw_rec, std_rec, source='delly'):
    """
    Copies basic record data and standardizes INFO/FORMAT fields.

    Arguments
    ---------
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
    std_rec.chrom = raw_rec.chrom
    std_rec.pos = raw_rec.pos
    std_rec.id = raw_rec.id
    std_rec.ref = raw_rec.ref
    std_rec.alts = raw_rec.alts

    # Strip filters
    std_rec.filter.add('PASS')

    # Copy defined INFO fields
    if source == 'delly':
        std_rec = standardize_delly(raw_rec, std_rec)
    else:
        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']
        std_rec.info['CHR2'] = raw_rec.chrom
        std_rec.info['END'] = raw_rec.pos + 1
        std_rec.info['SVLEN'] = 0
        std_rec.info['SOURCE'] = 'source'

    # Add per-sample formats
    for sample in raw_rec.samples:
        std_rec.samples[sample]['GT'] = raw_rec.samples[sample]['GT']

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


def standardize_delly(raw_rec, std_rec):
    """
    Standardize Delly record.

    1) Rename 'TRA' to 'BND'.
    2) Swap CHROM and CHR2 in translocations.
    3) Add END.
    4) Rename 'CT' to 'STRANDS' and convert notation.
    5) Compute SVLEN.
    6) Add SOURCE.
    7) Standardize ALT to VCF spec.
    """

    # Rename TRA to BND
    svtype = raw_rec.info['SVTYPE']
    if svtype == 'TRA':
        svtype = 'BND'
    std_rec.info['SVTYPE'] = svtype

    # Convert strandedness notation
    raw_strands = raw_rec.info['CT']
    if raw_strands == '5to3':
        strands = '-+'
    elif raw_strands == '3to3':
        strands = '++'
    elif raw_strands == '5to5':
        strands = '--'
    elif raw_strands == '3to5':
        strands = '+-'
    std_rec.info['STRANDS'] = strands

    pos, end = raw_rec.pos, raw_rec.info['END']

    # Swap CHR2/CHROM if necessary and update ALT
    if svtype == 'BND':
        chrom, chr2 = raw_rec.chrom, raw_rec.info['CHR2']

        # swap chr2/chrom, pos/end, and reverse strandedness
        if not is_smaller_chrom(chrom, chr2):
            std_rec.pos, end = end, pos
            std_rec.chrom, chr2 = chr2, chrom
            std_rec.info['STRANDS'] = strands[::-1]

        # make standard BND ALT
        alt = make_bnd_alt(chr2, end, strands)
        std_rec.alts = (alt, )
        pass
    else:
        chr2 = raw_rec.chrom

    # Add CHR2 and END
    std_rec.info['END'] = end
    std_rec.info['CHR2'] = chr2

    # Add SVLEN
    if std_rec.chrom == std_rec.info['CHR2']:
        std_rec.info['SVLEN'] = std_rec.info['END'] - std_rec.pos
    else:
        std_rec.info['SVLEN'] = -1

    return std_rec
