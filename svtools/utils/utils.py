# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.
"""
utils.py

Helper functions for svtools.
"""

from collections import deque
import pysam
import pybedtools as pbt


NULL_GT = [(0, 0), (None, None), (0, ), (None, )]


def is_smaller_chrom(chrA, chrB):
    """
    Test if chrA is naturally less than chrB

    Returns True if chrA == chrB so comparison will default to position
    """

    if chrA.startswith('chr'):
        chrA = chrA[3:]
    if chrB.startswith('chr'):
        chrB = chrB[3:]

    # Numeric comparison, if possible
    if chrA.isdigit() and chrB.isdigit():
        return int(chrA) <= int(chrB)

    # String comparison for X/Y
    elif not chrA.isdigit() and not chrB.isdigit():
        return chrA <= chrB

    # Numeric is always less than X/Y
    else:
        return chrA.isdigit()


def recip(startA, endA, startB, endB, frac):
    """
    Test if two intervals share a specified reciprocal overlap.
    """

    if frac == 0:
        return True

    start = max(startA, startB)
    end = min(endA, endB)
    olen = end - start
    lenA = endA - startA
    lenB = endB - startB

    try:
        lapA = olen / float(lenA)
        lapB = olen / float(lenB)
    except ZeroDivisionError:
        return False

    return (olen > 0) and (lapA >= frac) and (lapB >= frac)


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


def get_called_samples(record, include_null=False):
    """
    Return list of samples with variant call

    Parameters
    ----------
    record : pysam.VariantRecord
    include_null : bool
        Include samples without an explicit reference (0/0) call (i.e. ./.)

    Returns
    -------
    samples : list of str
        Sorted list of sample IDs with a variant call
    """

    samples = deque()
    for sample in record.samples.keys():
        if record.samples[sample]['GT'] not in NULL_GT:
            samples.append(sample)

    return sorted(samples)


# TODO: handle other end of interchromosomal BND
# TODO: check if record is CPX and make entry per complex interval
def vcf2bedtool(vcf):
    """
    Wrap VCF as a bedtool. Necessary as pybedtools does not support SV in VCF.

    Parameters
    ----------
    vcf : str or pysam.VariantFile

    Returns
    -------
    bt : pybedtools.BedTool
        SV converted to Bedtool. Ends of BND records are assigned as pos + 1.
        Included columns: chrom, start, end, name, svtype, strands
    """

    if not isinstance(vcf, pysam.VariantFile):
        vcf = pysam.VariantFile(vcf)

    # Convert each record in vcf to bed entry
    def _converter():
        bed = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'
        for record in vcf:
            if record.info['SVTYPE'] == 'BND':
                end = record.pos + 1
            else:
                end = record.info['END']
            yield bed.format(record.chrom, record.pos, end, record.id,
                             record.info['SVTYPE'], record.info['STRANDS'])

    return pbt.BedTool(_converter()).saveas()
