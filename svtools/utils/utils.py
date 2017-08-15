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


# TODO: check if record is CPX and make entry per complex interval
def vcf2bedtool(vcf, split_bnd=True, include_samples=False,
                include_strands=True, split_cpx=False, include_infos=None):
    """
    Wrap VCF as a bedtool. Necessary as pybedtools does not support SV in VCF.

    Parameters
    ----------
    vcf : str or pysam.VariantFile
    split_bnd : bool, optional
        Provide two records for each BND, one per breakend
    include_samples : bool, optional
        Provide comma-delimited list of called samples
    include_strands : bool, optional
        Provide breakpoint strandedness
    include_infos : list of str, optional
        INFO fields to add as columns in output. If "ALL" is present in the
        list, all INFO fields will be reported.

    Returns
    -------
    bt : pybedtools.BedTool
        SV converted to Bedtool. Ends of BND records are assigned as pos + 1.
        Included columns: chrom, start, end, name, svtype, strands
    """

    if not isinstance(vcf, pysam.VariantFile):
        vcf = pysam.VariantFile(vcf)

    entry = '{chrom}\t{start}\t{end}\t{name}\t{svtype}'
    if include_strands:
        entry += '\t{strands}'
    if include_samples:
        entry += '\t{samples}'
    if include_infos:
        if 'ALL' in include_infos:
            include_infos = vcf.header.info.keys()
        entry += '\t{infos}'
    entry += '\n'

    def _format_info(info):
        if info is None:
            return 'NA'
        elif isinstance(info, tuple) or isinstance(info, list):
            return ','.join([str(x) for x in info])
        else:
            return str(info)

    # Convert each record in vcf to bed entry
    def _converter():
        for record in vcf:
            chrom = record.chrom
            start = record.pos
            name = record.id
            svtype = record.info['SVTYPE']

            if include_strands:
                strands = record.info.get('STRANDS')
                strands = '.' if strands is None else strands
            if include_samples:
                samples = ','.join(get_called_samples(record))
            if include_infos:
                infos = [record.info.get(key) for key in include_infos]
                infos = [_format_info(v) for v in infos]
                infos = '\t'.join(infos)

            if record.info['SVTYPE'] == 'BND':
                # First end of breakpoint
                end = record.pos + 1
                yield entry.format(**locals())

                # Second end of breakpoint
                if split_bnd:
                    chrom = record.info['CHR2']
                    start = record.stop
                    end = record.stop + 1
                    yield entry.format(**locals())

            elif record.info['SVTYPE'] == 'INS':
                # Only yield insertion sinks for now
                # Treat them as deletions
                # TODO: rename CPX_INTERVALS to SOURCE for insertions
                svtype = 'DEL'
                end = record.stop

                # We permit start > end in insertions in cases of
                # microdup/microhomology
                # Reorder start/end so bedtools doesn't break
                start, end = sorted([start, end])
                yield entry.format(**locals())

            elif record.info['SVTYPE'] == 'CTX':
                end = record.pos + 1
                yield entry.format(**locals())

            elif 'CPX_INTERVALS' in record.info and split_cpx:
                # If complex, all constituent intervals are in CPX_INTERVALS
                for interval in record.info['CPX_INTERVALS']:
                    svtype, region = interval.split('_')
                    chrom, coords = region.split(':')
                    start, end = coords.split('-')
                    yield entry.format(**locals())

            else:
                end = record.stop
                yield entry.format(**locals())

    return pbt.BedTool(_converter()).saveas()
