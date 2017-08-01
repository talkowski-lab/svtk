#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Resolve complex SV from inversion/translocation breakpoints and CNV intervals.
"""

import argparse
import sys
from collections import deque
import pysam
import svtools.utils as svu
from svtools.cxsv import link_cpx, ComplexSV


CPX_INFO = [
    '##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">',
    '##ALT=<ID=CPX,Description="Complex SV">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=UNR,Description="Unresolved breakend or complex SV">',
    '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Class of complex variant.">',
    '##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Genomic intervals constituting complex variant.">',
    '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">'
]


def _merge_records(vcf, cpx_records, cpx_record_ids):
    """
    r1, r2 : iter of pysam.VariantRecord
    """

    SENTINEL = 'SENTINEL'

    def _next_record():
        try:
            return next(vcf)
        except StopIteration:
            return SENTINEL

    def _next_cpx():
        try:
            return cpx_records.popleft()
        except IndexError:
            return SENTINEL

    # Initialize merge
    curr_record = _next_record()
    curr_cpx = _next_cpx()

    while curr_record != SENTINEL and curr_cpx != SENTINEL:
        # Remove VCF records that were included in complex event
        if curr_record.id in cpx_record_ids:
            curr_record = _next_record()
            continue

        # Merge sort remaining
        if curr_record.chrom == curr_cpx.chrom:
            if curr_record.pos <= curr_cpx.pos:
                yield curr_record
                curr_record = _next_record()
            else:
                yield curr_cpx
                curr_cpx = _next_cpx()

        elif svu.is_smaller_chrom(curr_record.chrom, curr_cpx.chrom):
            yield curr_record
            curr_record = _next_record()
        else:
            yield curr_cpx
            curr_cpx = _next_cpx()

    # After one iterator is exhausted, return rest of other iterator
    if curr_record == SENTINEL:
        for cpx in cpx_records:
            yield cpx

    elif curr_cpx == SENTINEL:
        for record in vcf:
            if record.id not in cpx_record_ids:
                yield record


def resolve_complex_sv(vcf, variant_prefix='CPX_'):
    """
    Resolve complex SV from CNV intervals and BCA breakpoints.

    Yields all resolved events, simple or complex, in sorted order.

    Parameters
    ----------
    vcf : pysam.VariantFile
    variant_prefix : str
        Prefix to assign to resolved variants

    Yields
    ------
    sv : pysam.VariantRecord
    """

    clusters = link_cpx(vcf)

    resolved_idx = unresolved_idx = 1

    if not variant_prefix.endswith('_'):
        variant_prefix += '_'

    cpx_records = deque()
    cpx_record_ids = set()

    for cluster in clusters:
        cpx = ComplexSV(cluster)
        cpx_record_ids = cpx_record_ids.union(cpx.record_ids)

        if cpx.svtype == 'UNR':
            for i, record in enumerate(cpx.records):
                record.info['EVENT'] = 'UNRESOLVED_{0}'.format(unresolved_idx)
                cpx_records.append(record)
            unresolved_idx += 1

        else:
            cpx.vcf_record.id = variant_prefix + str(resolved_idx)
            cpx_records.append(cpx.vcf_record)
            resolved_idx += 1

    # Output all variants
    vcf.reset()

    for record in _merge_records(vcf, cpx_records, cpx_record_ids):
        yield record


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtools link-cpx',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Breakpoint VCFs.')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        help='Resolved complex variants.')
    parser.add_argument('unresolved', type=argparse.FileType('w'),
                        help='Unresolved complex breakpoints.')
    parser.add_argument('-p', '--prefix', default='CPX_',
                        help='Variant prefix [CPX_]')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)
    for line in CPX_INFO:
        vcf.header.add_line(line)

    resolved_f = pysam.VariantFile(args.fout, 'w', header=vcf.header)
    unresolved_f = pysam.VariantFile(args.unresolved, 'w', header=vcf.header)

    for record in resolve_complex_sv(vcf):
        if record.info['SVTYPE'] == 'UNR':
            unresolved_f.write(record)
        else:
            resolved_f.write(record)

    resolved_f.close()
    unresolved_f.close()


if __name__ == '__main__':
    main(sys.argv[1:])
