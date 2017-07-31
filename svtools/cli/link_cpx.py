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
import pysam
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


def resolve_complex_sv(vcf):
    """
    Resolve complex SV from CNV intervals and BCA breakpoints.

    Yields all resolved events, simple or complex, in sorted order.

    Parameters
    ----------
    vcf : pysam.VariantFile

    Yields
    ------
    sv : pysam.VariantRecord
    """

    #


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

    clusters = link_cpx(vcf)

    resolved_idx = unresolved_idx = 1

    for cluster in clusters:
        cpx = ComplexSV(cluster)

        if cpx.svtype == 'UNR':
            for i, record in enumerate(cpx.records):
                record.info['EVENT'] = 'UNRESOLVED_{0}'.format(unresolved_idx)
                unresolved_f.write(record)
            unresolved_idx += 1

        else:
            cpx.vcf_record.id = args.prefix + str(resolved_idx)
            resolved_f.write(cpx.vcf_record)
            resolved_idx += 1


if __name__ == '__main__':
    main(sys.argv[1:])
