#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Convert a VCF to a BED.
"""

import argparse
import sys
import pysam
import svtools.utils as svu


def vcf2bed(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='VCF to convert.')
    parser.add_argument('bed', help='Converted bed. Specify `-` or `stdout` '
                        'to write to stdout.')
    parser.add_argument('--no-samples', dest='include_samples',
                        action='store_false', default=True,
                        help='Don\'t include comma-delimited list of called '
                        'samples for each variant.')
    parser.add_argument('-i', '--info', action='append',
                        help='INFO field to include as column in output. '
                        'May be specified more than once. To include all INFO '
                        'fields, specify `--info ALL`. INFO fields are '
                        'reported in the order in which they are requested. '
                        'If ALL INFO fields are requested, they are reported '
                        'in the order in which they appear in the VCF header.')
    parser.add_argument('--split-bnd', action='store_true', default=False,
                        help='Report two entries in bed file for each BND.')
    parser.add_argument('--no-header', dest='header', action='store_false',
                        default=True, help='Suppress header.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)

    header = '#chrom start end name svtype'.split()
    if args.include_samples:
        header.append('samples')
    if args.info:
        if 'ALL' in args.info:
            header = header + vcf.header.info.keys()
        else:
            header = header + args.info
    header = '\t'.join(header)

    bt = svu.vcf2bedtool(vcf,
                         split_bnd=args.split_bnd,
                         include_samples=args.include_samples,
                         include_strands=False,
                         include_infos=args.info,
                         annotate_ins=False,
                         report_alt=True)

    if args.bed in 'stdout -'.split():
        if args.header:
            sys.stdout.write(header + '\n')
        sys.stdout.write(str(bt))
    else:
        if args.header:
            bt.saveas(args.bed, trackline=header)
        else:
            bt.saveas(args.bed)
