#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Calculate enrichment of clipped reads or discordant pairs at SV breakpoints.
"""

import argparse
import sys
import pysam
from svtools.pesr import SRTestRunner, PETestRunner


def sr_test(argv):
    parser = argparse.ArgumentParser(
        description="Calculate enrichment of clipped reads at SV breakpoints.",
        prog='svtools sr-test',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf',
                        help='VCF of variant calls. Standardized to include '
                        'CHR2, END, SVTYPE, STRANDS in INFO.')
    parser.add_argument('countfile', help='Tabix indexed file of split counts.'
                        ' Columns: chrom,pos,clip,count,sample')
    parser.add_argument('fout',
                        help='Output table of most significant start/end'
                        'positions.')

    parser.add_argument('-w', '--window', type=int, default=100,
                        help='Window around variant start/end to consider for '
                        'split read support. [100]')
    parser.add_argument('-b', '--background', type=int, default=160,
                        help='Number of background samples to choose for '
                        'comparison in t-test. [160]')
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None,
                        help='Whitelist of samples to restrict testing to.')
    # TODO: add normalization
    # parser.add_argument('--coverage-csv', default=None,
    #                     help='Median coverage statistics for each library '
    #                     '(optional). If provided, each sample\'s split '
    #                     'counts will be normalized accordingly. CSV: '
    #                     'sample,MEDIAN_COVERAGE,MAD_COVERAGE')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.vcf)
    countfile = pysam.TabixFile(args.counts)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    header = 'name coord pos log_pval called_median bg_median'.split()
    fout.write('\t'.join(header) + '\n')

    if args.samples is not None:
        whitelist = [s.strip() for s in args.samples.readlines()]
    else:
        whitelist = None

    runner = SRTestRunner(vcf, countfile, fout, args.background, args.window,
                          whitelist)
    runner.run()


def pe_test(argv):
    parser = argparse.ArgumentParser(
        description="Calculate enrichment of discordant pairs at SV breakpoints.",
        prog='svtools pe-test',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Variants.')
    parser.add_argument('disc', help='Table of discordant pair coordinates.')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        help='Output table of PE counts.')
    parser.add_argument('-o', '--window-out', type=int, default=500,
                        help='Window outside breakpoint to query for '
                        'discordant pairs. [500]')
    parser.add_argument('-i', '--window-in', type=int, default=50,
                        help='Window inside breakpoint to query for '
                        'discordant pairs. [50]')
    parser.add_argument('-b', '--background', type=int, default=160,
                        help='Number of background samples to sample for PE '
                        'evidence. [160]')
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None,
                        help='Whitelist of samples to restrict testing to.')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.variants in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.variants)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')
    header = 'name log_pval called_median bg_median'.split()
    args.fout.write('\t'.join(header) + '\n')

    if args.samples is not None:
        whitelist = [s.strip() for s in args.samples.readlines()]
    else:
        whitelist = None

    discfile = pysam.Tabixfile(args.disc)

    runner = PETestRunner(vcf, discfile, fout, args.background,
                          args.window_in, args.window_out, whitelist)
    runner.run()
