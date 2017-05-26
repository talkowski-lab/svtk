#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Count the instances of each SVTYPE observed in each sample in a VCF.

A record is only counted towards a sample's total if the sample received a
non-reference (0/0) or a non-null (./.) call. Records without an annotated
SVTYPE in the INFO field are counted as "NO_SVTYPE".

Counts are reported in a three column table (sample, svtype, count).
"""

import argparse
import sys
from collections import defaultdict
import pandas as pd
from pysam import VariantFile
from svtools.utils import NULL_GT


def count_svtypes(vcf):
    """
    Count instances of each SVTYPE in each sample in a VCF.

    Parameters
    ----------
    vcf : pysam.VariantFile

    Returns
    -------
    counts : pd.DataFrame
        Columns: sample, svtype, count
    """

    samples = list(vcf.header.samples)

    # Initialize counts per sample - each dict is keyed on svtype
    count_dict = {}
    for sample in samples:
        count_dict[sample] = defaultdict(int)

    for record in vcf:
        for sample in samples:
            # Only count called variants
            gt = record.samples[sample]['GT']
            if gt in NULL_GT:
                continue

            # Count the SVTYPE if it's present, otherwise increment NO_SVTYPE
            if 'SVTYPE' in record.info.keys():
                count_dict[sample][record.info['SVTYPE']] += 1
            else:
                count_dict[sample]['NO_SVTYPE'] += 1

    # Convert to dataframe, adding zeros to samples with no instances of a
    # given svtype
    counts = pd.DataFrame.from_dict(count_dict, orient='index')\
                         .fillna(0).astype(int)\
                         .reset_index().rename(columns={'index': 'sample'})

    # Tidy data from "column-per-svtype" format
    counts = pd.melt(counts, id_vars=['sample'],
                     var_name='svtype', value_name='count')

    return counts


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('--no-header', action='store_true', default=False,
                        help="Don't include header in output")
    parser.add_argument('fout', type=argparse.FileType('w'), nargs='?',
                        default=sys.stdout, help='Output file [stdout]')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = VariantFile(args.vcf)
    counts = count_svtypes(vcf)

    header = not args.no_header
    counts.to_csv(args.fout, sep='\t', index=False, header=header)


if __name__ == '__main__':
    main(sys.argv[1:])
