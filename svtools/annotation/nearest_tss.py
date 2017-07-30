#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import warnings
import pybedtools as pbt
import svtools.utils as svu


def annotate_nearest_tss(sv, gencode):
    def _make_tss(feature):
        feature.end = feature.start + 1
        return feature

    transcripts = gencode.filter(lambda r: r.fields[7] == 'transcript')
    tss = transcripts.each(_make_tss).saveas()

    nearest_tss = sv.sort().closest(tss.sort(), D='b', id=True).saveas()

    # Suppress warning about column names
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        nearest_tss = nearest_tss.to_dataframe()

    nearest_tss = nearest_tss[[3, 4, 9]].copy()
    nearest_tss.columns = 'name svtype gene_id'.split()
    nearest_tss['effect'] = 'NEAREST_TSS'

    return nearest_tss


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Structural variants.')
    parser.add_argument('gencode_annotation', help='Gencode annotation bed.')
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    sv = svu.vcf2bedtool(args.vcf)
    gencode = pbt.BedTool(args.gencode_annotation)

    nearest_tss = annotate_nearest_tss(sv, gencode)
    nearest_tss.to_csv(args.fout, index=False, sep='\t')


if __name__ == '__main__':
    main()
