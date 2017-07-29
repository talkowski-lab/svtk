#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import pandas as pd


def classify_del(disrupt_dict):
    """All deletion hits are DISRUPTING"""
    regions = list(disrupt_dict.keys())
    if 'exon' in regions:
        return 'LOF'
    if 'UTR' in regions:
        return 'UTR'
    if 'gene' in regions:
        return 'INTRONIC'
    if 'promoter' in regions:
        return 'promoter'

    return 'no_effect'


def classify_dup(disrupt_dict):
    elements = disrupt_dict.keys()
    if 'exon' in elements:
        # duplication internal to exon
        if 'BOTH-INSIDE' in disrupt_dict['exon']:
            return 'LOF'

        # duplication internal to gene, disrupting exon, is LoF
        # duplication spanning gene is copy gain
        # duplication overlapping gene is no effect (one good copy left)
        if 'BOTH-INSIDE' in disrupt_dict['gene']:
            return 'LOF'
        elif 'SPAN' in disrupt_dict['gene']:
            return 'COPY_GAIN'
        else:
            return 'no_effect'

    if 'UTR' in elements:
        if 'BOTH-INSIDE' in disrupt_dict['UTR']:
            return 'UTR'

    if 'gene' in elements:
        if 'BOTH-INSIDE' in disrupt_dict['gene']:
            return 'INTRONIC'
        else:
            return 'GENE_OTHER'

    if 'promoter' in elements:
        if 'BOTH-INSIDE' in disrupt_dict['promoter']:
            return 'promoter'

    return 'no_effect'


def classify_inv(disrupt_dict):
    elements = disrupt_dict.keys()
    if 'exon' in elements:
        # breakpoint disrupts exon -> LoF
        if ('BOTH-INSIDE' in disrupt_dict['exon'] or
                'ONE-INSIDE' in disrupt_dict['exon']):
            return 'LOF'

        # if breakpoint spans exon but gene is disrupted -> LoF
        elif ('BOTH-INSIDE' in disrupt_dict['gene'] or
                'ONE-INSIDE' in disrupt_dict['gene']):
            return 'LOF'

        # inversion spanning gene -> no effect
        else:
            return 'no_effect'

    if 'UTR' in elements:
        if ('BOTH-INSIDE' in disrupt_dict['UTR'] or
                'ONE-INSIDE' in disrupt_dict['UTR']):
            return 'UTR'

    if 'gene' in elements:
        if 'BOTH-INSIDE' in disrupt_dict['gene']:
            return 'INTRONIC'
        else:
            return 'GENE_OTHER'

    if 'promoter' in elements:
        if ('BOTH-INSIDE' in disrupt_dict['promoter'] or
                'ONE-INSIDE' in disrupt_dict['promoter']):
            return 'promoter'

    return 'no_effect'


def classify_ctx(disrupt_dict):
    elements = disrupt_dict.keys()
    # BND hitting any part of gene is LOF
    if 'exon' in elements:
        return 'LOF'
    if 'gene' in elements:
        return 'LOF'
    if 'UTR' in elements:
        return 'UTR'
    if 'promoter' in elements:
        return 'promoter'
    return 'no_effect'


def classify_disrupt(disrupt_dict, svtype):
    """

    Exonic disruptions count as LOF.

    If exonic span or copy gain, check gene:
        gene disruptions count as LOF
        CDS copy_gain/span count as themselves

    If no exonic hit, check UTR:
        UTR disruption counts as UTR

    If no UTR, check promoter:
        disruption counts
    """

    if svtype == 'DEL':
        return classify_del(disrupt_dict)
    if svtype == 'DUP':
        return classify_dup(disrupt_dict)
    if svtype == 'INV':
        return classify_inv(disrupt_dict)
    if svtype == 'BND':
        return classify_ctx(disrupt_dict)


def classify_effect(hits):
    hits['element_hit'] = hits['element_type'] + '__' + hits['hit_type']
    effects = hits.groupby('name svtype gene_id'.split())['element_hit']\
                  .agg(lambda s: ','.join(sorted(set(s))))\
                  .reset_index()

    def _apply_classify(row):
        svtype = row.svtype
        element_hit = row.element_hit.split(',')
        element_hit = [h.split('__') for h in element_hit]
        element_hit = {h[0]: h[1] for h in element_hit}

        return classify_disrupt(element_hit, svtype)

    effects['effect'] = effects.apply(_apply_classify, axis=1)

    return effects


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gencode_hits',
                        help='Overlap with gencode annotations')
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    hits = pd.read_table(args.gencode_hits)
    effects = classify_effect(hits)
    effects.to_csv(args.fout, index=False, sep='\t')


if __name__ == '__main__':
    main()
