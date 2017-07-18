#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 ec2-user <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
from collections import defaultdict
import numpy as np
import scipy.stats as ss
import pandas as pd
import pysam
import svtools.pesr as pesr


class _DiscPair:
    def __init__(self, chrA, posA, strandA, chrB, posB, strandB, sample):
        self.chrA = chrA
        self.posA = int(posA)
        self.strandA = strandA
        self.chrB = chrB
        self.posB = int(posB)
        self.strandB = strandB
        self.sample = sample


class PEBreakpoint(pesr.Breakpoint):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.strandA, self.strandB = self.strands

    def pe_test(self, samples, discfile, n_background, window_in, window_out):
        self.choose_background(samples, n_background)

        self.load_counts(discfile, window_in, window_out)
        if self.pair_counts.shape[0] == 0:
            self.null_score()
            return

        self.process_counts()
        self.test_counts()

    def load_counts(self, discfile, window_in, window_out):
        reg = '{0}:{1}-{2}'

        def _get_coords(pos, strand):
            if strand == '+':
                start, end = pos - window_out, pos + window_in
            else:
                start, end = pos - window_in, pos + window_out
            return start, end

        startA, endA = _get_coords(self.posA, self.strandA)
        startB, endB = _get_coords(self.posB, self.strandB)

        region = reg.format(self.chrA, startA, endA)

        counts = defaultdict(int)
        pairs = discfile.fetch(region=region, parser=pysam.asTuple())

        for pair in pairs:
            pair = _DiscPair(*pair)

            # Pairs were selected based on window around chrA; check chrB
            if pair.chrB != self.chrB:
                continue
            if not (startB <= pair.posB < endB):
                continue

            # Require pairs match breakpoint strand
            if pair.strandA != self.strandA or pair.strandB != self.strandB:
                continue

            counts[pair.sample] += 1

        self.pair_counts = pd.DataFrame.from_dict({'count': counts})

    def process_counts(self):
        counts = self.pair_counts

        # Restrict to called or background samples
        samples = self.samples + self.background
        counts = counts.reindex(samples).fillna(0).astype(int)
        counts = counts.reset_index().rename(columns={'index': 'sample'})

        # Label samples with background
        is_called = counts['sample'].isin(self.samples)
        if is_called.any():
            counts.loc[is_called, 'call_status'] = 'called'
        if (~is_called).any():
            counts.loc[~is_called, 'call_status'] = 'background'

        self.pair_counts = counts

    def test_counts(self):
        statuses = 'called background'.split()
        medians = self.pair_counts.groupby('call_status')['count'].median()
        medians = medians.reindex(statuses).fillna(0).astype(int)

        pval = ss.poisson.cdf(medians.background, medians.called)

        stats = [medians.called, medians.background, -np.log10(pval)]
        columns = 'called_median bg_median log_pval'.split()
        self.stats = pd.DataFrame([stats], columns=columns)
        self.stats.log_pval = self.stats.log_pval.abs()

    def null_score(self):
        columns = 'called_median bg_median log_pval'.split()
        self.stats = pd.DataFrame([[0, 0, 0]], columns=columns)



class PETest:
    def __init__(self, variants, discfile,
                 window_in=50, window_out=500, n_background=160):
        """
        variants : pysam.VariantFile
        discfile : pysam.TabixFile
            chrA, posA, strandA, chrB, posB, strandB, sample
        window : int, optional
            Window around variant start/end to query for discordant pairs
        """

        self.variants = variants
        self.discfile = discfile
        self.window_in = window_in
        self.window_out = window_out
        self.samples = list(self.variants.header.samples)
        self.n_background = n_background

    def run(self):
        def _strand_check(record):
            return ('STRANDS' in record.info.keys() and
                    record.info['STRANDS'] in '++ +- -+ --'.split())

        for record in self.variants:
            # Skip non-stranded variants (e.g. WHAM inversions)
            # TODO: log skipped records
            if not _strand_check(record):
                continue

            breakpoint = PEBreakpoint.from_vcf(record)

            breakpoint.pe_test(self.samples, self.discfile, self.n_background,
                               self.window_in, self.window_out)

            stats = breakpoint.stats
            stats['name'] = record.id

            cols = 'name log_pval called_median bg_median'.split()
            yield stats[cols]

    def count_pairs(self, record):
        pairs = self.fetch_pairs(record)

        counts = defaultdict(int)
        for pair in pairs:
            counts[pair.sample] += 1

        counts = pd.DataFrame.from_dict({'counts': counts})
        return counts


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('variants', help='Variants (default=VCF).')
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

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.variants in ['-', 'stdin']:
        variantfile = pysam.VariantFile(sys.stdin)
    else:
        variantfile = pysam.VariantFile(args.variants)

    discfile = pysam.Tabixfile(args.disc)

    petest = PETest(variantfile, discfile, args.window_in, args.window_out,
                    args.background)

    header = 'name log_pval called_median bg_median'.split()
    args.fout.write('\t'.join(header) + '\n')

    for stats in petest.run():
        stats.to_csv(args.fout, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main(sys.argv[1:])
