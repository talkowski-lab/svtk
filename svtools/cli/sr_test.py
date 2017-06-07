#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Calculate enrichment of clipped reads at SV breakpoints.
"""

import argparse
import sys
import io
from collections import deque
import numpy as np
import scipy.stats as ss
import pandas as pd
import pysam
import svtools.utils as svu


class _Breakpoint:
    def __init__(self, chrA, posA, chrB, posB, name, samples, strands):
        """
        Attributes
        ----------
        chrA : str
        posA : int
        chrB : str
        posB : str
        name : str
        samples : list of str
        strands : str
            +-,-+,++,--
        background : list of str
        split_counts : pd.DataFrame
        """
        self.chrA = chrA
        self.posA = posA
        self.chrB = chrB
        self.posB = posB
        self.name = name
        self.samples = samples
        self.strands = strands

        # Constructed later
        self.background = []
        self.split_counts = None

    @classmethod
    def from_vcf(cls, record):
        """
        Parameters
        ----------
        record : pysam.VariantRecord
        """

        chrA = record.chrom
        posA = record.pos
        chrB = record.info['CHR2']
        posB = record.info['END']

        name = record.id
        strands = record.info['STRANDS']

        samples = svu.get_called_samples(record)

        return cls(chrA, posA, chrB, posB, name, samples, strands)

    @classmethod
    def from_bed(cls, chrom, start, end, name, samples, svtype):
        chrA = chrom
        posA = int(start)
        chrB = chrom
        posB = int(end)

        name = name
        samples = samples.split(',')

        if svtype.upper() not in 'DEL DUP'.split():
            msg = 'Invalid SV type: {0} [expected del,dup,DEL,DUP]'
            msg = msg.format(svtype)
            raise Exception(msg)

        if svtype.upper() == 'DEL':
            strands = '+-'
        else:
            strands = '-+'

        return cls(chrA, posA, chrB, posB, name, samples, strands)

    def choose_background(self, samples, n_background):
        """
        Choose background samples

        Parameters
        ----------
        samples : list of str
            All available samples
        n_background : int
            Max number of background samples
        """

        self.background = [s for s in samples if s not in self.samples]
        if len(self.background) >= n_background:
            self.background = np.random.choice(self.background, n_background,
                                               replace=False).tolist()

    def load_counts(self, countfile, window):
        """
        Generate dataframe of split counts

        Parameters
        ----------
        countfile : pysam.TabixFile
            chrom pos clip count sample
        window : int
        """

        reg = '{0}:{1}-{2}'

        # Check if regions overlap so duplicate lines aren't fetched
        if (self.chrA == self.chrB) and (self.posB - self.posA <= 2 * window):
            regions = [
                reg.format(self.chrA, self.posA - window, self.posB + window)
            ]
        else:
            regions = [
                reg.format(self.chrA, self.posA - window, self.posA + window),
                reg.format(self.chrB, self.posB - window, self.posB + window),
            ]

        counts = deque()
        for region in regions:
            lines = countfile.fetch(region=region)
            lines = [l for l in lines]
            counts.append('\n'.join(lines))

        counts = io.StringIO('\n'.join(counts))

        cols = 'chrom pos clip count sample'.split()
        dtypes = dict(chrom=str, pos=int, clip=str, count=int, sample=str)
        self.split_counts = pd.read_table(counts, names=cols, dtype=dtypes)

    def process_counts(self, window):
        """
        Filter to called/background samples and assign coordinates

        Parameters
        ----------
        window : int
        """

        counts = self.split_counts

        # Restrict to called or background samples
        samples = self.samples + self.background
        counts = counts.loc[counts['sample'].isin(samples)].copy()

        # Determine whether clipped read supports start or end
        self.add_coords(counts)

        # Filter to within window
        # (excludes spurious left/right clips at other coord)
        counts = counts.loc[counts['dist'].abs() <= window].copy()

        self.split_counts = counts

        # Fill empty samples
        self.fill_counts()

        # Label samples with background
        is_called = self.split_counts['sample'].isin(self.samples)
        if is_called.any():
            self.split_counts.loc[is_called, 'call_status'] = 'called'
        if (~is_called).any():
            self.split_counts.loc[~is_called, 'call_status'] = 'background'

    def add_coords(self, df):
        df['posA'] = (df.pos - self.posA).abs()
        df['posB'] = (df.pos - self.posB).abs()

        strandA, strandB = self.strands

        # Map clip direction to strandedness
        clip_map = {'right': '+', 'left': '-'}
        df['strand'] = df['clip'].replace(clip_map)

        # If strands match, pick closest pos
        if strandA == strandB:
            cols = 'posA posB'.split()
            df['coord'] = df[cols].idxmin(axis=1)
        # Else match clip to position strand
        else:
            coord_map = {strandA: 'posA', strandB: 'posB'}
            df['coord'] = df['strand'].replace(coord_map)

        # Choose dist to identified coord
        df['dist'] = df.lookup(df.index, df.coord)

        if strandA == strandB:
            df.loc[df.strand != strandA, 'dist'] = np.inf

    def add_dists(self, df):
        # Get distance of clip position from variant start/end
        def _coord_dist(row):
            coord = row['coord']  # 'start' or 'end'
            pos = getattr(self, coord)
            return pos - row.pos
        df['dist'] = df.apply(_coord_dist, axis=1)

    def fill_counts(self):
        """
        Fill zeros in for samples with no observed splits
        """
        counts = self.split_counts
        samples = self.samples + self.background

        sub_dfs = []
        cols = 'sample pos count'.split()
        for coord in 'posA posB'.split():
            df = counts.loc[counts.coord == coord, cols]

            # Consider only positions found in a called sample
            pos = df.loc[df['sample'].isin(self.samples), 'pos']
            pos = pos.unique()

            idx = pd.MultiIndex.from_product(iterables=[samples, pos])

            df = df.set_index('sample pos'.split())
            df = df.reindex(idx).fillna(0).astype(int).reset_index()
            df = df.rename(columns=dict(level_0='sample', level_1='pos'))

            df['coord'] = coord
            sub_dfs.append(df)

        self.split_counts = pd.concat(sub_dfs)

    def ttest_counts(self):
        pvals = self.split_counts.groupby('coord pos'.split())\
                                 .apply(self.calc_ttest)
        cols = {
            0: 'called_mean',
            1: 'called_std',
            2: 'bg_mean',
            3: 'bg_std',
            4: 'called_n',
            5: 'bg_n',
            6: 'log_pval'
        }

        self.pvals = pvals.rename(columns=cols).reset_index()

    def choose_best_coords(self):
        # Pick coordinates with most significant enrichment
        max_pvals = self.pvals.groupby('coord')['log_pval'].max()\
                         .reset_index()\
                         .rename(columns={'log_pval': 'max_pval'})

        pvals = pd.merge(self.pvals, max_pvals, on='coord', how='left')
        pvals = pvals.loc[pvals.log_pval == pvals.max_pval].copy()
        pvals = pvals.drop('max_pval', axis=1)

        for coord in 'posA posB'.split():
            if coord not in pvals.coord.values:
                pvals = pd.concat([pvals, self.null_series(coord)])

        # Use distance as tiebreaker
        if pvals.shape[0] > 2:
            self.add_dists(pvals)
            closest = pvals.groupby('coord')['dist'].min().reset_index()\
                           .rename(columns={'dist': 'min_dist'})
            pvals = pd.merge(pvals, closest, on='coord', how='left')
            pvals = pvals.loc[pvals.dist == pvals.min_dist].copy()
            pvals = pvals.drop('dist min_dist'.split(), axis=1)

        pvals['name'] = self.name
        self.best_pvals = pvals

    def sr_test(self, samples, countfile, n_background, window):
        self.choose_background(samples, n_background)
        self.load_counts(countfile, window)
        self.process_counts(window)

        if self.split_counts.shape[0] == 0:
            self.null_score()
        else:
            self.ttest_counts()
            self.choose_best_coords()

    def null_score(self):
        cols = ('coord pos called_mean called_std bg_mean bg_std called_n '
                'bg_n log_pval name').split()
        self.best_pvals = pd.DataFrame([
            ['posA', 0, 0, 0, 0, 0, 0, 0, 0, self.name],
            ['posB', 0, 0, 0, 0, 0, 0, 0, 0, self.name]],
            columns=cols)

    def null_series(self, coord):
        cols = ('coord pos called_mean called_std bg_mean bg_std called_n '
                'bg_n log_pval name').split()
        return pd.DataFrame([[coord, 0, 0, 0, 0, 0, 0, 0, 0, self.name]],
                            columns=cols)

    @staticmethod
    def calc_ttest(df):
        def _summary_stats(series):
            return series.mean(), series.std(), series.shape[0]

        called = df.loc[df.call_status == 'called', 'count']
        called_mu, called_sigma, called_n = _summary_stats(called)

        background = df.loc[df.call_status == 'background', 'count']
        bg_mu, bg_sigma, bg_n = _summary_stats(background)

        if bg_n == 0:
            p = 1
        else:
            t, _p = ss.ttest_ind_from_stats(called_mu, called_sigma, called_n,
                                            bg_mu, bg_sigma, bg_n)

            # One-sided test
            p = _p / 2 if t > 0 else 1

        return pd.Series([called_mu, called_sigma, bg_mu, bg_sigma,
                          called_n, bg_n, -np.log10(p)])


def _BreakpointParser(variantfile, bed=False):
    """
    variantfile : pysam.Tabixfile or file
    bed : bool, optional
        variants is a bed
    """

    for record in variantfile:
        if bed:
            yield _Breakpoint.from_bed(*record.strip().split()[:6])
        else:
            yield _Breakpoint.from_vcf(record)


class SRTest():
    def __init__(self, variantfile, countfile, bed=False, samples=None,
                 window=100, n_background=160):
        """
        variantfile : file
            filepath of variants file
        countfile : pysam.TabixFile
            per-coordinate, per-sample split counts
            chrom pos clip count sample
        bed : bool, optional
            Variants file is a BED, not a VCF.
            BED columns: chrom start end name samples svtype
        samples : list of str, optional
            List of all samples to consider. By default, all samples in VCF
            header are considered. Required when specifying `bed=True`.
        window : int, optional
            Window around breakpoint to consider for split read enrichment
        n_background : int, optional
            Number of background samples to choose for comparison in t-test
        """

        if bed:
            if samples is None:
                msg = 'samples is required when providing calls in BED format.'
                raise ValueError(msg)
        else:
            variantfile = pysam.VariantFile(variantfile)
            samples = list(variantfile.header.samples)

        self.breakpoints = _BreakpointParser(variantfile, bed)

        self.countfile = countfile
        self.samples = sorted(samples)

        self.window = window
        self.n_background = n_background
        self.pvals = None

    def run(self):
        for breakpoint in self.breakpoints:
            breakpoint.sr_test(self.samples, self.countfile, self.n_background,
                               self.window)
            pvals = breakpoint.best_pvals
            cols = ('name coord pos log_pval called_mean called_std bg_mean '
                    'bg_std called_n bg_n').split()
            pvals = pvals[cols].fillna(0)

            for col in 'pos called_n bg_n'.split():
                pvals[col] = pvals[col].astype(int)

            yield pvals


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('variants',
                        help='VCF of variant calls. Standardized to include '
                        'CHR2, END, SVTYPE, STRANDS in INFO.')

    parser.add_argument('--bed', action='store_true', default=False,
                        help='Variants file is in bed format. First six cols: '
                        'chrom,start,end,name,samples,svtype')
    parser.add_argument('--samples', type=argparse.FileType('r'), default=None,
                        help='File listing all sample IDs. Parsed from VCF '
                        'header by default, required for --bed.')

    # TODO: permit direct querying around bams
    parser.add_argument('counts', help='Tabix indexed file of split counts. '
                        'Columns: chrom,pos,clip,count,sample')

    parser.add_argument('fout', type=argparse.FileType('w'),
                        help='Output table of most significant start/end'
                        'positions')

    parser.add_argument('-w', '--window', type=int, default=100,
                        help='Window around variant start/end to consider for '
                        'split read support.')
    parser.add_argument('-b', '--background', type=int, default=160,
                        help='Number of background samples to choose for '
                        'comparison in t-test [160]')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.samples is None:
        if args.bed:
            msg = '--samples is required when providing calls in BED format.'
            raise argparse.ArgumentError(msg)
    else:
        samples = [s.strip() for s in args.samples.readlines()]

    if args.variants in ['-', 'stdin']:
        variantfile = sys.stdin
    else:
        variantfile = open(args.variants)

    countfile = pysam.TabixFile(args.counts)

    srtest = SRTest(variantfile, countfile, args.bed, samples,
                    args.window, args.background)

    header = ('name coord pos log_pval called_mean called_std bg_mean '
              'bg_std called_n bg_n').split()
    args.fout.write('\t'.join(header) + '\n')

    for pvals in srtest.run():
        pvals.to_csv(args.fout, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main(sys.argv[1:])
