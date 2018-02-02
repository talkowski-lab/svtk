#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import pandas as pd
import numpy as np
from svtk.adjudicate import rf_classify, labelers


ALLOSOMES = 'X Y chrX chrY'.split()


def adjudicate_BAF(metrics, labeler, name):
    # Deletions
    testable = metrics.loc[(metrics.svtype == 'DEL') &
                           (metrics.svsize >= 5000)]
    trainable = testable.loc[(testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]
    features = 'BAF_snp_ratio BAF_del_loglik'.split()
    cutoffs = {'indep': ['BAF_snp_ratio'], 'dep': ['BAF_del_loglik']}

    del_cutoffs = rf_classify(metrics, trainable, testable, features,
                                       labeler, cutoffs, name)

    # Duplications
    testable = metrics.loc[(metrics.svtype == 'DUP') &
                           (metrics.svsize >= 5000)]
    trainable = testable.loc[(testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]
    features = 'BAF_KS_stat BAF_KS_log_pval'.split()
    cutoffs = {'indep': ['BAF_KS_stat'], 'dep': ['BAF_KS_log_pval']}

    dup_cutoffs = rf_classify(metrics, trainable, testable, features,
                                       labeler, cutoffs, name)

    # Combine cutoffs
    del_cutoffs['svtype'] = 'DEL'
    dup_cutoffs['svtype'] = 'DUP'
    cutoffs = pd.concat([del_cutoffs, dup_cutoffs]).reset_index()
    cutoffs['test'] = 'BAF'
    cutoffs['max_svsize'] = np.nan
    cutoffs['min_svsize'] = 5000
    cutoffs['algtype'] = 'any'

    return cutoffs


def adjudicate_BAF1(metrics):
    cutoffs = adjudicate_BAF(metrics, labelers.BAF1TrainingLabeler(), 'BAF1_prob')
    cutoffs['test'] = 'BAF1'

    return cutoffs


def adjudicate_BAF2(metrics):
    cutoffs = adjudicate_BAF(metrics, labelers.BAF2TrainingLabeler(), 'BAF2_prob')
    cutoffs['test'] = 'BAF2'
    return cutoffs


def adjudicate_SR1(metrics):
    testable = metrics.loc[~metrics.name.str.contains('_depth_')]
    trainable = testable.loc[(testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]
    features = ['SR_sum_log_pval', 'SR_sum_called_median', 'SR_sum_bg_median']
    cutoffs = {'indep': ['SR_sum_log_pval'], 'dep': []}
    labeler = labelers.SR1TrainingLabeler()

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoffs, 'SR1_prob')

    cutoffs['test'] = 'SR1'
    cutoffs['svtype'] = 'CNV'
    cutoffs['algtype'] = 'PESR'

    return cutoffs


def adjudicate_RD(metrics):
    features = ["RD_Median_Separation", "RD_log_pval", "RD_log_2ndMaxP"]
    cutoff_features = {'indep': ['RD_log_pval', 'RD_Median_Separation'], 
                       'dep': ['RD_log_2ndMaxP']}
    labeler = labelers.RDTrainingLabeler()
    cutoff_dfs = []

    # PE/SR >1 kb
    testable = metrics.loc[~metrics.name.str.contains('_depth_') &
                           (metrics.svsize >= 1000)]
    trainable = testable.loc[(testable.svsize >= 5000) &
                             (testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoff_features, 'RD_prob')

    cutoff_dfs.append(cutoffs)
    cutoff_dfs[0]['algtype'] = 'PESR'
    cutoff_dfs[0]['max_svsize'] = np.nan
    cutoff_dfs[0]['min_svsize'] = 1000
    cutoff_dfs[0]['svtype'] = 'CNV'

    # PE/SR <1 kb
    testable = metrics.loc[~metrics.name.str.contains('_depth_') &
                           (metrics.svsize < 1000)]
    trainable = testable.loc[(testable.svsize >= 100) &
                             (testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoff_features, 'RD_prob')

    cutoff_dfs.append(cutoffs)
    cutoff_dfs[1]['algtype'] = 'PESR'
    cutoff_dfs[1]['max_svsize'] = 1000
    cutoff_dfs[1]['min_svsize'] = 0
    cutoff_dfs[1]['svtype'] = 'CNV'

    # Depth dels
    cutoff_features = {'indep': ['RD_log_pval', 'RD_Median_Separation'], 
                       'dep': []}
    testable = metrics.loc[metrics.name.str.contains('_depth_') &
                           (metrics.svtype == 'DEL')]
    trainable = testable.loc[(testable.svsize >= 5000) &
                             (testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoff_features, 'RD_prob',
                                   clean_cutoffs=True)
    
    cutoff_dfs.append(cutoffs)
    cutoff_dfs[2]['algtype'] = 'Depth'
    cutoff_dfs[2]['max_svsize'] = np.nan
    cutoff_dfs[2]['min_svsize'] = 5000
    cutoff_dfs[2]['svtype'] = 'DEL'
    
    # Depth dups
    testable = metrics.loc[metrics.name.str.contains('_depth_') &
                           (metrics.svtype == 'DUP')]
    trainable = testable.loc[(testable.svsize >= 5000) &
                             (testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoff_features, 'RD_prob',
                                   clean_cutoffs=True)
    
    cutoff_dfs.append(cutoffs)
    cutoff_dfs[3]['algtype'] = 'Depth'
    cutoff_dfs[3]['max_svsize'] = np.nan
    cutoff_dfs[3]['min_svsize'] = 5000
    cutoff_dfs[3]['svtype'] = 'DUP'

    # Fail depth-only below 5 kb
    metrics.loc[metrics.name.str.contains('_depth_') &
                (metrics.svsize < 5000) &
                (metrics.RD_prob >= 0.5), 'RD_prob'] = 0.499

    cutoffs = pd.concat(cutoff_dfs)
    cutoffs['test'] = 'RD'

    return cutoffs


def adjudicate_PE(metrics):
    testable = metrics.loc[~metrics.name.str.contains('_depth_')]
    trainable = testable.loc[(testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]
    features = ['PE_log_pval', 'PE_called_median', 'PE_bg_median']
    cutoffs = {'indep': ['PE_log_pval'], 'dep': []}
    labeler = labelers.PETrainingLabeler()

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoffs, 'PE_prob')

    cutoffs['test'] = 'PE'
    cutoffs['svtype'] = 'CNV'
    cutoffs['algtype'] = 'PESR'

    return cutoffs


def adjudicate_SR2(metrics):
    testable = metrics.loc[~metrics.name.str.contains('_depth_')]
    trainable = testable.loc[(testable.svsize >= 5000) &
                             (testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]
    features = ['SR_sum_log_pval', 'SR_sum_called_median', 'SR_sum_bg_median']
    cutoffs = {'indep': ['SR_sum_log_pval'], 'dep': []}
    labeler = labelers.SR2TrainingLabeler()

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoffs, 'SR2_prob')

    cutoffs['test'] = 'SR2'
    cutoffs['svtype'] = 'CNV'
    cutoffs['algtype'] = 'PESR'

    return cutoffs

def adjudicate_PESR(metrics):
    testable = metrics.loc[~metrics.name.str.contains('_depth_')]
    trainable = testable.loc[(testable.poor_region_cov < 0.3) &
                             ~testable.chrom.isin(ALLOSOMES)]
    features = ['PESR_log_pval', 'PESR_called_median', 'PESR_bg_median']
    cutoffs = {'indep': ['PESR_log_pval'], 'dep': []}
    labeler = labelers.PESRTrainingLabeler()

    cutoffs = rf_classify(metrics, trainable, testable, features,
                                   labeler, cutoffs, 'PESR_prob')

    cutoffs['test'] = 'PESR'
    cutoffs['svtype'] = 'CNV'
    cutoffs['algtype'] = 'PESR'

    return cutoffs


def consolidate_score(metrics, cutoffs):
    """Assign final score to each variant"""

    # Score PE/SR <1 kb
    lt1kb = ~metrics.name.str.contains('depth') & (metrics.svsize < 1000)
    pesr_cols = 'PE_prob SR2_prob'.split()
    metrics.loc[lt1kb, 'score'] = metrics[pesr_cols].max(axis=1)

    # Score depth
    depth = metrics.name.str.contains('depth')
    metrics.loc[depth, 'score'] = metrics.RD_prob

    # Score PE/SR >1 kb
    gt1kb = ~metrics.name.str.contains('depth') & (metrics.svsize >= 1000)
    PESR_pass = (metrics.PE_prob >= 0.5) | (metrics.SR2_prob >= 0.5)
    RD_pass = (metrics.RD_prob >= 0.5)
    prob_cols = 'PE_prob SR2_prob RD_prob'.split()

    # If variants pass both RD + PE/SR or fail both, use max across the three
    metrics.loc[gt1kb & ~(PESR_pass ^ RD_pass), 'score'] = metrics[prob_cols].max(axis=1)

    # If variants pass PE/SR but not RD, use max PE/SR and reclassify as BND
    metrics.loc[gt1kb & PESR_pass & ~RD_pass, 'score'] = metrics[pesr_cols].max(axis=1)
    metrics.loc[gt1kb & PESR_pass & ~RD_pass, 'svtype'] = 'BND'

    # If variants pass RD but not PE/SR, pass if they pass depth-based cutoffs
    depth_dels = gt1kb & ~PESR_pass & RD_pass & (metrics.svtype == 'DEL')

    depth_cutoffs = cutoffs.loc[(cutoffs.algtype == 'Depth') &
                                (cutoffs.test == 'RD') &
                                (cutoffs.svtype == 'DEL')]
    depth_pass = metrics.svsize >= depth_cutoffs.min_svsize[0]
    for idx, row in depth_cutoffs.iterrows():
        depth_pass = depth_pass & metrics[row['metric']] >= row['cutoff']

    metrics.loc[depth_dels & depth_pass, 'score'] = metrics.RD_prob
    metrics.loc[depth_dels & ~depth_pass, 'score'] = 0.495
    
    # dups
    depth_dups = gt1kb & ~PESR_pass & RD_pass & (metrics.svtype == 'DUP')

    depth_cutoffs = cutoffs.loc[(cutoffs.algtype == 'Depth') &
                                (cutoffs.test == 'RD') &
                                (cutoffs.svtype == 'DUP')]
    depth_pass = metrics.svsize >= depth_cutoffs.min_svsize[0]
    for idx, row in depth_cutoffs.iterrows():
        depth_pass = depth_pass & metrics[row['metric']] >= row['cutoff']

    metrics.loc[depth_dups & depth_pass, 'score'] = metrics.RD_prob
    metrics.loc[depth_dups & ~depth_pass, 'score'] = 0.495

    probs = 'BAF2_prob SR2_prob RD_prob PE_prob PESR_prob'.split()
    return metrics['name svtype score'.split() + probs].copy()


def adjudicate_SV(metrics):
    if 'chrom' not in metrics.columns:
        metrics['chrom'] = metrics.name.str.split('_').str[-2]
    cutoffs = np.empty(7, dtype=object)
    cutoffs[0] = adjudicate_BAF1(metrics)
    cutoffs[1] = adjudicate_SR1(metrics)
    cutoffs[2] = adjudicate_RD(metrics)
    cutoffs[3] = adjudicate_PE(metrics)
    cutoffs[4] = adjudicate_BAF2(metrics)
    cutoffs[5] = adjudicate_SR2(metrics)
    cutoffs[6] = adjudicate_PESR(metrics)

    cutoffs = pd.concat(cutoffs)

    scores = consolidate_score(metrics, cutoffs)

    return scores, cutoffs
