#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Link complex records
"""

import itertools
from collections import deque
import numpy as np
import scipy.sparse as sps
import natsort
import svtk.utils as svu


def samples_overlap(recA, recB, upper_thresh=0.5, lower_thresh=0.5):
    """
    Report if the samples called in two VCF records overlap sufficiently.

    The fraction of each record's samples which are shared with the other
    record is calculated. The record with a greater fraction of shared samples
    must exceed the upper threshold AND the record with a lesser fraction of
    shared samples must exceed the lower threshold. This is intended to
    maximize sensitivity in rare variants with a false negative in one
    breakpoint.

    Parameters
    ----------
    recA : pysam.VariantRecord
    recB : pysam.VariantRecord
    upper_thresh : float, optional
        Minimum sample overlap in record with greater overlap
    lower_thresh : float, optional
        Minimum sample overlap in record with lesser overlap

    Returns
    -------
    samples_overlap : bool
        Samples shared between records meet required thresholds.
    """

    # Get lists of called samples for each record
    samplesA = set(svu.get_called_samples(recA))
    samplesB = set(svu.get_called_samples(recB))

    # Compute fraction of each record's samples which are shared
    shared = samplesA & samplesB
    fracA = len(shared) / len(samplesA)
    fracB = len(shared) / len(samplesB)

    min_frac, max_frac = sorted([fracA, fracB])

    return min_frac >= lower_thresh and max_frac >= upper_thresh


def extract_breakpoints(vcf, bkpt_idxs):
    """
    Extract all VCF records in list of IDs
    (Assumes VCF is sorted by variant ID)

    Parameters
    ----------
    vcfpath : str
        Path to VCF
    bkpt_idxs : dict of {str : int}
        Mapping of variant IDs to array index

    Returns
    -------
    bkpts : list of pysam.VariantRecord
    """

    #  vcf = pysam.VariantFile(vcfpath)
    n_bkpts = len(bkpt_idxs)
    bkpts = np.empty(n_bkpts, dtype=object)

    for record in vcf:
        idx = bkpt_idxs.get(record.id)
        if idx is not None:
            bkpts[idx] = record

    return bkpts


def link_cpx(vcf, bkpt_window=300, cpx_dist=20000):
    """
    Parameters
    ----------
    vcfpath : str
        Path to breakpoint VCF
    """

    bt = svu.vcf2bedtool(vcf.filename, annotate_ins=False)

    # Identify breakpoints which overlap within specified window
    overlap = bt.window(bt, w=bkpt_window).saveas()

    # Exclude self-hits
    #  overlap = overlap.filter(lambda b: b.fields[3] != b.fields[9]).saveas()

    # # Restrict to overlaps involving a BCA breakpoint
    # cnvtypes = 'DEL DUP'.split()
    # overlap = overlap.filter(lambda b: b.fields[4] not in cnvtypes).saveas()

    # Get linked variant IDs
    links = [(b[3], b[9]) for b in overlap.intervals]
    linked_IDs = natsort.natsorted(set(itertools.chain.from_iterable(links)))
    linked_IDs = np.array(linked_IDs)

    # Map variant IDs to indices
    bkpt_idxs = {ID: i for i, ID in enumerate(linked_IDs)}
    indexed_links = np.array([(bkpt_idxs[a], bkpt_idxs[b]) for a, b in links])

    # Extract VariantRecords corresponding to breakpoints
    n_bkpts = len(linked_IDs)
    bkpts = extract_breakpoints(vcf, bkpt_idxs)

    # Exclude wildly disparate overlaps
    def close_enough(r1, r2):
        distA = np.abs(r1.pos - r2.pos)
        distB = np.abs(r1.stop - r2.stop)
        return distA < cpx_dist or distB < cpx_dist

    # Build sparse graph from links
    G = sps.eye(n_bkpts, dtype=np.uint16, format='lil')
    for i, j in indexed_links:
        if (samples_overlap(bkpts[i], bkpts[j]) and
                close_enough(bkpts[i], bkpts[j])):
            G[i, j] = 1

    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])

    # # Remove clusters of only CNV - leftover from shared sample filtering
    # def _ok_cluster(cluster):
    #     ok = any([record.info['SVTYPE'] not in cnvtypes for record in cluster])
    #     return ok

    clusters = [c for c in clusters if _ok_cluster(c)]
    #  clusters = [c for c in clusters if len(c) > 1]

    return clusters
