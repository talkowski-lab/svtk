#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import itertools
from collections import deque
import numpy as np
import scipy.sparse as sps
import pysam
import pybedtools as pbt
import natsort
import svtools.utils as svu
from .cpx_inv import classify_complex_inversion


def samples_overlap(recA, recB, upper_thresh=0.8, lower_thresh=0.5):
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


def extract_breakpoints(vcfpath, bkpt_idxs):
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

    vcf = pysam.VariantFile(vcfpath)
    n_bkpts = len(bkpt_idxs)
    bkpts = np.empty(n_bkpts, dtype=object)

    for record in vcf:
        idx = bkpt_idxs.get(record.id)
        if idx is not None:
            bkpts[idx] = record

    return bkpts


def vcf2bedtool(vcfpath):
    """
    Wrap VCF as a bedtool. Necessary as pybedtools does not support SV in VCF.

    Parameters
    ----------
    vcfpath : str
        File path to VCF

    Returns
    -------
    bt : pybedtools.BedTool
        SV converted to Bedtool. Ends of BND records are assigned as pos + 1.
        Included columns: chrom, start, end, name, svtype, strands
    """

    vcf = pysam.VariantFile(vcfpath)

    # Convert each record in vcf to bed entry
    def _converter():
        bed = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'
        for record in vcf:
            if record.info['SVTYPE'] == 'BND':
                end = record.pos + 1
            else:
                end = record.info['END']
            yield bed.format(record.chrom, record.pos, end, record.id,
                             record.info['SVTYPE'], record.info['STRANDS'])

    return pbt.BedTool(_converter()).saveas()


def resolve_cpx(cluster):
    """
    Resolve complex variant structure from a cluster of breakpoints.

    Parameters
    ----------
    cluster : list of pysam.VariantRecord

    Returns
    -------
    variant : pysam.VariantRecord
    """

    inversions = [rec for rec in cluster if rec.info['SVTYPE'] == 'INV']
    tlocs = [rec for rec in cluster if rec.info['SVTYPE'] == 'BND']

    cnvtypes = 'DEL DUP'.split()
    cnvs = [rec for rec in cluster if rec.info['SVTYPE'] in cnvtypes]

    # Restrict to double-ended inversion events with appropriate strand pairing
    if len(inversions) == 0:
        if len(tlocs) > 0:
            cluster_type = 'INTERCHROMOSOMAL'
        else:
            cluster_type = 'ERROR_CNV_ONLY'
    elif len(inversions) == 1:
        cluster_type = 'SINGLE_ENDER'
    elif len(inversions) > 2:
        cluster_type = 'COMPLEX_3plus'
    elif inversions[0].info['STRANDS'] == inversions[1].info['STRANDS']:
        cluster_type = 'MATCHED_STRANDS'
    elif len(tlocs) > 0:
        cluster_type = 'MIXED_INV_TLOC'
    else:
        cluster_type = 'CANDIDATE'

    if cluster_type != 'CANDIDATE':
        return cluster_type, make_unresolved_entry(cluster, cluster_type)

    # Assign stranded breakpoints
    if inversions[0].info['STRANDS'] == '++':
        FF, RR = inversions
    else:
        RR, FF = inversions

    svtype = classify_complex_inversion(FF, RR, cnvs)

    if svtype != 'UNK':
        return svtype, make_bed_entry(FF, RR, cnvs, svtype)
    else:
        return svtype, make_unresolved_entry(cluster, svtype)


def make_unresolved_entry(cluster, cluster_type):
    """
    Parameters
    ----------
    cluster : list of pysam.VariantRecord
    """

    entry = ('{chrom}\t{start}\t{end}\t{ID}\t{svtype}\t{strands}\t'
             '{{name}}\t{cluster_type}\t{samples}\n')
    entries = []

    for record in cluster:
        chrom = record.chrom
        start, end = record.pos, record.info['END']
        ID = record.id
        svtype = record.info['SVTYPE']
        strands = record.info['STRANDS']
        samples = ','.join(svu.get_called_samples(record))

        entries.append(entry.format(**locals()))

    return ''.join(entries)


def make_bed_entry(FF, RR, cnvs, svtype):
    start = min(FF.pos, RR.pos)
    end = max(FF.info['END'], RR.info['END'])
    chrom = FF.chrom

    FF_samples = svu.get_called_samples(FF)
    RR_samples = svu.get_called_samples(RR)
    samples = sorted(set().union(FF_samples, RR_samples))

    FF_name = FF.id
    RR_name = RR.id
    if len(cnvs) == 0:
        CNV_names = '.'
    else:
        CNV_names = ','.join([cnv.id for cnv in cnvs])

    intervals = []
    inv_start = RR.pos
    inv_end = FF.info['END']
    intervals = ['INV_{chrom}:{inv_start}-{inv_end}'.format(**locals())]

    interval = '{cnv_type}_{chrom}:{cnv_start}-{cnv_end}'
    if svtype.startswith('del'):
        cnv_start = FF.pos
        cnv_end = RR.pos
        cnv_type = 'DEL5'
        intervals.append(interval.format(**locals()))

    if svtype.endswith('del'):
        cnv_start = FF.info['END']
        cnv_end = RR.info['END']
        cnv_type = 'DEL3'
        intervals.append(interval.format(**locals()))

    if svtype.startswith('dup'):
        cnv_start = RR.pos
        cnv_end = FF.pos
        cnv_type = 'DUP5'
        intervals.append(interval.format(**locals()))

    if svtype.endswith('dup'):
        cnv_start = RR.info['END']
        cnv_end = FF.info['END']
        cnv_type = 'DUP3'
        intervals.append(interval.format(**locals()))

    intervals = ';'.join(intervals)
    samples = ','.join(samples)

    entry = ('{chrom}\t{start}\t{end}\t{{name}}\t{svtype}\t'
             '{intervals}\t{FF_name}\t{RR_name}\t{CNV_names}\t{samples}\n')
    entry = entry.format(**locals())

    return entry


def link_cpx(vcfpath, bkpt_window=100):
    """
    Parameters
    ----------
    vcfpath : str
        Path to breakpoint VCF
    """

    bt = vcf2bedtool(vcfpath)

    # Identify breakpoints which overlap within specified window
    overlap = bt.window(bt, w=bkpt_window).saveas()

    # Exclude self-hits
    overlap = overlap.filter(lambda b: b.fields[3] != b.fields[9]).saveas()

    # Restrict to overlaps involving a BCA breakpoint
    cnvtypes = 'DEL DUP'.split()
    overlap = overlap.filter(lambda b: b.fields[4] not in cnvtypes).saveas()

    # Get linked variant IDs
    links = [(b[3], b[9]) for b in overlap.intervals]
    linked_IDs = natsort.natsorted(set(itertools.chain.from_iterable(links)))
    linked_IDs = np.array(linked_IDs)

    # Map variant IDs to indices
    bkpt_idxs = {ID: i for i, ID in enumerate(linked_IDs)}
    indexed_links = np.array([(bkpt_idxs[a], bkpt_idxs[b]) for a, b in links])

    # Extract VariantRecords corresponding to breakpoints
    n_bkpts = len(linked_IDs)
    bkpts = extract_breakpoints(vcfpath, bkpt_idxs)

    # Build sparse graph from links
    G = sps.eye(n_bkpts, dtype=np.uint16, format='lil')
    for i, j in indexed_links:
        if samples_overlap(bkpts[i], bkpts[j]):
            G[i, j] = 1

    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])

    # Remove clusters of one variant - leftover from shared sample filtering
    clusters = [c for c in clusters if len(c) > 1]

    return clusters
