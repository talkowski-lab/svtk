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


def make_unresolved_entry(cluster, cluster_type):
    """
    Parameters
    ----------
    cluster : list of pysam.VariantRecord
    """

    entry = ('{chrom}\t{start}\t{end}\t{chr2}\t{ID}\t{svtype}\t{strands}\t'
             '{{name}}\t{cluster_type}\t{samples}\n')
    entries = []

    for record in cluster:
        chrom = record.chrom
        start, end = record.pos, record.info['END']
        ID = record.id
        svtype = record.info['SVTYPE']
        strands = record.info['STRANDS']
        chr2 = record.info['CHR2']
        samples = ','.join(svu.get_called_samples(record))

        entries.append(entry.format(**locals()))

    return ''.join(entries)


class ComplexSV:
    def __init__(self, records):
        """
        Parameters
        ----------
        records : list of pysam.VariantRecord
            Clustered records to resolve
        """

        self.records = records

        self.inversions = [r for r in records if r.info['SVTYPE'] == 'INV']
        self.tlocs = [r for r in records if r.info['SVTYPE'] == 'BND']

        cnvtypes = 'DEL DUP'.split()
        self.cnvs = [r for r in records if r.info['SVTYPE'] in cnvtypes]

        self.resolve()
        self.make_record()

    def resolve(self):
        self._set_cluster_type()

        if self.cluster_type == 'CANDIDATE_INVERSION':
            self.resolve_inversion()
        elif self.cluster_type == 'CANDIDATE_TRANSLOCATION':
            pass
        else:
            self.svtype = 'UNR'
            self.cpx_type = self.cluster_type

        pass

    def resolve_inversion(self):
        if self.inversions[0].info['STRANDS'] == '++':
            FF, RR = self.inversions
        else:
            RR, FF = self.inversions

        self.cpx_type = classify_complex_inversion(FF, RR, self.cnvs)
        if self.cpx_type == 'INV':
            self.svtype = 'INV'
        elif self.cpx_type == 'UNK':
            self.svtype = 'UNR'
        elif 'INS' in self.cpx_type:
            self.svtype = 'INS'
        else:
            self.svtype = 'CPX'

    def _set_cluster_type(self):
        # Restrict to double-ended inversion events with appropriate
        # strand pairing
        if len(self.inversions) == 0 and len(self.tlocs) == 0:
            self.cluster_type = 'ERROR_CNV_ONLY'
        elif len(self.inversions) == 1 and len(self.tlocs) == 0:
            self.cluster_type = 'SINGLE_ENDER_INV'
        elif len(self.inversions) == 0 and len(self.tlocs) == 1:
            self.cluster_type = 'SINGLE_ENDER_TLOC'
        elif len(self.inversions) == 1 and len(self.tlocs) == 1:
            self.cluster_type = 'MIXED_INV_TLOC'
        elif len(self.inversions) == 2 and len(self.tlocs) == 0:
            if (self.inversions[0].info['STRANDS'] ==
                    self.inversions[1].info['STRANDS']):
                self.cluster_type = 'MATCHED_STRANDS'
            else:
                self.cluster_type = 'CANDIDATE_INVERSION'
        elif len(self.inversions) == 0 and len(self.tlocs) == 2:
            self.cluster_type = 'CANDIDATE_TRANSLOCATION'
        elif len(self.inversions) + len(self.tlocs) >= 3:
            self.cluster_type = 'COMPLEX_3plus'
        else:
            self.cluster_type = 'ERROR_UNCLASSIFIED'

    def make_record(self):
        self.vcf_record = cluster[0].copy()
        called_samples = set(svu.get_called_samples(self.vcf_record))
    
        # Take union of called samples
        for record in cluster[1:]:
            cs = svu.get_called_samples(record)
            for sample in cs:
                if sample not in called_samples:
                    new_record.samples[sample]['GT'] = (0, 1)
                    called_samples.add(sample)
    
        new_record.alts = ('<{0}>'.format(svtype), )
        new_record.info['SVTYPE'] = svtype
    
        new_record.info['CPX_TYPE'] = cpx_type

    #  if cluster_type != 'CANDIDATE':
        #  return cluster_type, make_unresolved_entry(cluster, cluster_type)

    #  svtype = classify_complex_inversion(FF, RR, cnvs)

    #  if svtype != 'UNK':
        #  return svtype, make_bed_entry(FF, RR, cnvs, svtype)
    #  else:
        #  return svtype, make_unresolved_entry(cluster, svtype)




def make_vcf_record(cluster, svtype, cpx_type):
    """
    """






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

    bt = svu.vcf2bedtool(vcfpath)

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
