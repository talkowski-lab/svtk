#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import itertools
from collections import deque
import numpy as np
import scipy.sparse as sps
import pysam
import pybedtools as pbt
import natsort
import svtools.utils as svu


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


def extract_breakpoints(vcfpath, IDs):
    """
    Extract all VCF records in list of IDs
    (Assumes VCF is sorted by variant ID)

    Parameters
    ----------
    vcfpath : str
        Path to VCF
    IDs : list of str
        Variant IDs to extract

    Returns
    -------
    bkpts : list of pysam.VariantRecord
    """

    vcf = pysam.VariantFile(vcfpath)
    n_bkpts = len(IDs)
    bkpts = np.empty(n_bkpts, dtype=object)
    idx = 0

    for record in vcf:
        if record.id in IDs:
            bkpts[idx] = record
            idx += 1
            if idx == n_bkpts:
                break

    # TODO: fix upstream VCF output so input IDs are sorted
    #  for record in vcf:
    #      if record.id == IDs[idx]:
    #          bkpts[idx] = record
    #          idx += 1
    #          if idx == n_bkpts:
    #              break

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


def resolve_cx(cluster):
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
    cnvtypes = 'DEL DUP'.split()
    cnvs = [rec for rec in cluster if rec.info['SVTYPE'] in cnvtypes]

    # Restrict to double-ended inversion events with appropriate strand pairing
    if len(inversions) == 1:
        cluster_type = 'SINGLE_ENDER'
    elif len(inversions) > 2:
        cluster_type = 'COMPLEX_3plus'
    elif inversions[0].info['STRANDS'] == inversions[1].info['STRANDS']:
        cluster_type = 'MATCHED_STRANDS'
    else:
        cluster_type = 'CANDIDATE'

    if cluster_type != 'CANDIDATE':
        return cluster_type, cluster


def cx_link(vcfpath, bkpt_window=100):
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
    link_key = {ID: i for i, ID in enumerate(linked_IDs)}
    keyed_links = np.array([(link_key[a], link_key[b]) for a, b in links])

    # Extract VariantRecords corresponding to breakpoints
    n_bkpts = len(linked_IDs)
    bkpts = extract_breakpoints(vcfpath, linked_IDs)

    # Build sparse graph from links
    G = sps.eye(n_bkpts, dtype=np.uint16, format='lil')
    for i, j in keyed_links:
        if samples_overlap(bkpts[i], bkpts[j]):
            G[i, j] = 1

    # Generate lists of clustered breakpoints
    n_comp, comp_list = sps.csgraph.connected_components(G)
    clusters = [deque() for x in range(n_comp)]
    for i, c_label in enumerate(comp_list):
        clusters[c_label].append(bkpts[i])

    # Remove clusters of one variant - leftover from shared sample filtering
    clusters = [c for c in clusters if len(c) > 1]

    for cluster in clusters:
        variant = resolve_cx(cluster)
    return clusters


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Breakpoint VCFs.')
    args = parser.parse_args()

    cx_link(args.vcf)


if __name__ == '__main__':
    main()
