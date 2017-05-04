#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Cluster a bed produced by a bedtools intersection of a bed with itself.

A self-intersected bed has two sets of columns in each row, corresponding to
two entries in the original bed that shared sufficient overlap. The clustering
is performed by linking each such pair and then returning all clusters of
linked calls.
"""

import argparse
import sys
from collections import namedtuple, deque
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
from collections import defaultdict


class FieldNumberError(Exception):
    """Number of fields detected don't match number specified"""


BedCall = namedtuple('BedCall', 'chrom start end ID sample svtype'.split())


def rmsstd(calls):
    starts = np.array([call.start for call in calls])
    ends = np.array([call.end for call in calls])

    def _meanSS(X):
        mu = np.mean(X)
        return np.sum((X - mu) ** 2) / len(X)

    SS = _meanSS(starts) + _meanSS(ends)
    return np.sqrt(SS)


def bedcluster(bed, variant_list, preserve_links=False):
    """
    Cluster a self-intersected bed.

    Parameters
    ----------
    bed : file
        Product of a bedtools intersect.
        Columns: chrA, startA, endA, nameA, sampleA, svtypeA, chrB, startB,
        endB, nameB, sampleB, svtypeB
    variant_list : file
        Unique variant IDs in bed. Necessary for sparse graph

    Returns
    -------
    clusters : list of BedCall
        (chrom, start, end, ID, sample, svtype)
    samples : list of str
        List of unique samples found (used in population frequency check)
    """

    def _is_null(bedcall):
        return bedcall.chrom == '.'

    variant_indexes = {}
    variant_IDs = [l.strip() for l in variant_list.readlines()]
    for i, variant in enumerate(variant_IDs):
        variant_indexes[variant.strip()] = i
    num_variants = len(variant_indexes)

    #  G = nx.Graph()
    G = sparse.eye(num_variants, dtype=np.uint16, format='lil')

    samples = set()
    prev_c1 = None

    variants = {}

    for line in bed:
        # Skip header
        if line.startswith('#'):
            continue

        data = line.strip().split()

        try:
            c1 = data[:6]
            c2 = data[6:]

            # Cast start/end to ints
            c1[1], c1[2] = int(c1[1]), int(c1[2])
            c2[1], c2[2] = int(c2[1]), int(c2[2])

            c1 = BedCall(*c1)
            c2 = BedCall(*c2)

            if c1.ID not in variants:
                idx = variant_indexes[c1.ID]
                variants[idx] = c1
            if c2.ID not in variants:
                idx = variant_indexes[c2.ID]
                variants[idx] = c2

        except:
            c = len(data) / 2
            b = '%s: ' % bed.name
            msg = b + '%d fields per region specified, %d found' % (6, c)
            raise FieldNumberError(msg)

        samples.add(c1.sample)

        # Link the two calls from the current line
        if not _is_null(c2) and c1.svtype == c2.svtype:
            idx1 = variant_indexes[c1.ID]
            idx2 = variant_indexes[c2.ID]
            G[idx1, idx2] = 1
            samples.add(c2.sample)

        if preserve_links and prev_c1 is not None and c1.ID == prev_c1.ID:
            idx1 = variant_indexes[c1.ID]
            idx2 = variant_indexes[prev_c1.ID]
            G[idx1, idx2] = 1

        prev_c1 = c1

    #  clusters = list(nx.connected_components(G))
    n_comp, labels = csgraph.connected_components(G, connection='weak')

    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[label].append(variants[idx])

    clusters = list(clusters.values())

    clusters = [sorted(cluster, key=lambda c: (c.chrom, int(c.start)))
                for cluster in clusters]
    clusters = sorted(clusters, key=lambda c: (c[0].chrom, int(c[0].start)))

    return clusters, samples


def collapse_sample_calls(cluster):
    """
    Merges multiple variants in same sample

    Parameters
    ----------
    cluster : list of BedCall

    Returns
    -------
    cluster : list of BedCall
    """

    calldict = defaultdict(list)
    variants = deque()

    # Get all calls in each sample
    for call in cluster:
        calldict[call.sample].append(call)

    # If a sample has only one call, keep it, otherwise merge
    for sample, calls in calldict.items():
        if len(calls) == 1:
            variants.append(calls[0])
            continue

        # Track IDs of merged variants
        ID = ','.join([call.ID for call in calls])

        # To merge variants in a sample, take broadest range
        start = np.min([call.start for call in calls])
        end = np.max([call.end for call in calls])

        chrom = calls[0].chrom
        sample = calls[0].sample
        svtype = calls[0].svtype

        merged = BedCall(chrom, start, end, ID, sample, svtype)
        variants.append(merged)

    return list(variants)


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtools bedcluster',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='SV calls to cluster. Columns: #chr, '
                        'start, end, name, sample, svtype')
    parser.add_argument('-r', '--region', help='Region to cluster '
                        '(chrom:start-end). Requires tabixed bed.')
    parser.add_argument('-p', '--prefix', default='prefix',
                        help='Cluster ID prefix')
    parser.add_argument('-f', '--max-freq',
                        type=float, default=1.0,
                        help='Max population frequency to allow [1.0]')
    parser.add_argument('-l', '--preserve-links',
                        action='store_true', default=False,
                        help='Link two consecutive bed entries that share ID. '
                        'Used if clustering the intersection of two clustered '
                        'beds. Requires all entries from same cluster be '
                        'adjacent.')
    parser.add_argument('-m', '--merge-coordinates',
                        action='store_true', default=False,
                        help='Report median of start and end positions in '
                        'each cluster as final coordinates of cluster. '
                        'Recommended to turn off with preserve-links.')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout,
                        help='Clustered bed.')
    args = parser.parse_args(argv)

    bed = pbt.BedTool(args.bed)
    if args.region:
        bed = bed.tabix_intervals(args.region)

    # Drop any columns beyond those required
    bed = bed.cut(range(6))

    


    header = ('#chrom start end name svtype sample call_name vaf vac '
              'pre_rmsstd post_rmsstd')
    header = '\t'.join(header.split()) + '\n'
    args.fout.write(header)

    clusters, samples = bedcluster(args.bed, args.variant_list,
                                   args.preserve_links)
    num_samples = float(len(samples))

    for i, cluster in enumerate(clusters):
        # Calculate RMSSTD before merging per-sample variants
        pre_RMSSTD = rmsstd(cluster)

        # Make a single variant for each sample
        cluster = collapse_sample_calls(cluster)

        # Re-calculate RMSSTD after merging per-sample variants
        post_RMSSTD = rmsstd(cluster)

        # Merge coordinates AFTER getting min/max per sample
        if args.merge_coordinates:
            # Report median region of overlap
            start = str(int(np.median([int(call.start) for call in cluster])))
            end = str(int(np.median([int(call.end) for call in cluster])))
            cluster = [BedCall(c.chrom, start, end, *c[3:]) for c in cluster]

        # Get variant frequency info
        vac = len(set([call.sample for call in cluster]))
        vaf = vac / num_samples

        # Filter clusters exceeding max population freq
        if vaf > args.max_freq:
            continue

        # Assign cluster ID
        cid = args.prefix + ('_%d' % i)

        for call in cluster:
            entry = ('{chrom}\t{start}\t{end}\t{{cid}}\t{svtype}\t'
                     '{sample}\t{ID}').format(**call._asdict())
            entry = (entry + '\t{vaf:.3f}\t{vac}\t{pre_RMSSTD:.3f}\t'
                     '{post_RMSSTD:.3f}\n').format(**locals())

            args.fout.write(entry)


if __name__ == '__main__':
    main(sys.argv[1:])
