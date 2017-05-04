#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Intersect SV called by PE/SR-based algorithms.

Paired-end and split-read callers provide a reasonably precise estimation of
an SV breakpoint. This program identifies variant calls that fall within
the expected margin of error made by these programs and clusters them together.
The cluster distance defaults to 500 bp but it is recommended to use the
maximum individual clustering distance across the libraries being analyzed.
(Generally median + 7 * MAD)
"""

import argparse
import os
import sys
from collections import deque
from pysam import VariantFile, TabixFile

from svtools.vcfcluster import VCFCluster


def flatten_pos(records, name, fout):
    for record in records:
        chrom = record.CHROM
        start = record.POS
        end = record.INFO['END']
        quad = record.samples[0].sample.split('.')[0]
        source = record.source
        ID = record.ID
        svtype = record.INFO['SVTYPE']

        entry = [str(x) for x in [chrom, start, end, quad, source, ID, name,
                                  svtype]]
        entry = '\t'.join(entry) + '\n'
        fout.write(entry)


def parse_filepaths(filepaths):
    """
    Parameters
    ----------
    filepaths : list of str
        List of paths to standardized VCFs

    Returns
    -------
    vcfs : list of pysam.VariantFile
    """

    vcfs = deque()
    for path in filepaths:
        if len(path.split()) != 1:
            raise ValueError('File list must be single column')
        if not os.path.isfile(path):
            raise FileNotFoundError('VCF {0} not found'.format(path))

        vcf = VariantFile(path)
        vcfs.append(vcf)

    return vcfs


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtools vcfcluster',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filelist', type=argparse.FileType('r'),
                        help='List of paths to standardized VCFS')
    parser.add_argument('fout', help='Clustered VCF.')
    parser.add_argument('-r', '--region', default=None,
                        help='Restrict clustering to genomic region.')
    parser.add_argument('-d', '--dist',
                        type=int, default=500,
                        help='Maximum clustering distance. Suggested to use '
                        'max of median + 7*MAD over samples. [500]')
    parser.add_argument('-f', '--frac',
                        type=float, default=0.1,
                        help='Minimum reciprocal overlap between variants. '
                        '[0.1]')
    parser.add_argument('-m', '--match-strands', default=False,
                        action='store_true')
    parser.add_argument('-x', '--blacklist', metavar='BED.GZ',
                        type=TabixFile, default=None,
                        help='Tabix indexed bed of blacklisted regions. Any '
                        'SV with a breakpoint falling inside one of these '
                        'regions is filtered from output.')
    parser.add_argument('-z', '--svsize', type=int, default=500,
                        help='Minimum SV size to report for intrachromosomal '
                        'events. [0]')
    parser.add_argument('-p', '--prefix',
                        default='MERGED',
                        help='Prefix for merged variant IDs. [MERGED]')
    parser.add_argument('-t', '--svtypes', default='DEL,DUP,INV,BND',
                        help='Comma delimited list of svtypes to cluster '
                        '[DEL,DUP,INV,BND]')
    #  parser.add_argument('--cluster-bed', type=argparse.FileType('w'),
    #                      help='Bed of constituent calls in each cluster')
    args = parser.parse_args(argv)

    # Parse SV files and lists of samples and sources
    filepaths = [line.strip() for line in args.filelist.readlines()]
    vcfs = parse_filepaths(filepaths)

    #  if args.svtypes is not None:
    #      svtypes = args.svtypes.split(',')
    #      permitted_types = 'del dup inv tloc'.split()
    #      for svtype in svtypes:
    #          if svtype not in permitted_types:
    #              raise Exception('Unpermitted svtype: {0}'.format(svtype))
    #  else:
    svtypes = args.svtypes.split(',')

    svc = VCFCluster(vcfs, dist=args.dist, blacklist=args.blacklist,
                     frac=args.frac, match_strands=args.match_strands,
                     svtypes=svtypes, region=args.region)

    # Open new file
    fout = VariantFile(args.fout, mode='w', header=svc.header)

    for i, (record, cluster) in enumerate(svc.cluster()):
        # Name record
        if args.prefix:
            name = [args.prefix]
        else:
            name = ['SV']
        if args.region:
            chrom = args.region.split(':')[0]
            name.append(chrom)
        name.append(str(i + 1))
        record.id = '_'.join(name)

        fout.write(record)

        # Size filter (CTX have size -1)
        if -1 < record.info['SVLEN'] < args.svsize:
            continue

        #  if args.cluster_bed is not None:
            #  flatten_pos(cluster, record.ID, args.cluster_bed)

    fout.close()

if __name__ == '__main__':
    main(sys.argv[1:])
