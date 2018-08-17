#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Resolve complex SV from inversion/translocation breakpoints and CNV intervals.
"""

import argparse
import sys
import subprocess
import random
import string
from collections import deque
import itertools
import pysam
import pandas as pd
import pybedtools as pbt
import svtk.utils as svu
from svtk.cxsv import link_cpx, ComplexSV, rescan_single_ender,link_cpx_V2


CPX_INFO = [
    '##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">',
    '##ALT=<ID=CPX,Description="Complex SV">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=UNR,Description="Unresolved breakend or complex SV">',
    '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of inserted sequence.">',
    '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Class of complex variant.">',
    '##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Genomic intervals constituting complex variant.">',
    '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
    '##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="Variant is unresolved.">',
    '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">'
]


def _merge_records(vcf, cpx_records, cpx_record_ids):
    """
    r1, r2 : iter of pysam.VariantRecord
    """
    def _next_record():
        try:
            return next(vcf)
        except StopIteration:
            return None
    def _next_cpx():
        try:
            return cpx_records.popleft()
        except IndexError:
            return None
    # Initialize merge
    curr_record = _next_record()
    curr_cpx = _next_cpx()
    while curr_record is not None and curr_cpx is not None:
        # Remove VCF records that were included in complex event
        if curr_record.id in cpx_record_ids:
            curr_record = _next_record()
            continue
        # Merge sort remaining
        if curr_record.chrom == curr_cpx.chrom:
            if curr_record.pos <= curr_cpx.pos:
                yield curr_record
                curr_record = _next_record()
            else:
                yield curr_cpx
                curr_cpx = _next_cpx()
        elif svu.is_smaller_chrom(curr_record.chrom, curr_cpx.chrom):
            yield curr_record
            curr_record = _next_record()
        else:
            yield curr_cpx
            curr_cpx = _next_cpx()
    # After one iterator is exhausted, return rest of other iterator
    if curr_record is None:
        for cpx in itertools.chain([curr_cpx], cpx_records):
            yield cpx
    elif curr_cpx is None:
        for record in itertools.chain([curr_record], vcf):
            if record.id not in cpx_record_ids:
                yield record

def remove_CPX_from_INV(resolve_CPX, resolve_INV):
    cpx_interval = [[i.chrom, i.pos, i.stop, i] for i in resolve_CPX]
    inv_interval = [[i.chrom, i.pos, i.stop, i] for i in resolve_INV]
    out=[]
    for i in inv_interval:
        rec = 0
        for j in cpx_interval:
            if i[0]==j[0]:
                if i[2]< j[1] or i[1]> j[2]:
                    continue
                else:
                    rec=rec+1
        if rec ==0:
            out.append(i[3])
    return out

def cluster_INV(independent_INV):
    inv_hash={}
    for i in independent_INV:
        if not i.chrom in inv_hash.keys():
            inv_hash[i.chrom]={}
        if not i.pos in inv_hash[i.chrom].keys():
            inv_hash[i.chrom][i.pos] ={}
        if not i.stop in inv_hash[i.chrom][i.pos].keys():
            inv_hash[i.chrom][i.pos][i.stop]=i
    list_INV = {}
    for i in inv_hash.keys():
        list_INV[i]=[]
        for j in sorted(inv_hash[i].keys()):
            for k in sorted(inv_hash[i][j].keys()):
                list_INV[i].append(inv_hash[i][j][k])   
    out = []
    for i in list_INV.keys():
        out+=cluster_INV_list(list_INV[i])
    return out

def cluster_INV_list(independent_INV):
    out = []
    rec = 0
    for i in independent_INV:
        if out ==[]:
            out.append([i])
            rec = i.stop
        else:
            if i.pos > rec:
                out.append([i])
                rec = i.stop
            else:
                out[-1].append(i)
                if i.stop > rec:
                    rec = i.stop
    return out


def cluster_single_cleanup(cluster):
    #for clusters including both FF and RR inversions, only keep inversions for CPX resolution
    cluster_index = [i for i in range(len(cluster))]
    cluster_svtype = [[i.info['SVTYPE'], i.info['STRANDS']] for i in cluster]
    out=deque()
    if ['INV', '--'] in cluster_svtype and ['INV', '++'] in cluster_svtype:
        out.append(cluster[cluster_svtype.index(['INV', '--'])])
        out.append(cluster[cluster_svtype.index(['INV', '++'])])
    else:
        out=cluster
    return out

def clusters_cleanup(clusters):
    out=deque()
    for cluster in clusters:
        out.append(cluster_single_cleanup(cluster))
    return out


def resolve_complex_sv(vcf, cytobands, disc_pairs, mei_bed,variant_prefix='CPX_', min_rescan_support=4, pe_blacklist=None):
    """
    Resolve complex SV from CNV intervals and BCA breakpoints.
    Yields all resolved events, simple or complex, in sorted order.
    Parameters
    ----------
    vcf : pysam.VariantFile
    cytobands : pysam.TabixFile
    disc_pairs : pysam.TabixFile
    mei_bed : pybedtools.BedTool
    variant_prefix : str
        Prefix to assign to resolved variants
    min_rescan_support : int
        Number of pairs required to count a sample as 
        supported during PE rescan
    pe_blacklist : pysam.TabixFile, optional
        Blacklisted genomic regions. Anomalous pairs in these regions will be
        removed prior to clustering.
    Yields
    ------
    sv : pysam.VariantRecord
    """

    clusters = link_cpx(vcf)
    clusters = clusters_cleanup(clusters)

    # resolved_idx = unresolved_idx = 1

    if not variant_prefix.endswith('_'):
        variant_prefix += '_'

    cpx_records = deque()
    cpx_record_ids = set()

    for cluster in clusters:
        # Try finding opposite strand support for single ender inversions
        if len(cluster) == 1 and cluster[0].info['SVTYPE'] == 'INV':
            rec, opp = rescan_single_ender(cluster[0], disc_pairs, 
                                           min_rescan_support, 
                                           pe_blacklist=pe_blacklist)
            if opp is not None:
                cluster = deque([rec, opp])
        # if cxsv overlap pulled in unrelated insertions, keep them separate
        if all([r.info['SVTYPE'] == 'INS' for r in cluster]):
            for record in cluster:
                cpx = ComplexSV([record], cytobands, mei_bed)
                cpx_record_ids = cpx_record_ids.union(cpx.record_ids)
                
                # Assign random string as resolved ID to handle sharding
                cpx.vcf_record.id = variant_prefix + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
                cpx_records.append(cpx.vcf_record)
                # resolved_idx += 1
        else:
            cpx = ComplexSV(cluster, cytobands, mei_bed)
            cpx_record_ids = cpx_record_ids.union(cpx.record_ids)
            if cpx.svtype == 'UNR':
                # Assign random string as unresolved ID to handle sharding
                unresolved_vid = 'UNRESOLVED_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
                for i, record in enumerate(cpx.records):
                    record.info['EVENT'] = unresolved_vid
                    record.info['CPX_TYPE'] = cpx.cpx_type
                    record.info['UNRESOLVED'] = True
                    cpx_records.append(record)
                # unresolved_idx += 1
            else:
                cpx.vcf_record.id = variant_prefix + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
                cpx_records.append(cpx.vcf_record)
                # resolved_idx += 1

    # Output all variants
    vcf.reset()

    for record in _merge_records(vcf, cpx_records, cpx_record_ids):
        if 'CPX_TYPE' in record.info.keys():
            if 'UNRESOLVED' in record.info.keys():
                record.info['UNRESOLVED_TYPE'] = record.info['CPX_TYPE']
                record.info.pop('CPX_TYPE')
            else:
                record.info.pop('STRANDS')
        if 'CIPOS' in record.info.keys():
            record.info.pop('CIPOS')
        if 'CIEND' in record.info.keys():
            record.info.pop('CIEND')
        if 'RMSSTD' in record.info.keys():
            record.info.pop('RMSSTD')
        yield record

def cluster_cleanup(clusters_v2):
    rank  =  [i for i in range(len(clusters_v2))]
    cluster_pos = []
    cluster_info = []    
    for i in range(len(clusters_v2)):
        info = clusters_v2[i]
        info_pos = '_'.join(sorted([','.join([str(k) for k in [j.pos, j.stop, j.info['SVTYPE']]]) for j in info]))
        if not info_pos in cluster_info:
            cluster_info.append(info_pos)
            cluster_pos.append(i)
    return [clusters_v2[i] for i in cluster_pos]

def resolve_complex_sv_v2(resolve_CPX, resolve_INV, resolve_CNV, cytobands,disc_pairs, mei_bed,variant_prefix='CPX_', min_rescan_support=4, pe_blacklist=None):
    #resolve_CPX = [i for i in out_rec if i.info['SVTYPE']=='CPX']
    #resolve_INV = [i for i in out_rec if i.info['SVTYPE']=='INV']
    independent_INV = remove_CPX_from_INV(resolve_CPX, resolve_INV)
    linked_INV = cluster_INV(independent_INV)
    clusters_v2=link_cpx_V2(linked_INV,resolve_CNV,cpx_dist=2000 )
    clusters_v2=cluster_cleanup(clusters_v2)
    cpx_records_v2 = deque()
    cpx_record_ids_v2 = set()
    for cluster in clusters_v2:
        # Try finding opposite strand support for single ender inversions
        if len(cluster) == 1 and cluster[0].info['SVTYPE'] == 'INV':
            rec, opp = rescan_single_ender(cluster[0], disc_pairs, min_rescan_support, pe_blacklist=pe_blacklist)
            if opp is not None:
                cluster = deque([rec, opp])
        # if cxsv overlap pulled in unrelated insertions, keep them separate
        if all([r.info['SVTYPE'] == 'INS' for r in cluster]):
            for record in cluster:
                cpx = ComplexSV([record], cytobands, mei_bed)
                cpx_record_ids = cpx_record_ids.union(cpx.record_ids)
                
                # Assign random string as resolved ID to handle sharding
                cpx.vcf_record.id = variant_prefix + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
                cpx_records.append(cpx.vcf_record)
                # resolved_idx += 1
        else:
            cpx = ComplexSV(cluster, cytobands, mei_bed)
            cpx_record_ids_v2 = cpx_record_ids_v2.union(cpx.record_ids)
            if cpx.svtype == 'UNR':
                # Assign random string as unresolved ID to handle sharding
                unresolved_vid = 'UNRESOLVED_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
                for i, record in enumerate(cpx.records):
                    record.info['EVENT'] = unresolved_vid
                    record.info['CPX_TYPE'] = cpx.cpx_type
                    record.info['UNRESOLVED'] = True
                    cpx_records_v2.append(record)
                # unresolved_idx += 1
            else:
                cpx.vcf_record.id = variant_prefix + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
                cpx_records_v2.append(cpx.vcf_record)
    for i in cpx_records_v2:
        if i.info['SVTYPE'] == 'CPX':
            for info in 'UNRESOLVED EVENT UNRESOLVED_TYPE'.split():
                if info in i.info.keys():
                    i.info.pop(info)
    return cpx_records_v2

def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk resolve',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('raw', help='Filtered breakpoints and CNV intervals.')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--discfile', default=None,
                       help='Scraped discordant pairs. Required '
                       'to attempt to resolve single-ender inversions.')
    group.add_argument('--discfile-list', default=None,
                       type=argparse.FileType('r'),
                       help='Tab-delimited list of discordant pair files '
                       'and indices')
    parser.add_argument('resolved', type=argparse.FileType('w'),
                        help='Resolved simple and complex variants.')
    parser.add_argument('--mei-bed', help='Mobile element insertion bed. '
                        'Required to classify inverted insertions.',
                        required=True)
    parser.add_argument('--cytobands', help='Cytoband file. Required to '
                        'correctly classify interchromosomal events.',
                        required=True)
    #  parser.add_argument('--bincov', help='Bincov file.', required=True)
    #  parser.add_argument('--medianfile', help='Medianfile', required=True)
    #  parser.add_argument('--famfile', help='Fam file', required=True)
    #  parser.add_argument('--cutoffs', help='Random forest cutoffs',
                        #  required=True)
    parser.add_argument('--min-rescan-pe-support', type=int, default=4, 
                        help='Minumum discordant pairs required during '
                        'single-ender rescan.')
    parser.add_argument('-x', '--pe-blacklist', metavar='BED.GZ',
                        default=None, help='Tabix indexed bed of blacklisted '
                        'regions. Any anomalous pair falling inside one '
                        'of these regions is excluded from PE rescanning.')
    parser.add_argument('-u', '--unresolved', type=argparse.FileType('w'),
                        help='Unresolved complex breakpoints and CNV.')
    parser.add_argument('-p', '--prefix', default='CPX_',
                        help='Variant prefix [CPX_]')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = pysam.VariantFile(args.raw)
    for line in CPX_INFO:
        vcf.header.add_line(line)

    resolved_pipe = subprocess.Popen(['vcf-sort', '-c'],
                                     stdin=subprocess.PIPE,
                                     stdout=args.resolved)

    resolved_f = pysam.VariantFile(resolved_pipe.stdin, 'w', header=vcf.header)
    unresolved_f = pysam.VariantFile(args.unresolved, 'w', header=vcf.header)

    cytobands = pysam.TabixFile(args.cytobands)

    mei_bed = pbt.BedTool(args.mei_bed)
    if args.pe_blacklist is not None:
        blacklist = pysam.TabixFile(args.pe_blacklist)
    else:
        blacklist = None
    #  cutoffs = pd.read_table(args.cutoffs)
    #  rdtest = svu.RdTest(args.bincov, args.medianfile, args.famfile, 
                        #  list(vcf.header.samples), cutoffs)

    if args.discfile is not None:
        disc_pairs = pysam.TabixFile(args.discfile)
    else:
        tabixfiles = []
        for line in args.discfile_list:
            fname, idx = line.strip().split()
            tabixfiles.append(pysam.TabixFile(fname, index=idx))
        disc_pairs = svu.MultiTabixFile(tabixfiles)

    resolve_CPX=[]
    resolve_INV=[]
    resolve_CNV=[]
    cpx_dist=20000
    for record in resolve_complex_sv(vcf, cytobands, disc_pairs, mei_bed, args.prefix, args.min_rescan_pe_support, blacklist):
        #Move members to existing variant IDs unless variant is complex
        if record.info['SVTYPE'] != 'CPX' and 'CPX' not in record.id.split('_'):
            record.info['MEMBERS'] = record.id
        #Treat 
        if record.info['UNRESOLVED']:
            unresolved_f.write(record)
        else:
            resolved_f.write(record)
        if record.info['SVTYPE']=='CPX': 
            resolve_CPX.append(record)
        if record.info['SVTYPE']=='INV' and record.stop-record.pos > cpx_dist: 
            resolve_INV.append(record)
        if record.info['SVTYPE'] in ['DEL','DUP'] and record.stop-record.pos > cpx_dist/2:
            resolve_CNV.append(record) 

    #out_rec = resolve_complex_sv(vcf, cytobands, disc_pairs, mei_bed, args.prefix, args.min_rescan_pe_support, blacklist)
    cpx_records_v2 = resolve_complex_sv_v2(resolve_CPX, resolve_INV, resolve_CNV ,cytobands, disc_pairs, mei_bed, args.prefix, args.min_rescan_pe_support, blacklist)

    for record in cpx_records_v2:
        #Move members to existing variant IDs unless variant is complex
        if record.info['SVTYPE'] != 'CPX' and 'CPX' not in record.id.split('_'):
            record.info['MEMBERS'] = record.id
        if record.info['UNRESOLVED']:
            unresolved_f.write(record)
        else:
            resolved_f.write(record)
    resolved_f.close()
    unresolved_f.close()

    stdout, stderr = resolved_pipe.communicate()


if __name__ == '__main__':
    main(sys.argv[1:])
