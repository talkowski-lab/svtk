#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import os
import tempfile
import pkg_resources
import subprocess as sp
import svtk.utils as svu


def _record_to_bed(record):
    entry = '{chrom}\t{start}\t{end}\t{name}\t{samples}\t{svtype}\n'

    if record.info['SVTYPE'] in 'DEL DUP'.split():
        return entry.format(chrom=record.chrom,
                            start=record.pos,
                            end=record.stop,
                            name=record.id,
                            samples=','.join(svu.get_called_samples(record)),
                            svtype=record.info['SVTYPE'])
    else:
        return ''


def _make_rdtest_bed(variants):
    """
    Make temporary bed file for RdTest

    Parameters
    ----------
    variants : list of pysam.VariantRecord

    Returns
    -------
    bed : tempfile.NamedTemporaryFile
    """

    bed = tempfile.NamedTemporaryFile(dir=os.getcwd())

    for variant in variants:
        bed.write(_record_to_bed(variants))

    return bed


def call_rdtest(variants, bincov_file, medianfile, famfile, whitelist):
    """
    Utility wrapper around a basic RdTest call

    Parameters
    ----------
    variants : list of pysam.VariantRecord
    bincov_file : str
    medianfile : str
    famfile : str
    whitelist : str or list of str
        Filepath to sample whitelist or list of whitelisted sample IDs

    Returns
    -------
    metrics : pd.DataFrame
        RdTest metrics for the provided variants
    """
    
    if not os.path.exists(bincov_file):
        raise Exception('Bincov file does not exist: {0}'.format(bincov_file))

    if not os.path.exists(medianfile):
        raise Exception('Medianfile does not exist: {0}'.format(medianfile))
    
    if not os.path.exists(famfile):
        raise Exception('Famfile does not exist: {0}'.format(famfile))

    if isinstance(whitelist, str):
        whitelist_filename = whitelist
    elif isinstance(whitelist, list):
        whitelist_file = tempfile.NamedTemporaryFile(dir=os.getcwd())
        for sample in whitelist:
            whitelist_file.write(sample + '\n')
        whitelist_filename = whitelist_file.name
    else:
        msg = 'Invalid type for whitelist: {0}\n'.format(str(type(whitelist)))
        msg += 'Must be str or list of str'
        raise Exception(msg)

    output_dir = tempfile.TemporaryDirectory(dir=os.getcwd())

    bed = _make_rdtest_bed(variants)

    RdTest = pkg_resources.resource_filename('svtk', 'RdTest/RdTest.R')

    sp.run(['Rscript', RdTest,
            '-b', bed.name,
            '-o', output_dir.name,
            '-n', 'tmp',
            '-c', bincov_file,
            '-m', medianfile,
            '-f', famfile,
            '-w', whitelist_filename,
            flags])

    metrics = pd.read_table(os.path.join(output_dir.name, 'tmp.metrics'))

    return metrics
