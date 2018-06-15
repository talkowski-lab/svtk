==============================================================================
SVTK CLI - Command line tools for processing and analyzing structural variants
==============================================================================

Introduction
============

vcfcluster
==========

.. argparse::
   :module: svtk.cli.vcfcluster
   :func: make_argparse
   :prog: svtk vcfcluster

annotate
========

.. argparse::
   :module: svtk.cli.annotate
   :func: make_argparse
   :prog: svtk annotate

bedcluster
==========

.. argparse::
   :module: svtk.cli.bedcluster
   :func: make_argparse
   :prog: svtk bedcluster

bincov
======

.. argparse::
   :module: svtk.cli.bincov
   :func: make_argparse
   :prog: svtk bincov

collect-pesr
============

.. argparse::
   :module: svtk.cli.collect_pesr
   :func: make_argparse
   :prog: svtk collect-pesr


count-svtypes
=============

.. argparse::
   :module: svtk.cli.count_svtypes
   :func: make_argparse
   :prog: svtk count-svtypes

resolve
=======

.. argparse::
   :module: svtk.cli.resolve
   :func: make_argparse
   :prog: svtk resolve

standardize
===========

.. argparse::
   :module: svtk.cli.standardize_vcf
   :func: make_argparse
   :prog: svtk standardize

vcf2bed
=======

.. argparse::
   :module: svtk.cli.utils
   :func: make_vcf2bed_argparse
   :prog: svtk vcf2bed

remote_tabix
============

.. argparse::
   :module: svtk.cli.utils
   :func: make_remote_tabix_argparse
   :prog: svtk remote-tabix

pe-test
=======

.. argparse::
   :module: svtk.cli.pesr_test
   :func: make_pe_test_argparse
   :prog: svtk pe-test

sr-test
=======

.. argparse::
   :module: svtk.cli.pesr_test
   :func: make_sr_test_argparse
   :prog: svtk sr-test

count-pe
========

.. argparse::
   :module: svtk.cli.pesr_test
   :func: make_count_pe_argparse
   :prog: svtk count-pe

count-sr
========

.. argparse::
   :module: svtk.cli.pesr_test
   :func: make_count_sr_argparse
   :prog: svtk count-sr
