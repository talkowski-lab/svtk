=======================================================================
svtk.utils - Utility function and classes for working with SV in Python
=======================================================================

Introduction
============

The svtk.utils module includes a variety of utility functions and classes that
are used throughout the other components of SVTK and are intended for general
purpose use in SV analysis.

SVTK relies heavily on pysam and the utils module is frequently imported in
tandem::

   import pysam
   import svtk.utils as svu

API
===

Converting a VCF to a BedTool
-----------------------------

.. autofunction:: svtk.utils.vcf2bedtool

Manipulating and querying the samples called in a VariantRecord
---------------------------------------------------------------

.. autofunction:: svtk.utils.get_called_samples

.. autofunction:: svtk.utils.samples_overlap

.. autofunction:: svtk.utils.set_null

Chromosome and genomic interval comparisons and logic
-----------------------------------------------------

.. autofunction:: svtk.utils.is_smaller_chrom

.. autofunction:: svtk.utils.recip

.. autofunction:: svtk.utils.reciprocal_overlap

.. autofunction:: svtk.utils.overlap_frac

Filtering alignments
--------------------

.. autofunction:: svtk.utils.is_excluded

.. autofunction:: svtk.utils.is_soft_clipped

Writing to a bgzipped file
--------------------------

.. autoclass:: svtk.utils.BgzipFile

Miscellaneous helpers
---------------------

.. autofunction:: svtk.utils.make_bnd_alt
