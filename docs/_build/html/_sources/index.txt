.. bcflash documentation master file, created by
   sphinx-quickstart on Tue Sep 25 15:02:08 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Bound-Constrained Flash (BCFLASH)
=================================

BCFLASH is a Steighaug Conjugate-Gradient type solver for optimization problems of the form:

.. math::

   \min_{x \in \mathbb{R}^n}\enspace f(x) \enspace\mbox{subject to}\enspace \ell \le x \le u

where :math:`f` has two continuous derivatives; it is essentially a MATLAB implementation of TRON_.

.. _TRON: https://www.mcs.anl.gov/~more/tron/

Contents
========

.. toctree::
   :maxdepth: 2

   installation-usage
   input-output

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

