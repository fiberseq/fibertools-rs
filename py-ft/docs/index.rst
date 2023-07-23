.. pyft documentation master file, created by
   sphinx-quickstart on Wed Jul 19 20:00:42 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. toctree::
   :maxdepth: 2
   :caption: Contents:

pyft: python bindings for fibertools-rs
=======================================
.. image:: https://readthedocs.org/projects/py-ft/badge/?version=latest
    :target: https://py-ft.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://badge.fury.io/py/pyft.svg
    :target: https://badge.fury.io/py/pyft


**pyft** provides a python API for the rust library `fibertools-rs <https://github.com/fiberseq/fibertools-rs>`_. The inspiration for this API is to make analysis in python easier and faster; therefore, only extraction of data from a fiberseq bam is supported and not writing. 

Install
=======
.. code-block:: bash

    pip install pyft


Example
=======
.. highlight:: python
.. include:: ../example.py
   :literal:
.. highlight:: none

API Reference
==================

.. automodule:: pyft
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
