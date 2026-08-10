Getting started
===============

pyft: a fiber-seq inspection toolkit for Python
-----------------------------------------------

.. image:: https://readthedocs.org/projects/py-ft/badge/?version=latest
    :target: https://py-ft.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://badge.fury.io/py/pyft.svg
    :target: https://badge.fury.io/py/pyft


**pyft** is a lightweight toolkit for inspecting and visualizing `fiber-seq
<https://github.com/fiberseq/fibertools-rs>`_ BAMs from Python. It reads a
record with `pysam <https://pysam.readthedocs.io>`_ and decodes its fiber-seq
annotations through the `molecular-annotation
<https://pypi.org/project/molecular-annotation/>`_ library, then exposes them
under the fiber-seq vocabulary: ``m6a``, ``cpg`` (5mC), ``nuc``, ``msp`` and
``fire``. It is read/inspection-focused — for producing or editing fiber-seq
annotations, use `fibertools-rs <https://github.com/fiberseq/fibertools-rs>`_.


Install
=======
.. code-block:: bash

    pip install pyft

The plotting helpers (``pyft.plot`` / ``pyft.utils`` and the ``fiberplot`` CLI)
require extra dependencies:

.. code-block:: bash

    pip install "pyft[viz]"


Quick example
=============
.. code-block:: python

    from pyft import fetch

    for fiber in fetch("fiberseq.bam"):
        print(fiber.name, fiber.strand, "m6A:", len(fiber.m6a), "MSPs:", len(fiber.msp))


Vignettes
=========
The :doc:`vignettes <vignettes/index>` are a good place to start.


Indices and tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: pyft

    self
    api.rst
    vignettes/index.rst
