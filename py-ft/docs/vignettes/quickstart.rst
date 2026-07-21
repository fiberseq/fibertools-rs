Quickstart
==========

``pyft`` turns each record of a fiber-seq BAM into a :class:`~pyft.Fiberdata`
object that exposes the fiber-seq features under familiar names.

Iterating a BAM
---------------

:func:`~pyft.fetch` yields one :class:`~pyft.Fiberdata` per read. Without a
region it walks the whole file, including **unmapped** reads, so it works on raw
(unaligned) fiber-seq BAMs as well as aligned ones:

.. code-block:: python

    from pyft import fetch

    for fiber in fetch("fiberseq.bam"):
        print(
            fiber.name,
            f"{fiber.chrom}:{fiber.start}-{fiber.end}",
            fiber.strand,
            "m6A:", len(fiber.m6a),
            "5mC:", len(fiber.cpg),
            "nuc:", len(fiber.nuc),
            "MSP:", len(fiber.msp),
        )

To restrict to a region of an indexed BAM, pass a ``(chrom, start, end)`` tuple
(this returns the mapped reads overlapping the region):

.. code-block:: python

    for fiber in fetch("fiberseq.bam", region=("chr1", 100_000, 101_000)):
        ...

Annotation accessors
--------------------

Each accessor returns a list of :class:`~pyft.Feature`. Coordinates are 0-based
half-open; ``start``/``end`` are in molecular (read) orientation, while
``ref_start``/``ref_end`` are lifted to the reference (``None`` where the
annotation falls outside an aligned block). ``qual`` is the per-feature quality
(the ML probability for base mods, MSP precision for ``msp``/``fire``):

.. code-block:: python

    fiber = next(fetch("fiberseq.bam"))

    # m6A calls that lifted to the reference
    m6a_ref = [f.ref_start for f in fiber.m6a if f.ref_start is not None]

    # nucleosomes as reference intervals
    nucs = [(f.ref_start, f.ref_end) for f in fiber.nuc]

    # high-confidence m6A only
    strong = [f for f in fiber.m6a if f.qual >= 200]

The available accessors are ``m6a``, ``cpg`` (5mC), ``nuc``, ``msp`` and
``fire`` (a distinct type if present, otherwise the high-precision MSP subset).

Mapped vs. unmapped reads
-------------------------

Unmapped reads still carry their molecular annotations; only the reference
coordinates are unavailable. Filter with :attr:`~pyft.Fiberdata.is_mapped`:

.. code-block:: python

    mapped = [f for f in fetch("fiberseq.bam") if f.is_mapped]

Plotting
--------

The accessors feed directly into any plotting library. This draws a browser-style
row per fiber — nucleosomes as blocks, m6A as ticks — using reference
coordinates (requires ``matplotlib``):

.. code-block:: python

    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from pyft import fetch

    fig, ax = plt.subplots(figsize=(12, 6))
    for y, fiber in enumerate(fetch("fiberseq.bam", region=("chr1", 100_000, 101_000))):
        for f in fiber.nuc:
            if f.ref_start is not None:
                ax.add_patch(Rectangle((f.ref_start, y - 0.2),
                                       f.ref_end - f.ref_start, 0.4, color="0.7"))
        xs = [f.ref_start for f in fiber.m6a if f.ref_start is not None]
        ax.plot(xs, [y] * len(xs), "|", color="#800080", markersize=5)
    ax.set(xlabel="chr1 position", ylabel="fiber")
    fig.savefig("fibers.png", dpi=150)
