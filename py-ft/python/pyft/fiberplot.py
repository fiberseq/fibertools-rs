#!/usr/bin/env python3
"""
fiberplot - Generate an Altair HTML visualization of a specific genomic region from an indexed BAM file

Usage:
    fiberplot <bam_file> <region> [options]

Example:
    fiberplot my_data.bam chr22:26354169-26354500 -o output.html
    fiberplot my_data.bam chr22:26354169-26354500 --width 1000 --height 600
"""

import argparse
import sys
import re
from pathlib import Path

import pyft
from pyft import plot, utils


def parse_region(region_str):
    """
    Parse a region string in the format 'chr:start-end' or 'chr:start'

    Args:
        region_str: String in format like "chr22:26354169-26354500" or "chr22:26354169"

    Returns:
        Tuple of (chrom, start, end) where end may be start+1 if not specified
    """
    # first remove any commas and whitespace
    region_str = region_str.replace(",", "").strip()
    # Match patterns like chr22:26354169-26354500 or chr22:26354169
    match = re.match(r"([^:]+):(\d+)(?:-(\d+))?", region_str)

    if not match:
        raise ValueError(
            f"Invalid region format: '{region_str}'. "
            "Expected format: 'chr:start-end' or 'chr:start' (e.g., 'chr22:26354169-26354500')"
        )

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3)) if match.group(3) else start + 1

    if start >= end:
        raise ValueError(
            f"Start position ({start}) must be less than end position ({end})"
        )

    return (chrom, start, end)


def create_fiberplot(
    bam_file,
    region,
    output_file=None,
    width=600,
    height=800,
    centered=False,
    strand="+",
    max_flank=None,
    min_basemod_qual=125,
):
    """
    Create an Altair visualization of fibers in a genomic region

    Args:
        bam_file: Path to indexed BAM file
        region: Tuple of (chrom, start, end) or string "chr:start-end"
        output_file: Path to output HTML file (default: fiberplot.html)
        width: Chart width in pixels
        height: Chart height in pixels
        centered: If True, center positions on the region
        strand: Strand for centering ('+' or '-')
        max_flank: Maximum flanking distance for centered plots
        min_basemod_qual: Minimum base modification quality score

    Returns:
        Path to the output HTML file
    """
    # Parse region if it's a string
    if isinstance(region, str):
        region = parse_region(region)

    # Set default output file
    if output_file is None:
        output_file = "fiberplot.html"

    # Load the BAM file
    print(f"Loading BAM file: {bam_file}", file=sys.stderr)
    fiberbam = pyft.Fiberbam(str(bam_file))

    # Fetch and convert to dataframe
    print(f"Fetching region: {region[0]}:{region[1]}-{region[2]}", file=sys.stderr)

    if centered:
        print(f"Creating centered plot (strand={strand})", file=sys.stderr)
        df = utils.region_to_centered_df(
            fiberbam,
            region,
            strand=strand,
            max_flank=max_flank,
            min_basemod_qual=min_basemod_qual,
        )
        chart = plot.centered_chart(df, width=width, height=height)
    else:
        print("Creating region plot", file=sys.stderr)
        df = utils.region_to_df(fiberbam, region, min_basemod_qual=min_basemod_qual)
        # Add centering columns required by centered_chart
        df["centering_position"] = region[1]
        df["centering_strand"] = strand
        # Make sure numeric columns are properly typed (not object/list types)
        for col in ["start", "end", "qual"]:
            if col in df.columns:
                df[col] = df[col].astype("float64")
        chart = plot.centered_chart(
            df,
            width=width,
            height=height,
            auto_sort=True,
            default_start=region[1],
            default_end=region[2],
        )

    # Save the chart
    print(f"Saving chart to: {output_file}", file=sys.stderr)
    # Use inline data transformer to avoid vegafusion issues with certain data types
    import altair as alt

    with alt.data_transformers.enable("default"):
        chart.save(str(output_file))
    print(f"Successfully created: {output_file}", file=sys.stderr)

    return output_file


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description="Generate an Altair HTML visualization of fibers in a genomic region",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - plot a specific region
  fiberplot data.bam chr22:26354169-26354500 -o my_plot.html

  # Create a centered plot on the positive strand
  fiberplot data.bam chr22:26354169-26354170 --centered --strand + -o centered.html

  # Create a centered plot with custom dimensions and flanking region
  fiberplot data.bam chr22:26354169-26354170 --centered --max-flank 250 --width 1000 --height 600

  # Adjust minimum base modification quality threshold
  fiberplot data.bam chr22:26354169-26354500 --min-basemod-qual 200
        """,
    )

    parser.add_argument("bam", type=str, help="Path to indexed BAM file")

    parser.add_argument(
        "region",
        type=str,
        help="Genomic region in format chr:start-end (e.g., chr22:26354169-26354500)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="fiberplot.html",
        help="Output HTML file path (default: fiberplot.html)",
    )

    parser.add_argument(
        "--width", type=int, default=600, help="Chart width in pixels (default: 600)"
    )

    parser.add_argument(
        "--height", type=int, default=800, help="Chart height in pixels (default: 800)"
    )

    parser.add_argument(
        "--centered",
        action="store_true",
        help="Create a centered plot (positions relative to region start/end)",
    )

    parser.add_argument(
        "--strand",
        type=str,
        choices=["+", "-"],
        default="+",
        help="Strand for centering (default: +)",
    )

    parser.add_argument(
        "--max-flank",
        type=int,
        default=None,
        help="Maximum flanking distance for centered plots (optional)",
    )

    parser.add_argument(
        "--min-basemod-qual",
        type=int,
        default=125,
        help="Minimum base modification quality score (default: 125)",
    )

    args = parser.parse_args()

    # Validate BAM file exists
    if not Path(args.bam).exists():
        print(f"Error: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    # Check for BAM index
    bam_path = Path(args.bam)
    if not bam_path.with_suffix(bam_path.suffix + ".bai").exists():
        print(f"Warning: BAM index (.bai) not found for {args.bam}", file=sys.stderr)
        print("The BAM file should be indexed for optimal performance", file=sys.stderr)

    try:
        # Parse and validate region
        region = parse_region(args.region)

        # Create the plot
        output_file = create_fiberplot(
            bam_file=args.bam,
            region=region,
            output_file=args.output,
            width=args.width,
            height=args.height,
            centered=args.centered,
            strand=args.strand,
            max_flank=args.max_flank,
            min_basemod_qual=args.min_basemod_qual,
        )

        print(f"\n✓ Plot successfully created: {output_file}")

    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error creating plot: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
