import pandas as pd
import altair as alt
import polars as pl
import numpy as np
import polars.selectors as cs

alt.data_transformers.enable("vegafusion")

WIDE_COLUMNS = ["start", "end", "qual"]
M6A_COLOR = "#800080"
M5C_COLOR = "#8B4513"
NUC_COLOR = "#A9A9A9"
MSP_COLOR = "#9370db"


def empty_data_dict():
    return {
        "chrom": [],
        "fiber_start": [],
        "fiber_end": [],
        "fiber_name": [],
        "strand": [],
        "type": [],
        "start": [],
        "end": [],
        "qual": [],
    }
    

def split_to_ints(df, col, sep=",", trim=True):
    """Split a columns with list of ints separated by
    "sep" into a numpy array quickly.

    Args:
        df (dataframe): dataframe that is like a bed12 file.
        col (str): column name within the dataframe to split up.
        sep (str, optional): Defaults to ",".
        trim (bool, optional): Remove the first and last call from the bed12 (removes bookendings). Defaults to True.

    Returns:
        column: New column that is a list of numpy array of ints.
    """
    if trim:
        return df[col].map_elements(
            lambda x: np.fromstring(x, sep=sep, dtype=np.int32)[1:-1],
            return_dtype=pl.datatypes.List(pl.datatypes.Int32),
        )
    return df[col].map_elements(
        lambda x: np.fromstring(x, sep=sep, dtype=np.int32),
        return_dtype=pl.datatypes.List(pl.datatypes.Int32),
    )



def read_extract_all(f: str, pandas=False, n_rows=None, long=False):
    """Read a table made with fibertools-rs. Specifically `ft extract --all`.
    Args:
        f (str): File path to the table. Can be compressed.
    Returns:
        pl.DataFrame: Dataframe of the table.
    """
    cols_with_lists = [
        "nuc_starts",
        "nuc_lengths",
        "ref_nuc_starts",
        "ref_nuc_lengths",
        "msp_starts",
        "msp_lengths",
        "fire",
        "ref_msp_starts",
        "ref_msp_lengths",
        "m6a",
        "ref_m6a",
        "m6a_qual",
        "5mC",
        "ref_5mC",
        "5mC_qual",
    ]
    df = pl.read_csv(
        f,
        separator="\t",
        n_rows=n_rows,
        null_values=["."],
        comment_prefix=None,
    )
    # clean up comment char
    df.columns = list(map(lambda x: x.strip("#"), df.columns))
    if df.shape[0] > 0:
        for col in cols_with_lists:
            col_index = df.columns.index(col)
            df.replace_column(col_index, split_to_ints(df, col, trim=False))
    
    if long:
        explode_sets = {
            "m6a":["m6a", "ref_m6a", "m6a_qual"],
            "5mC": ["5mC", "ref_5mC", "5mC_qual"],
            "msp": ["msp_starts", "msp_lengths", "ref_msp_starts", "ref_msp_lengths", "fire"],
            "nuc": ["nuc_starts", "nuc_lengths", "ref_nuc_starts", "ref_nuc_lengths"],
        }
        def my_rename(col):
            if col == "m6a" or col == "5mC" or col == "msp_starts" or col == "nuc_starts":
                return "long_start"
            if col == "ref_m6a" or col == "ref_5mC" or col == "ref_msp_starts" or col == "ref_nuc_starts":
                return "long_ref_start"
            if col == "m6a_qual" or col == "5mC_qual" or col == "fire":
                return "long_qual"
            if col == "msp_lengths" or col == "nuc_lengths":
                return "zlength"
            if col == "ref_msp_lengths" or col == "ref_nuc_lengths":
                return "ref_zlength"
            else:
                return col
            
        
        all_explode_cols = [s for _t, sublist in explode_sets.items() for s in sublist] 
        dfs = []
        for cur_type, s in explode_sets.items():
            not_in_s = [x for x in all_explode_cols if x not in s]
            tdf = df.drop(not_in_s).explode(s).with_columns(
                pl.lit(cur_type).alias("type")
            ).rename(my_rename)
            if cur_type == "m6a" or cur_type == "5mC":
                tdf = tdf.with_columns(
                    zlength=1,
                    ref_zlength=1,
                )
            if cur_type == "nuc":
                tdf = tdf.with_columns(
                    long_qual=0,
                )
            tdf = tdf.with_columns(
                long_end=pl.col("long_start") + pl.col("zlength"),
                long_ref_end=pl.col("long_ref_start") + pl.col("ref_zlength"),
            ).drop(
                cs.contains("zlength")
            )
            tdf=tdf.select(sorted(tdf.columns))
            dfs.append(tdf)
        # concat all of the dataframes
        df = pl.concat(dfs)
    if pandas:
        df = pd.DataFrame(df.to_dicts())
    return df



def footprint_code_to_vec(footprint_code, n_modules):
    """
    Convert a footprint code to a boolean vector of length n_modules
    """
    rtn = []
    for i in range(1, n_modules + 1):
        val = (footprint_code & (1 << i)) > 0
        rtn.append(val)
    return rtn


def read_footprint_table(f, long=False):
    """
    Read a ft-footprint bed into a pandas dataframe
    """
    wide_cols = ["footprint_codes", "fire_quals", "fiber_names"]
    df = pd.read_csv(f, sep="\t")
    # filter out columns that do not contain spanning fibers
    df = df[df["n_spanning_fibers"] > 0]
    # split the wide columns
    for col in wide_cols:
        df[col] = df[col].str.split(",")
        if col != "fiber_names":
            df[col] = df[col].apply(lambda x: [int(i) for i in x])

    # infer the number of modules form the column names
    n_modules = df.columns.str.contains("module:").sum()
    # save the names of the module columns to a list
    module_names = df.columns[df.columns.str.contains("module:")]
    df["n_modules"] = n_modules

    if long:
        df = df.explode(wide_cols)
        df["has_spanning_msp"] = df["footprint_codes"].apply(
            lambda x: (x & (1 << 0)) > 0
        )
        df["footprinted_modules"] = df["footprint_codes"].apply(
            lambda x: footprint_code_to_vec(x, n_modules)
        )
        # drop the module columns
        df.drop(
            df.columns[df.columns.str.contains("module_")],
            axis=1,
            inplace=True,
        )
        # split footprinted_modules into separate columns
        df[module_names] = pd.DataFrame(
            df["footprinted_modules"].tolist(), index=df.index
        )
        df.drop("footprinted_modules", axis=1, inplace=True)
        # rename fiber_names to fiber_name
        df.rename(
            columns={"fiber_names": "fiber_name", "fire_quals": "fire_qual"},
            inplace=True,
        )

    # get the right int columns
    df = df.infer_objects()
    df.rename(
        columns={"#chrom": "chrom", "start": "motif_start", "end": "motif_end"},
        inplace=True,
    )
    return df


def read_and_center_footprint_table(f):
    """
    Read a ft-footprint bed into a pandas dataframe and reformat the data to reflect ft-center output
    """
    df = read_footprint_table(f, long=True)
    # add how many modules are present in each row
    module_columns = [c for c in df.columns if c.startswith("module:")]
    df["n_footprints"] = df[module_columns].sum(axis=1)
    df["first_footprint"] = df[module_columns].idxmax(axis=1)

    # example of reading in a footprinting table
    not_module_columns = [col for col in df.columns if not col.startswith("module:")]

    dfm = df.melt(
        ignore_index=False, id_vars=not_module_columns, value_name="footprinted"
    ).reset_index()
    dfm[["module", "start", "end"]] = dfm["variable"].str.split(":|-", n=2, expand=True)
    dfm["start"] = dfm["start"].astype(int)
    dfm["end"] = dfm["end"].astype(int)
    dfm = dfm.infer_objects()

    # rename the columns to match the centering output style
    dfm["centering_position"] = dfm["motif_start"]
    # swap start and end if the footprint is on the reverse strand
    dfm.loc[dfm.strand == "-", "centering_position"] = (
        dfm.loc[dfm.strand == "-", "motif_end"] - 1
    )
    dfm["centering_strand"] = dfm["strand"]
    dfm["type"] = "not-footprinted"
    dfm.loc[dfm.has_spanning_msp & dfm.footprinted, "type"] = "footprinted"

    dfm.drop(
        columns=[
            "module",
            "variable",
            "n_modules",
            "index",
            "first_footprint",
            "n_footprints",
            "n_spanning_fibers",
            "n_spanning_msps",
            "n_overlapping_nucs",
        ],
        inplace=True,
        axis=1,
    )
    return dfm


def read_center_table(f):
    """
    Read a ft-center bed into a pandas dataframe
    """
    df = pd.read_csv(f, sep="\t")
    df = df.infer_objects()
    return df


def _add_fiber_to_data_dict(
    fiber,
    data_dict,
    min_fire_qual=200,
):
    # sub function to add fiber data to a dictionary
    def add_standard_columns():
        data_dict["chrom"].append(fiber.chrom)
        data_dict["fiber_start"].append(fiber.start)
        data_dict["fiber_end"].append(fiber.end)
        data_dict["strand"].append(fiber.strand)
        data_dict["fiber_name"].append(fiber.qname)

    def base_mod_end(list_of_starts):
        # add one to the vaules as long as they are not None
        return [x + 1 if x is not None else None for x in list_of_starts]

    # for msp
    add_standard_columns()
    data_dict["type"].append("msp")
    data_dict["start"].append(fiber.msp.reference_starts)
    data_dict["end"].append(fiber.msp.reference_ends)
    data_dict["qual"].append(fiber.msp.qual)

    # for FIRE (high-quality MSPs)
    fire_starts = []
    fire_ends = []
    fire_quals = []
    for start, end, qual in zip(fiber.msp.reference_starts, fiber.msp.reference_ends, fiber.msp.qual):
        if qual >= min_fire_qual:
            fire_starts.append(start)
            fire_ends.append(end)
            fire_quals.append(qual)

    add_standard_columns()
    data_dict["type"].append("fire")
    data_dict["start"].append(fire_starts)
    data_dict["end"].append(fire_ends)
    data_dict["qual"].append(fire_quals)

    # for nuc
    add_standard_columns()
    data_dict["type"].append("nuc")
    data_dict["start"].append(fiber.nuc.reference_starts)
    data_dict["end"].append(fiber.nuc.reference_ends)
    data_dict["qual"].append(fiber.nuc.qual)

    # for m6a
    add_standard_columns()
    data_dict["type"].append("m6a")
    data_dict["start"].append(fiber.m6a.reference_starts)
    data_dict["end"].append(fiber.m6a.get_reference_ends())
    data_dict["qual"].append(fiber.m6a.ml)

    # for 5mC
    add_standard_columns()
    data_dict["type"].append("5mC")
    data_dict["start"].append(fiber.cpg.reference_starts)
    data_dict["end"].append(fiber.cpg.get_reference_ends())
    data_dict["qual"].append(fiber.cpg.ml)


def region_to_centered_df(
    fiberbam, region, strand="+", max_flank=None, min_basemod_qual=125
):
    """
    Takes a fiberbam and a region and returns a pandas dataframe with reference centered positions
    """
    data_dict = empty_data_dict()
    for fiber in fiberbam.center(
        region[0], start=region[1], end=region[2], strand=strand
    ):
        _add_fiber_to_data_dict(fiber, data_dict)
    df = pd.DataFrame.from_dict(data_dict).explode(WIDE_COLUMNS)
    df["centering_position"] = region[1]
    if strand == "-":
        df["centering_position"] = region[2] - 1
    df["centering_strand"] = strand

    # trim the dataframe to only include fibers that overlap
    if max_flank is not None:
        df = df[(df.start < +max_flank) & (df.end > -max_flank)]

    is_basemod = df.type.isin(["m6a", "5mC"])
    df = df[~(is_basemod & (df.qual < min_basemod_qual))]
    return df


def region_to_df(fiberbam, region, min_basemod_qual=125):
    """
    Takes a fiberbam and a region and returns a pandas dataframe with fibers
    that overlap the region
    """
    data_dict = empty_data_dict()
    for fiber in fiberbam.fetch(
        region[0],
        start=region[1],
        end=region[2],
    ):
        _add_fiber_to_data_dict(fiber, data_dict)
    df = pd.DataFrame.from_dict(data_dict).explode(WIDE_COLUMNS)
    is_basemod = df.type.isin(["m6a", "5mC"])
    df = df[~(is_basemod & (df.qual < min_basemod_qual))]
    return df
