import pandas as pd
import altair as alt

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
            df.columns[df.columns.str.contains("module_|footprint_codes")],
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
