import altair as alt
import pandas as pd

alt.data_transformers.enable("vegafusion")

WIDE_COLUMNS = ["start", "end", "qual"]
M6A_COLOR = "#800080"
M5C_COLOR = "#8B4513"
NUC_COLOR = "#A9A9A9"
MSP_COLOR = "#9370db"
FIRE_COLOR = "#FF0000"  # Red for FIRE elements (matches fibertools-rs)


def centered_chart(
    in_df,
    height=800,
    width=800,
    m6a_color=M6A_COLOR,
    m5c_color=M5C_COLOR,
    nuc_color=NUC_COLOR,
    msp_color=MSP_COLOR,
    fire_color=FIRE_COLOR,
    auto_sort=True,
    default_start=None,
    default_end=None,
):
    """Make an altair chart from a dataframe of centered positions
    Args:
        in_df (_type_): _description_
        m6a_color (_type_, optional): _description_. Defaults to M6A_COLOR.
        m5c_color (_type_, optional): _description_. Defaults to M5C_COLOR.
        nuc_color (_type_, optional): _description_. Defaults to NUC_COLOR.
        msp_color (_type_, optional): _description_. Defaults to MSP_COLOR.
        auto_sort (bool, optional): _description_. Defaults to True.

    Returns:
        altair.Chart
    """
    dfm = (
        in_df.copy()[
            [
                "chrom",
                "centering_position",
                "centering_strand",
                "start",
                "end",
                "type",
                "fiber_name",
            ]
        ]
        .dropna()
        .reset_index(drop=True)
        .infer_objects()
    )  # .query("centered_position_type != '5mC'")

    if auto_sort:
        dfm = _auto_sort_centered_chart(dfm)

    # get the fiber names in order they appear in the dfm
    # fiber_names = dfm["fiber_name"].unique()

    dfm["region"] = (
        dfm["chrom"]
        + ":"
        + dfm["centering_position"].astype(str)
        + " "
        + dfm["centering_strand"]
    )

    # set the colors
    color_m6a = alt.param(value=m6a_color, bind=alt.binding(input="color", name="m6a"))
    color_5mc = alt.param(value=m5c_color, bind=alt.binding(input="color", name="5mC"))
    color_nuc = alt.param(value=nuc_color, bind=alt.binding(input="color", name="nuc"))
    color_msp = alt.param(value=msp_color, bind=alt.binding(input="color", name="msp"))
    color_fire = alt.param(
        value=fire_color, bind=alt.binding(input="color", name="fire")
    )
    color_fp = alt.param(
        value="green", bind=alt.binding(input="color", name="footprinted")
    )
    color_not_fp = alt.param(
        value="lightgray", bind=alt.binding(input="color", name="not-footprinted")
    )

    domain = ["5mC", "m6a", "nuc", "msp", "fire", "footprinted", "not-footprinted"]
    range_ = [
        color_5mc,
        color_m6a,
        color_nuc,
        color_msp,
        color_fire,
        color_fp,
        color_not_fp,
    ]
    opacity = dict(zip(domain, [1.0, 1.0, 0.05, 0.10, 0.8, 0.1, 0.01]))

    # add opacity column to the dataframe
    dfm = dfm.assign(opacity=dfm["type"].map(opacity))

    input_dropdown = alt.binding_select(
        # Add the empty selection which shows all when clicked
        options=dfm.region.unique(),
        name="Region: ",
    )

    selection = alt.selection_point(
        fields=["region"],
        bind=input_dropdown,
        value=dfm.region[0],
    )

    bind_range_w = alt.binding_range(min=200, max=1600, name="Chart width: ")
    param_width = alt.param("width", bind=bind_range_w)
    bind_range_h = alt.binding_range(min=200, max=1600, name="Chart height: ")
    param_height = alt.param("height", bind=bind_range_h)

    # Add legend selection for toggling visibility (multi-select with all selected by default)
    type_selection = alt.selection_point(
        fields=["type"],
        bind="legend",
        toggle="true",
        value=[{"type": t} for t in domain],  # All types selected by default
    )

    base = (
        alt.Chart(dfm)
        .encode(
            x="start:Q",
            x2="end:Q",
            color=alt.Color("type:O").scale(domain=domain, range=range_),
            y=alt.Y("fiber_name:O", sort=None),
            opacity=alt.condition(
                type_selection, alt.Opacity("opacity:Q", legend=None), alt.value(0.0)
            ),
        )
        .transform_filter(
            selection,
        )
        .add_params(type_selection)
    )

    if "footprinted" in dfm["type"].unique():
        footprint_lines = _add_footprint_lines_to_centered_chart(dfm)
        chart = base.mark_rect() + footprint_lines.transform_filter(selection)
    else:
        chart = base.mark_rect()

    # Configure legend to show opaque symbols
    chart = chart.configure_legend(
        symbolOpacity=1.0,  # Make legend symbols fully opaque
        symbolStrokeWidth=0,
    )

    # if there is a default start and end, set the x axis to that range
    if default_start is not None and default_end is not None:
        chart = chart.encode(
            x=alt.X(
                "start:Q",
                scale=alt.Scale(domain=(default_start, default_end)),
            )
        )

    chart = (
        chart.properties(width=width, height=height)
        .add_params(
            selection,
            param_width,
            param_height,
            color_m6a,
            color_5mc,
            color_nuc,
            color_msp,
            color_fire,
            color_fp,
            color_not_fp,
        )
        .interactive()
    )

    return chart


def _auto_sort_centered_chart(dfm):
    z = (
        dfm.groupby(
            [
                "chrom",
                "fiber_name",
                "centering_position",
                "centering_strand",
                "type",
            ]
        )
        .size()
        .reset_index(name="count")
    )

    # combine the centered_position_type and count into a set of wide columns
    z = (
        z.pivot_table(
            index=["fiber_name", "centering_position", "centering_strand"],
            columns="type",
            values="count",
        )
        .reset_index()
        .fillna(0)
    )

    # join z with both_dfs
    dfm = dfm.merge(
        z, on=["fiber_name", "centering_position", "centering_strand"], how="left"
    )

    # sort on footprinted, fire, m6a, and msp if they exist
    sort_cols = ["footprinted", "fire", "m6a", "msp"]
    # filter sort_cols to only include columns that exist in the dataframe
    sort_cols = [x for x in sort_cols if x in dfm.columns]
    # sort
    base_ascending = [
        True,
        True,
        True,
    ]  # for chrom, centering_position, centering_strand
    sort_ascending = [False] * len(sort_cols)  # descending for all sort columns
    dfm.sort_values(
        ["chrom", "centering_position", "centering_strand"] + sort_cols,
        inplace=True,
        ascending=base_ascending + sort_ascending,
    )

    return dfm


def _add_footprint_lines_to_centered_chart(dfm):
    # draw vertical lines at the start and end of the footprints
    fp_df = dfm[["start", "end", "type", "region"]].query("type == 'footprinted'")
    # pivot the start and end columns into a long format
    fp_df = (
        fp_df.melt(
            id_vars=["type", "region"],
            value_vars=["start", "end"],
            value_name="position",
        )
        .reset_index(drop=True)
        .drop_duplicates()
    )

    footprint_lines = (
        alt.Chart(fp_df)
        .encode(
            x="position:Q",
            # y=alt.Y('fiber_name:O', sort=None),
        )
        .mark_rule(
            color="black",
            strokeWidth=0.75,
            strokeDash=[3, 3],
            opacity=0.75,
        )
    )
    return footprint_lines


def extract_chart(
    in_df,
    height=600,
    width=800,
    m6a_color=M6A_COLOR,
    m5c_color=M5C_COLOR,
    nuc_color=NUC_COLOR,
    msp_color=MSP_COLOR,
    max_fibers=None,
    ref=True,
):
    """Make an altair chart from a dataframe of centered positions
    Args:
        in_df (_type_): Pandas DataFrame from utils.read_extract_all with long=True.
        m6a_color (_type_, optional): _description_. Defaults to M6A_COLOR.
        m5c_color (_type_, optional): _description_. Defaults to M5C_COLOR.
        nuc_color (_type_, optional): _description_. Defaults to NUC_COLOR.
        msp_color (_type_, optional): _description_. Defaults to MSP_COLOR.
        max_fibers (_type_, optional): maximum number of fibers to plot. Defaults to None.
        ref (_type_, optional): Whether to use reference (ref) coordinates. Defaults to True.
    Returns:
        altair.Chart
    """
    if ref:
        start = "long_ref_start"
        end = "long_ref_end"
    else:
        start = "long_start"
        end = "long_end"
    dfm = (
        in_df.copy()[
            [
                start,
                end,
                "type",
                "fiber",
            ]
        ]
        .dropna()
        .reset_index(drop=True)
        .infer_objects()
    )  # .query("centered_position_type != '5mC'")

    # set the colors
    color_m6a = alt.param(value=m6a_color, bind=alt.binding(input="color", name="m6a"))
    color_5mc = alt.param(value=m5c_color, bind=alt.binding(input="color", name="5mC"))
    color_nuc = alt.param(value=nuc_color, bind=alt.binding(input="color", name="nuc"))
    color_msp = alt.param(value=msp_color, bind=alt.binding(input="color", name="msp"))

    domain = ["5mC", "m6a", "nuc", "msp"]
    range_ = [color_5mc, color_m6a, color_nuc, color_msp]
    opacity = dict(zip(domain, [1.0, 1.0, 0.5, 0.20]))

    # add opacity column to the dataframe
    dfm = dfm.assign(opacity=dfm["type"].map(opacity))

    dfm["long_end"] = dfm["long_end"]

    # convert fiber column to integer
    dfm["height"] = 0.4
    dfm.loc[dfm.type == "nuc", "height"] = 0.2
    dfm.loc[dfm.type == "msp", "height"] = 0.3
    dfm["y"] = pd.factorize(dfm["fiber"])[0] - dfm["height"]
    dfm["y2"] = pd.factorize(dfm["fiber"])[0] + dfm["height"]
    if max_fibers is not None:
        dfm = dfm.query("y < @max_fibers")

    bind_range_w = alt.binding_range(min=200, max=1600, name="Chart width: ")
    param_width = alt.param("width", bind=bind_range_w)
    bind_range_h = alt.binding_range(min=200, max=1600, name="Chart height: ")
    param_height = alt.param("height", bind=bind_range_h)

    base = alt.Chart(dfm).encode(
        x=alt.X(f"{start}:Q", scale=alt.Scale(domain=(500, 2500))),
        x2=f"{end}:Q",
        color=alt.Color("type:O").scale(domain=domain, range=range_),
        y="y:Q",
        y2="y2:Q",
        opacity=alt.Opacity("opacity"),
        # opacity=alt.condition(type_selection, "opacity", alt.value(0))
    )

    chart = base.mark_rect()
    chart = (
        chart.properties(width=width, height=height)
        .add_params(
            param_width,
            param_height,
            color_m6a,
            color_5mc,
            color_nuc,
            color_msp,
        )
        .interactive(bind_y=False)
    )

    return chart
