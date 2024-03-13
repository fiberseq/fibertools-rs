import altair as alt

alt.data_transformers.enable("vegafusion")

WIDE_COLUMNS = ["start", "end", "qual"]
M6A_COLOR = "#800080"
M5C_COLOR = "#8B4513"
NUC_COLOR = "#A9A9A9"
MSP_COLOR = "#9370db"


def centered_chart(
    in_df,
    height=800,
    width=800,
    m6a_color=M6A_COLOR,
    m5c_color=M5C_COLOR,
    nuc_color=NUC_COLOR,
    msp_color=MSP_COLOR,
    auto_sort=True,
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
    color_fp = alt.param(
        value="green", bind=alt.binding(input="color", name="footprinted")
    )
    color_not_fp = alt.param(
        value="lightgray", bind=alt.binding(input="color", name="not-footprinted")
    )

    domain = ["5mC", "m6a", "nuc", "msp", "footprinted", "not-footprinted"]
    range_ = [color_5mc, color_m6a, color_nuc, color_msp, color_fp, color_not_fp]
    opacity = dict(zip(domain, [1.0, 1.0, 0.05, 0.10, 0.1, 0.01]))

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

    # type_selection = alt.selection_point(fields=["type"], bind="legend")

    base = (
        alt.Chart(dfm)
        .encode(
            x="start:Q",
            x2="end:Q",
            color=alt.Color("type:O").scale(domain=domain, range=range_),
            y=alt.Y("fiber_name:O", sort=None),
            opacity=alt.Opacity("opacity"),
            # opacity=alt.condition(type_selection, "opacity", alt.value(0))
        )
        .transform_filter(
            selection,
        )
    )

    if "footprinted" in dfm["type"].unique():
        footprint_lines = _add_footprint_lines_to_centered_chart(dfm)
        chart = base.mark_rect() + footprint_lines.transform_filter(selection)
    else:
        chart = base.mark_rect()

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

    # sort on footprinted, m6a, and msp if they exist
    sort_cols = ["footprinted", "m6a", "msp"]
    # filter sort_cols to only include columns that exist in the dataframe
    sort_cols = [x for x in sort_cols if x in dfm.columns]
    # sort
    dfm.sort_values(
        ["chrom", "centering_position", "centering_strand"] + sort_cols,
        inplace=True,
        ascending=[True, True, True, False, False, False],
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
