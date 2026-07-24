"""microplotpy Shiny app: upload -> validate -> plot.

Dev use: from the microplotpy repo root, `pip install -e ".[app,plots,diversity]"`
then `shiny run app/app.py`.
"""
from __future__ import annotations

import io

import matplotlib

matplotlib.use("Agg")

import pandas as pd
from shiny import App, reactive, render, ui

from microplotpy import MicrobiomeData, validate, validate_assembly_summary, validate_contig_lengths, validate_mag
from microplotpy.plots import (
    alpha_diversity_plot,
    assembly_nx_plot,
    assembly_summary_barplot,
    beta_diversity_plot,
    function_barplot,
    function_heatmap,
    function_treemap,
    mag_quality_distribution,
    mag_quality_plot,
    taxa_barplot,
    taxa_bubbleplot,
    taxa_gradient_plot,
    taxa_heatmap,
    taxa_treemap,
)

NONE_CHOICE = "(none)"


def read_tsv(fileinfo) -> pd.DataFrame | None:
    if not fileinfo:
        return None
    return pd.read_csv(fileinfo[0]["datapath"], sep="\t")


def validation_html(report) -> ui.TagChild:
    if report is None:
        return ui.div("Upload the required file(s) to validate.", class_="alert alert-info")

    errs = report.errors
    warns = report.warnings
    parts = []
    if errs:
        parts.append(
            ui.div(
                ui.strong("Errors (fix before plotting):"),
                ui.tags.ul(*[ui.tags.li(f"[{f.file}] {f.field}: {f.message}") for f in errs]),
                class_="alert alert-danger",
            )
        )
    if warns:
        parts.append(
            ui.div(
                ui.strong("Warnings (plotting still allowed):"),
                ui.tags.ul(*[ui.tags.li(f"[{f.file}] {f.field}: {f.message}") for f in warns]),
                class_="alert alert-warning",
            )
        )
    if not errs and not warns:
        parts.append(ui.div("Input looks good.", class_="alert alert-success"))
    return ui.TagList(*parts)


def numeric_cols(df: pd.DataFrame, exclude: str) -> list[str]:
    cols = [c for c in df.columns if c != exclude]
    return [c for c in cols if pd.to_numeric(df[c], errors="coerce").notna().all()]


def none_or(value):
    return None if value in (None, NONE_CHOICE) else value


def download_controls(prefix: str):
    """Width/height (inches) + format (PNG/SVG) inputs + download button,
    shared by all four tabs."""
    return ui.TagList(
        ui.row(
            ui.column(6, ui.input_numeric(f"{prefix}_dl_width", "Width (in)", value=8, min=2, step=0.5)),
            ui.column(6, ui.input_numeric(f"{prefix}_dl_height", "Height (in)", value=6, min=2, step=0.5)),
        ),
        ui.input_radio_buttons(f"{prefix}_dl_format", "Format", {"png": "PNG", "svg": "SVG"}, inline=True),
        ui.download_button(f"{prefix}_download", "Download plot"),
    )


def save_fig_bytes(fig, width: float, height: float, fmt: str) -> bytes:
    """Resize a matplotlib Figure to user-specified dimensions and return it
    encoded in the requested format -- the same call the download handlers
    of every tab make."""
    fig.set_size_inches(width, height)
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf.read()


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

def taxonomy_panel():
    return ui.nav_panel(
        "Taxonomy",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_file("tax_feature_table", "Feature table (.tsv)", accept=[".tsv"]),
                ui.input_file("tax_taxonomy", "Taxonomy (.tsv)", accept=[".tsv"]),
                ui.input_file("tax_metadata", "Metadata (.tsv)", accept=[".tsv"]),
                ui.output_ui("tax_validation"),
                ui.hr(),
                ui.output_ui("tax_controls"),
            ),
            ui.output_plot("tax_plot", height="600px"),
            download_controls("tax"),
        ),
    )


def function_panel():
    return ui.nav_panel(
        "Functional profile",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_file("func_gene_counts", "Gene/KO count table (.tsv)", accept=[".tsv"]),
                ui.input_file("func_annotation", "Function annotation (.tsv)", accept=[".tsv"]),
                ui.output_ui("func_validation"),
                ui.hr(),
                ui.output_ui("func_controls"),
            ),
            ui.output_plot("func_plot", height="600px"),
            download_controls("func"),
        ),
    )


def mag_panel():
    return ui.nav_panel(
        "MAG quality",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_file("mag_table_file", "MAG quality table (.tsv, CheckM/CheckM2)", accept=[".tsv"]),
                ui.output_ui("mag_validation"),
                ui.hr(),
                ui.output_ui("mag_controls"),
            ),
            ui.output_plot("mag_plot", height="600px"),
            download_controls("mag"),
        ),
    )


def assembly_panel():
    return ui.nav_panel(
        "Assembly QC",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_select("asm_plot_type", "Plot type", {"nx": "Nx curve", "summary": "Summary statistic"}),
                ui.panel_conditional(
                    "input.asm_plot_type === 'nx'",
                    ui.input_file("asm_contig_lengths", "Contig lengths (.tsv)", accept=[".tsv"]),
                ),
                ui.panel_conditional(
                    "input.asm_plot_type === 'summary'",
                    ui.input_file("asm_summary_file", "Assembly summary (.tsv, QUAST report)", accept=[".tsv"]),
                ),
                ui.output_ui("asm_validation"),
                ui.hr(),
                ui.output_ui("asm_controls"),
            ),
            ui.output_plot("asm_plot", height="600px"),
            download_controls("asm"),
        ),
    )


app_ui = ui.page_navbar(
    taxonomy_panel(),
    function_panel(),
    mag_panel(),
    assembly_panel(),
    title="microplotpy",
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def server(input, output, session):
    # ---- Taxonomy ----

    @reactive.calc
    def tax_data():
        ft = read_tsv(input.tax_feature_table())
        tax = read_tsv(input.tax_taxonomy())
        meta = read_tsv(input.tax_metadata())
        if ft is None or tax is None or meta is None:
            return None
        return MicrobiomeData(feature_table=ft, taxonomy=tax, metadata=meta)

    @reactive.calc
    def tax_report():
        d = tax_data()
        return validate(d) if d is not None else None

    @render.ui
    def tax_validation():
        return validation_html(tax_report())

    @render.ui
    def tax_controls():
        d = tax_data()
        r = tax_report()
        if d is None or r is None or not r.is_valid:
            return None
        tax_cols = [c for c in d.taxonomy.columns if c != "Feature_ID"]
        meta_cols = [c for c in d.metadata.columns if c != "Sample_ID"]
        nmc = numeric_cols(d.metadata, "Sample_ID")

        return ui.TagList(
            ui.input_select("tax_plot_type", "Plot type", {
                "barplot": "Barplot", "heatmap": "Heatmap", "bubbleplot": "Bubbleplot",
                "gradient": "Gradient", "alpha": "Alpha diversity", "beta": "Beta diversity",
                "treemap": "Treemap",
            }),
            ui.panel_conditional(
                "input.tax_plot_type !== 'alpha' && input.tax_plot_type !== 'beta'",
                ui.input_select("tax_rank", "Rank", tax_cols, selected="Genus" if "Genus" in tax_cols else tax_cols[0]),
                ui.input_select("tax_group_rank", "Group rank (upper level)", [NONE_CHOICE] + tax_cols,
                                 selected="Phylum" if "Phylum" in tax_cols else NONE_CHOICE),
                ui.input_numeric("tax_top_n", "Top N taxa", value=10, min=1),
            ),
            ui.panel_conditional(
                "input.tax_plot_type === 'gradient'",
                ui.input_select("tax_gradient_var", "Gradient variable (numeric)", nmc),
            ),
            ui.panel_conditional(
                "input.tax_plot_type === 'alpha' || input.tax_plot_type === 'beta'",
                ui.input_select("tax_group_var", "Group variable", meta_cols),
            ),
            ui.panel_conditional(
                "input.tax_plot_type === 'beta'",
                ui.input_select("tax_method", "Dissimilarity", ["bray", "jaccard"]),
                ui.input_select("tax_ordination", "Ordination", ["pcoa", "nmds"]),
            ),
        )

    @reactive.calc
    def tax_plot_obj():
        d = tax_data()
        r = tax_report()
        if d is None or r is None or not r.is_valid:
            return None
        pt = input.tax_plot_type()
        gr = none_or(input.tax_group_rank()) if pt not in ("alpha", "beta") else None

        if pt == "barplot":
            return taxa_barplot(d, rank=input.tax_rank(), group_rank=gr, top_n=input.tax_top_n())
        if pt == "heatmap":
            return taxa_heatmap(d, rank=input.tax_rank(), top_n=input.tax_top_n())
        if pt == "bubbleplot":
            return taxa_bubbleplot(d, rank=input.tax_rank(), group_rank=gr, top_n=input.tax_top_n())
        if pt == "gradient":
            return taxa_gradient_plot(d, gradient_var=input.tax_gradient_var(), rank=input.tax_rank(),
                                       group_rank=gr, top_n=input.tax_top_n())
        if pt == "alpha":
            return alpha_diversity_plot(d, group_var=input.tax_group_var())
        if pt == "beta":
            return beta_diversity_plot(d, group_var=input.tax_group_var(), method=input.tax_method(),
                                        ordination=input.tax_ordination())
        if pt == "treemap":
            return taxa_treemap(d, rank=input.tax_rank(), group_rank=gr, top_n=input.tax_top_n())
        return None

    @render.plot(bbox_inches="tight")
    def tax_plot():
        return tax_plot_obj()

    @render.download(filename=lambda: f"microplotpy_{input.tax_plot_type()}.{input.tax_dl_format()}")
    def tax_download():
        fig = tax_plot_obj()
        if fig is not None:
            yield save_fig_bytes(fig, input.tax_dl_width(), input.tax_dl_height(), input.tax_dl_format())

    # ---- Functional profile ----

    @reactive.calc
    def func_tables():
        ft = read_tsv(input.func_gene_counts())
        fa = read_tsv(input.func_annotation())
        if ft is None or fa is None:
            return None
        return ft, fa

    @reactive.calc
    def func_report():
        t = func_tables()
        if t is None:
            return None
        ft, fa = t
        stub_meta = pd.DataFrame({"Sample_ID": [c for c in ft.columns if c != "Feature_ID"]})
        return validate(MicrobiomeData(feature_table=ft, taxonomy=fa, metadata=stub_meta), required_ranks=())

    @render.ui
    def func_validation():
        return validation_html(func_report())

    @render.ui
    def func_controls():
        t = func_tables()
        r = func_report()
        if t is None or r is None or not r.is_valid:
            return None
        _, fa = t
        ac = [c for c in fa.columns if c != "Feature_ID"]
        return ui.TagList(
            ui.input_select("func_plot_type", "Plot type",
                             {"barplot": "Barplot", "heatmap": "Heatmap", "treemap": "Treemap"}),
            ui.input_select("func_rank", "Rank", ac, selected="KEGG_ko" if "KEGG_ko" in ac else ac[0]),
            ui.input_select("func_group_rank", "Group rank (upper level)", [NONE_CHOICE] + ac,
                             selected="COG_category" if "COG_category" in ac else NONE_CHOICE),
            ui.input_numeric("func_top_n", "Top N", value=10, min=1),
        )

    @reactive.calc
    def func_plot_obj():
        t = func_tables()
        r = func_report()
        if t is None or r is None or not r.is_valid:
            return None
        ft, fa = t
        pt = input.func_plot_type()
        gr = none_or(input.func_group_rank())

        if pt == "barplot":
            return function_barplot(ft, fa, rank=input.func_rank(), group_rank=gr,
                                     top_n=input.func_top_n(), nested_legend=gr is not None)
        if pt == "heatmap":
            return function_heatmap(ft, fa, rank=input.func_rank(), top_n=input.func_top_n())
        if pt == "treemap":
            return function_treemap(ft, fa, rank=input.func_rank(), group_rank=gr, top_n=input.func_top_n())
        return None

    @render.plot(bbox_inches="tight")
    def func_plot():
        return func_plot_obj()

    @render.download(filename=lambda: f"microplotpy_function_{input.func_plot_type()}.{input.func_dl_format()}")
    def func_download():
        fig = func_plot_obj()
        if fig is not None:
            yield save_fig_bytes(fig, input.func_dl_width(), input.func_dl_height(), input.func_dl_format())

    # ---- MAG quality ----

    @reactive.calc
    def mag_table():
        return read_tsv(input.mag_table_file())

    @reactive.calc
    def mag_report():
        mt = mag_table()
        return validate_mag(mt) if mt is not None else None

    @render.ui
    def mag_validation():
        return validation_html(mag_report())

    @render.ui
    def mag_controls():
        mt = mag_table()
        r = mag_report()
        if mt is None or r is None or not r.is_valid:
            return None
        cols = [c for c in mt.columns if c not in ("Name", "Completeness", "Contamination")]
        return ui.TagList(
            ui.input_select("mag_plot_type", "Plot type",
                             {"scatter": "Completeness vs contamination", "distribution": "Quality distribution"}),
            ui.panel_conditional(
                "input.mag_plot_type === 'scatter'",
                ui.input_select("mag_size_col", "Size by", [NONE_CHOICE] + cols,
                                 selected="Genome_Size" if "Genome_Size" in cols else NONE_CHOICE),
                ui.input_select("mag_color_col", "Color by", [NONE_CHOICE] + cols),
            ),
        )

    @reactive.calc
    def mag_plot_obj():
        mt = mag_table()
        r = mag_report()
        if mt is None or r is None or not r.is_valid:
            return None
        if input.mag_plot_type() == "scatter":
            return mag_quality_plot(mt, size_col=none_or(input.mag_size_col()), color_col=none_or(input.mag_color_col()))
        return mag_quality_distribution(mt)

    @render.plot(bbox_inches="tight")
    def mag_plot():
        return mag_plot_obj()

    @render.download(filename=lambda: f"microplotpy_mag_{input.mag_plot_type()}.{input.mag_dl_format()}")
    def mag_download():
        fig = mag_plot_obj()
        if fig is not None:
            yield save_fig_bytes(fig, input.mag_dl_width(), input.mag_dl_height(), input.mag_dl_format())

    # ---- Assembly QC ----

    @reactive.calc
    def asm_contigs():
        return read_tsv(input.asm_contig_lengths())

    @reactive.calc
    def asm_summary():
        return read_tsv(input.asm_summary_file())

    @reactive.calc
    def asm_report():
        if input.asm_plot_type() == "nx":
            c = asm_contigs()
            return validate_contig_lengths(c) if c is not None else None
        s = asm_summary()
        return validate_assembly_summary(s) if s is not None else None

    @render.ui
    def asm_validation():
        return validation_html(asm_report())

    @render.ui
    def asm_controls():
        r = asm_report()
        if r is None or not r.is_valid or input.asm_plot_type() != "summary":
            return None
        s = asm_summary()
        cols = [c for c in s.columns if c != "Assembly_ID"]
        return ui.input_select("asm_stat_col", "Statistic", cols, selected="N50" if "N50" in cols else cols[0])

    @reactive.calc
    def asm_plot_obj():
        r = asm_report()
        if r is None or not r.is_valid:
            return None
        if input.asm_plot_type() == "nx":
            return assembly_nx_plot(asm_contigs())
        return assembly_summary_barplot(asm_summary(), stat_col=input.asm_stat_col())

    @render.plot(bbox_inches="tight")
    def asm_plot():
        return asm_plot_obj()

    @render.download(filename=lambda: f"microplotpy_assembly_{input.asm_plot_type()}.{input.asm_dl_format()}")
    def asm_download():
        fig = asm_plot_obj()
        if fig is not None:
            yield save_fig_bytes(fig, input.asm_dl_width(), input.asm_dl_height(), input.asm_dl_format())


app = App(app_ui, server)
