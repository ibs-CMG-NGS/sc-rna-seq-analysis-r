"""
scRNA-seq 3rd Analysis Pipeline
Usage: snakemake --profile profiles/default --configfile configs/config_<experiment>.yml --cores <N>
"""

import os
from pathlib import Path

configfile: "configs/template/config.yml"


def join_cfg(value, default="null"):
    """Render a list/string config value for a comma-separated CLI argument."""
    if value is None:
        return default
    if isinstance(value, (list, tuple)):
        return ",".join(map(str, value)) if value else default
    value = str(value).strip()
    return value if value else default


# ── Derived variables ──────────────────────────────────────────────────────────
OUTPUT     = config["output_dir"]
START_FROM = config.get("start_from", "input")
if START_FROM not in {"input", "loaded"}:
    raise ValueError("start_from must be 'input' or 'loaded'")
LOADED_RDS = f"{OUTPUT}/loaded/seurat_validated.rds"
LOAD_FLAG  = f"{OUTPUT}/loaded/load_complete.flag"
OPTIONAL   = config.get("optional_modules", {})
RUN_TRAJ   = OPTIONAL.get("trajectory", False)
RUN_COMM   = OPTIONAL.get("cell_communication", True)
RUN_LBL    = OPTIONAL.get("label_transfer", False)
# Core module toggles (default: keep previous behaviour = all on)
MODULES    = config.get("modules", {})
RUN_DE     = MODULES.get("pseudobulk_de", True)
RUN_COMP   = MODULES.get("composition", True)
RUN_ENR    = MODULES.get("enrichment", True)
# Exploratory gene-panel expression (works without replicates).
# Either optional_modules.exploration or exploration.enabled enables it.
RUN_EXPLORE = bool(OPTIONAL.get("exploration", False) or
                   config.get("exploration", {}).get("enabled", False))

EXPR_CFG   = config.get("expression_export", {})
VIEWER_CFG = config.get("viewer_export")
RUN_VIEWER = bool(VIEWER_CFG and VIEWER_CFG.get("enabled", False))

# Unified viewer export (preferred). One switch + one Import Folder; internally
# analysis (DE/GO/composition/MG) and expression summary are built separately.
if RUN_VIEWER:
    RUN_VIEWER_ANALYSIS   = bool(VIEWER_CFG.get("analysis", True))
    RUN_VIEWER_EXPRESSION = bool(VIEWER_CFG.get("expression", False))
    RUN_EXPORT     = False  # legacy analysis-only rule is superseded
    RUN_EXPR_EXPORT = RUN_VIEWER_EXPRESSION
    if not (RUN_VIEWER_ANALYSIS or RUN_VIEWER_EXPRESSION):
        raise ValueError("viewer_export.enabled is true but both analysis and expression are false")
    if RUN_VIEWER_ANALYSIS and not RUN_DE:
        raise ValueError("viewer_export.analysis requires modules.pseudobulk_de")
else:
    RUN_VIEWER_ANALYSIS   = False
    RUN_VIEWER_EXPRESSION = False
    # Backward compatibility: legacy separate exports still work when the
    # viewer_export block is absent.
    RUN_EXPORT      = config.get("export_seqviewer", {}).get("enabled", True)
    RUN_EXPR_EXPORT = EXPR_CFG.get("enabled", False)

RUN_EXPR_PLOTS = RUN_EXPR_EXPORT and EXPR_CFG.get("plots", True)
DIST_CFG = config.get("expression_distribution", {})
RUN_DISTRIBUTION = bool(RUN_VIEWER and DIST_CFG.get("enabled", False))

# Canonical CSV/Excel table layer. CSV tables are produced whenever the
# analysis/expression exporters run; these switches control merge + Excel.
TABLE_CFG        = config.get("table_export", {})
RUN_TABLE_EXPORT = bool(TABLE_CFG.get("enabled", True))
RUN_TABLE_ANALYSIS   = RUN_TABLE_EXPORT and RUN_DE
RUN_TABLE_EXPRESSION = RUN_TABLE_EXPORT and RUN_EXPR_EXPORT
RUN_TABLE_MERGE      = RUN_TABLE_EXPORT and (RUN_TABLE_ANALYSIS or RUN_TABLE_EXPRESSION)
RUN_TABLE_EXCEL      = RUN_TABLE_MERGE and bool(TABLE_CFG.get("excel", True))
RUN_ANALYSIS_EXPORT  = RUN_DE and (RUN_VIEWER_ANALYSIS or RUN_EXPORT or RUN_TABLE_EXPORT)

if RUN_ENR and not RUN_DE:
    raise ValueError("modules.enrichment requires modules.pseudobulk_de; "
                     "set modules.enrichment: false")

# Phase 0a: when enabled, load_input consumes the annotated object produced
# by label transfer instead of the raw input.
QUERY_INPUT = f"{OUTPUT}/label_transfer/annotated.rds" if RUN_LBL \
              else config["single_cell"]["input_path"]


# ── Helper: save config snapshot ──────────────────────────────────────────────
import shutil, yaml, datetime

def save_config_snapshot(config, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    out = os.path.join(output_dir, "config_used.yml")
    with open(out, "w") as f:
        yaml.dump(dict(config), f, default_flow_style=False)


onstart:
    if not workflow.dryrun:
        save_config_snapshot(config, OUTPUT)


# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand("{out}/label_transfer/label_transfer.flag", out=OUTPUT) if RUN_LBL else [],
        expand("{out}/pseudobulk_de/de_complete.flag",   out=OUTPUT) if RUN_DE else [],
        expand("{out}/composition/composition_complete.flag", out=OUTPUT) if RUN_COMP else [],
        expand("{out}/enrichment/enrichment_complete.flag",   out=OUTPUT) if RUN_ENR else [],
        expand("{out}/exploration/exploration_complete.flag", out=OUTPUT) if RUN_EXPLORE else [],
        expand("{out}/export/expression/expression_complete.flag", out=OUTPUT) if RUN_EXPR_EXPORT else [],
        expand("{out}/export/expression/plot_complete.flag", out=OUTPUT) if RUN_EXPR_PLOTS else [],
        expand("{out}/trajectory/trajectory_complete.flag",   out=OUTPUT) if RUN_TRAJ else [],
        expand("{out}/cell_communication/comm_complete.flag", out=OUTPUT) if RUN_COMM else [],
        expand("{out}/export/seqviewer/export_complete.flag",  out=OUTPUT) if (RUN_VIEWER or RUN_EXPORT) else [],
        expand("{out}/tables/tables_complete.flag",            out=OUTPUT) if RUN_TABLE_MERGE else [],
        expand("{out}/export/excel/excel_complete.flag",       out=OUTPUT) if RUN_TABLE_EXCEL else [],
        expand("{out}/report/report.html",                    out=OUTPUT),


# ── Phase 0a: Label transfer (optional, v2.2) ────────────────────────────────
if RUN_LBL:
    rule label_transfer:
        input:
            data = config["single_cell"]["input_path"],
        output:
            rds  = f"{OUTPUT}/label_transfer/annotated.rds",
            flag = f"{OUTPUT}/label_transfer/label_transfer.flag",
        log:
            f"logs/label_transfer.log",
        params:
            fmt        = config["single_cell"]["input_format"],
            method     = config["label_transfer"]["method"],
            label_col  = config["label_transfer"].get("label_column", "C66_named"),
            max_epochs = config["label_transfer"].get("max_epochs", 20),
            suffix     = config["label_transfer"].get("suffix", "query"),
            ct_col     = config["single_cell"]["cell_type_column"],
            prov_out   = f"{OUTPUT}/label_transfer/label_provenance.csv",
        shell:
            """
            Rscript scripts/00_label_transfer.R \
                --input          {input.data} \
                --format         {params.fmt} \
                --method         {params.method} \
                --label_column   "{params.label_col}" \
                --max_epochs     {params.max_epochs} \
                --suffix         "{params.suffix}" \
                --cell_type_col  "{params.ct_col}" \
                --output_rds     {output.rds} \
                --provenance_out {params.prov_out} \
                2>&1 | tee {log}
            touch {output.flag}
            """
else:
    rule label_transfer:
        output:
            flag = f"{OUTPUT}/label_transfer/label_transfer.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 1: Load & validate input ───────────────────────────────────────────
# start_from=loaded skips input_loader and uses the existing validated RDS.
if START_FROM != "loaded":
    rule load_input:
        input:
            data = QUERY_INPUT,
        output:
            rds  = LOADED_RDS,
            flag = LOAD_FLAG,
        log:
            f"logs/load_input.log",
        params:
            format         = config["single_cell"]["input_format"],
            cell_type_col  = config["single_cell"]["cell_type_column"],
            sample_col     = config["single_cell"]["sample_column"],
            condition_col  = config["single_cell"]["condition_column"],
            cell_types     = join_cfg(config["single_cell"].get("cell_types")),
            min_cells      = config["single_cell"]["min_cells_per_group"],
            counts_layer   = config["single_cell"].get("counts_layer", "counts"),
        shell:
            """
            Rscript scripts/input_loader.R \
                --input         {input.data} \
                --format        {params.format} \
                --cell_type_col {params.cell_type_col} \
                --sample_col    {params.sample_col} \
                --condition_col {params.condition_col} \
                --cell_types    "{params.cell_types}" \
                --min_cells     {params.min_cells} \
                --counts_layer  {params.counts_layer} \
                --output_rds    {output.rds} \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 2a: Pseudobulk DE ───────────────────────────────────────────────────
if RUN_DE:
    rule pseudobulk_de:
        input:
            rds       = LOADED_RDS,
        output:
            flag = f"{OUTPUT}/pseudobulk_de/de_complete.flag",
        log:
            f"logs/pseudobulk_de.log",
        params:
            method      = config["de_analysis"]["method"],
            ref_level   = config["de_analysis"].get("reference_level", ""),
            covariates  = ",".join(config["de_analysis"].get("covariates", [])),
            padj        = config["padj_cutoff"],
            log2fc      = config["log2fc_cutoff"],
            output_dir  = f"{OUTPUT}/pseudobulk_de",
        shell:
            """
            Rscript scripts/01_pseudobulk_de.R \
                --input        {input.rds} \
                --method       {params.method} \
                --ref_level    "{params.ref_level}" \
                --covariates   "{params.covariates}" \
                --padj_cutoff  {params.padj} \
                --log2fc_cutoff {params.log2fc} \
                --output_dir   {params.output_dir} \
                2>&1 | tee {log}
            touch {output.flag}
            """
else:
    rule pseudobulk_de:
        output:
            flag = f"{OUTPUT}/pseudobulk_de/de_complete.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 2b: Compositional analysis ─────────────────────────────────────────
if RUN_COMP:
    rule composition:
        input:
            rds       = LOADED_RDS,
        output:
            flag = f"{OUTPUT}/composition/composition_complete.flag",
        log:
            f"logs/composition.log",
        params:
            output_dir = f"{OUTPUT}/composition",
        shell:
            """
            Rscript scripts/01b_composition.R \
                --input      {input.rds} \
                --output_dir {params.output_dir} \
                2>&1 | tee {log}
            touch {output.flag}
            """
else:
    rule composition:
        output:
            flag = f"{OUTPUT}/composition/composition_complete.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 3: Enrichment ───────────────────────────────────────────────────────
if RUN_ENR:
    rule enrichment:
        input:
            flag = rules.pseudobulk_de.output.flag,
        output:
            flag = f"{OUTPUT}/enrichment/enrichment_complete.flag",
        log:
            f"logs/enrichment.log",
        params:
            species     = config["species"],
            go_onto     = ",".join(config["enrichment"]["go_ontologies"]),
            run_kegg    = str(config["enrichment"].get("run_kegg", True)).upper(),
            pval        = config["enrichment"]["pvalue_cutoff"],
            qval        = config["enrichment"]["qvalue_cutoff"],
            min_gs      = config["enrichment"].get("min_gs_size", 10),
            max_gs      = config["enrichment"].get("max_gs_size", 500),
            de_dir      = f"{OUTPUT}/pseudobulk_de",
            output_dir  = f"{OUTPUT}/enrichment",
        shell:
            """
            Rscript scripts/02_enrichment.R \
                --de_dir       {params.de_dir} \
                --species      {params.species} \
                --go_ontologies {params.go_onto} \
                --run_kegg     {params.run_kegg} \
                --pvalue_cutoff {params.pval} \
                --qvalue_cutoff {params.qval} \
                --min_gs_size  {params.min_gs} \
                --max_gs_size  {params.max_gs} \
                --output_dir   {params.output_dir} \
                2>&1 | tee {log}
            touch {output.flag}
            """
else:
    rule enrichment:
        output:
            flag = f"{OUTPUT}/enrichment/enrichment_complete.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 4: Trajectory (optional) ───────────────────────────────────────────
if RUN_TRAJ:
    rule trajectory:
        input:
            rds       = LOADED_RDS,
        output:
            flag = f"{OUTPUT}/trajectory/trajectory_complete.flag",
        log:
            f"logs/trajectory.log",
        params:
            method         = config["trajectory"]["method"],
            root_cell_type = config["trajectory"].get("root_cell_type", ""),
            root_cluster   = config["trajectory"].get("root_cluster", ""),
            output_dir     = f"{OUTPUT}/trajectory",
        shell:
            """
            Rscript scripts/03_trajectory.R \
                --input          {input.rds} \
                --method         {params.method} \
                --root_cell_type "{params.root_cell_type}" \
                --root_cluster   "{params.root_cluster}" \
                --output_dir     {params.output_dir} \
                2>&1 | tee {log}
            touch {output.flag}
            """
else:
    rule trajectory:
        output:
            flag = f"{OUTPUT}/trajectory/trajectory_complete.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 5: Cell-Cell Communication (optional) ───────────────────────────────
if RUN_COMM:
    rule cell_communication:
        input:
            rds       = LOADED_RDS,
        output:
            flag = f"{OUTPUT}/cell_communication/comm_complete.flag",
        log:
            f"logs/cell_communication.log",
        params:
            method      = config["cell_communication"]["method"],
            conditions  = ",".join(config["cell_communication"].get("conditions_to_compare", [])),
            species     = config["species"],
            output_dir  = f"{OUTPUT}/cell_communication",
        shell:
            """
            Rscript scripts/04_cell_communication.R \
                --input      {input.rds} \
                --method     {params.method} \
                --conditions "{params.conditions}" \
                --species    {params.species} \
                --output_dir {params.output_dir} \
                2>&1 | tee {log}
            touch {output.flag}
            """
else:
    rule cell_communication:
        output:
            flag = f"{OUTPUT}/cell_communication/comm_complete.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 7: Analysis tables + cmg-seqviewer staging ─────────────────────────
# 06 writes both the viewer parquet staging bundle and canonical CSV tables.
if RUN_ANALYSIS_EXPORT:
    rule export_analysis:
        input:
            de_flag   = rules.pseudobulk_de.output.flag,
            comp_flag = rules.composition.output.flag if RUN_COMP else [],
            enr_flag  = rules.enrichment.output.flag if RUN_ENR else [],
        output:
            manifest      = f"{OUTPUT}/export/analysis/seqviewer_manifest.json",
            index         = f"{OUTPUT}/export/analysis/export_index.tsv",
            flag          = f"{OUTPUT}/export/analysis/analysis_complete.flag",
            table_manifest = f"{OUTPUT}/tables/analysis/manifest.json" if RUN_TABLE_ANALYSIS else [],
        log:
            "logs/export_analysis.log",
        params:
            de_dir          = f"{OUTPUT}/pseudobulk_de",
            enrichment_dir  = f"{OUTPUT}/enrichment" if RUN_ENR else "",
            composition_dir = f"{OUTPUT}/composition" if RUN_COMP else "",
            output_dir      = f"{OUTPUT}/export/analysis",
            csv_dir         = f"{OUTPUT}/tables/analysis" if RUN_TABLE_ANALYSIS else "",
            schema_version  = config.get("export_seqviewer", {}).get("schema_version", "1.0"),
            table_schema_version = "1.0",
            species         = config["species"],
            padj            = config["padj_cutoff"],
            log2fc          = config["log2fc_cutoff"],
        shell:
            """
            Rscript scripts/06_export_seqviewer.R \
                --de_dir          "{params.de_dir}" \
                --enrichment_dir  "{params.enrichment_dir}" \
                --composition_dir "{params.composition_dir}" \
                --output_dir      "{params.output_dir}" \
                --csv_dir         "{params.csv_dir}" \
                --schema_version  "{params.schema_version}" \
                --table_schema_version "{params.table_schema_version}" \
                --species         "{params.species}" \
                --padj_cutoff     {params.padj} \
                --log2fc_cutoff   {params.log2fc} \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 7: cmg-seqviewer handoff export (Track B, E1) ───────────────────────
if RUN_EXPORT and RUN_DE:
    rule export_seqviewer:
        input:
            flag = rules.export_analysis.output.flag,
        output:
            manifest = f"{OUTPUT}/export/seqviewer/seqviewer_manifest.json",
            index    = f"{OUTPUT}/export/seqviewer/export_index.tsv",
            flag     = f"{OUTPUT}/export/seqviewer/export_complete.flag",
        log:
            "logs/export_seqviewer.log",
        params:
            analysis_dir   = f"{OUTPUT}/export/analysis",
            output_dir     = f"{OUTPUT}/export/seqviewer",
            schema_version = config.get("export_seqviewer", {}).get("schema_version", "1.0"),
        shell:
            """
            Rscript scripts/10_merge_viewer_bundle.R \
                --analysis_dir   "{params.analysis_dir}" \
                --expression_dir "" \
                --output_dir     "{params.output_dir}" \
                --schema_version "{params.schema_version}" \
                2>&1 | tee {log}
            touch {output.flag}
            """
elif not RUN_VIEWER:
    rule export_seqviewer:
        output:
            flag = f"{OUTPUT}/export/seqviewer/export_complete.flag",
        shell:
            "mkdir -p $(dirname {output.flag}) && touch {output.flag}"


# ── Phase 5b: Exploratory gene-panel expression (optional) ────────────────────
if RUN_EXPLORE:
    EXP = config.get("exploration", {})

    rule explore_expression:
        input:
            data = config["single_cell"]["input_path"],
        output:
            csv  = f"{OUTPUT}/exploration/gene_expression_by_celltype.csv",
            flag = f"{OUTPUT}/exploration/exploration_complete.flag",
        log:
            "logs/explore_expression.log",
        params:
            genes      = ",".join(EXP.get("genes", [])),
            group_by   = EXP.get("group_by", config["single_cell"]["cell_type_column"]),
            split_by   = EXP.get("split_by", ""),
            min_cells  = EXP.get("min_cells", 10),
            excl_col   = EXP.get("exclude_column", ""),
            excl_vals  = ",".join(EXP.get("exclude_values", [])),
            assay      = EXP.get("assay", "RNA"),
            top_n      = EXP.get("top_n", 20),
            output_dir = f"{OUTPUT}/exploration",
        shell:
            """
            Rscript scripts/05_explore_expression.R \
                --input        {input.data} \
                --genes        "{params.genes}" \
                --group_by     "{params.group_by}" \
                --split_by     "{params.split_by}" \
                --min_cells    {params.min_cells} \
                --exclude_col  "{params.excl_col}" \
                --exclude_vals "{params.excl_vals}" \
                --assay        "{params.assay}" \
                --top_n        {params.top_n} \
                --output_dir   {params.output_dir} \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 5c: Full expression export (optional) ──────────────────────────────
if RUN_EXPR_EXPORT:
    EXPR = config.get("expression_export", {})

    rule export_expression:
        input:
            data = QUERY_INPUT,
        output:
            flag     = f"{OUTPUT}/export/expression/expression_complete.flag",
            manifest = f"{OUTPUT}/tables/expression/manifest.json" if RUN_TABLE_EXPRESSION else [],
        log:
            "logs/export_expression.log",
        params:
            cell_type_column = EXPR.get("cell_type_column",
                                        config["single_cell"]["cell_type_column"]),
            cell_types       = ",".join(EXPR.get("cell_types") or []),
            genes            = ",".join(EXPR.get("genes") or []),
            by_dataset       = str(EXPR.get("by_dataset", False)).upper(),
            dataset_column   = EXPR.get("dataset_column", "Dataset"),
            datasets         = ",".join(EXPR.get("datasets") or []),
            min_cells        = EXPR.get("min_cells", 10),
            assay            = EXPR.get("assay", "RNA"),
            value_layer      = EXPR.get("value_layer", "data"),
            counts_layer     = EXPR.get("counts_layer", "counts"),
            exclude_col      = EXPR.get("exclude_column", ""),
            exclude_vals     = ",".join(EXPR.get("exclude_values") or []),
            drop_zero        = str(EXPR.get("drop_zero_rows", True)).upper(),
            gene_id_type     = EXPR.get("gene_id_type", "symbol"),
            fmt              = EXPR.get("format", "parquet"),
            species          = EXPR.get("species") or config["species"],
            output_dir       = f"{OUTPUT}/export/expression",
            csv_dir          = f"{OUTPUT}/tables/expression" if RUN_TABLE_EXPRESSION else "",
            table_schema_version = "1.0",
        shell:
            """
            Rscript scripts/07_export_expression.R \
                --input            {input.data} \
                --output_dir       {params.output_dir} \
                --cell_type_column "{params.cell_type_column}" \
                --cell_types       "{params.cell_types}" \
                --genes            "{params.genes}" \
                --by_dataset       {params.by_dataset} \
                --dataset_column   "{params.dataset_column}" \
                --datasets         "{params.datasets}" \
                --min_cells        {params.min_cells} \
                --assay            "{params.assay}" \
                --value_layer      "{params.value_layer}" \
                --counts_layer     "{params.counts_layer}" \
                --exclude_col      "{params.exclude_col}" \
                --exclude_vals     "{params.exclude_vals}" \
                --drop_zero_rows   {params.drop_zero} \
                --gene_id_type     "{params.gene_id_type}" \
                --format           "{params.fmt}" \
                --csv_dir          "{params.csv_dir}" \
                --table_schema_version "{params.table_schema_version}" \
                --species          "{params.species}" \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 5c-2: Expression distribution (optional) ────────────────────────────
if RUN_DISTRIBUTION:
    rule export_distribution:
        input:
            data = QUERY_INPUT,
        output:
            manifest = f"{OUTPUT}/export/distribution/seqviewer_manifest.json",
            flag     = f"{OUTPUT}/export/distribution/distribution_complete.flag",
        log:
            "logs/export_distribution.log",
        params:
            output_dir       = f"{OUTPUT}/export/distribution",
            prefix           = DIST_CFG.get("output_prefix", "distribution"),
            alias            = DIST_CFG.get("alias", "Expression distribution"),
            genes            = ",".join(DIST_CFG.get("genes") or []),
            cell_type_column = DIST_CFG.get("cell_type_column",
                                            config["single_cell"]["cell_type_column"]),
            cell_types       = ",".join(DIST_CFG.get("cell_types") or []),
            sample_column    = DIST_CFG.get("sample_column",
                                            config["single_cell"]["sample_column"]),
            dataset_column   = DIST_CFG.get("dataset_column", "Dataset"),
            condition_column = DIST_CFG.get("condition_column",
                                            config["single_cell"]["condition_column"]),
            min_cells        = DIST_CFG.get("min_cells", 10),
            assay            = DIST_CFG.get("assay", "RNA"),
            value_layer      = DIST_CFG.get("value_layer", "data"),
            counts_layer     = DIST_CFG.get("counts_layer", "counts"),
            exclude_col      = DIST_CFG.get("exclude_column", ""),
            exclude_vals     = ",".join(DIST_CFG.get("exclude_values") or []),
            gene_id_type     = DIST_CFG.get("gene_id_type", "symbol"),
            species          = DIST_CFG.get("species") or config["species"],
            source           = DIST_CFG.get("source", "HypoMap"),
            license          = DIST_CFG.get("license", "CC BY 4.0"),
            citation         = DIST_CFG.get("citation", ""),
            notes            = DIST_CFG.get("notes", ""),
        shell:
            """
            Rscript scripts/13_export_distribution.R \
                --input            {input.data} \
                --output_dir       "{params.output_dir}" \
                --prefix           "{params.prefix}" \
                --alias            "{params.alias}" \
                --genes            "{params.genes}" \
                --cell_type_column "{params.cell_type_column}" \
                --cell_types       "{params.cell_types}" \
                --sample_column    "{params.sample_column}" \
                --dataset_column   "{params.dataset_column}" \
                --condition_column "{params.condition_column}" \
                --min_cells        {params.min_cells} \
                --assay            "{params.assay}" \
                --value_layer      "{params.value_layer}" \
                --counts_layer     "{params.counts_layer}" \
                --exclude_col      "{params.exclude_col}" \
                --exclude_vals     "{params.exclude_vals}" \
                --gene_id_type     "{params.gene_id_type}" \
                --species          "{params.species}" \
                --source           "{params.source}" \
                --license          "{params.license}" \
                --citation         "{params.citation}" \
                --notes            "{params.notes}" \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 5d: Static summary figures (optional) ──────────────────────────────
if RUN_EXPR_PLOTS:
    EXPR = config.get("expression_export", {})

    rule plot_expression:
        input:
            flag = f"{OUTPUT}/export/expression/expression_complete.flag",
        output:
            flag = f"{OUTPUT}/export/expression/plot_complete.flag",
        log:
            "logs/plot_expression.log",
        params:
            input_dir       = f"{OUTPUT}/export/expression",
            output_dir      = f"{OUTPUT}/export/expression/plots",
            top_genes       = EXPR.get("plots_top_genes", 50),
            marker_top_n    = EXPR.get("plots_marker_top_n", 40),
            genes           = ",".join(EXPR.get("plots_genes") or []),
            max_cell_types  = EXPR.get("plots_max_cell_types", 0),
            min_pct         = EXPR.get("plots_min_pct", 10),
            min_mean        = EXPR.get("plots_min_mean", 0.05),
            dpi             = EXPR.get("plots_dpi", 150),
            pdf             = str(EXPR.get("plots_pdf", True)).upper(),
        shell:
            """
            Rscript scripts/09_plot_expression.R \
                --input_dir       "{params.input_dir}" \
                --output_dir      "{params.output_dir}" \
                --top_genes       {params.top_genes} \
                --marker_top_n    {params.marker_top_n} \
                --genes           "{params.genes}" \
                --max_cell_types  {params.max_cell_types} \
                --min_pct         {params.min_pct} \
                --min_mean        {params.min_mean} \
                --dpi             {params.dpi} \
                --pdf             "{params.pdf}" \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 5e: Unified viewer export ──────────────────────────────────────────
# Preferred entry point when config has a `viewer_export` block. Builds the
# analysis and expression pieces in separate staging dirs, then merges them into
# one Import Folder. Legacy export_seqviewer / expression_export keep working
# when `viewer_export` is absent.
if RUN_VIEWER:
    VIEWER_SCHEMA = VIEWER_CFG.get(
        "schema_version",
        config.get("export_seqviewer", {}).get("schema_version", "1.0"))
    VIEWER_ORGANISM = "Mus musculus" if config.get("species") == "mouse" else "Homo sapiens"
    VIEWER_EXPR_PREFIX = VIEWER_CFG.get("expression_prefix", "expression")
    VIEWER_EXPR_ALIAS = VIEWER_CFG.get("expression_alias", VIEWER_EXPR_PREFIX)
    VIEWER_EXPR_CELL_TYPE_LEVEL = VIEWER_CFG.get("expression_cell_type_level", "")
    VIEWER_EXPR_PCT_SCALE = VIEWER_CFG.get("expression_pct_scale", "0-100")
    VIEWER_EXPR_GENE_INDEX = VIEWER_CFG.get(
        "expression_gene_index_name", f"{VIEWER_EXPR_PREFIX}_gene_index.parquet")
    VIEWER_EXPR_SUMMARY = VIEWER_CFG.get(
        "expression_summary_name", f"{VIEWER_EXPR_PREFIX}_cell_type_summary.parquet")
    VIEWER_EXPR_DATASET_TYPE = VIEWER_CFG.get("expression_dataset_type", "expression_summary")
    VIEWER_EXPR_PRIMARY = VIEWER_CFG.get("expression_primary", "wide")
    VIEWER_EXPR_INCLUDE_LONG = str(VIEWER_CFG.get("expression_include_long", False)).upper()
    VIEWER_EXPR_SOURCE = VIEWER_CFG.get("expression_source", "scRNA-seq pipeline")
    VIEWER_EXPR_LICENSE = VIEWER_CFG.get("expression_license", "CC BY 4.0")
    VIEWER_EXPR_CITATION = VIEWER_CFG.get("expression_citation", "")

    def viewer_export_inputs(wildcards):
        d = {}
        if RUN_VIEWER_ANALYSIS:
            d["analysis_flag"] = rules.export_analysis.output.flag
        if RUN_VIEWER_EXPRESSION:
            d["expr_flag"] = rules.export_expression.output.flag
        if RUN_DISTRIBUTION:
            d["distribution_flag"] = rules.export_distribution.output.flag
        return d

    rule viewer_export:
        input:
            unpack(viewer_export_inputs)
        output:
            manifest = f"{OUTPUT}/export/seqviewer/seqviewer_manifest.json",
            index    = f"{OUTPUT}/export/seqviewer/export_index.tsv",
            flag     = f"{OUTPUT}/export/seqviewer/export_complete.flag",
        log:
            "logs/viewer_export.log",
        params:
            expression            = str(RUN_VIEWER_EXPRESSION).upper(),
            work_dir              = f"{OUTPUT}/export/seqviewer_work",
            analysis_dir          = f"{OUTPUT}/export/analysis" if RUN_VIEWER_ANALYSIS else "",
            expression_dir        = f"{OUTPUT}/export/seqviewer_work/expression",
            distribution_dir      = f"{OUTPUT}/export/distribution" if RUN_DISTRIBUTION else "",
            output_dir            = f"{OUTPUT}/export/seqviewer",
            expression_export_dir = f"{OUTPUT}/export/expression",
            schema_version        = VIEWER_SCHEMA,
            organism              = VIEWER_ORGANISM,
            expr_prefix           = VIEWER_EXPR_PREFIX,
            expr_alias            = VIEWER_EXPR_ALIAS,
            expr_cell_type_level  = VIEWER_EXPR_CELL_TYPE_LEVEL,
            expr_pct_scale        = VIEWER_EXPR_PCT_SCALE,
            expr_gene_index       = VIEWER_EXPR_GENE_INDEX,
            expr_summary          = VIEWER_EXPR_SUMMARY,
            expr_dataset_type     = VIEWER_EXPR_DATASET_TYPE,
            expr_primary          = VIEWER_EXPR_PRIMARY,
            expr_include_long     = VIEWER_EXPR_INCLUDE_LONG,
            expr_source           = VIEWER_EXPR_SOURCE,
            expr_license          = VIEWER_EXPR_LICENSE,
            expr_citation         = VIEWER_EXPR_CITATION,
            expr_gene_id_type     = EXPR_CFG.get("gene_id_type", "symbol"),
            expr_normalization    = EXPR_CFG.get("normalization", "log1p_normalized"),
            expr_detection        = EXPR_CFG.get("detection", "count > 0"),
        shell:
            r"""
            set -euo pipefail
            rm -rf "{params.expression_dir}"
            mkdir -p "{params.expression_dir}"

            if [ "{params.expression}" = "TRUE" ]; then
              Rscript scripts/08_build_viewer_bundle.R \
                --export_dir      "{params.expression_export_dir}" \
                --bundle_dir      "{params.expression_dir}" \
                --prefix          "{params.expr_prefix}" \
                --alias           "{params.expr_alias}" \
                --dataset_type    "{params.expr_dataset_type}" \
                --primary         "{params.expr_primary}" \
                --cell_type_level "{params.expr_cell_type_level}" \
                --organism        "{params.organism}" \
                --source          "{params.expr_source}" \
                --license         "{params.expr_license}" \
                --citation        "{params.expr_citation}" \
                --pct_scale       "{params.expr_pct_scale}" \
                --gene_id_type    "{params.expr_gene_id_type}" \
                --normalization   "{params.expr_normalization}" \
                --detection       "{params.expr_detection}" \
                --gene_index_name "{params.expr_gene_index}" \
                --summary_name    "{params.expr_summary}" \
                --include_long    "{params.expr_include_long}" \
                --fallback_multi_group FALSE
            fi

            Rscript scripts/10_merge_viewer_bundle.R \
              --analysis_dir     "{params.analysis_dir}" \
              --expression_dir   "{params.expression_dir}" \
              --distribution_dir "{params.distribution_dir}" \
              --output_dir       "{params.output_dir}" \
              --schema_version   "{params.schema_version}" \
              2>&1 | tee {log}

            rm -rf "{params.expression_dir}"
            rmdir "{params.work_dir}" 2>/dev/null || true
            """


# ── Phase 5f: Canonical table merge + Excel ──────────────────────────────────
if RUN_TABLE_MERGE:
    rule merge_tables:
        input:
            analysis   = f"{OUTPUT}/tables/analysis/manifest.json" if RUN_TABLE_ANALYSIS else [],
            expression = f"{OUTPUT}/tables/expression/manifest.json" if RUN_TABLE_EXPRESSION else [],
        output:
            manifest = f"{OUTPUT}/tables/manifest.json",
            index    = f"{OUTPUT}/tables/index.tsv",
            readme   = f"{OUTPUT}/tables/README.md",
            flag     = f"{OUTPUT}/tables/tables_complete.flag",
        log:
            "logs/merge_tables.log",
        params:
            tables_dir          = f"{OUTPUT}/tables",
            analysis_manifest   = f"{OUTPUT}/tables/analysis/manifest.json" if RUN_TABLE_ANALYSIS else "",
            expression_manifest = f"{OUTPUT}/tables/expression/manifest.json" if RUN_TABLE_EXPRESSION else "",
            table_schema_version = "1.0",
        shell:
            """
            Rscript scripts/11_merge_tables.R \
                --tables_dir          "{params.tables_dir}" \
                --analysis_manifest   "{params.analysis_manifest}" \
                --expression_manifest "{params.expression_manifest}" \
                --table_schema_version "{params.table_schema_version}" \
                2>&1 | tee {log}
            touch {output.flag}
            """


if RUN_TABLE_EXCEL:
    rule export_excel:
        input:
            manifest = f"{OUTPUT}/tables/manifest.json",
        output:
            manifest = f"{OUTPUT}/export/excel/excel_manifest.json",
            readme   = f"{OUTPUT}/export/excel/README.md",
            flag     = f"{OUTPUT}/export/excel/excel_complete.flag",
        log:
            "logs/export_excel.log",
        params:
            tables_dir = f"{OUTPUT}/tables",
            output_dir = f"{OUTPUT}/export/excel",
        shell:
            """
            Rscript scripts/12_export_excel.R \
                --tables_dir "{params.tables_dir}" \
                --output_dir "{params.output_dir}" \
                2>&1 | tee {log}
            touch {output.flag}
            """


# ── Phase 6: Report ───────────────────────────────────────────────────────────
def report_inputs(wildcards):
    d = {}
    if RUN_DE:    d["de_flag"]   = rules.pseudobulk_de.output.flag
    if RUN_COMP:  d["comp_flag"] = rules.composition.output.flag
    if RUN_ENR:   d["enr_flag"]  = rules.enrichment.output.flag
    if RUN_TRAJ:  d["traj_flag"] = rules.trajectory.output.flag
    if RUN_COMM:  d["comm_flag"] = rules.cell_communication.output.flag
    # unpack() requires a non-empty dict; depend on the template so the
    # report also re-renders when report.Rmd changes.
    if not d:
        d["report_rmd"] = "scripts/report.Rmd"
    return d

rule report:
    input:
        unpack(report_inputs)
    output:
        html = f"{OUTPUT}/report/report.html",
    log:
        f"logs/report.log",
    params:
        output_dir       = f"{OUTPUT}",
        report_dir       = f"{OUTPUT}/report",
        report_title     = config["report"]["title"],
        author           = config["report"].get("author", ""),
        methods_draft    = str(config["report"].get("generate_methods_draft", True)).upper(),
        padj             = config["padj_cutoff"],
        log2fc           = config["log2fc_cutoff"],
        min_cells        = config["single_cell"]["min_cells_per_group"],
        species          = config["species"],
        run_traj         = str(RUN_TRAJ).upper(),
        run_comm         = str(RUN_COMM).upper(),
        run_de           = str(RUN_DE).upper(),
        run_comp         = str(RUN_COMP).upper(),
        run_enr          = str(RUN_ENR).upper(),
        run_explore      = str(RUN_EXPLORE).upper(),
    shell:
        """
        Rscript -e "rmarkdown::render(
            'scripts/report.Rmd',
            output_dir   = '{params.report_dir}',
            output_file  = 'report.html',
            params = list(
                output_dir    = '{params.output_dir}',
                title         = '{params.report_title}',
                author        = '{params.author}',
                padj_cutoff   = {params.padj},
                log2fc_cutoff = {params.log2fc},
                min_cells     = {params.min_cells},
                species       = '{params.species}',
                run_traj      = {params.run_traj},
                run_comm      = {params.run_comm},
                run_de        = {params.run_de},
                run_comp      = {params.run_comp},
                run_enr       = {params.run_enr},
                run_explore   = {params.run_explore},
                methods_draft = {params.methods_draft}
            )
        )" 2>&1 | tee {log}
        """
