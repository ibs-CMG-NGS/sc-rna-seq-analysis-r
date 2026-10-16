"""
scRNA-seq 3rd Analysis Pipeline
Usage: snakemake --configfile configs/config_<experiment>.yml --cores <N>
"""

import os
from pathlib import Path

configfile: "configs/template/config.yml"


# ── Derived variables ──────────────────────────────────────────────────────────
OUTPUT     = config["output_dir"]
OPTIONAL   = config.get("optional_modules", {})
RUN_TRAJ   = OPTIONAL.get("trajectory", False)
RUN_COMM   = OPTIONAL.get("cell_communication", True)


# ── Helper: save config snapshot ──────────────────────────────────────────────
import shutil, yaml, datetime

def save_config_snapshot(config, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    out = os.path.join(output_dir, "config_used.yml")
    with open(out, "w") as f:
        yaml.dump(dict(config), f, default_flow_style=False)


# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        expand("{out}/pseudobulk_de/de_complete.flag",   out=OUTPUT),
        expand("{out}/composition/composition_complete.flag", out=OUTPUT),
        expand("{out}/enrichment/enrichment_complete.flag",   out=OUTPUT),
        expand("{out}/trajectory/trajectory_complete.flag",   out=OUTPUT) if RUN_TRAJ else [],
        expand("{out}/cell_communication/comm_complete.flag", out=OUTPUT) if RUN_COMM else [],
        expand("{out}/report/report.html",                    out=OUTPUT),


# ── Phase 1: Load & validate input ───────────────────────────────────────────
rule load_input:
    input:
        data = config["single_cell"]["input_path"],
    output:
        rds  = f"{OUTPUT}/loaded/seurat_validated.rds",
        flag = f"{OUTPUT}/loaded/load_complete.flag",
    log:
        f"logs/load_input.log",
    params:
        format         = config["single_cell"]["input_format"],
        cell_type_col  = config["single_cell"]["cell_type_column"],
        sample_col     = config["single_cell"]["sample_column"],
        condition_col  = config["single_cell"]["condition_column"],
        cell_types     = config["single_cell"].get("cell_types") or "null",
        min_cells      = config["single_cell"]["min_cells_per_group"],
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
            --output_rds    {output.rds} \
            2>&1 | tee {log}
        touch {output.flag}
        """


# ── Phase 2a: Pseudobulk DE ───────────────────────────────────────────────────
rule pseudobulk_de:
    input:
        rds  = rules.load_input.output.rds,
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


# ── Phase 2b: Compositional analysis ─────────────────────────────────────────
rule composition:
    input:
        rds  = rules.load_input.output.rds,
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


# ── Phase 3: Enrichment ───────────────────────────────────────────────────────
rule enrichment:
    input:
        flag = rules.pseudobulk_de.output.flag,
        de_dir = f"{OUTPUT}/pseudobulk_de",
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
        output_dir  = f"{OUTPUT}/enrichment",
    shell:
        """
        Rscript scripts/02_enrichment.R \
            --de_dir       {input.de_dir} \
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


# ── Phase 4: Trajectory (optional) ───────────────────────────────────────────
if RUN_TRAJ:
    rule trajectory:
        input:
            rds  = rules.load_input.output.rds,
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
            rds  = rules.load_input.output.rds,
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


# ── Phase 6: Report ───────────────────────────────────────────────────────────
rule report:
    input:
        de_flag   = rules.pseudobulk_de.output.flag,
        comp_flag = rules.composition.output.flag,
        enr_flag  = rules.enrichment.output.flag,
        traj_flag = rules.trajectory.output.flag,
        comm_flag = rules.cell_communication.output.flag,
    output:
        html = f"{OUTPUT}/report/report.html",
    log:
        f"logs/report.log",
    params:
        output_dir       = f"{OUTPUT}",
        report_title     = config["report"]["title"],
        author           = config["report"].get("author", ""),
        methods_draft    = str(config["report"].get("generate_methods_draft", True)).upper(),
        padj             = config["padj_cutoff"],
        log2fc           = config["log2fc_cutoff"],
        species          = config["species"],
        run_traj         = str(RUN_TRAJ).upper(),
        run_comm         = str(RUN_COMM).upper(),
    shell:
        """
        Rscript -e "rmarkdown::render(
            'scripts/report.Rmd',
            output_file  = '../{output.html}',
            params = list(
                output_dir    = '{params.output_dir}',
                title         = '{params.report_title}',
                author        = '{params.author}',
                padj_cutoff   = {params.padj},
                log2fc_cutoff = {params.log2fc},
                species       = '{params.species}',
                run_traj      = {params.run_traj},
                run_comm      = {params.run_comm},
                methods_draft = {params.methods_draft}
            )
        )" 2>&1 | tee {log}
        """
