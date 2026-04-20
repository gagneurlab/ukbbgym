#!/usr/bin/env python3

import os
import yaml
import logging
import dxpy
import polars as pl
import pandas as pd
import numpy as np
import itertools
from scipy import stats
from plotnine import *

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

WORK_DIR = "PATH_TO_FILE"
IN_DIR = f"{WORK_DIR}/in"
OUT_DIR = f"{WORK_DIR}/out"

def resolve_input(dxlink, local_name):
    local_path = os.path.join(IN_DIR, local_name)
    logger.info(f"Downloading {local_name}...")
    dxpy.download_dxfile(dxlink, local_path)
    return local_path

def load_annotation_config(override_path=None, bundled_name="config_correlations.yaml"):
    path = override_path or os.path.join(os.path.dirname(__file__), bundled_name)
    logger.info(f"Loading annotation config from {path}")
    with open(path) as f:
        return yaml.safe_load(f)

def load_variant_class_config(variant_class, override_yaml_path=None):
    if variant_class == "other":
        if override_yaml_path is None:
            raise ValueError("variant_class='other' requires variant_class_config_yaml input")
        config_path = override_yaml_path
    else:
        config_path = os.path.join(os.path.dirname(__file__), "config_variant_classes.yaml")
    logger.info(f"Loading variant class config from {config_path}")
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    if variant_class != "other" and variant_class not in cfg:
        raise ValueError(f"Unknown variant_class '{variant_class}'. Options: {list(cfg.keys())}")
    return cfg[variant_class]

def load_and_merge_custom_scores(ann_lazy, custom_path, model_labels_str, model_directions_str, fill_missing):
    custom = pl.read_parquet(custom_path)
    has_region = "region" in custom.columns
    join_on = ["id", "region"] if has_region else ["id"]

    meta_cols = {"id", "region"}
    model_cols = [c for c in custom.columns if c not in meta_cols]

    labels = [s.strip() for s in model_labels_str.split(",")] if model_labels_str else model_cols
    directions = [int(s.strip()) for s in model_directions_str.split(",")] if model_directions_str else [1] * len(model_cols)

    assert len(labels) == len(model_cols), "model_labels count must match score columns"
    assert len(directions) == len(model_cols), "model_directions count must match score columns"

    logger.info(f"Merging custom scores: {len(model_cols)} models")
    ann_lazy = ann_lazy.join(custom.lazy(), on=join_on, how="left")

    if fill_missing:
        ann_lazy = ann_lazy.with_columns([pl.col(c).fill_null(0.0).cast(pl.Float32) for c in model_cols])

    return ann_lazy, model_cols, labels, directions

def build_annotation_config_df(config, custom_models, custom_labels, custom_directions):
    records = [
        {
            "category": category,
            "annotation": anno,
            "color": props["color"],
            "label": props["label"],
            "annotation_dir": props.get("direction", 1),
        }
        for category, annos in config["rare_variant_annotations"].items()
        for anno, props in annos.items()
    ]

    # Add custom models
    for model, label, direction in zip(custom_models, custom_labels, custom_directions):
        records.append({
            "category": "custom_model",
            "annotation": model,
            "color": "#808080",  # Gray for custom models
            "label": label,
            "annotation_dir": direction,
        })

    return pl.DataFrame(records).with_columns(pl.col("annotation_dir").cast(pl.Int8))

@dxpy.entry_point("main")
def main(
    associations_parquet,
    annotations_parquet,
    appv_parquet,
    custom_scores_parquet=None,
    model_labels=None,
    model_directions=None,
    fill_missing_scores=False,
    annotation_categories="plof,missense,conservation,splicing,regulatory",
    variant_class="all_variants",
    variant_class_config_yaml=None,
    annotation_config_yaml=None,
    mac_threshold=20,
    only_snps=True,
    only_clinvar=False,
    exclude_clinvar=False,
):
    os.makedirs(IN_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    # Download inputs
    assoc_path = resolve_input(associations_parquet, "associations.parquet")
    anno_path = resolve_input(annotations_parquet, "annotations.parquet")
    appv_path = resolve_input(appv_parquet, "appv.parquet")

    custom_path = None
    if custom_scores_parquet is not None:
        custom_path = resolve_input(custom_scores_parquet, "custom_scores.parquet")

    variant_class_yaml_path = None
    if variant_class_config_yaml is not None:
        variant_class_yaml_path = resolve_input(variant_class_config_yaml, "variant_class.yaml")

    annotation_yaml_path = None
    if annotation_config_yaml is not None:
        annotation_yaml_path = resolve_input(annotation_config_yaml, "annotation_config.yaml")

    # Load configurations
    anno_config = load_annotation_config(annotation_yaml_path, "config_correlations.yaml")
    vc = load_variant_class_config(variant_class, variant_class_yaml_path)

    # Load associations with LOFTEE correlations
    gene_trait_df = (
        pl.read_parquet(assoc_path)
        .select(["region", "phenotype", "correlation"])
        .with_columns(
            loftee_corr=pl.col("correlation"),
            loftee_corr_abs=pl.col("correlation").abs(),
            loftee_corr_dir=pl.col("correlation")/pl.col("correlation").abs(),
        )
        .select(["region", "phenotype", "loftee_corr", "loftee_corr_abs", "loftee_corr_dir"])
    )
    logger.info(f"Loaded {gene_trait_df.shape[0]} gene-trait associations")

    # Load annotations and apply filters
    ann = pl.scan_parquet(anno_path)

    # Build dynamic filters from variant_class.yaml
    dynamic_filters = []
    vc_filters = vc.get("variant_filtering", [])
    if vc_filters:
        for f in vc_filters:
            dynamic_filters.append(eval(f))

    if exclude_clinvar:
        dynamic_filters.append(pl.col("clinical_significance").is_null())
    elif only_clinvar:
        dynamic_filters.append(pl.col("clinical_significance").is_not_null())

    if only_snps:
        dynamic_filters.append((pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1))

    ann = (
        ann
        .filter(pl.col("region").is_in(gene_trait_df["region"].unique()))
        .filter(*dynamic_filters) if dynamic_filters else ann
    )

    # Merge custom scores if provided
    custom_models = []
    custom_labels = []
    custom_directions = []
    if custom_path:
        ann, custom_models, custom_labels, custom_directions = load_and_merge_custom_scores(
            ann, custom_path, model_labels, model_directions, fill_missing_scores
        )

    # Build annotation config dataframe (built-in + custom)
    anno_config_df = build_annotation_config_df(anno_config, custom_models, custom_labels, custom_directions)

    # Select annotation categories
    selected_categories = [c.strip() for c in annotation_categories.split(",")]
    logger.info(f"Selected annotation categories: {selected_categories}")

    # Get all possible annotations for selected categories
    all_annotation_list = anno_config_df.filter(
        pl.col("category").is_in(selected_categories)
    )["annotation"].to_list()

    # Filter to annotations that actually exist in the annotation parquet
    anno_schema = ann.collect_schema()
    existing_annos = [a for a in all_annotation_list if a in anno_schema.names()]

    logger.info(f"Using {len(existing_annos)} annotations from selected categories")

    # Select and collect annotation data
    anno = (
        ann
        .select(["id", "region"] + existing_annos)
        .collect(engine="streaming")
    )
    logger.info(f"Loaded annotations: {anno.shape}")

    # Melt annotations to long format
    melted_anno = (
        anno.lazy()
        .unpivot(
            index=["id", "region"],
            on=existing_annos,
            variable_name="annotation",
            value_name="annotation_score"
        )
        .with_columns(pl.col("annotation_score").cast(pl.Float32))
        .collect(engine="streaming")
    )
    logger.info(f"Melted annotations: {melted_anno.shape}")

    # Load APPV
    appv = (
        pl.scan_parquet(appv_path)
        .filter(pl.col("n_individuals") <= mac_threshold)
        .select(["id", "phenotype", "mean_pheno_value", "n_individuals"])
        .collect(engine="streaming")
    )
    logger.info(f"Loaded APPV: {appv.shape}")

    # Build id->region mapping
    id_region = anno.select(["id", "region"]).unique().lazy()

    # Main join and correlation computation
    final_lazy_plan = (
        appv.lazy()
        .join(id_region, on="id", how="inner")
        .join(
            gene_trait_df[["region", "phenotype"]].lazy(),
            on=["region", "phenotype"],
            how="inner"
        )
        .join(
            melted_anno.lazy(),
            on=["id", "region"],
            how="inner"
        )
        # Spearman correlation via ranking
        .with_columns(
            annotation_score_rank=pl.col("annotation_score")
            .rank("average")
            .over(["region", "phenotype", "annotation"]),
            mean_pheno_value_rank=pl.col("mean_pheno_value")
            .rank("average")
            .over(["region", "phenotype", "annotation"]),
        )
        .group_by(["region", "phenotype", "annotation"])
        .agg(
            n_variants=pl.col("id").count(),
            correlation=pl.when(
                (pl.col("annotation_score_rank").n_unique() > 1) &
                (pl.col("mean_pheno_value_rank").n_unique() > 1)
            )
            .then(pl.corr("annotation_score_rank", "mean_pheno_value_rank", propagate_nans=True))
            .otherwise(None)
        )
        .drop_nulls()
    )

    correlation_df = final_lazy_plan.collect(engine="streaming")
    logger.info(f"Computed correlations: {correlation_df.shape}")

    # Join with metadata and direction-correct
    correlation_df = (
        correlation_df
        .join(anno_config_df, on="annotation", how="left")
        .join(gene_trait_df, on=["region", "phenotype"], how="left")
        .with_columns(
            corr_beta=pl.col("correlation") * pl.col("loftee_corr_dir") * pl.col("annotation_dir")
        )
        .with_columns(
            corr_beta_rescaled=(pl.col("corr_beta") / pl.col("loftee_corr_abs")) * pl.col("loftee_corr_dir")
        )
        .drop_nans()
        .drop_nulls()
    )

    # Filter: n_variants > 100
    filt_corr_df = correlation_df.filter(pl.col("n_variants") > 100)
    logger.info(f"After filtering (n_variants > 100): {filt_corr_df.shape}")

    # Output results parquet
    out_results_path = os.path.join(OUT_DIR, "correlations_olink_results.parquet")
    filt_corr_df.write_parquet(out_results_path)
    logger.info(f"Results written to {out_results_path}")

    # Plotting
    # Main boxplot
    plotting_col = "corr_beta"
    dashed_line_value = 0

    plot_corr_pl = (
        filt_corr_df
        .drop_nans()
        .with_columns(
            median_corr_beta=pl.col(plotting_col).median().over("annotation")
        )
    )

    # Order annotations by median
    ordered_labels = (
        plot_corr_pl
        .sort("median_corr_beta", descending=False)
        .select("label")
        .unique(maintain_order=True)
        .to_series()
        .to_list()
    )

    plot_corr_pl = plot_corr_pl.with_columns(
        pl.col("label").cast(pl.Enum(ordered_labels))
    )

    color_dict = dict(plot_corr_pl.select("annotation", "color").unique().iter_rows())

    # Convert to pandas for plotnine
    plot_df = plot_corr_pl.to_pandas()

    boxplot = (
        ggplot(plot_df, aes(x="label", y=plotting_col, fill="annotation"))
        + geom_hline(aes(yintercept=dashed_line_value), color="black", linetype="dotted")
        + geom_boxplot(alpha=1, outlier_shape=None)
        + theme_minimal()
        + scale_fill_manual(values=color_dict)
        + labs(
            x="",
            y="Spearman correlation",
            title=f"Variant Annotation Correlations ({plot_corr_pl.select('region').n_unique()} genes)",
        )
        + coord_flip()
        + theme(
            figure_size=(8, max(5, plot_corr_pl.select("annotation").n_unique() / 3 + 0.5)),
            legend_position="none",
            axis_text=element_text(size=13),
            axis_title=element_text(size=13),
            plot_background=element_rect(fill="white", color="white"),
        )
    )

    boxplot_path = os.path.join(OUT_DIR, "boxplot.svg")
    boxplot.save(boxplot_path, verbose=False)
    logger.info(f"Boxplot saved to {boxplot_path}")

    # Pairwise Wilcoxon heatmap
    plot_annotations = plot_corr_pl.select("annotation").unique().to_series().to_list()
    anno_to_label = dict(
        plot_corr_pl.select(["annotation", "label"]).unique().iter_rows()
    )

    heatmap_data = []
    for ann_x, ann_y in itertools.product(plot_annotations, repeat=2):
        label_x = anno_to_label.get(ann_x, ann_x)
        label_y = anno_to_label.get(ann_y, ann_y)

        if ann_x == ann_y:
            heatmap_data.append({
                "Tool_X": label_x,
                "Tool_Y": label_y,
                "mean_diff": 0.0,
                "sig": ""
            })
            continue

        # Inner join on (region, phenotype) for paired values
        df_x = plot_corr_pl.filter(pl.col("annotation") == ann_x).select(["region", "phenotype", "corr_beta"])
        df_y = plot_corr_pl.filter(pl.col("annotation") == ann_y).select(["region", "phenotype", "corr_beta"])

        paired = df_x.join(df_y, on=["region", "phenotype"], suffix="_y")
        c_x = paired["corr_beta"].to_numpy()
        c_y = paired["corr_beta_y"].to_numpy()

        mean_diff = float(c_y.mean() - c_x.mean())

        if len(c_x) >= 10 and not np.allclose(c_x, c_y):
            _, p_val = stats.wilcoxon(c_x, c_y, alternative="two-sided")
        else:
            p_val = 1.0

        if p_val < 0.001:
            sig = "***"
        elif p_val < 0.01:
            sig = "**"
        elif p_val < 0.05:
            sig = "*"
        else:
            sig = ""

        heatmap_data.append({
            "Tool_X": label_x,
            "Tool_Y": label_y,
            "mean_diff": mean_diff,
            "sig": sig
        })

    df_heat = pd.DataFrame(heatmap_data)
    df_heat["Tool_X"] = pd.Categorical(df_heat["Tool_X"], categories=ordered_labels, ordered=True)
    df_heat["Tool_Y"] = pd.Categorical(df_heat["Tool_Y"], categories=ordered_labels[::-1], ordered=True)

    n_tools = len(plot_annotations)

    heatmap = (
        ggplot(df_heat, aes(x="Tool_X", y="Tool_Y", fill="mean_diff"))
        + geom_tile(color="#333333", size=0.5)
        + geom_text(aes(label="sig"), color="black", size=12, va="center", nudge_y=-0.1)
        + scale_fill_gradient2(low="#2C7BB6", mid="#FFFFFF", high="#D7191C", midpoint=0)
        + labs(
            title=f"Pairwise Wilcoxon Signed-Rank Test ({plot_corr_pl.select('region').n_unique()} genes)",
            subtitle="*** p<0.001, ** p<0.01, * p<0.05",
            x="Tool X",
            y="Tool Y",
            fill="Tool Y − X\n(avg correlation)"
        )
        + theme_minimal()
        + theme(
            figure_size=(n_tools * 0.6 + 2.5, n_tools * 0.6 + 1.5),
            aspect_ratio=1,
            axis_title=element_text(size=13),
            axis_text=element_text(size=13),
            axis_text_x=element_text(rotation=45, hjust=1),
            panel_grid=element_blank(),
            legend_background=element_rect(fill="white", color="white", alpha=0.8),
            plot_background=element_rect(fill="white", color="white"),
        )
    )

    heatmap_path = os.path.join(OUT_DIR, "heatmap.svg")
    heatmap.save(heatmap_path, verbose=False)
    logger.info(f"Heatmap saved to {heatmap_path}")

    # Upload outputs
    uploaded_results = dxpy.upload_local_file(out_results_path)
    uploaded_boxplot = dxpy.upload_local_file(boxplot_path)
    uploaded_heatmap = dxpy.upload_local_file(heatmap_path)

    logger.info("Upload complete")

    return {
        "results_parquet": dxpy.dxlink(uploaded_results),
        "boxplot_svg": dxpy.dxlink(uploaded_boxplot),
        "heatmap_svg": dxpy.dxlink(uploaded_heatmap),
    }

dxpy.run()
