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

WORK_DIR = "/home/dnanexus"
IN_DIR = f"{WORK_DIR}/in"
OUT_DIR = f"{WORK_DIR}/out"

def resolve_input(dxlink, local_name):
    local_path = os.path.join(IN_DIR, local_name)
    logger.info(f"Downloading {local_name}...")
    dxpy.download_dxfile(dxlink, local_path)
    return local_path

def load_annotation_config(override_path=None, bundled_name="config_odds.yaml"):
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
    for model, label, direction in zip(custom_models, custom_labels, custom_directions):
        records.append({
            "category": "custom_model",
            "annotation": model,
            "color": "#808080",
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
    target_k="5,10",
    max_gene_rank=50,
):
    os.makedirs(IN_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    # Parse target_k
    target_k = np.array([int(x.strip()) for x in target_k.split(",")], dtype=np.int64)
    logger.info(f"Target k values: {target_k}")

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
    anno_config = load_annotation_config(annotation_yaml_path, "config_odds.yaml")
    vc = load_variant_class_config(variant_class, variant_class_yaml_path)

    # Load associations
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

    # Load annotations
    ann = pl.scan_parquet(anno_path)
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
    ann = ann.filter(pl.col("region").is_in(gene_trait_df["region"].unique()))
    if dynamic_filters:
        ann = ann.filter(*dynamic_filters)

    # Merge custom scores if provided
    custom_models = []
    custom_labels = []
    custom_directions = []
    if custom_path:
        ann, custom_models, custom_labels, custom_directions = load_and_merge_custom_scores(
            ann, custom_path, model_labels, model_directions, fill_missing_scores
        )

    # Build annotation config
    anno_config_df = build_annotation_config_df(anno_config, custom_models, custom_labels, custom_directions)
    selected_categories = [c.strip() for c in annotation_categories.split(",")]
    all_annotation_list = anno_config_df.filter(
        pl.col("category").is_in(selected_categories)
    )["annotation"].to_list()
    anno_schema = ann.collect_schema()
    existing_annos = [a for a in all_annotation_list if a in anno_schema.names()]
    logger.info(f"Using {len(existing_annos)} annotations")

    anno = (
        ann
        .select(["id", "region"] + existing_annos)
        .collect(engine="streaming")
    )
    logger.info(f"Loaded annotations: {anno.shape}")

    # Melt annotations
    melted_anno = (
        anno.lazy()
        .unpivot(
            index=["id", "region"],
            on=existing_annos,
            variable_name="annotation",
            value_name="annotation_score"
        )
        .with_columns(pl.col("annotation_score").cast(pl.Float32))
        .join(
            anno_config_df.select(["annotation", "category", "annotation_dir"]).lazy(),
            on="annotation",
            how="left"
        )
        .filter(pl.col("category").is_in(selected_categories))
        .with_columns(
            annotation_score_dircor=pl.col("annotation_score") * pl.col("annotation_dir").cast(pl.Float32)
        )
        .drop_nulls("annotation_score")
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

    # Direction-correct phenotype values
    variant_zscores = (
        appv.lazy()
        .join(id_region, on="id", how="inner")
        .join(gene_trait_df.lazy(), on=["region", "phenotype"], how="inner")
        .with_columns(
            mean_pheno_value=pl.col("mean_pheno_value") * pl.col("loftee_corr_dir").cast(pl.Float32)
        )
        .select(["id", "region", "mean_pheno_value"])
        .collect(engine="streaming")
    )
    logger.info(f"Variant z-scores: {variant_zscores.shape[0]} variant-gene pairs")

    # Rank variants and compute per-gene cumulative stats
    ranked_zscores_genewise = (
        variant_zscores.lazy()
        .join(
            melted_anno.lazy().select(["id", "region", "annotation", "annotation_score_dircor"]),
            on=["id", "region"],
            how="inner"
        )
        .with_columns(
            gene_rank=pl.col("annotation_score_dircor")
            .rank(method="max", descending=True)
            .over(["annotation", "region"])
            .cast(pl.Int32)
        )
        .with_columns(
            mean_pheno_value_by_rank=pl.col("mean_pheno_value").mean().over(["annotation", "region", "gene_rank"])
        )
        .unique(subset=["annotation", "region", "gene_rank"], keep="any")
        .sort(["annotation", "region", "gene_rank"])
        .with_columns(
            ordinal_rank=(pl.int_range(pl.len()).over(["annotation", "region"]) + 1).cast(pl.Int32)
        )
        .collect(engine="streaming")
    )
    logger.info(f"Ranked variants: {ranked_zscores_genewise.shape[0]}")

    # Cumulative stats
    cum_genewise_stats = (
        ranked_zscores_genewise.lazy()
        .filter(pl.col("ordinal_rank") <= max_gene_rank)
        .sort(["annotation", "region", "ordinal_rank"])
        .with_columns(
            cum_mean_pheno=(
                pl.col("mean_pheno_value_by_rank").cum_sum().over(["annotation", "region"])
                / pl.col("ordinal_rank")
            )
        )
        .group_by(["annotation", "ordinal_rank"])
        .agg(
            mean_cum_zscore=pl.col("cum_mean_pheno").mean(),
            std_cum_zscore=pl.col("cum_mean_pheno").std(),
            n_genes=pl.col("region").n_unique(),
        )
        .with_columns(
            se_cum_zscore=pl.col("std_cum_zscore") / pl.col("n_genes").cast(pl.Float64).sqrt(),
        )
        .with_columns(
            ci_lower=pl.col("mean_cum_zscore") - pl.col("se_cum_zscore"),
            ci_upper=pl.col("mean_cum_zscore") + pl.col("se_cum_zscore"),
        )
        .join(anno_config_df.select(["annotation", "label", "color"]).lazy(), on="annotation", how="left")
        .sort(["annotation", "ordinal_rank"])
        .collect()
    )
    logger.info(f"Cumulative stats: {cum_genewise_stats.shape}")

    # Output cumulative stats
    cum_out_path = os.path.join(OUT_DIR, "cum_stats.parquet")
    cum_genewise_stats.write_parquet(cum_out_path)

    # Per-gene top-N means
    gene_topn_list = []
    for n_val in target_k:
        gdf = (
            ranked_zscores_genewise.lazy()
            .filter(pl.col("ordinal_rank") <= n_val)
            .group_by(["region", "annotation"])
            .agg(
                gene_topn_zscore=pl.col("mean_pheno_value_by_rank").mean(),
                n_variants_used=pl.len(),
            )
            .with_columns(top_n=pl.lit(n_val, dtype=pl.Int64))
            .collect(engine="streaming")
        )
        gene_topn_list.append(gdf)
    gene_topn = pl.concat(gene_topn_list)

    gene_topn_out_path = os.path.join(OUT_DIR, "gene_topn.parquet")
    gene_topn.write_parquet(gene_topn_out_path)
    logger.info(f"Per-gene top-N: {gene_topn.shape}")

    # Plotting
    color_dict = dict(cum_genewise_stats.select(["label", "color"]).unique().iter_rows())

    cum_plot_df = cum_genewise_stats.to_pandas()

    lineplot = (
        ggplot(cum_plot_df, aes(x="ordinal_rank", y="mean_cum_zscore", color="label", fill="label"))
        + geom_ribbon(aes(ymin="ci_lower", ymax="ci_upper"), alpha=0.1, color="none")
        + geom_line()
        + scale_color_manual(values=color_dict)
        + scale_fill_manual(values=color_dict)
        + labs(
            title=f"Cumulative Mean Z-Score by Rank",
            subtitle=f"Ribbon = 1× SEM across {ranked_zscores_genewise.select('region').n_unique()} genes",
            x="Rank within gene (N)",
            y="Cumulative mean z-score",
            color="Annotation",
            fill="Annotation",
        )
        + theme_minimal()
        + theme(
            figure_size=(10, 6),
            axis_text=element_text(size=11),
            axis_title=element_text(size=12),
            legend_text=element_text(size=11),
            legend_background=element_rect(fill="white", color="white", alpha=0.8),
            plot_background=element_rect(fill="white", color="white"),
        )
    )

    lineplot_path = os.path.join(OUT_DIR, "lineplot.svg")
    lineplot.save(lineplot_path, verbose=False)
    logger.info(f"Lineplot saved to {lineplot_path}")

    # Pairwise Wilcoxon heatmap
    plot_annotations = gene_topn.select("annotation").unique().to_series().to_list()
    anno_to_label = dict(anno_config_df.select(["annotation", "label"]).unique().iter_rows())
    tool_means = {ann: float(gene_topn.filter(pl.col("annotation") == ann)["gene_topn_zscore"].mean())
                  for ann in plot_annotations}
    ordered_tools = sorted(tool_means.keys(), key=lambda x: tool_means[x], reverse=True)
    ordered_labels = [anno_to_label.get(t, t) for t in ordered_tools]

    heatmap_data = []
    for n_val in target_k:
        for ann_x, ann_y in itertools.product(ordered_tools, repeat=2):
            label_x = anno_to_label.get(ann_x, ann_x)
            label_y = anno_to_label.get(ann_y, ann_y)

            if ann_x == ann_y:
                heatmap_data.append({
                    "Window": f"Top {n_val:,}",
                    "Window_Int": n_val,
                    "Tool_X": label_x,
                    "Tool_Y": label_y,
                    "mean_diff": 0.0,
                    "sig": ""
                })
                continue

            df_x = gene_topn.filter((pl.col("annotation") == ann_x) & (pl.col("top_n") == n_val)).select(["region", "gene_topn_zscore"])
            df_y = gene_topn.filter((pl.col("annotation") == ann_y) & (pl.col("top_n") == n_val)).select(["region", "gene_topn_zscore"])

            paired = df_x.join(df_y, on="region", suffix="_y")
            z_x = paired["gene_topn_zscore"].to_numpy()
            z_y = paired["gene_topn_zscore_y"].to_numpy()

            mean_diff = float(z_y.mean() - z_x.mean())
            if len(z_x) >= 10 and not np.allclose(z_x, z_y):
                _, p_val = stats.wilcoxon(z_x, z_y, alternative="two-sided")
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
                "Window": f"Top {n_val:,}",
                "Window_Int": n_val,
                "Tool_X": label_x,
                "Tool_Y": label_y,
                "mean_diff": mean_diff,
                "sig": sig
            })

    df_heat = pd.DataFrame(heatmap_data)
    df_heat["Window"] = pd.Categorical(df_heat["Window"], categories=df_heat["Window"].unique(), ordered=True)
    df_heat["Tool_X"] = pd.Categorical(df_heat["Tool_X"], categories=ordered_labels, ordered=True)
    df_heat["Tool_Y"] = pd.Categorical(df_heat["Tool_Y"], categories=ordered_labels[::-1], ordered=True)

    n_facets = len(target_k)
    n_tools = len(ordered_tools)

    heatmap = (
        ggplot(df_heat, aes(x="Tool_X", y="Tool_Y", fill="mean_diff"))
        + geom_tile(color="#333333", size=0.5)
        + geom_text(aes(label="sig"), color="black", size=10, va="center", nudge_y=-0.1)
        + facet_wrap("~Window", ncol=n_facets)
        + scale_fill_gradient2(low="#2C7BB6", mid="#FFFFFF", high="#D7191C", midpoint=0)
        + labs(
            title=f"Pairwise Wilcoxon Signed-Rank Test",
            subtitle="*** p<0.001, ** p<0.01, * p<0.05",
            x="Tool X",
            y="Tool Y",
            fill="Tool Y - X"
        )
        + theme_minimal()
        + theme(
            figure_size=(n_facets * (n_tools * 0.5 + 0.8) + 1, n_tools * 0.5 + 1),
            aspect_ratio=1,
            axis_title=element_text(size=11),
            axis_text=element_text(size=10),
            axis_text_x=element_text(rotation=45, hjust=1),
            strip_text=element_text(size=10),
            panel_grid=element_blank(),
            legend_background=element_rect(fill="white", color="white", alpha=0.8),
            plot_background=element_rect(fill="white", color="white"),
        )
    )

    heatmap_path = os.path.join(OUT_DIR, "heatmap.svg")
    heatmap.save(heatmap_path, verbose=False)
    logger.info(f"Heatmap saved to {heatmap_path}")

    # Upload outputs
    uploaded_cum = dxpy.upload_local_file(cum_out_path)
    uploaded_topn = dxpy.upload_local_file(gene_topn_out_path)
    uploaded_lineplot = dxpy.upload_local_file(lineplot_path)
    uploaded_heatmap = dxpy.upload_local_file(heatmap_path)

    logger.info("Upload complete")

    return {
        "cum_stats_parquet": dxpy.dxlink(uploaded_cum),
        "gene_topn_parquet": dxpy.dxlink(uploaded_topn),
        "lineplot_svg": dxpy.dxlink(uploaded_lineplot),
        "heatmap_svg": dxpy.dxlink(uploaded_heatmap),
    }

dxpy.run()
