"""Shared variant-filtering helpers for genebass_sumstats/analysis notebooks.

Consolidates the load_config / add_derived / scan_variants / appv_of / pick_annos /
coverage-filter boilerplate that used to be copy-pasted (with small, undocumented drift)
into every analysis notebook. See genebass_sumstats/configs/config_variant_classes.yaml
for the variant-class filter definitions these helpers apply.
"""

import os
from pathlib import Path

import yaml
import polars as pl
from huggingface_hub import hf_hub_download


def fetch_hf_data(filename, repo_root=None, repo_id='gagneurlab/ukbbgym'):
    """data/<filename>, downloading it from the private gagneurlab/ukbbgym HF dataset first
    if missing (needs `huggingface-cli login` or HF_TOKEN -- see genebass/README.md)."""
    root = Path(repo_root) if repo_root else next(
        p for p in [Path.cwd(), *Path.cwd().parents] if (p / 'utils' / 'variant_filtering.py').exists())
    local = root / 'data' / filename
    if local.exists():
        return str(local)
    return hf_hub_download(repo_id=repo_id, repo_type='dataset', filename=filename,
                            local_dir=str(root / 'data'))


# merge_new_score()'s auto-detected join keys, in try-order, as master-column-name -> the
# aliases a score file's own column might use for it. Genomic goes first (variant-level,
# works for any consequence); protein position is the fallback for scores that only make
# sense per amino-acid substitution (e.g. a per-protein LLR table with no genomic coordinates
# at all) -- both are genuinely used by the pipeline's own annotation files (see
# utils/annotations/add_more_annotations.py), not a guess at what a score file might look like.
_GENOMIC_KEY = {'chrom': ('chrom', 'chromosome', 'chr'), 'pos': ('pos', 'position'),
                'ref': ('ref',), 'alt': ('alt',)}
_PROTEIN_KEY = {'uniprot_id': ('uniprot_id',), 'position': ('position', 'protein_position'),
                'aa_ref': ('aa_ref',), 'aa_alt': ('aa_alt',)}


def _resolve_key(score_cols, key_aliases):
    """key_aliases values checked in order against score_cols; {master_col: score_col} if
    every key resolves, else None."""
    resolved = {}
    for master_col, aliases in key_aliases.items():
        hit = next((a for a in aliases if a in score_cols), None)
        if hit is None:
            return None
        resolved[master_col] = hit
    return resolved


def merge_new_score(master_path, score_path, out_path=None, on=None):
    """Left-join score_path's new columns onto master_path's variants, writing a new parquet.
    Returns (out_path, new_numeric_column_names) -- the latter for register_score_in_config():
    every new column is merged in, but non-numeric ones (transcript/uniprot ids, gene symbols,
    ...) aren't predictors, so they're excluded from what gets registered as one.
    score_path can be .parquet or .csv(.gz).

    `on` picks the join key: pass an explicit column name (or list, for a compound key) to
    force it. Left as None (the default), it's auto-detected -- genomic (chrom, pos, ref, alt,
    tried first) if score_path has those columns (under common aliases), else protein
    (uniprot_id, position, aa_ref, aa_alt) if it has those instead. `uniprot_id` is folded
    into an already-genomic key too, when both sides have it, to disambiguate a score given
    per (variant, isoform) rather than per variant.

    Rows in score_path that repeat a join key are collapsed to their first occurrence (a
    warning is printed with the count) so the join can't inflate master_path's row count --
    e.g. a genomic key against a file scored per isoform, where only the join key's own
    columns actually needed the isoform to be found, not to stay distinct in the output.

    Columns score_path shares with the master table (beyond the join key) are dropped, so a
    rerun with an updated score file is safe. out_path defaults to
    '<master>_plus_<score file stem>.parquet' next to master_path."""
    master_lf = pl.scan_parquet(master_path)
    if str(score_path).endswith('.csv') or str(score_path).endswith('.csv.gz'):
        score_lf = pl.scan_csv(score_path)
    else:
        score_lf = pl.scan_parquet(score_path)
    master_cols = set(master_lf.collect_schema().names())
    score_schema = score_lf.collect_schema()
    score_cols = score_schema.names()

    if on is not None:
        key = {k: k for k in (on if isinstance(on, list) else [on])}
    else:
        key = _resolve_key(score_cols, _GENOMIC_KEY) or _resolve_key(score_cols, _PROTEIN_KEY)
        if key is None:
            raise ValueError(
                f"{score_path}: found neither a genomic key {list(_GENOMIC_KEY)} nor a "
                f"protein key {list(_PROTEIN_KEY)} (checked common aliases) -- pass `on` "
                f"explicitly if it uses different column names")
        if set(key) == set(_GENOMIC_KEY) and 'uniprot_id' in master_cols and 'uniprot_id' in score_cols:
            key['uniprot_id'] = 'uniprot_id'

    score_lf = score_lf.rename({v: k for k, v in key.items() if v != k})
    join_keys = list(key)

    value_cols = [c for c in score_cols if c not in join_keys and c not in master_cols]
    if not value_cols:
        raise ValueError(f"{score_path} has no columns beyond the join key {join_keys} that "
                          f"aren't already in {master_path}")
    numeric_cols = [c for c in value_cols if score_schema[c].is_numeric()]

    # A score file (e.g. a whole-proteome table) can be far bigger than what's relevant to
    # master_path's ~611 genes -- semi-join down to master's own key values first. A left
    # join only keeps rows that match one of these anyway, so this doesn't change the
    # result, just how much of score_path everything below has to look at. Collected once,
    # eagerly, right here: score_lf is still a lazy plan rooted at the full score file, so
    # each further .collect() on it would otherwise re-run this scan+semi-join from scratch.
    score_df = (score_lf.select(join_keys + value_cols)
                .join(master_lf.select(join_keys).unique(), on=join_keys, how='semi')
                .collect())

    n_before = score_df.height
    n_after = score_df.select(join_keys).unique().height
    if n_after < n_before:
        print(f"merge_new_score: {score_path} has {n_before - n_after} rows sharing a join "
              f"key {join_keys} with another row -- keeping the first of each")
        score_df = score_df.unique(subset=join_keys, keep='first')

    if len(numeric_cols) < len(value_cols):
        print(f"merge_new_score: {score_path} -- merged in all of {value_cols}, but only the "
              f"numeric columns {numeric_cols} will be registered as predictors")

    out_path = out_path or f"{Path(master_path).with_suffix('')}_plus_{Path(score_path).stem}.parquet"
    master_lf.join(score_df.lazy(), on=join_keys, how='left').sink_parquet(out_path)
    return str(out_path), numeric_cols


def register_score_in_config(config_dir, config_file, columns, category='missense', color='#e6194b'):
    """Copy config_file with each name in `columns` added as a predictor under `category`
    (creating the category if it's new), written to config_dir/_custom_score.yaml
    (gitignored) so a merged-in score gets picked up by pick_annos() without touching the
    checked-in config. Returns '_custom_score.yaml' -- pass that as CONFIG_FILE."""
    cfg = yaml.safe_load(open(f'{config_dir}/{config_file}'))
    section = cfg.setdefault('rare_variant_annotations', {}).setdefault(category, {})
    for name in columns:
        section[name] = {'color': color, 'label': name}
    out_name = '_custom_score.yaml'
    yaml.safe_dump(cfg, open(f'{config_dir}/{out_name}', 'w'), sort_keys=False)
    return out_name


def env_override(name, default, cast=str):
    """Read UKBBGYM_<name> from the environment if set, else return `default` unchanged.

    Lets a notebook's parameter cell double as both its own default (when run interactively)
    and something a driver script can override without editing the notebook (e.g. run_all.sh
    exporting UKBBGYM_MAC=10 before invoking nbconvert).

    cast: str (default) | int | bool (case-insensitive '1'/'true'/'yes' -> True, else False) |
    'list' (comma-split, whitespace-trimmed, empty items dropped).
    """
    val = os.environ.get(f'UKBBGYM_{name}')
    if val is None:
        return default
    if cast is bool:
        return val.strip().lower() in ('1', 'true', 'yes')
    if cast == 'list':
        return [x.strip() for x in val.split(',') if x.strip()]
    return cast(val)


def load_config(config_dir, config_file):
    """config_*.yaml -> (anno_config_df, all_annotation_list)."""
    cfg = yaml.safe_load(open(f'{config_dir}/{config_file}'))
    df = pl.DataFrame([
        {'category': c, 'annotation': a, 'color': p['color'], 'label': p['label'],
         'annotation_dir': p.get('direction', 1)}
        for c, annos in cfg['rare_variant_annotations'].items() for a, p in annos.items()
    ]).with_columns(pl.col('annotation_dir').cast(pl.Int8))
    return df, df['annotation'].to_list()


def load_variant_class(config_dir, variant_class):
    """Entry from config_variant_classes.yaml."""
    return yaml.safe_load(open(f'{config_dir}/config_variant_classes.yaml'))[variant_class]


def add_derived(lf):
    """Annotations the configs name but the master table does not store -- cheap to
    recompute, so they are derived on the fly rather than materialised on 20.8M rows."""
    return (
        lf
        .with_columns(
            is_ins = pl.col('ref').str.len_chars() < pl.col('alt').str.len_chars(),
            is_del = pl.col('ref').str.len_chars() > pl.col('alt').str.len_chars(),
            non_ted_domain = pl.col('ted_domain') == False,
            high_plddt = (pl.col('plddt') > 70).cast(pl.Int8),
            low_plddt  = (pl.col('plddt') <= 70).cast(pl.Int8),
            loftee_disorder       = pl.all_horizontal(pl.col('loftee_hc') == True, pl.col('mobi_curated_disorder_priority') == True),
            loftee_lip            = pl.all_horizontal(pl.col('loftee_hc') == True, pl.col('mobi_full_lip_priority') == True),
            loftee_ted            = pl.all_horizontal(pl.col('loftee_hc') == True, pl.col('ted_domain') == True),
            loftee_low_complexity = pl.all_horizontal(pl.col('loftee_hc') == True, pl.col('low_complexity_domain') == True),
            has_inter_chain_hydrogen_bond_pdb          = pl.col('inter_chain_hydrogen_bond_pdb_count') > 0,
            has_inter_chain_non_bonded_interaction_pdb = pl.col('inter_chain_non_bonded_interaction_pdb_count') > 0,
            has_inter_chain_disulfide_bond_pdb         = pl.col('inter_chain_disulfide_bond_pdb_count') > 0,
            has_inter_chain_salt_bridge_pdb            = pl.col('inter_chain_salt_bridge_pdb_count') > 0,
            encode_any_tf    = pl.any_horizontal(pl.col('encode_tf') == 1, pl.col('encode_ca-tf') == 1),
            promoterai_abs   = pl.col('promoterai').abs(),
            promoterai_under = pl.col('promoterai'),
            promoterai_over  = pl.col('promoterai'),
            # NOTE: these regexes are non-exclusive by design (Pathogenic also matches
            # Likely_pathogenic), matching avg_zscore_categories.ipynb.
            clinvar_patho        = pl.col('clinical_significance').str.contains('(?i)Pathogenic').fill_null(False),
            clinvar_likely_patho = pl.col('clinical_significance').str.contains('(?i)Likely_pathogenic').fill_null(False),
            clinvar_benign       = pl.col('clinical_significance').str.contains('(?i)Benign').fill_null(False),
            clinvar_likely_benign= pl.col('clinical_significance').str.contains('(?i)Likely_benign').fill_null(False),
        )
        .with_columns(
            inframe_deletion  = (pl.col('variant_length') % 3 == 1) & (pl.col('is_del') == True),
            inframe_insertion = (pl.col('variant_length') % 3 == 1) & (pl.col('is_ins') == True),
        )
    )


def build_variant_filters(vc, only_snps=False, only_clinvar=False, exclude_clinvar=False,
                           max_variant_length=None):
    """variant-class filter expressions + the common ad-hoc extras, as a list of pl.Expr.

    Split out from scan_variants so it can be applied directly to a LazyFrame that hasn't
    gone through add_derived (e.g. an annotation file that isn't the master table)."""
    f = [eval(e) for e in (vc.get('variant_filtering') or [])]
    if exclude_clinvar:
        f.append(pl.col('clinical_significance').is_null())
    elif only_clinvar:
        f.append(pl.col('clinical_significance').is_not_null())
    if only_snps:
        f.append((pl.col('ref').str.len_chars() == 1) & (pl.col('alt').str.len_chars() == 1))
    if max_variant_length is not None:
        f.append(pl.col('variant_length') <= max_variant_length)
    return f


def scan_variants(source, vc, only_snps=False, only_clinvar=False, exclude_clinvar=False,
                   max_variant_length=None):
    """Master table (path or LazyFrame) + derived columns, with the variant-class
    filters applied."""
    lf = source if isinstance(source, pl.LazyFrame) else pl.scan_parquet(source)
    lf = add_derived(lf)
    f = build_variant_filters(vc, only_snps=only_snps, only_clinvar=only_clinvar,
                               exclude_clinvar=exclude_clinvar,
                               max_variant_length=max_variant_length)
    return lf.filter(*f) if f else lf


def derived_schema(master_path):
    """Column names after add_derived() -- the schema pick_annos() needs, without every
    notebook re-deriving the whole LazyFrame by hand just to read off its column names."""
    return add_derived(pl.scan_parquet(master_path)).collect_schema().names()


def appv_of(lf, mac):
    """Per-variant genebass beta at the MAC<=mac proxy. Rows without a beta are dropped here."""
    return (lf.filter(pl.col('mean_pheno_value').is_not_null())
              .filter(pl.col('AF') <= mac / (2 * pl.col('n_cases')))
              .select(['id', 'region', 'phenotype', 'mean_pheno_value', 'AF'])
              .unique())


def pick_annos(anno_config_df, all_annotation_list, selected_categories, schema_names):
    """Config annotations for the chosen categories that actually exist after add_derived."""
    want = anno_config_df.filter(pl.col('category').is_in(selected_categories))['annotation'].to_list()
    return list(set(want) & {c for c in all_annotation_list if c in schema_names})


def build_gene_trait_tool_correlations(lf, mac, annos, anno_config_df, selected_categories,
                                        extra_group_cols=()):
    """Per (gene, trait, tool[, extra...]) direction-corrected Spearman rho at MAC<=mac.
    No coverage/min-variants filter applied -- pass the result through filter_covered()
    yourself, or use gene_trait_tool_correlations() below for the common case where you want
    that applied immediately.

    `extra_group_cols` adds further grouping keys (e.g. `('in_domain',)`) to both the rank
    `.over()` and the `.group_by()`, for analyses that split the correlation by more than
    just (gene, trait, tool).
    """
    group_cols = ['region', 'phenotype', 'annotation'] + list(extra_group_cols)
    melt_index = ['id', 'region'] + list(extra_group_cols)
    melted = (lf
        .select(set(melt_index) | set(annos))
        .unpivot(index=melt_index, on=annos, variable_name='annotation', value_name='annotation_score')
        .with_columns(pl.col('annotation_score').cast(pl.Float32), pl.col('region').cast(pl.Utf8)))

    return (
        appv_of(lf, mac).select(['id', 'region', 'phenotype', 'mean_pheno_value'])
        .join(melted, on=['id', 'region'], how='inner')
        .with_columns(pl.col(c).rank('average').over(group_cols).alias(f'{c}_rank')
                      for c in ['annotation_score', 'mean_pheno_value'])
        .group_by(group_cols)
        .agg(n_variants=pl.col('id').count(),
             correlation=pl.when((pl.col('annotation_score_rank').n_unique() > 1) &
                                 (pl.col('mean_pheno_value_rank').n_unique() > 1))
               .then(pl.corr('annotation_score_rank', 'mean_pheno_value_rank', propagate_nans=True))
               .otherwise(None))
        .drop_nans().drop_nulls()
        .collect(engine='streaming')
        .join(anno_config_df.filter(pl.col('category').is_in(selected_categories)), on='annotation')
        .join(lf.select(['region', 'phenotype', 'loftee_corr', 'loftee_corr_dir'])
                .unique().collect(engine='streaming'), on=['region', 'phenotype'])
        .with_columns(corr_beta=pl.col('correlation') * pl.col('loftee_corr_dir') * pl.col('annotation_dir'))
        .with_columns(corr_beta_rescaled=(pl.col('corr_beta') / pl.col('loftee_corr')) * pl.col('loftee_corr_dir'))
    )


def gene_trait_tool_correlations(lf, mac, annos, anno_config_df, selected_categories, min_variants,
                                  coverage_group_col=('region', 'phenotype'), extra_group_cols=()):
    """build_gene_trait_tool_correlations(...) + filter_covered(...) -- the single canonical
    "which (gene, trait) pairs have enough, sufficiently-covered data" computation. Every
    analysis that needs to know which UKBBGym genes qualify under a given tool set/mac/
    min_variants should call this (or build_gene_trait_tool_correlations directly, if the
    coverage filter needs to be deferred -- e.g. until a cross-dataset common tool set is
    known) instead of re-deriving its own version, so they agree on the exact same universe
    by construction rather than by coincidence."""
    correlation_df = build_gene_trait_tool_correlations(
        lf, mac, annos, anno_config_df, selected_categories, extra_group_cols=extra_group_cols)
    return filter_covered(correlation_df, group_col=list(coverage_group_col),
                          n_variants_col='n_variants', min_variants=min_variants)


def pick_best_trait_per_gene(df, group_col='region'):
    """Among a gene's surviving (region, phenotype) rows, keep only the phenotype with the
    strongest |loftee_corr| -- the same "most confidently associated trait" choice notebooks
    made before filtering, now applied after, so a gene is never dropped just because its
    single best trait happened to lack coverage."""
    best_pheno = (df.select([group_col, 'phenotype', 'loftee_corr']).unique()
                    .with_columns(loftee_corr_abs=pl.col('loftee_corr').abs())
                    .sort('loftee_corr_abs', descending=True)
                    .unique(subset=[group_col], keep='first', maintain_order=True)
                    .select([group_col, 'phenotype']))
    return df.join(best_pheno, on=[group_col, 'phenotype'], how='semi')


def filter_covered(df, group_col, annotation_col='annotation', n_variants_col=None,
                    min_variants=None, enabled=True):
    """Keep only rows for groups (genes, protein-units, ...) scored by every annotation
    present in `df` (the "coverage filter"), optionally also dropping groups with too few
    variants.

    `enabled=False` skips the coverage semi-join (still applies `min_variants` if given) --
    used where the coverage filter is intentionally turned off for a given notebook.

    Coverage is computed against the full (unfiltered) `df` -- i.e. before `min_variants`
    is applied -- matching the original per-notebook implementations this replaces."""
    covered = None
    if enabled:
        k = df[annotation_col].n_unique()
        covered = (df.group_by(group_col).agg(n=pl.col(annotation_col).n_unique())
                     .filter(pl.col('n') == k).select(group_col))
    out = df
    if n_variants_col is not None and min_variants is not None:
        out = out.filter(pl.col(n_variants_col) > min_variants)
    if covered is not None:
        out = out.join(covered, on=group_col, how='semi')
    return out
