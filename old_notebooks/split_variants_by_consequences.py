import os
import polars as pl

# define consequence groups
variant_groups = {
    "intronic": ['consequence_intron_variant'],
    "missense": ['consequence_missense_variant'],
    "plof": [
        'consequence_frameshift_variant',
        'consequence_start_lost',
        'consequence_stop_gained',
        'consequence_stop_lost',
        'consequence_splice_donor_variant',
        'consequence_splice_acceptor_variant'
    ],
    "5utr_upstream10kb": [
        'consequence_5_prime_utr_variant',
        'consequence_upstream_gene_variant'
    ],
    "3utr_downstream10kb": [
        'consequence_3_prime_utr_variant',
        'consequence_downstream_gene_variant'
    ],
    "other_low_impact": [
        'consequence_splice_donor_5th_base_variant',
        'consequence_splice_donor_region_variant',
        'consequence_splice_polypyrimidine_tract_variant',
        'consequence_splice_region_variant',
        'consequence_inframe_deletion',
        'consequence_inframe_insertion',
        'consequence_intergenic_variant',
        'consequence_non_coding_transcript_exon_variant',
        'consequence_non_coding_transcript_variant',
        'consequence_protein_altering_variant',
        'consequence_start_retained_variant',
        'consequence_stop_retained_variant',
        'consequence_synonymous_variant',
        'consequence_coding_sequence_variant',
        'consequence_NMD_transcript_variant'
    ],
    # NOTE: Splice variants also in plof and other_low_impact. All other groups are mutually exclusive and cover all variants
    "splicing": [
        'consequence_splice_donor_variant',
        'consequence_splice_acceptor_variant',
        'consequence_splice_donor_5th_base_variant',
        'consequence_splice_donor_region_variant',
        'consequence_splice_polypyrimidine_tract_variant',
        'consequence_splice_region_variant',
    ],
}

out_dir = "PATH_TO_FILE"
os.makedirs(out_dir, exist_ok=True)

# available columns in the anno LazyFrame
# anno = pl.scan_parquet("PATH_TO_FILE")
anno = pl.scan_parquet("PATH_TO_FILE")
available_cols = anno.collect_schema().names()

for grp_name, cols in variant_groups.items():
    # keep only cols that actually exist in the anno schema
    cols_present = [c for c in cols if c in available_cols]
    if not cols_present:
        print(f"Skipping {grp_name} (no matching columns)")
        continue

    # build boolean mask: any of the consequence columns == 1
    mask = None
    for c in cols_present:
        expr = (pl.col(c) == 1)
        mask = expr if mask is None else (mask | expr)

    # select a few identifying cols + the consequence cols
    select_cols = ['id', 'region'] + cols_present
    lf = anno.filter(mask)#.select([c for c in select_cols if c in available_cols])

    out_path = os.path.join(out_dir, f"genebass1e6_genes_10kb_{grp_name}_variants_annotated_250925.parquet")

    # schedule sink and trigger execution with streaming engine
    lf.sink_parquet(out_path, engine='streaming')

# After writing group-specific files, also write variants with NONE of the consequence columns == 1
cons_annos = sorted({col for cols_list in variant_groups.values() for col in cols_list})
cols_cons_present = [c for c in cons_annos if c in available_cols]
mask_none = None
for c in cols_cons_present:
    expr = (pl.col(c) == 1)
    mask_none = (~expr) if mask_none is None else (mask_none & ~expr)

out_path_none = os.path.join(out_dir, "genebass1e6_genes_10kb_NOconsequence_variants_250925.parquet")
lf_none = anno.filter(mask_none)#.select([c for c in ['id', 'region'] + cols_cons_present if c in available_cols])
lf_none.sink_parquet(out_path_none, engine='streaming')
print(f"Wrote unannotated variants -> {out_path_none}")
