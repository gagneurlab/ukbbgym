import os
import polars as pl

# define consequence groups
variant_groups = {
    "intronic": ['Consequence_intron_variant'],
    "missense": ['Consequence_missense_variant'],
    "plof": [
        'Consequence_frameshift_variant',
        'Consequence_start_lost',
        'Consequence_stop_gained',
        'Consequence_stop_lost',
        'Consequence_splice_donor_variant',
        'Consequence_splice_acceptor_variant'
    ],
    # Splice variants also in plof and other_low_impact. All other groups are mutually exclusive and cover all variants
    "splice": [
        'Consequence_splice_donor_variant',
        'Consequence_splice_acceptor_variant',
        'Consequence_splice_donor_5th_base_variant',
        'Consequence_splice_donor_region_variant',
        'Consequence_splice_polypyrimidine_tract_variant',
        'Consequence_splice_region_variant'
    ],
    "utr": [
        'Consequence_3_prime_UTR_variant',
        'Consequence_5_prime_UTR_variant'
    ],
    "upstream_downstream": [
        'Consequence_upstream_gene_variant',
        'Consequence_downstream_gene_variant'
    ],
    "other_low_impact": [
        'Consequence_splice_donor_5th_base_variant',
        'Consequence_splice_donor_region_variant',
        'Consequence_splice_polypyrimidine_tract_variant',
        'Consequence_splice_region_variant',
        'Consequence_inframe_deletion',
        'Consequence_inframe_insertion',
        'Consequence_intergenic_variant',
        'Consequence_non_coding_transcript_exon_variant',
        'Consequence_non_coding_transcript_variant',
        'Consequence_protein_altering_variant',
        'Consequence_start_retained_variant',
        'Consequence_stop_retained_variant',
        'Consequence_synonymous_variant',
        'Consequence_coding_sequence_variant',
        'Consequence_NMD_transcript_variant'
    ]
}

out_dir = "PATH_TO_FILE"
os.makedirs(out_dir, exist_ok=True)

# available columns in the anno LazyFrame
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

    out_path = os.path.join(out_dir, f"genebass1e6_{grp_name}_variants.parquet")

    # schedule sink and trigger execution with streaming engine
    lf.sink_parquet(out_path, engine='streaming')


# After writing group-specific files, also write variants with NONE of the consequence columns == 1
cols_cons_present = [c for c in cons_annos if c in available_cols]
mask_none = None
for c in cols_cons_present:
    expr = (pl.col(c) == 1)
    mask_none = (~expr) if mask_none is None else (mask_none & ~expr)

out_path_none = os.path.join(out_dir, "genebass1e6_NOconsequence_variants.parquet")
lf_none = anno.filter(mask_none)#.select([c for c in ['id', 'region'] + cols_cons_present if c in available_cols])
lf_none.sink_parquet(out_path_none, engine='streaming')
print(f"Wrote unannotated variants -> {out_path_none}")
