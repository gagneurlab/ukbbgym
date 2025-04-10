import os
import sys
import yaml
import zarr
import pandas as pd
import numpy as np
from tqdm import tqdm
from anngeno import AnnGeno

def get_gene_burdens(
    anngeno_obj, gene_num, annotation_list, max_burden=False
):
    try:
        gene = anngeno_obj.get_region(gene_num)
    except Exception as e:
        print(f"Cannot find {gene_num} in AnnGeno, Returning NaNs. Error: {e}")
        nan_gis = np.zeros((len(anngeno_obj.samples), len(annotation_list))) * np.nan
        return nan_gis, nan_gis
    
    geno = gene["genotypes"]
    no_variant_mask = geno.sum(axis = 1) == 0

    try:
        var_scores = gene["annotations"][annotation_list].fillna(0).to_numpy()
    except Exception as e:
        print(f"Cannot find {annotation_list} in AnnGeno, Returning NaNs. Error: {e}")
        nan_gis = np.zeros((len(anngeno_obj.samples), len(annotation_list))) * np.nan
        return nan_gis, nan_gis
    
    # Calculate sum burden directly
    gis_sum = np.dot(geno, var_scores)
    gis_sum[no_variant_mask, :] = np.nan

    # If max_burden is False, return sum burden
    if not max_burden:
        return gis_sum, gis_sum

    # For max calculation, process in chunks of annotations
    n_samples, n_variants = geno.shape
    n_annotations = var_scores.shape[1]
    gis_max = np.zeros((n_samples, n_annotations))
    
    # Determine chunk size based on available memory
    # You can adjust this parameter based on your system
    chunk_size = max(1, n_annotations // 5)  # Process ~10% of annotations at a time
    
    for start_idx in range(0, n_annotations, chunk_size):
        end_idx = min(start_idx + chunk_size, n_annotations)
        chunk_scores = var_scores[:, start_idx:end_idx]

        # Process a chunk of annotations at once using broadcasting on a smaller scale
        # Create a view of geno that can be broadcasted with chunk_scores
        geno_view = geno[:, :, np.newaxis]  # Shape: (samples, variants, 1)
        scores_view = chunk_scores[np.newaxis, :, :]  # Shape: (1, variants, chunk_size)
        
        # Compute the product and max for this chunk
        # This uses less memory than processing all annotations at once
        chunk_result = np.max(geno_view * scores_view, axis=1)  # Shape: (samples, chunk_size)
        
        # Store the result for this chunk
        gis_max[:, start_idx:end_idx] = chunk_result
        
        # Explicitly delete temporary arrays to free memory
        del geno_view, scores_view, chunk_result
    
    # set gene scores for sample with no variant to none
    gis_max[no_variant_mask, :] = np.nan
    return gis_sum, gis_max


def get_burdens_array(
    config, 
    associations_df_path, 
    annotation_list,
    anno_scores_df=None, 
    max_burden=False,
    only_snps=False,
):
    maf = config.get("association_testing_maf")

    associations_df = pd.read_parquet(associations_df_path)
    genes = associations_df.gene.unique()

    print("Loading AnnGeno file")
    anngeno_file = config.get("anngeno_file")
    ag = AnnGeno(filename=anngeno_file, filemode="r")

    print(f"Filtering for variants with MAF < {maf}")
    variants_to_keep = set(ag.annotations.query("MAF < @maf & region in @genes")["id"])
    ag.subset_variants(variants_to_keep)

    if only_snps:
        print(f"Filtering for SNPs only")
        snp_variants = set(ag.annotations.query("(ref.str.len()==1) & (alt.str.len()==1)")["id"])
        ag.subset_variants(snp_variants)

    if anno_scores_df is not None:
        # TODO improve this (subset anngeno)
        ## CAUTION very hacky and not stable
        print("Adding annotations to AnnGeno file")
        # new_cols = list(set(annos.columns) - set(ag.annotations.columns))
        new_cols = list(set(anno_scores_df.columns) - set(ag.annotations.columns))
        print(f"Adding new annotations to anngeno: {new_cols}")
        merged = ag.annotations\
            .merge(anno_scores_df[["chrom", "pos", "ref", "alt", "region", *new_cols]], how = "left", on = ["chrom", "pos", "ref", "alt", "region"])
        ag._set_annotations(merged)

    gene_id_list = list(genes)
    gene_burdens_sum = []
    gene_burdens_max = []
    print(f"Starting to compute gene burdens for {len(genes)} genes")
    for gene in tqdm(gene_id_list):
        gis_sum, gis_max = get_gene_burdens(ag, gene, annotation_list, max_burden)
        gene_burdens_sum.append(gis_sum)
        if max_burden:
            gene_burdens_max.append(gis_max)

    gene_burdens_sum_df = np.stack(gene_burdens_sum, axis=1)  # (n_samples, n_genes, n_annotations)
    if max_burden:
        gene_burdens_max_df = np.stack(gene_burdens_max, axis=1)  # (n_samples, n_genes, n_annotations)
        print("Returning sum and max burden.")
        return gene_burdens_sum_df, gene_burdens_max_df, ag.samples, gene_id_list

    print("Max burden is False, returning only sum burden.")
    return gene_burdens_sum_df, gene_burdens_sum_df, ag.samples, gene_id_list


def compute_burdens(
    config_path,
    output_zarr,
    max_burden=False,
    only_snps=False,
    overwrite=False, # TODO add function to overwrite zarr file
):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    associations_df_path = config.get("associations_df_path")
    rare_variant_annotations_dict = config.get('rare_variant_annotations') # Get the nested dict
    all_annotation_list = [] # Initialize empty list
    if rare_variant_annotations_dict: # Flatten the nested dict into a single list
        for category in rare_variant_annotations_dict.values():
            all_annotation_list.extend(category)

    # --- Zarr Writing and Metadata ---
    # Define the path to the zarr file
    zarr_file_path = output_zarr
    anno_chunk_size = 1

    if os.path.exists(zarr_file_path):  # Check if zarr_file_path exists
        print(
            f"Zarr file exists at {zarr_file_path}, checking for new annotations."
        )

        if only_snps:
            print(f"only_snps is True, adding annotations only for SNPs.")

        # Open the zarr group in read/write mode ('r+')
        root = zarr.group(zarr_file_path)
        sum_burdens = root["sum_burdens"]

        existing_annotations = list(root["annotations"][:])
        new_annotation_list = [
            ann for ann in all_annotation_list if ann not in existing_annotations
        ]  # Find annotations not already in zarr

        if new_annotation_list:
            print(f"New annotations found: {new_annotation_list}")
            new_gene_burdens_sum_df, new_gene_burdens_max_df, _, _ = get_burdens_array(config, associations_df_path, new_annotation_list, max_burden=max_burden, only_snps=only_snps)

            current_shape = sum_burdens.shape
            new_shape = (
                current_shape[0],
                current_shape[1],
                current_shape[2] + len(new_annotation_list),
            )
            sum_burdens.resize(new_shape)  # Resize along axis 2 to accommodate new annotations
            sum_burdens[:, :, current_shape[2] :] = (new_gene_burdens_sum_df)

            if max_burden:
                max_burdens = root["max_burdens"]
                print("Max burden is True, appending max_burdens array.")
                max_burdens.resize(new_shape)  # Resize along axis 2 to accommodate new annotations
                max_burdens[:, :, current_shape[2] :] = (new_gene_burdens_max_df)

            print(f"Appended data for new annotations: {new_annotation_list}")

            # Update zarr_annotations array
            zarr_annotations = root["annotations"]  # Assuming 'annotations' array was created in the initial run
            zarr_annotations_current_shape = zarr_annotations.shape
            zarr_annotations_new_shape = (
                len(all_annotation_list),
            )  # Reshape to the new length of annotations list
            
            if zarr_annotations_new_shape != zarr_annotations_current_shape:
                zarr_annotations.resize(zarr_annotations_new_shape)

            zarr_annotations[:] = (all_annotation_list)
            print("Updated 'zarr_annotations' array.")

        else:
            print("No new annotations to add.\nExiting.")

    else:
        print(f"Zarr file does not exist at {zarr_file_path}, creating a new one.")

        gene_burdens_sum_df, gene_burdens_max_df, sample_id_arr, gene_id_list = get_burdens_array(config, associations_df_path, all_annotation_list, max_burden=max_burden, only_snps=only_snps)

        root = zarr.group(zarr_file_path)
        sum_burdens = root.create_array(
            "sum_burdens",
            shape=gene_burdens_sum_df.shape,
            dtype=gene_burdens_sum_df.dtype,
            chunks=(
                gene_burdens_sum_df.shape[0],
                gene_burdens_sum_df.shape[1],
                anno_chunk_size,
            ),  # Chunk only along axis=2
        )
        sum_burdens[:] = gene_burdens_sum_df

        if max_burden:
            print("Max burden is True, creating max_burdens array.")
            max_burdens = root.create_array(
                "max_burdens",
                shape=gene_burdens_max_df.shape,
                dtype=gene_burdens_max_df.dtype,
                chunks=(
                    gene_burdens_max_df.shape[0],
                    gene_burdens_max_df.shape[1],
                    anno_chunk_size,
                ),  # Chunk only along axis=2
            )
            max_burdens[:] = gene_burdens_max_df

        zarr_samples = root.create_array("samples", shape=sample_id_arr.shape, dtype="str")
        zarr_samples[:] = sample_id_arr

        zarr_genes = root.create_array("genes", shape=len(gene_id_list), dtype="str")
        zarr_genes[:] = gene_id_list

        zarr_annotations = root.create_array(
            "annotations",
            shape=len(all_annotation_list),  # Use all_annotation_list here
            dtype="str",
        )
        zarr_annotations[:] = all_annotation_list  # Use all_annotation_list here
        print("Created new zarr array 'gene_burdens', 'samples', 'genes', 'annotations' and metadata.")


def add_new_anno_burdens(
    config_path,
    output_zarr,
    anno_scores_path,
    max_burden=False,
    overwrite=False, # TODO add function to overwrite exiting burdens file
):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    associations_df_path = config.get("associations_df_path")

    anno_scores_df = pd.read_parquet(anno_scores_path)
    # TODO add assertions to check if anno_scores_df is valid
    all_annotation_list = list(set(anno_scores_df.columns) - set(["chrom", "pos", "ref", "alt", "region"]))

    # --- Zarr Writing and Metadata ---
    # Define the path to the zarr file
    zarr_file_path = output_zarr
    anno_chunk_size = 1

    if os.path.exists(zarr_file_path):  # Check if zarr_file_path exists
        print(
            f"Zarr file exists at {zarr_file_path}, checking for new annotations."
        )
        # Open the zarr group in read/write mode ('r+')
        root = zarr.group(zarr_file_path)
        sum_burdens = root["sum_burdens"]

        existing_annotation_list = list(root["annotations"][:])
        new_annotation_list = list(set(all_annotation_list) - set(existing_annotation_list))  # Find annotations not already in zarr

        union_annotation_list = existing_annotation_list + new_annotation_list  # Find annotations already in zarr
        # new_annotation_list = [
        #     ann for ann in all_annotation_list if ann not in existing_annotations
        # ]  # Find annotations not already in zarr

        if new_annotation_list:
            print(f"New annotations found: {new_annotation_list}")
            new_gene_burdens_sum_df, new_gene_burdens_max_df, _, _ = get_burdens_array(config, associations_df_path, new_annotation_list, anno_scores_df=anno_scores_df, max_burden=max_burden)

            current_shape = sum_burdens.shape
            new_shape = (
                current_shape[0],
                current_shape[1],
                current_shape[2] + len(new_annotation_list),
            )
            sum_burdens.resize(new_shape)  # Resize along axis 2 to accommodate new annotations
            sum_burdens[:, :, current_shape[2] :] = (new_gene_burdens_sum_df)

            if max_burden:
                print("Max burden is True, appending max_burdens array.")
                max_burdens = root["max_burdens"]
                max_burdens.resize(new_shape)  # Resize along axis 2 to accommodate new annotations
                max_burdens[:, :, current_shape[2] :] = (new_gene_burdens_max_df)

            print(f"Appended data for new annotations: {new_annotation_list}")

            # Update zarr_annotations array
            zarr_annotations = root["annotations"]  # Assuming 'annotations' array was created in the initial run
            zarr_annotations_current_shape = zarr_annotations.shape
            zarr_annotations_new_shape = (
                len(union_annotation_list),
            )  # Reshape to the new length of annotations list
            
            if zarr_annotations_new_shape != zarr_annotations_current_shape:
                zarr_annotations.resize(zarr_annotations_new_shape)

            zarr_annotations[:] = (union_annotation_list)
            print("Updated 'zarr_annotations' array.")

        else:
            print("No new annotations to add.\nExiting.")

    else:
        print(f"Zarr file does not exist at {zarr_file_path}.")
        print(f"To add new annotations, please create a new zarr file first.")
        print(f"Exiting.")
        sys.exit(1)
