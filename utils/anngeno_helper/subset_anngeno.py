import polars as pl
from pathlib import Path
import numpy as np
import dask.array as da
import pandas as pd

import time
import dask.array as da
import zarr
import os
import numpy as np
from dask.distributed import Client, LocalCluster
import os
import time
from tqdm import tqdm
from anngeno import AnnGeno

import polars as pl
from pathlib import Path
import numpy as np
import dask.array as da
import time
import zarr
import os
from dask.distributed import Client, LocalCluster
from anngeno import AnnGeno

if __name__ == "__main__":
    test_out_dir = Path("/home/dnanexus/data_dir/")
    ag = AnnGeno(
        test_out_dir / "dms_anngeno.ag", low_mem=True, filemode="r", sanity_check=False
    )

    ids_to_keep = (
        ag.annotations.filter(
            (pl.col("relative_cds_position_is_nan") == 0)
            | (pl.col("loftee_hc_is_nan") == 0)
        )
        .select("id")
        .collect()["id"]
        .unique()
        .to_list()
    )
    print(len(ids_to_keep))

    variant_mask = ag.all_variant_metadata["id"].is_in(ids_to_keep).to_numpy()
    variant_mask_concrete = variant_mask
    final_new_chunk_shape = (10_000, 1000, 2)
    # target_zarr_path = f"/home/dnanexus/data_dir/seed_genes_{final_new_chunk_shape[0]}_{final_new_chunk_shape[1]}.zarr"
    target_zarr_path = "/home/dnanexus/data_dir/new_genotypes.zarr"
    worker_memory_limit = "6GiB"  # Or '16GiB', '20GiB', etc.
    source_zarr_path = "/home/dnanexus/data_dir/dms_anngeno.ag/zarr_store/genotypes/"

    n_workers = os.cpu_count()
    n_workers = 16
    cluster = LocalCluster(
        n_workers=n_workers,
        processes=True,
        threads_per_worker=1,
        memory_limit=worker_memory_limit,
    )
    client = Client(cluster)
    print(f"Dask Dashboard: {client.dashboard_link}")
    print("Waiting for workers to start...")
    client.wait_for_workers(n_workers=n_workers)
    print("Dask cluster ready.")

    print("Subsetting variants and samples .")
    source_array = da.from_zarr(source_zarr_path)
    masked_array = source_array[variant_mask_concrete, :, :]

    print(f"Masked array shape: {masked_array.shape}")

    print("rechunking")
    rechunked_masked_array = masked_array.rechunk(final_new_chunk_shape)
    print(f"Writing rechunked and masked array to: {target_zarr_path}")
    print("This will take a significant amount of time and disk space...")

    rechunked_masked_array.to_zarr(
        target_zarr_path, overwrite=False, zarr_format=3
    )  # overwrite=False is a safety measure

    print("done")
