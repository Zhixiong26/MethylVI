#!/usr/bin/env python3
"""Plot sample and IR/NR labels on the original ALLCools embedding."""

from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import pandas as pd

from mvi_14_utils_pipeline import env_path, load_annotations
from mvi_15_utils_plot import categorical_embedding_plot


def main() -> None:
    h5ad = env_path("MVI_H5AD")
    root = env_path("MVI_RESULTS")
    annotation_string = os.environ.get("MVI_ANNOTATION")
    annotation = Path(annotation_string).expanduser().resolve() if annotation_string else None
    sample_metadata = env_path("MVI_SAMPLE_METADATA")
    sample_id_regex = os.environ.get("MVI_SAMPLE_ID_REGEX", r"^([^_]+_[^_]+)_")
    adata = ad.read_h5ad(h5ad, backed="r")
    if "X_umap" in adata.obsm:
        coordinates = adata.obsm["X_umap"]
        x, y, name = "UMAP1", "UMAP2", "UMAP"
    elif "X_tsne" in adata.obsm:
        coordinates = adata.obsm["X_tsne"]
        x, y, name = "tSNE1", "tSNE2", "t-SNE"
    else:
        raise ValueError("Original H5AD contains neither X_umap nor X_tsne")
    annotations, _ = load_annotations(adata.obs_names, annotation, sample_metadata, sample_id_regex)
    table = pd.DataFrame({x: coordinates[:, 0], y: coordinates[:, 1]}, index=adata.obs_names)
    for column, label in (("sample_id", "sample"), ("condition", "condition (IR/NR)")):
        table[column] = annotations[column].to_numpy()
        categorical_embedding_plot(
            table,
            x,
            y,
            column,
            root / f"allcools_original_embedding_{column}.pdf",
            f"Original ALLCools {name} — {label}",
        )


if __name__ == "__main__":
    main()
