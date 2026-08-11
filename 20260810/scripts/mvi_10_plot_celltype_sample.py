#!/usr/bin/env python3
"""Plot the MethylVI UMAP by cell type and sample."""

from __future__ import annotations

import pandas as pd

from mvi_13_utils_pipeline import env_path
from mvi_14_utils_plot import categorical_embedding_plot


def main() -> None:
    results = env_path("MVI_RESULTS")
    figures = env_path("MVI_FIGURES_AFTER_DIR")
    table = pd.read_csv(results / "cell_annotations_umap.tsv.gz", sep="\t", index_col=0)
    for column, label in (("cell_type", "cell type"), ("sample_id", "sample")):
        categorical_embedding_plot(
            table,
            "UMAP1",
            "UMAP2",
            column,
            figures / f"methylvi_umap_{column}.pdf",
            f"MethylVI UMAP — {label}",
        )


if __name__ == "__main__":
    main()
