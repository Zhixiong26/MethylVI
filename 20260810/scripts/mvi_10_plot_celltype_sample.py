#!/usr/bin/env python3
"""Plot the MethylVI UMAP by cell type and sample."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from mvi_13_utils_pipeline import env_path, load_annotations
from mvi_14_utils_plot import categorical_embedding_plot


def main() -> None:
    results = env_path("MVI_RESULTS")
    figures = env_path("MVI_FIGURES_AFTER_DIR")
    table = pd.read_csv(results / "cell_annotations_umap.tsv.gz", sep="\t", index_col=0)
    annotation_string = os.environ.get("MVI_ANNOTATION")
    annotation = Path(annotation_string).expanduser().resolve() if annotation_string else None
    sample_metadata = env_path("MVI_SAMPLE_METADATA")
    sample_id_regex = os.environ.get("MVI_SAMPLE_ID_REGEX", r"^([^_]+_[^_]+)_")
    annotations, stats = load_annotations(
        table.index,
        annotation,
        sample_metadata,
        sample_id_regex,
    )
    for column in annotations.columns:
        table[column] = annotations[column].to_numpy()
    print(
        f"已从当前注释表刷新cell type: "
        f"matched={stats.get('fully_annotated_selected_cells', 0):,}, "
        f"unmatched={stats.get('annotation_unmatched_selected_cells', 0):,}",
        flush=True,
    )
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
