#!/usr/bin/env python3
"""Small plotting helpers shared by the workflow's figure stages."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def categorical_embedding_plot(
    table: pd.DataFrame,
    x: str,
    y: str,
    color: str,
    output: Path,
    title: str,
    seed: int = 0,
) -> None:
    if color not in table:
        raise ValueError(f"Missing plotting annotation: {color}")
    data = table[[x, y, color]].copy()
    data[color] = data[color].fillna("Unknown").astype(str)
    categories = sorted(data[color].unique())
    cmap = plt.get_cmap("tab20" if len(categories) <= 20 else "gist_ncar")
    colors = {category: cmap(index / max(1, len(categories) - 1)) for index, category in enumerate(categories)}
    order = np.random.default_rng(seed).permutation(len(data))
    data = data.iloc[order]

    width = 11 if len(categories) > 12 else 8
    figure, axis = plt.subplots(figsize=(width, 7))
    for category in categories:
        subset = data[data[color] == category]
        axis.scatter(
            subset[x],
            subset[y],
            s=4,
            alpha=0.75,
            linewidths=0,
            color=colors[category],
            label=category,
        )
    axis.set(xlabel=x, ylabel=y, title=title)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.spines[:].set_visible(False)
    columns = 2 if len(categories) > 20 else 1
    axis.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        markerscale=2.5,
        fontsize=7,
        ncol=columns,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
