#!/usr/bin/env python3
"""Fast synthetic check of the pinned MethylVI/scvi-tools API."""

from __future__ import annotations

import anndata as ad
import mudata
import numpy as np
import pandas as pd
from methyl_vi.model import MethylVI


def main() -> None:
    rng = np.random.default_rng(0)
    n_cells, n_features = 32, 64
    cov = rng.integers(0, 8, size=(n_cells, n_features), dtype=np.uint16)
    probabilities = rng.beta(2, 5, size=(n_cells, n_features))
    mc = rng.binomial(cov, probabilities).astype(np.uint16)
    obs = pd.DataFrame(
        {
            "sample_id": pd.Categorical(["IR_01"] * 16 + ["NR_01"] * 16),
            "condition": pd.Categorical(["IR"] * 16 + ["NR"] * 16),
        },
        index=[f"cell_{index}" for index in range(n_cells)],
    )
    adata = ad.AnnData(X=None, obs=obs, var=pd.DataFrame(index=[f"bin_{i}" for i in range(n_features)]))
    adata.layers["mc"] = mc
    adata.layers["cov"] = cov
    mdata = mudata.MuData({"mCG": adata})
    MethylVI.setup_mudata(
        mdata,
        mc_layer="mc",
        cov_layer="cov",
        batch_key="sample_id",
        methylation_modalities={"mCG": "mCG"},
        covariate_modalities={"batch_key": "mCG"},
    )
    model = MethylVI(mdata, n_latent=4, n_hidden=16, n_layers=1)
    model.train(
        max_epochs=2,
        batch_size=16,
        early_stopping=False,
        accelerator="cpu",
        devices=1,
    )
    latent = model.get_latent_representation()
    if latent.shape != (n_cells, 4) or not np.isfinite(latent).all():
        raise RuntimeError(f"Unexpected synthetic latent representation: {latent.shape}")
    print("MethylVI smoke test passed: latent shape", latent.shape)


if __name__ == "__main__":
    main()
