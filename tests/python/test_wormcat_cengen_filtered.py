import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "wormcat_cengen_filtered", ROOT / "scripts" / "analysis" / "05_wormcat_cengen_filtered.py"
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_bh_matches_r_p_adjust():
    # R: p.adjust(c(0.01, 0.04, 0.03, 0.2), method = "BH")
    adjusted = module.bh(pd.Series([0.01, 0.04, 0.03, 0.2]))
    np.testing.assert_allclose(adjusted, [0.04, 0.16 / 3, 0.16 / 3, 0.2])


def test_filter_and_enrichment():
    pairs = pd.DataFrame({
        "cell_type": ["A"] * 4 + ["B"],
        "stable_id": ["g1", "g2", "g3", "g4", "g1"],
        "final_hvg": [True, True, False, False, True],
        "cat1": ["X", "X", "Y", None, "X"],
        "cat2": [None] * 5,
        "cat3": [None] * 5,
        # g2 not expressed in A; g4 has no call; type B has no calls at all.
        "expressed": [True, False, True, np.nan, np.nan],
    })
    filtered = module.cengen_filter(pairs)
    assert list(filtered.stable_id) == ["g1", "g3", "g4"]

    result = module.wormcat_enrichment(filtered).set_index("wormcat_family")
    # 1 HVG of 3 background genes; category X holds the HVG and 1 gene.
    assert result.loc["X", "RGS"] == 1 and result.loc["X", "AC"] == 1
    assert result.loc["X", "fold_enrichment"] == (1 / 1) / (1 / 3)
    assert "Y" not in result.index  # categories without HVGs are not tested
