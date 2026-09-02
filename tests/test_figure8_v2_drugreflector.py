import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "figure8_v2_drugreflector_inference.py"
if not SCRIPT.exists():
    raise RuntimeError("Expected RED failure: Figure 8 v2 DrugReflector adapter is not implemented")

spec = importlib.util.spec_from_file_location("figure8_v2_dr", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_checkpoint_hashes_are_frozen_and_distinct() -> None:
    records = module.validate_checkpoints(ROOT / "metadata" / "driver" / "drugreflector_checkpoints")
    assert len(records) == 3
    assert len({row["md5"] for row in records}) == 3
    assert [row["md5"] for row in records] == [
        "0a27e253713c37f4874318b5ba0c27a9",
        "0e785196fd046d946f84e4480c81ff53",
        "d8e36f6a8f9fa7a22feda7acdd0bee86",
    ]


def test_model_order_requires_all_978_columns_in_exact_order(tmp_path: Path) -> None:
    model = pd.read_csv(
        ROOT / "metadata" / "driver" / "figure8_transcriptomic_reversal" / "figure8_drugreflector_model_genes.tsv",
        sep="\t",
    )
    frame = pd.DataFrame([[0.0] * 978], columns=model["gene"], index=["sig"])
    module.validate_model_order(frame, model["gene"].tolist())
    with pytest.raises(ValueError, match="exact frozen 978-gene order"):
        module.validate_model_order(frame.loc[:, list(reversed(frame.columns))], model["gene"].tolist())


def test_rank_descending_is_one_based_and_deterministic() -> None:
    scores = np.asarray([[2.0, 1.0, 1.0, -1.0]])
    ranks = module.rank_descending(scores)
    assert ranks.dtype == np.int32
    assert ranks.tolist() == [[1, 2, 2, 4]]

