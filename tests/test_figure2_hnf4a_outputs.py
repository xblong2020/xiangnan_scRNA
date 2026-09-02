from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def test_figure2b_baseline_is_exactly_equivalent_to_sox4() -> None:
    audit = read_tsv(ROOT / "metadata/driver/figure2b_hnf4a/figure2b_baseline_equivalence_audit.tsv")
    assert audit["values_exactly_equal"].astype(bool).all()
    assert (pd.to_numeric(audit["max_absolute_numeric_difference"]) == 0).all()


def test_figure2c_matched_vectors_are_hnf4a_not_sox4() -> None:
    matched = read_tsv(ROOT / "metadata/driver/figure2c_hnf4a/figure2c_hnf4a_matched_cells.tsv.gz")
    shifts = read_tsv(ROOT / "metadata/driver/celloracle_module6_8_cell_shift_summary.tsv.gz")
    hnf4a = shifts.loc[shifts["tf"].eq("HNF4A"), ["cell_id", "delta_embedding_1", "delta_embedding_2"]]
    sox4 = shifts.loc[shifts["tf"].eq("SOX4"), ["cell_id", "delta_embedding_1", "delta_embedding_2"]]
    actual = matched[["cell_id", "delta_embedding_1", "delta_embedding_2"]]
    h = actual.merge(hnf4a, on="cell_id", suffixes=("_actual", "_expected"), validate="one_to_one")
    s = actual.merge(sox4, on="cell_id", suffixes=("_actual", "_sox4"), validate="one_to_one")
    assert np.allclose(h[["delta_embedding_1_actual", "delta_embedding_2_actual"]],
                       h[["delta_embedding_1_expected", "delta_embedding_2_expected"]])
    assert not np.allclose(s[["delta_embedding_1_actual", "delta_embedding_2_actual"]],
                           s[["delta_embedding_1_sox4", "delta_embedding_2_sox4"]])


def test_figure2e_contains_only_fdr_significant_non_target_genes() -> None:
    dat = read_tsv(ROOT / "metadata/driver/figure2e_hnf4a/figure2e_hnf4a_significant_perturbed_genes.tsv")
    assert not dat["gene"].eq("HNF4A").any()
    assert (pd.to_numeric(dat["p.adj"]) < 0.05).all()
    assert 1 <= len(dat) <= 20


def test_figure2f_never_promotes_nominal_pathways() -> None:
    report = read_json(ROOT / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_report.json")
    plot_data = read_tsv(ROOT / "metadata/driver/figure2f_hnf4a/figure2f_hnf4a_plot_data.tsv")
    assert report["target_tf"] == "HNF4A"
    assert len(plot_data) == report["n_plotted"]
    if len(plot_data):
        assert (pd.to_numeric(plot_data["p.adjust"]) < 0.05).all()
    else:
        assert report["n_significant_pathways"] == 0
        assert report["figure_generated"] is False


def test_validation_report_is_complete_and_sox4_is_preserved() -> None:
    report = read_json(ROOT / "metadata/driver/figure2_hnf4a_b_to_f_validation_report.json")
    assert report["target_tf"] == "HNF4A"
    assert report["all_required_checks_passed"] is True
    assert report["n_failed"] == 0
    assert report["n_sox4_hashes_checked"] > 0
