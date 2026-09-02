from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "scripts" / "figure8_v2_10_plot_mainfigure.R"
EXT = ROOT / "scripts" / "figure8_v2_11_plot_extended_data.R"


def test_r_only_plot_contract() -> None:
    assert MAIN.exists(), "Expected RED failure: main-figure R script is absent"
    assert EXT.exists(), "Expected RED failure: Extended Data R script is absent"
    text = MAIN.read_text(encoding="utf-8") + "\n" + EXT.read_text(encoding="utf-8")
    assert "figure8_v2_theme.R" in text
    assert "plot_layout" in text
    assert "tag_levels = \"a\"" in text
    assert "figure8_v2_save_plot" in text
    forbidden = ("matplotlib", "seaborn", "plotly", "viridis", "rainbow")
    assert not any(token.lower() in text.lower() for token in forbidden)


def test_main_plot_declares_all_seven_panels() -> None:
    text = MAIN.read_text(encoding="utf-8")
    for panel in "abcdefg":
        assert f"panel_{panel}" in text
    assert "aes(label = node)" in text
    assert text.count("wrap_elements(full") >= 3
    assert 'guides = "collect"' not in text
    assert 'top20[1:12]' in text
    assert 'aes(x = malignant, label' not in text
