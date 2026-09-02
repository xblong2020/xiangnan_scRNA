#!/usr/bin/env python3
"""Create a review-only HNF4A Figure 2B-F montage from generated PNG panels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def run(root: Path, out_dir: Path) -> dict:
    candidates = [
        ("B", root / "figures/driver/figure2b_hnf4a/figure2b_hnf4a_baseline_umap.png"),
        ("C", root / "figures/driver/figure2c_hnf4a/figure2c_hnf4a_inner_product_umap.png"),
        ("D", root / "figures/driver/figure2d_hnf4a/figure2d_hnf4a_pseudotime_inner_product_umap.png"),
        ("E", root / "figures/driver/figure2e_hnf4a/figure2e_hnf4a_significant_perturbed_genes.png"),
        ("F", root / "figures/driver/figure2f_hnf4a/figure2f_hnf4a_pathway_enrichment.png"),
    ]
    available = [(label, path) for label, path in candidates if path.exists()]
    if len(available) < 4:
        raise ValueError(f"At least panels B-E are required for montage; found {len(available)}")
    tile_w, tile_h, margin = 1700, 1300, 80
    cols = 2
    rows = (len(available) + cols - 1) // cols
    canvas = Image.new("RGB", (cols * tile_w + (cols + 1) * margin,
                               rows * tile_h + (rows + 1) * margin), "white")
    for i, (label, path) in enumerate(available):
        with Image.open(path) as image:
            image = image.convert("RGB")
            image.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
            x = margin + (i % cols) * (tile_w + margin)
            y = margin + (i // cols) * (tile_h + margin)
            canvas.paste(image, (x, y))
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "figure2_hnf4a_b_to_f_preview.png"
    pdf = out_dir / "figure2_hnf4a_b_to_f_preview.pdf"
    canvas.save(png, dpi=(300, 300))
    canvas.save(pdf, resolution=300)
    report = {"target_tf": "HNF4A", "purpose": "review-only montage, not final assembly",
              "panels_included": [x[0] for x in available],
              "png": str(png.resolve()), "pdf": str(pdf.resolve())}
    (out_dir / "figure2_hnf4a_b_to_f_preview_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "figures/driver/figure2_hnf4a_preview")
    args = parser.parse_args()
    print(json.dumps(run(args.project_root, args.out_dir), ensure_ascii=False))


if __name__ == "__main__":
    main()
