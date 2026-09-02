#!/usr/bin/env python3
"""Create a review-only Figure 3A-F montage without inventing absent panels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, JpegImagePlugin

try:
    from figure3_egr1_common import PROJECT_ROOT, TARGET_TF, write_json
except ModuleNotFoundError:
    from scripts.figure3_egr1_common import PROJECT_ROOT, TARGET_TF, write_json


DEFAULT_OUT_DIR = PROJECT_ROOT / "figures/driver/figure3_egr1_preview"


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
    ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


def fit_image(path: Path, width: int, height: int) -> Image.Image:
    with Image.open(path) as source:
        image = source.convert("RGB")
    image.thumbnail((width, height), Image.Resampling.LANCZOS)
    tile = Image.new("RGB", (width, height), "white")
    tile.paste(image, ((width - image.width) // 2, (height - image.height) // 2))
    return tile


def placeholder(width: int, height: int, message: str) -> Image.Image:
    tile = Image.new("RGB", (width, height), "#F7F7F7")
    draw = ImageDraw.Draw(tile)
    draw.rectangle((1, 1, width - 2, height - 2), outline="#A0A0A0", width=3)
    font = load_font(31)
    lines = []
    words = message.split()
    current = ""
    for word in words:
        trial = f"{current} {word}".strip()
        if draw.textlength(trial, font=font) > width - 100 and current:
            lines.append(current)
            current = word
        else:
            current = trial
    if current:
        lines.append(current)
    line_height = 42
    y = (height - line_height * len(lines)) // 2
    for line in lines:
        box = draw.textbbox((0, 0), line, font=font)
        draw.text(((width - (box[2] - box[0])) // 2, y), line, fill="#555555", font=font)
        y += line_height
    return tile


def panel_tile(
    label: str,
    path: Path | None,
    width: int,
    height: int,
    missing_message: str,
) -> tuple[Image.Image, dict]:
    header = 72
    body_height = height - header
    if path is not None and path.exists() and path.stat().st_size:
        body = fit_image(path, width, body_height)
        status = "included"
    else:
        body = placeholder(width, body_height, missing_message)
        status = "fdr_suppressed_or_missing"
    tile = Image.new("RGB", (width, height), "white")
    tile.paste(body, (0, header))
    draw = ImageDraw.Draw(tile)
    draw.text((18, 12), f"Figure 3{label}", font=load_font(40, bold=True), fill="black")
    draw.rectangle((0, 0, width - 1, height - 1), outline="#D0D0D0", width=2)
    return tile, {"panel": label, "path": str(path.resolve()) if path else None, "status": status}


def run(project_root: Path, out_dir: Path) -> dict:
    panel_paths = {
        "A": project_root / "figures/driver/figure3a_stress_transition/figure3a_stress_transition_selection.png",
        "B": project_root / "figures/driver/figure3b_egr1/figure3b_egr1_baseline_umap.png",
        "C": project_root / "figures/driver/figure3c_egr1/figure3c_egr1_inner_product_umap.png",
        "D": project_root / "figures/driver/figure3d_egr1/figure3d_egr1_pseudotime_inner_product_umap.png",
        "E": project_root / "figures/driver/figure3e_egr1/figure3e_egr1_significant_perturbed_genes.png",
        "F": project_root / "figures/driver/figure3f_egr1/figure3f_egr1_pathway_enrichment.png",
    }
    required = ["A", "B", "C", "D"]
    missing_required = [label for label in required if not panel_paths[label].exists()]
    if missing_required:
        raise FileNotFoundError(f"Required montage panels are missing: {missing_required}")

    margin, gap, full_width = 70, 34, 3600
    inner_width = full_width - 2 * margin
    half_width = (inner_width - gap) // 2
    heights = {"A": 1500, "B": 1300, "C": 1300, "D": 1300, "E": 1300, "F": 1400}
    records = []
    tiles = {}
    missing_messages = {
        "E": "Formal panel suppressed when no FDR-significant perturbed genes are available. See Figure 3E report.",
        "F": "Formal panel suppressed when no globally BH FDR-significant pathways are available. See Figure 3F report.",
    }
    for label in ["A", "F"]:
        tiles[label], record = panel_tile(
            label,
            panel_paths[label],
            inner_width,
            heights[label],
            missing_messages.get(label, "Panel unavailable"),
        )
        records.append(record)
    for label in ["B", "C", "D", "E"]:
        tiles[label], record = panel_tile(
            label,
            panel_paths[label],
            half_width,
            heights[label],
            missing_messages.get(label, "Panel unavailable"),
        )
        records.append(record)

    total_height = (
        2 * margin
        + heights["A"]
        + heights["B"]
        + heights["D"]
        + heights["F"]
        + 3 * gap
    )
    canvas = Image.new("RGB", (full_width, total_height), "white")
    y = margin
    canvas.paste(tiles["A"], (margin, y))
    y += heights["A"] + gap
    canvas.paste(tiles["B"], (margin, y))
    canvas.paste(tiles["C"], (margin + half_width + gap, y))
    y += heights["B"] + gap
    canvas.paste(tiles["D"], (margin, y))
    canvas.paste(tiles["E"], (margin + half_width + gap, y))
    y += heights["D"] + gap
    canvas.paste(tiles["F"], (margin, y))

    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "figure3_egr1_a_to_f_preview.png"
    pdf_path = out_dir / "figure3_egr1_a_to_f_preview.pdf"
    canvas.save(png_path, dpi=(300, 300), optimize=True)
    canvas.save(pdf_path, resolution=300)
    report = {
        "module": "Figure 3 EGR1 review montage",
        "target_tf": TARGET_TF,
        "purpose": "Review-only montage; not final SCI assembly",
        "panels": records,
        "outputs": {"png": str(png_path.resolve()), "pdf": str(pdf_path.resolve())},
        "caveat": "Absent FDR-supported E/F panels are represented by explicit review placeholders, never fabricated plots.",
    }
    write_json(report, out_dir / "figure3_egr1_a_to_f_preview_report.json")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    print(json.dumps(run(args.project_root, args.out_dir), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
