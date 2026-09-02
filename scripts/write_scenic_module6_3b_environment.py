from __future__ import annotations

import argparse
import importlib.metadata as metadata
import platform
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_text(command: list[str]) -> str:
    try:
        result = subprocess.run(command, capture_output=True, timeout=30, check=False)
        chunks = []
        for raw in [result.stdout or b"", result.stderr or b""]:
            text = raw.decode("utf-8", "replace")
            if text.count("\x00") > 2:
                text = raw.decode("utf-16-le", "replace")
            chunks.append(text.replace("\x00", ""))
        return "".join(chunks).strip()
    except Exception as exc:
        return f"{type(exc).__name__}: {exc}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Write the Module 6.3b runtime environment summary.")
    parser.add_argument("--output", type=Path, default=ROOT / "metadata/driver/scenic_module6_3b_environment.txt")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    packages = ["numpy", "pandas", "scipy", "pyarrow", "pyscenic", "arboreto", "ctxcore", "dask", "distributed", "loompy", "anndata", "scanpy", "cellrank"]
    lines = [
        "Module 6.3b environment",
        f"python={sys.version}",
        f"platform={platform.platform()}",
        f"processor={platform.processor()}",
        f"wsl_status={run_text(['wsl', '--status'])}",
        f"wsl_distributions={run_text(['wsl', '-l', '-v'])}",
        f"docker_version={run_text(['docker', '--version'])}",
    ]
    for package in packages:
        try:
            value = metadata.version(package)
        except metadata.PackageNotFoundError:
            value = "not_installed"
        lines.append(f"{package}={value}")
    lines.extend(
        [
            "formal_grn_method=GRNBoost2",
            "formal_grn_seed=777",
            "formal_expression_cells=9512",
            "formal_expression_genes_after_all_zero_cleanup=11923",
            "formal_expressed_tfs=1767",
            "formal_grn_workers=32",
            "formal_ctx_workers=8",
        ]
    )
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
