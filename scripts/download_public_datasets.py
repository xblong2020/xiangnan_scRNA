from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import requests


ROOT = Path(__file__).resolve().parents[1]
GEO_ROOT = ROOT / "data" / "public" / "geo"
FIGSHARE_ROOT = ROOT / "data" / "public" / "figshare" / "HCC_atlas"
META_ROOT = ROOT / "metadata"


@dataclass(frozen=True)
class DownloadItem:
    dataset: str
    category: str
    filename: str
    url: str
    expected_size: int
    out_dir: Path

    @property
    def path(self) -> Path:
        return self.out_dir / self.filename


def geo_url(gse: str, filename: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/{filename}"


def build_manifest() -> list[DownloadItem]:
    geo_files: dict[str, list[tuple[str, int]]] = {
        "GSE202379": [
            ("GSE202379_RAW.tar", 595_322_880),
            ("GSE202379_SeuratObject_AllCells.rds.gz", 2_840_561_837),
            ("filelist.txt", 8_962),
        ],
        "GSE174748": [
            ("GSE174748_RAW.tar", 415_252_480),
            ("GSE174748_hl_nuclei.rds.gz", 207_972_747),
            ("filelist.txt", 954),
        ],
        "GSE185477": [
            ("GSE185477_Final_Metadata.txt.gz", 5_340_345),
            ("GSE185477_GSM3178784_C41_SC_raw_counts.zip", 39_286_693),
            ("GSE185477_RAW.tar", 1_274_798_080),
            ("filelist.txt", 883),
        ],
        "GSE212046": [
            ("GSE212046_RAW.tar", 260_433_920),
            ("filelist.txt", 402),
        ],
        "GSE149614": [
            ("GSE149614_HCC.metadata.updated.txt.gz", 489_882),
            ("GSE149614_HCC.scRNAseq.S71915.count.txt.gz", 165_349_783),
            ("GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz", 1_330_022_464),
        ],
        "GSE151530": [
            ("GSE151530_Info.txt.gz", 330_527),
            ("GSE151530_barcodes.tsv.gz", 256_634),
            ("GSE151530_genes.tsv.gz", 146_052),
            ("GSE151530_matrix.mtx.gz", 304_902_525),
        ],
    }

    items: list[DownloadItem] = []
    for gse, files in geo_files.items():
        out_dir = GEO_ROOT / gse
        for filename, expected_size in files:
            items.append(
                DownloadItem(
                    dataset=gse,
                    category="GEO",
                    filename=filename,
                    url=geo_url(gse, filename),
                    expected_size=expected_size,
                    out_dir=out_dir,
                )
            )

    figshare_files = [
        ("HCC_atlas_myeloid_release.rds", "https://ndownloader.figshare.com/files/41624958", 445_374_991),
        ("HCC_atlas_Bcell_release.rds", "https://ndownloader.figshare.com/files/41624988", 84_692_238),
        ("HCC_atlas_endothelium_release.rds", "https://ndownloader.figshare.com/files/41624991", 94_343_876),
        ("HCC_atlas_fibroblast_release.rds", "https://ndownloader.figshare.com/files/41624994", 47_610_570),
        ("HCC_atlas_TNK_release.rds", "https://ndownloader.figshare.com/files/41625000", 823_737_307),
        ("HCC_atlas_metadata_batch_effect.csv", "https://ndownloader.figshare.com/files/41645349", 98_696_380),
        ("HCC_atlas_all_release.rds", "https://ndownloader.figshare.com/files/41655405", 2_323_723_464),
    ]
    for filename, url, expected_size in figshare_files:
        items.append(
            DownloadItem(
                dataset="HCC_atlas",
                category="figshare",
                filename=filename,
                url=url,
                expected_size=expected_size,
                out_dir=FIGSHARE_ROOT,
            )
        )
    return items


def fmt_size(size: int) -> str:
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024:
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} PB"


def file_size(path: Path) -> int:
    try:
        return path.stat().st_size
    except FileNotFoundError:
        return 0


def is_complete(item: DownloadItem) -> bool:
    size = file_size(item.path)
    return size == item.expected_size if item.expected_size else size > 0


def remote_size(url: str) -> int:
    try:
        response = requests.head(url, allow_redirects=True, timeout=30)
        response.raise_for_status()
        value = response.headers.get("content-length")
        return int(value) if value else 0
    except requests.RequestException:
        return 0


def download(item: DownloadItem) -> None:
    item.out_dir.mkdir(parents=True, exist_ok=True)
    expected_size = item.expected_size if item.filename == "filelist.txt" else remote_size(item.url) or item.expected_size
    item_for_check = DownloadItem(
        item.dataset, item.category, item.filename, item.url, expected_size, item.out_dir
    )
    if is_complete(item_for_check):
        print(f"SKIP complete: {item.path} ({fmt_size(expected_size)})", flush=True)
        return

    part = item.path.with_suffix(item.path.suffix + ".part")
    if item.path.exists() and not part.exists():
        item.path.rename(part)
    if part.exists() and file_size(part) == expected_size:
        part.replace(item.path)
        print(f"DONE {item.path} ({fmt_size(expected_size)})", flush=True)
        return

    existing = file_size(part)
    headers = {"User-Agent": "Mozilla/5.0"}
    mode = "wb"
    if existing:
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"

    print(
        f"DOWNLOAD {item.dataset}: {item.filename} "
        f"({fmt_size(expected_size)}, resume {fmt_size(existing)})",
        flush=True,
    )
    with requests.get(item.url, headers=headers, stream=True, timeout=(30, 120)) as response:
        if existing and response.status_code != 206:
            mode = "wb"
            existing = 0
        response.raise_for_status()
        downloaded = existing
        last_report = 0
        with part.open(mode + ("" if "b" in mode else "b")) as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if not chunk:
                    continue
                handle.write(chunk)
                downloaded += len(chunk)
                if downloaded - last_report >= 128 * 1024 * 1024:
                    pct = (downloaded / expected_size * 100) if expected_size else 0
                    print(
                        f"  {item.filename}: {fmt_size(downloaded)} / "
                        f"{fmt_size(expected_size)} ({pct:.1f}%)",
                        flush=True,
                    )
                    last_report = downloaded

    actual = file_size(part)
    if expected_size and actual != expected_size:
        raise RuntimeError(f"Size mismatch for {item.filename}: got {actual}, expected {expected_size}")
    part.replace(item.path)
    print(f"DONE {item.path} ({fmt_size(actual)})", flush=True)


def write_manifest(items: list[DownloadItem], name: str) -> None:
    META_ROOT.mkdir(parents=True, exist_ok=True)
    path = META_ROOT / name
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["category", "dataset", "filename", "expected_size", "local_size", "status", "url", "local_path"])
        for item in items:
            local_size = file_size(item.path)
            status = "complete" if is_complete(item) else "missing_or_partial"
            writer.writerow(
                [
                    item.category,
                    item.dataset,
                    item.filename,
                    item.expected_size,
                    local_size,
                    status,
                    item.url,
                    str(item.path),
                ]
            )
    print(f"WROTE {path}", flush=True)


def main() -> int:
    items = build_manifest()
    write_manifest(items, "public_dataset_download_manifest.before.tsv")
    for item in items:
        download(item)
    write_manifest(items, "public_dataset_download_manifest.after.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
