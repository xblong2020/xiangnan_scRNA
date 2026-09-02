from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def normalize_value(value: object) -> object:
    if isinstance(value, dict):
        return {key: normalize_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_value(item) for item in value]
    if isinstance(value, str) and value.startswith(str(ROOT)):
        try:
            return str(Path(value).resolve().relative_to(ROOT.resolve())).replace("\\", "/")
        except ValueError:
            return value
    return value


def main() -> None:
    json_paths = [
        *sorted((ROOT / "metadata/driver/scenic_module6_3b").glob("*.json")),
        ROOT / "metadata/driver/scenic_resources_v10/resource_validation.json",
        ROOT / "metadata/driver/driver_module6_3b_canonical_scenic_status.json",
    ]
    for path in json_paths:
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        path.write_text(json.dumps(normalize_value(data), indent=2), encoding="utf-8")
    root_text = str(ROOT)
    root_text_forward = root_text.replace("\\", "/")
    for path in sorted((ROOT / "reports").glob("module6_3b*.md")):
        text = path.read_text(encoding="utf-8")
        text = text.replace(root_text + "\\", "").replace(root_text_forward + "/", "")
        text = text.replace("\\", "/")
        path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
