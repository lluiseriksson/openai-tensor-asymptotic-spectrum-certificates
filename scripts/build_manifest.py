"""Generate or verify MANIFEST.sha256 for Paper 68."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "MANIFEST.sha256"

TRACKED_GLOBS = [
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "lean/lean-toolchain",
    "lean/lakefile.toml",
    "lean/TensorSpectrumCert.lean",
    "src/__init__.py",
    "src/tensor_spectrum_verifier.py",
    "scripts/generate_certificates.py",
    "scripts/generate_figures.py",
    "scripts/build_manifest.py",
    "tests/test_verifier.py",
    "certificates/*.json",
    "paper/main.tex",
    "paper/references.bib",
    "paper/main.pdf",
    "paper/figures/*.pdf",
    "paper/figures/*.png",
]


def collect_files() -> list[Path]:
    found: set[Path] = set()
    for pat in TRACKED_GLOBS:
        for p in ROOT.glob(pat):
            if p.is_file():
                found.add(p)
    return sorted(found, key=lambda x: x.relative_to(ROOT).as_posix())


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        lines = [ln.strip() for ln in MANIFEST.read_text(encoding="utf-8").splitlines() if ln.strip()]
        ok = 0
        for ln in lines:
            expected, rel = ln.split("  ", 1)
            actual = sha256_file(ROOT / rel)
            if actual != expected:
                print(f"[MISMATCH] {rel}: {actual} != {expected}", file=sys.stderr)
                return 1
            ok += 1
        print(f"MANIFEST.sha256 OK: {ok}/{len(lines)} files match.")
        return 0

    files = collect_files()
    out_lines = [f"{sha256_file(p)}  {p.relative_to(ROOT).as_posix()}" for p in files]
    MANIFEST.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"Wrote MANIFEST.sha256 with {len(out_lines)} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
