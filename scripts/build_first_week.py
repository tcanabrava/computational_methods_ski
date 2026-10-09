"""Generate the first-week figures and compile Tomaz's report.

Run from the repository root: python -m scripts.build_first_week
Requires matplotlib and typst.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, nargs="?",
                        default=ROOT / "tomaz-first-week.pdf",
                        help="Output PDF (default: tomaz-first-week.pdf)")
    args = parser.parse_args()
    typst = shutil.which("typst")
    if typst is None:
        parser.error("typst is not installed or is missing from PATH")

    output = args.output.resolve()
    try:
        for module in ("scripts.tomaz_aerodynamics", "scripts.aerodynamic_report",
                       "scripts.coefficient_precision_report"):
            subprocess.run([sys.executable, "-m", module], cwd=ROOT, check=True)
        subprocess.run([typst, "compile", "tomaz-first-week.typ", str(output)],
                       cwd=ROOT, check=True)
    except subprocess.CalledProcessError as error:
        raise SystemExit(error.returncode)
    print(f"PDF: {output}")


if __name__ == "__main__":
    main()
