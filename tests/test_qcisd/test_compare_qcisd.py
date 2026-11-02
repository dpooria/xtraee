import subprocess
import sys
from pathlib import Path

import pytest


def test_qcisd_compare():
    """``compareall`` between an EOM-CCSD and a CISD log, through ``python -m``."""
    cwd = Path(__file__).parent
    ccsd = cwd / "../test_qccsd/input2.log"
    cisd = cwd / "input_CISD1.log"
    missing = [str(path) for path in (ccsd, cisd) if not path.exists()]
    if missing:
        pytest.skip(f"missing fixtures: {', '.join(missing)}")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "xtraee",
            "compareall",
            str(ccsd),
            str(cisd),
            "--acc-method=inner-prod",
            f"--outdir={cwd.as_posix()}",
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"Script failed with error: {result.stderr}"
    print(result.stdout)


if __name__ == "__main__":
    test_qcisd_compare()
