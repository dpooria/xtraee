import pytest

from pathlib import Path

from xtraee.cli import main


def test_compareall():
    """``compareall`` on a CIS and an EOM-CCSD log writes both matrices."""
    cwd = Path(__file__).parent
    cis = cwd / "../test_qcis/input1.log"
    ccsd = cwd / "../test_qccsd/input2.log"
    missing = [str(path) for path in (cis, ccsd) if not path.exists()]
    if missing:
        pytest.skip(f"missing fixtures: {', '.join(missing)}")

    outdir = cwd / "_compare_out"
    exit_code = main(
        [
            "compareall",
            str(cis),
            str(ccsd),
            "--acc-method=inner-prod",
            f"--outdir={outdir}",
        ]
    )
    try:
        assert exit_code == 0
        produced = sorted(path.name for path in outdir.glob("*.csv"))
        assert any(name.startswith("singlet-") for name in produced)
        assert any(name.startswith("CIS_") for name in produced)
    finally:
        for path in outdir.glob("*"):
            path.unlink()
        if outdir.exists():
            outdir.rmdir()


if __name__ == "__main__":
    test_compareall()
