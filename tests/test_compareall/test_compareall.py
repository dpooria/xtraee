import subprocess
from pathlib import Path


def test_compareall():
    # Run the compareall.py script
    cwd = Path(__file__).parent
    result = subprocess.run(
        [
            "xtraee",
            "compareall",
            cwd / "../test_qcis/input1.log",
            cwd / "../test_qccsd/input2.log",
            "--acc-method=inner-prod",
            f"--outdir={cwd.as_posix()}",
        ],
        capture_output=True,
        text=True,
    )

    # Check if the script ran successfully
    assert result.returncode == 0, f"Script failed with error: {result.stderr}"

    # Optionally, print the output for debugging
    print(result.stdout)


if __name__ == "__main__":
    test_compareall()
