import subprocess


def test_compareall():
    # Run the compareall.py script
    result = subprocess.run(
        [
            "xtraee",
            "compareall",
            "../test_qcis/input1.log",
            "../test_qccsd/input2.log",
            "--acc-method=2"
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
