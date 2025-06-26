import filecmp
import sys
from pathlib import Path

import pytest
from xtraee.parser import Parser


@pytest.mark.parametrize(
    "inputfile, expected_output",
    [
        ("input1_qcc2.log", "happy_family1.txt"),
    ],
)
def test_happy_family(inputfile, expected_output, tmp_path):
    test_dir = Path(__file__).parent
    inputfile_path = test_dir / inputfile
    expected_output_path = test_dir / expected_output
    outputfile = tmp_path / ("test_" + expected_output)
    qccsd = Parser(inputfile_path, 0.2)
    qccsd.process_file()
    qccsd.write_dataset(outputfile)
    assert filecmp.cmp(
        outputfile, expected_output_path, shallow=False
    ), f"Files {outputfile} and {expected_output_path} do not match."


def usage_and_exit():
    print("Usage: %s {test,update-test}" % sys.argv[0])
    exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        usage_and_exit()
    mode = sys.argv[1]
    if mode not in ["test", "update-test"]:
        usage_and_exit()
    for i in [1]:
        inputfile = f"input{i}_qcc2.log"
        expected_output = f"happy_family{i}.txt"
        test_dir = Path(__file__).parent
        inputfile_path = test_dir / inputfile
        if mode == "test":
            test_happy_family(inputfile, expected_output, test_dir)
        elif mode == "update-test":
            outputfile = test_dir / expected_output
            qccsd = Parser(inputfile_path, 0.2)
            qccsd.process_file()
            qccsd.write_dataset(outputfile)
    exit(0)
