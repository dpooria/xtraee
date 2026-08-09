import sys
import filecmp
import pytest
from xtraee.parser import Parser
from pathlib import Path


@pytest.mark.parametrize(
    "inputfile, expected_output",
    [
        ("input1.log", "output1.txt"),
    ],
)
def test_qcisd(inputfile, expected_output):
    test_dir = Path(__file__).parent
    inputfile_path = test_dir / inputfile
    expected_output_path = test_dir / expected_output
    outputfile = test_dir / ("test_" + expected_output)
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
    for inputfile, output in [("input1.log", "output1.txt")]:
        if mode == "test":
            test_qcisd(inputfile, output)
        else:
            test_dir = Path(__file__).parent
            inputfile_path = test_dir / inputfile
            output_path = test_dir / output
            qccsd = Parser(inputfile_path, 0.2)
            qccsd.process_file()
            breakpoint()
            qccsd.write_dataset(output_path)
    exit(0)
