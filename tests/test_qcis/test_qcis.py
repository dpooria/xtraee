import filecmp
import pytest
from xtraee.qc import QCISParser
from pathlib import Path


@pytest.mark.xfail(raises=NotImplementedError)
@pytest.mark.parametrize(
    "inputfile, expected_output",
    [
        ("input1.log", "happy_family1.txt"),
    ],
)
def test_happy_family(inputfile, expected_output, tmp_path):
    test_dir = Path(__file__).parent
    inputfile_path = test_dir / inputfile
    expected_output_path = test_dir / expected_output
    outputfile = tmp_path / ("test_" + expected_output)
    qccsd = QCISParser(0.2, inputfile_path, ".tmp", outputfile)
    qccsd.process_file()
    qccsd.write_happy_family()
    assert filecmp.cmp(
        outputfile, expected_output_path, shallow=False
    ), f"Files {outputfile} and {expected_output_path} do not match."


if __name__ == "__main__":
    for i in range(1, 2):
        inputfile = f"input{i}.log"
        output = f"happy_family{i}.txt"
        test_dir = Path(__file__).parent
        inputfile_path = test_dir / inputfile
        outputfile = test_dir / output
        qccsd = QCISParser(0.2, inputfile_path, ".tmp", outputfile)
        qccsd.process_file()
        qccsd.write_happy_family()
