import filecmp
import pytest
from xtraee.parser import QCISParser
from pathlib import Path


@pytest.mark.parametrize(
    "inputfile, expected_output",
    [
        ("input1.log", "output1.txt"),
    ],
)
def test_qcis(inputfile, expected_output):
    test_dir = Path(__file__).parent
    inputfile_path = test_dir / inputfile
    expected_output_path = test_dir / expected_output
    outputfile = test_dir / ("test_" + expected_output)
    qccsd = QCISParser(inputfile_path, 0.2)
    qccsd.process_file()
    qccsd.write_dataset(outputfile)
    assert filecmp.cmp(
        outputfile, expected_output_path, shallow=False
    ), f"Files {outputfile} and {expected_output_path} do not match."


# if __name__ == "__main__":
#     for i in range(1, 2):
#         inputfile = f"input{i}.log"
#         output = f"happy_family{i}.txt"
#         test_dir = Path(__file__).parent
#         inputfile_path = test_dir / inputfile
#         outputfile = test_dir / output
#         qccsd = QCISParser(inputfile_path, 0.2)
#         qccsd.process_file()
#         qccsd.write_dataset(outputfile)
if __name__ == "__main__":
    exit(test_qcis("input1.log", "happy_family1.txt"))
