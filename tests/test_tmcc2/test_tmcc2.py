import filecmp
import pytest
from xtraee.parser import TMCC2Parser
from pathlib import Path


@pytest.mark.parametrize(
    "inputfile, expected_output",
    [
        ("ricc2.out", "happy_family_ricc2.out"),
    ],
)
def test_happy_family(inputfile, expected_output, tmp_path):
    test_dir = Path(__file__).parent
    inputfile_path = test_dir / inputfile
    expected_output_path = test_dir / expected_output
    outputfile = tmp_path / ("test_" + expected_output)
    tmcc2 = TMCC2Parser(inputfile_path, 0.0)
    tmcc2.process_file()
    tmcc2.write_dataset(outputfile)
    assert filecmp.cmp(
        outputfile, expected_output_path, shallow=False
    ), f"Files {outputfile} and {expected_output_path} do not match."


if __name__ == "__main__":
    for inputfile in ["ricc2.out"]:
        output = f"happy_family_{inputfile}"
        test_dir = Path(__file__).parent
        inputfile_path = test_dir / inputfile
        outputfile = test_dir / output
        tmcc2 = TMCC2Parser(inputfile_path, 0.0)
        tmcc2.process_file()
        tmcc2.write_dataset(outputfile)
