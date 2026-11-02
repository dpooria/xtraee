import filecmp
import pathlib

import pytest

from xtraee.parser import QCCSDParser, QCISParser


def test_sad_family(tmp_path):
    cwd = pathlib.Path(__file__).parent
    input_ccsd = cwd / "test_qccsd/input2.log"
    input_cis = cwd / "test_qcis/input1.log"
    expected = cwd / "sad_family_cis.txt"
    missing = [
        str(path) for path in (input_ccsd, input_cis, expected) if not path.exists()
    ]
    if missing:
        pytest.skip(f"missing fixtures: {', '.join(missing)}")

    ccsd_parser = QCCSDParser(input_ccsd, 0.2)
    ccsd_parser.process_file()
    ccsd_parser.write_dataset(tmp_path / "first_kid_ccsd.txt")

    cis_parser = QCISParser(input_cis, 0.2)
    cis_parser.process_file()
    cis_parser.write_full(tmp_path / "first_kid_cis.txt")
    cis_parser.write_dataset(tmp_path / "happy_family_cis.txt")
    produced = tmp_path / "new.sad_family_cis.txt"
    cis_parser.write_vs_std(produced, ccsd_parser.irreps)
    assert filecmp.cmp(produced, expected, shallow=False), "Files do not match."


if __name__ == "__main__":
    test_sad_family(pathlib.Path("."))
