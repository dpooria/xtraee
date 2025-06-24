import pathlib
import filecmp

from xtraee.parser import QCCSDParser, QCISParser


def test_sad_family():
    cwd = pathlib.Path(__file__).parent
    input_ccsd = cwd / "test_qccsd/input2.log"
    ccsd_parser = QCCSDParser(input_ccsd, 0.2)
    ccsd_parser.process_file()
    ccsd_parser.write_dataset(cwd / "first_kid_ccsd.txt")
    input_cis = cwd / "test_qcis/input1.log"
    cis_parser = QCISParser(input_cis, 0.2)
    cis_parser.process_file()
    cis_parser.write_full(cwd / "first_kid_cis.txt")
    cis_parser.write_dataset(cwd / "happy_family_cis.txt")
    cis_parser.write_vs_std(cwd / "new.sad_family_cis.txt", ccsd_parser.irreps_dict)
    assert filecmp.cmp(
        cwd / "new.sad_family_cis.txt",
        cwd / "sad_family_cis.txt",
        shallow=False,
    ), "Files do not match."


if __name__ == "__main__":
    test_sad_family()
