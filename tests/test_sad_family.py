import pathlib
from xtraee.qc import QCCSDParser, QCISParser


def test_sad_family():
    cwd = pathlib.Path(__file__).parent
    input_ccsd = cwd / "test_qccsd/input2.log"
    ccsd_parser = QCCSDParser(
        0.2, input_ccsd, cwd / "first_kid_ccsd.txt", cwd / "happy_family_ccsd.txt"
    )
    ccsd_parser.process_file()
    ccsd_parser.write_happy_family()
    input_cis = cwd / "test_qcis/input1.log"
    cis_parser = QCISParser(
        0.2, input_cis, cwd / "first_kid_cis.txt", cwd / "happy_family_cis.txt"
    )
    cis_parser.process_file()
    cis_parser.write_first_kid()
    cis_parser.write_happy_family()
    cis_parser.write_sad_family(ccsd_parser.irreps_dict, "sad_family_cis.txt")

    # ccsd_parser.
