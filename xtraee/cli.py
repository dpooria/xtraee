import argparse
from xtraee.qc import QCCSDParser, QCISParser


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract data from QChem CCSD output file"
    )
    parser.add_argument(
        "--input-ccsd", type=str, help="Input file for CCSD", default=None
    )
    parser.add_argument(
        "--firstkid-ccsd",
        type=str,
        default="first_kid_ccsd.txt",
        help="Output file for the complete information",
    )
    parser.add_argument(
        "--happyfamily-ccsd",
        type=str,
        default="happy_family_ccsd.txt",
        help="Output file for the happy family",
    )
    parser.add_argument(
        "--input-cis", type=str, help="Input file for CIS", default=None
    )
    parser.add_argument(
        "--firstkid-cis",
        type=str,
        default="first_kid_cis.txt",
        help="Output file for the complete information",
    )
    parser.add_argument(
        "--happyfamily-cis",
        type=str,
        default="happy_family_cis.txt",
        help="Output file for the happy family",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Threshold for filterning out the matching transitions",
    )
    parser.add_argument(
        "--sadfamily",
        type=str,
        default="sad-family.txt",
        help="Output file for the complete information",
    )
    subparsers = parser.add_subparsers(dest="command", default="compare")
    compare_parser = subparsers.add_parser("compare", help="Compare two different states")
    compare_parser.add_argument("state1", type=str, help="The first state")
    compare_parser.add_argument("state2", type=str, help="The second state")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.input_ccsd is None and args.input_cis is None:
        print(
            "Please provide input file for CCSD or CIS or both (type --help for help)"
        )
        return 1
    if args.input_ccsd is not None:
        qccsd = QCCSDParser(
            args.threshold, args.input_ccsd, args.firstkid_ccsd, args.happyfamily_ccsd
        )
        qccsd.process_file()
        qccsd.write_first_kid()
        qccsd.write_happy_family()
        if args.command == "compare":
            qccsd.compare_states(args.state1, args.state2)
    if args.input_cis is not None:
        qcis = QCISParser(
            args.threshold, args.input_cis, args.firstkid_cis, args.happyfamily_cis
        )
        qcis.process_file()
        qcis.write_first_kid()
        qcis.write_happy_family()
        if args.input_ccsd is not None:
            qcis.write_sad_family(qccsd.irreps_dict, args.sadfamily)
    return 0
