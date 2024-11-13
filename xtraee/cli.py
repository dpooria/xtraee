import argparse
import logging

from rich import print

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
    parser.add_argument("--debug", action="store_true", help="Debug mode")
    subparsers = parser.add_subparsers(dest="command")
    compare_parser = subparsers.add_parser(
        "compare", help="Compare two different states"
    )
    compare_parser.add_argument(
        "state1",
        type=str,
        help="The first state with format {method}:{excitation}-{state-id}/{irrep-id}, e.g. CIS:singlet-1/A",
    )
    compare_parser.add_argument("state2", type=str, help="The second state")
    compare_parser.add_argument("--acc-method", type=str, default="1")
    compareall_parser = subparsers.add_parser(
        "compareall", help="Compare all different states. "
    )
    compareall_parser.add_argument(
        "--output-ccsd",
        type=str,
        default="compare_ccsd.csv",
        help="output file format for comparison of all of the CCCSD states",
    )
    compareall_parser.add_argument(
        "--output-cis",
        type=str,
        default="compare_cis.csv",
        help="output file format for comparison of all of the CIS states",
    )
    compareall_parser.add_argument(
        "--output-mix",
        type=str,
        default="compare_ccsd_vs_cis.csv",
        help="output file format for comparison of all of the CIS and CCSD states",
    )
    compareall_parser.add_argument("--acc-method", type=str, default="1")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
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
    else:
        qccsd = None

    if args.input_cis is not None:
        qcis = QCISParser(
            args.threshold, args.input_cis, args.firstkid_cis, args.happyfamily_cis
        )
        qcis.process_file()
        qcis.write_first_kid()
        qcis.write_happy_family()
        if args.input_ccsd is not None:
            qcis.write_sad_family(qccsd.irreps_dict, args.sadfamily)
    else:
        qcis = None
    if args.command == "compare":
        state1 = args.state1
        state2 = args.state2
        try:
            method1, state1 = state1.split(":")
            method2, state2 = state2.split(":")
        except Exception:
            print(
                f"Could not parse states {state1} and {state2} "
                "the format required is {method}:{excitation}-{state-id}/{irrep-id}, e.g. CIS:singlet-1/A",
            )
            return 1
        state1_trblock = None
        state2_trblock = None
        if qccsd is not None:
            for irrep in qccsd.irreps_dict.values():
                if method1.lower() == "ccsd" and state1 in irrep.transitions_dict:
                    state1_trblock = irrep.transitions_dict[state1]
                if method2.lower() == "ccsd" and state2 in irrep.transitions_dict:
                    state2_trblock = irrep.transitions_dict[state2]
        if qcis is not None:
            for irrep in qcis.irreps_dict.values():
                if method1.lower() == "cis" and state1 in irrep.transitions_dict:
                    state1_trblock = irrep.transitions_dict[state1]
                if method2.lower() == "cis" and state2 in irrep.transitions_dict:
                    state2_trblock = irrep.transitions_dict[state2]
        if state1_trblock is None or state2_trblock is None:
            print("Couldn't find the states, sorry :(")
            return 1
        print(f"Comparing {state1} and {state2}")
        print("accuracy  | error | fraction matched")
        print(
            "|\t".join(
                map(str, state1_trblock.compare(state2_trblock, args.acc_method))
            )
        )
    elif args.command == "compareall":
        # ccsd
        if qccsd is not None:
            data = qccsd.compare_all(args.acc_method)
            for key, df in data.items():
                df.to_csv(f"{key}_{args.output_ccsd}")
        if qcis is not None:
            data = qcis.compare_all(args.acc_method)
            for key, df in data.items():
                df.to_csv(f"{key}_{args.output_cis}")

        if qccsd is not None and qcis is not None:
            data = qcis.compare_eomee(qccsd.irreps_dict, args.acc_method)
            for key, df in data.items():
                df.to_csv(f"{key}_{args.output_mix}")

    return 0
