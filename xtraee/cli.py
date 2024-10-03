import argparse
from xtraee.parser.ccsd import CCSDParser


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract data from QChem CCSD output file"
    )
    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument(
        "--first-kid",
        "-fk",
        type=str,
        default="first_kid.txt",
        help="Output file for the complete information",
    )
    parser.add_argument(
        "--happy-family",
        "-hf",
        type=str,
        default="happy_family.txt",
        help="Output file for the happy family",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    qccsd = CCSDParser(args.input, args.first_kid, args.happy_family)
    qccsd.process_file()
    qccsd.write_first_kid()
    qccsd.write_happy_family()
    return 0
