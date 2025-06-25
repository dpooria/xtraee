#!/usr/bin/env python3
import sys
import argparse
import logging
from contextlib import chdir
from pathlib import Path

import pandas as pd
from rich import print

from xtraee.parser import BaseParser, Parser


def extract_and_write(args) -> dict[str, BaseParser]:
    """Process inputs and write full/overview/vsstd outputs."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    parsers: dict[str, BaseParser] = {}
    for infile in args.input:
        p = Parser(infile, args.threshold)
        p.process_file()
        parsers[p.name] = p

    with chdir(outdir):
        for name, p in parsers.items():
            if args.full:
                p.write_full(f"{args.out_fdata}_{name}.txt")
            p.write_dataset(f"{args.out_data}_{name}.txt")

    ref = parsers.get(args.ref)
    if ref:
        with chdir(outdir):
            for name, p in parsers.items():
                if name != args.ref:
                    p.write_vs_std(f"{args.out_vsstd}_{name}.txt", ref.irreps_dict)

    return parsers


def do_compare(args, parsers: dict[str, BaseParser]) -> int:
    """Compare two individual states."""
    try:
        m1, s1 = args.state1.split(":", 1)
        m2, s2 = args.state2.split(":", 1)
        m1 = "EOM-CCSD" if m1.upper() == "CCSD" else m1
        m2 = "EOM-CCSD" if m2.upper() == "CCSD" else m2
    except ValueError:
        logging.error(
            "Could not parse %r and %r; format is "
            "{method}:{excitation}-{state-id}/{irrep-id}, e.g. CIS:singlet-1/A",
            args.state1,
            args.state2,
        )
        return 1

    p1, p2 = parsers.get(m1), parsers.get(m2)
    if not p1 or not p2:
        logging.error("Parser for %s or %s not found!", m1, m2)
        return 1

    tb1 = p1.get_transition_block(s1)
    tb2 = p2.get_transition_block(s2)
    if tb1 is None or tb2 is None:
        logging.error("State %r or %r not found!", s1, s2)
        return 1

    print(f"Comparing {s1} and {s2}")
    print("accuracy  | error | fraction matched")
    print("|\t".join(map(str, tb1.compare(tb2, args.acc_method))))
    return 0


def do_compare_all(args, parsers: dict[str, BaseParser]) -> None:
    """Compare every state (and vs CCSD if available)."""
    outdir = Path(args.outdir)
    refname = args.ref
    ref = parsers.get(refname)

    for name, p in parsers.items():
        data = p.compare_all(args.acc_method)
        for key, df in data.items():
            with chdir(outdir):
                df.to_csv(f"{key}_{name}_{args.output}.csv")

        if ref and name != refname:
            data_vs = p.compare_std(ref.irreps_dict, args.acc_method)
            for key, df in data_vs.items():
                with chdir(outdir):
                    df.to_csv(f"{key}_{args.output}.csv")


def do_descriptors(args, parsers: dict[str, BaseParser]) -> int:
    """Gather descriptors for EOM-CCSD and dump to CSV."""
    eom = parsers.get("EOM-CCSD")
    if not eom:
        logging.error("Descriptors are only implemented for EOM-CCSD")
        return 1

    df = pd.DataFrame(eom.gather_descriptors())
    df.to_csv(args.descriptors_file, index=False)
    return 0


def parse_args() -> argparse.Namespace:
    # If the first non-option token isn't a known sub-command, inject "extract"
    commands = ("extract", "compare", "compareall", "descriptors")
    if not any(arg in ("-h", "--help") for arg in sys.argv[1:]):
        for i, tok in enumerate(sys.argv[1:], start=1):
            if tok.startswith("-"):
                continue
            if tok in commands:
                break
            # insert default command at this position
            sys.argv.insert(i, "extract")
            break
        else:
            # no non-option token at all
            sys.argv.append("extract")

    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument("input", nargs="+", help="Input file(s)")
    parent.add_argument(
        "--full", action="store_true", help="Write the full data to file"
    )
    parent.add_argument(
        "--out-fdata",
        type=str,
        default="fulldata",
        help="Filename format for the full data",
    )
    parent.add_argument("--outdir", type=str, default=".", help="Output directory")
    parent.add_argument(
        "--out-data",
        type=str,
        default="data",
        help="Filename format for the overview data",
    )
    parent.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Threshold for filtering out matching transitions",
    )
    parent.add_argument(
        "--out-vsstd",
        type=str,
        default="vsstd",
        help="Filename format for the comparison with the reference calculation (usually EOM-CCSD)",
    )
    parent.add_argument(
        "--ref",
        type=str,
        default="EOM-CCSD",
        help="Reference parser name for comparison (usually EOM-CCSD)",
    )
    parent.add_argument("--debug", action="store_true", help="Debug mode")

    parser = argparse.ArgumentParser(
        description="Extract data from excited-state calculations output."
        "supported methods: QChem(CIS, CISD, EOM-CCSD), Turbomole(CC2)"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # default extract mode (hidden in help)
    _ext = sub.add_parser(
        "extract",
        parents=[parent],
        help="Extract data from excited-state calculations output file (default behaviour)",
        description="Extract data from excited-state calculations output file",
    )

    # compare
    cmp = sub.add_parser(
        "compare", parents=[parent], help="Compare two different states"
    )
    cmp.add_argument("state1", type=str, help="The first state, e.g. CIS:singlet-1/A")
    cmp.add_argument("state2", type=str, help="The second state")
    cmp.add_argument("--acc-method", type=str, default="inner-prod")
    # compareall
    cpa = sub.add_parser(
        "compareall", parents=[parent], help="Compare all different states."
    )
    cpa.add_argument(
        "--output",
        type=str,
        default="compare",
        help="Output file format for all-state comparison",
    )
    cpa.add_argument("--acc-method", type=str, default="inner-prod")
    # descriptors
    dsc = sub.add_parser(
        "descriptors", parents=[parent], help="Extract descriptors and write to CSV"
    )
    dsc.add_argument(
        "--descriptors-file",
        type=str,
        default="descriptors.csv",
        help="CSV filename for descriptors",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    # always do the extraction step first
    parsers = extract_and_write(args)

    # then dispatch
    if args.command == "compare":
        return do_compare(args, parsers)
    if args.command == "compareall":
        do_compare_all(args, parsers)
        return 0
    if args.command == "descriptors":
        return do_descriptors(args, parsers)
    # "extract" → nothing more
    return 0


if __name__ == "__main__":
    sys.exit(main())
