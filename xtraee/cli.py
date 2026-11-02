#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import TYPE_CHECKING

import pandas as pd

from xtraee import __version__
from xtraee.config import configure, get_logger
from xtraee.parser import Parser, QCCSDParser

if TYPE_CHECKING:
    from xtraee.parser import BaseParser

COMMANDS = ("extract", "compare", "compareall", "descriptors")

#: Options that consume the following token; the default-command injection has
#: to skip their values, otherwise ``--threshold 0.2 file.log`` becomes
#: ``--threshold extract``.
VALUE_OPTIONS = frozenset(
    {
        "--outdir",
        "--out-data",
        "--out-fdata",
        "--out-vsstd",
        "--ref",
        "--threshold",
        "--descriptors-file",
        "--acc-method",
        "--output",
    }
)


def read_inputs(args) -> dict[str, BaseParser]:
    """Parse every input file, keyed by the parser (method) name."""
    log = get_logger("reader")
    parsers: dict[str, BaseParser] = {}
    sources: dict[str, str] = {}
    for infile in args.input:
        parser = Parser(infile, args.threshold)
        if parser is None:  # pragma: no cover - Parser() raises instead
            raise ValueError(f"Could not find a suitable parser for {infile}")
        parser.process_file()
        if parser.name in parsers:
            log.warning(
                "%s already holds %s; %s replaces it",
                parser.name,
                sources[parser.name],
                infile,
            )
        parsers[parser.name] = parser
        sources[parser.name] = str(infile)
    return parsers


def extract_and_write(args) -> dict[str, BaseParser]:
    """Process inputs and write full/overview/vsstd outputs."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    parsers = read_inputs(args)

    for name, p in parsers.items():
        if args.full:
            p.write_full(outdir / f"{args.out_fdata}_{name}.txt")
        p.write_dataset(outdir / f"{args.out_data}_{name}.txt")

    ref = parsers.get(args.ref)
    if ref:
        for name, p in parsers.items():
            if name != args.ref:
                p.write_vs_std(outdir / f"{args.out_vsstd}_{name}.txt", ref.irreps)

    return parsers


def do_compare(args, parsers: dict[str, BaseParser]) -> int:
    """Compare two individual states and print the result."""
    log = get_logger("compare")
    try:
        m1, s1 = args.state1.split(":", 1)
        m2, s2 = args.state2.split(":", 1)
        m1 = "EOM-CCSD" if m1.upper() == "CCSD" else m1
        m2 = "EOM-CCSD" if m2.upper() == "CCSD" else m2
    except ValueError:
        log.error(
            "Could not parse %r and %r; format is "
            "{method}:{excitation}-{state-id}/{irrep-id}, e.g. CIS:singlet-1/A",
            args.state1,
            args.state2,
        )
        return 1

    p1, p2 = parsers.get(m1), parsers.get(m2)
    if p1 is None or p2 is None:
        available = ", ".join(sorted(parsers)) or "none"
        log.error(
            "Parser for %s or %s not found! Parsed methods: %s", m1, m2, available
        )
        return 1

    tb1 = p1.find_trblock(s1)
    tb2 = p2.find_trblock(s2)
    if tb1 is None or tb2 is None:
        missing = s1 if tb1 is None else s2
        log.error("State %r not found in %s", missing, m1 if tb1 is None else m2)
        return 1

    accuracy, error, matched = tb1.compare(tb2, args.acc_method)
    # print, not log: this is the answer the command exists for, and it must
    # not depend on the log level.
    print(f"Comparing {s1} and {s2}")
    print()
    print("accuracy          | error             | fraction matched")
    print(f"{accuracy:<17.10f} | {error:<17.10f} | {matched:.4f}")
    return 0


def do_compare_all(args, parsers: dict[str, BaseParser]) -> None:
    """Compare every state (and vs the reference, if available)."""
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    refname = args.ref
    ref = parsers.get(refname)

    for name, p in parsers.items():
        data = p.compare_all(args.acc_method)
        for key, df in data.items():
            df.to_csv(outdir / f"{key}_{name}_{args.output}.csv")

        if ref and name != refname:
            data_vs = p.compare_std(ref.irreps, args.acc_method)
            for key, df in data_vs.items():
                df.to_csv(outdir / f"{key}_{args.output}.csv")


def resolve_descriptors_path(args) -> Path:
    """``--descriptors-file`` is a name inside ``--outdir`` unless it is a path."""
    path = Path(args.descriptors_file)
    if path.is_absolute() or path.parent != Path("."):
        return path
    return Path(args.outdir) / path


def do_descriptors(args, parsers: dict[str, BaseParser]) -> int:
    """Gather descriptors for EOM-CCSD and dump to CSV."""
    log = get_logger("descriptor")
    not_ccsd = sorted(
        name for name, parser in parsers.items() if not isinstance(parser, QCCSDParser)
    )
    if not_ccsd:
        log.error(
            "Descriptors are only implemented for EOM-CCSD/CC2 outputs, got: %s",
            ", ".join(not_ccsd),
        )
        return 1

    eom = parsers.get(args.ref)
    if not isinstance(eom, QCCSDParser):
        available = ", ".join(sorted(parsers)) or "none"
        log.error(
            "No %s input to take descriptors from (parsed methods: %s)",
            args.ref,
            available,
        )
        return 1

    rows = eom.gather_descriptors()
    if not rows:
        log.warning("No excited states found in %s", eom.input_file)
        return 1
    path = resolve_descriptors_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    print(f"Wrote {len(rows)} states with descriptors to {path}")
    return 0


def inject_default_command(argv: list[str]) -> list[str]:
    """Make ``xtraee file.log`` mean ``xtraee extract file.log``.

    argparse only accepts the shared options *after* the sub-command, so a
    command found further right is hoisted to the front: ``xtraee --outdir out
    compare a.log b.log ...`` becomes ``compare --outdir out a.log b.log ...``.
    The tokens are walked the way argparse walks them, so the value of an
    option is never mistaken for the sub-command.
    """
    args = list(argv)
    if any(arg in ("-h", "--help", "--version") for arg in args):
        return args
    skip_next = False
    for i, tok in enumerate(args):
        if skip_next:
            skip_next = False
            continue
        if tok.startswith("-"):
            if "=" not in tok and tok in VALUE_OPTIONS:
                skip_next = True
            continue
        if tok in COMMANDS:
            return [tok] + args[:i] + args[i + 1:]
        break
    return ["extract"] + args


def build_parser() -> argparse.ArgumentParser:
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
    parent.add_argument("--outdir", type=str, default=".",
                        help="Output directory")
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
    parser.add_argument(
        "--version", action="version", version=f"xtraee {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # default extract mode (hidden in help)
    sub.add_parser(
        "extract",
        parents=[parent],
        help="Extract data from excited-state calculations output file (default behaviour)",
        description="Extract data from excited-state calculations output file",
    )

    # compare
    cmp = sub.add_parser(
        "compare", parents=[parent], help="Compare two different states"
    )
    cmp.add_argument("state1", type=str,
                     help="The first state, e.g. CIS:singlet-1/A")
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
        help="CSV filename for descriptors (inside --outdir unless it is a path)",
    )

    return parser


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    argv = list(sys.argv[1:] if argv is None else argv)
    return build_parser().parse_args(inject_default_command(argv))


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    configure(args.debug)

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
    return 0


if __name__ == "__main__":
    sys.exit(main())
