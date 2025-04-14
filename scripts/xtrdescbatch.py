#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd
from xtraee.qc import QCCSDParser


def find_log_files(start_dir: Path):
    """
    Recursively search for 'CCSD/tda_opt.log' within all subdirectories of start_dir.

    Args:
        start_dir (Path): The base directory to start the search.

    Returns:
        list[Path]: A list of paths where the file was found.
    """
    # Use the '**' pattern to search recursively.
    # This finds files that exactly match 'CCSD/tda_opt.log' in any depth.
    return list(start_dir.glob("**/CCSD/test/tda_opt.log"))


def process_file(file_path: Path):
    """
    Placeholder for file-specific processing.

    Args:
        file_path (Path): Path to the 'tda_opt.log' file.

    This function is where you can add your code to process the file.
    """
    print(f"Processing file: {file_path}")
    qccsd = QCCSDParser(file_path)
    qccsd.process_file()
    return qccsd.gather_descriptors(
        extra_id=file_path.parent.parent.relative_to(Path.cwd()).as_posix()
    )


def main():
    # Set up the command line argument parser.
    parser = argparse.ArgumentParser(
        description="Recursively gather descriptors from the directory"
    )
    parser.add_argument(
        "directory",
        type=str,
        nargs="?",
        default=".",
        help="The directory to search in (defaults to current working directory).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="descriptors.csv",
        help="The output file (defaults to descriptors.csv)",
    )

    args = parser.parse_args()

    # Resolve the base directory path.
    start_dir = Path(args.directory).resolve()
    if not start_dir.is_dir():
        print(f"Error: '{start_dir}' is not a valid directory.")
        return

    # Find all matching 'tda_opt.log' files.
    log_files = find_log_files(start_dir)

    if log_files:
        print(f"Found {len(log_files)} file(s):")
        data = []
        for log_file in log_files:
            print(f"  {log_file}")
            data.append(process_file(log_file))
        pd.DataFrame(data).to_csv(args.output, index=False)
    else:
        print("No 'CCSD/tda_opt.log' files found.")


if __name__ == "__main__":
    raise SystemExit(main())
