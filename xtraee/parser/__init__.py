import re
from pathlib import Path

from xtraee.parser.base import BaseParser
from xtraee.parser.qccsd import QCCSDParser
from xtraee.parser.qcis import QCISParser
from xtraee.parser.tmcc2 import TMCC2Parser


def Parser(input_file: str | Path, *args, **kwargs) -> BaseParser:
    """
    Detects the appropriate parser based on the content of the input file.
    """
    with open(input_file, "r") as f:
        for line in f:
            ls = re.sub(r"\s+", " ", line.strip().lower())
            if "method eom-ccsd" in ls:
                parser = QCCSDParser(input_file, *args, **kwargs)
                break
            elif "cis excitation energies" in ls:
                parser = QCISParser(input_file, *args, **kwargs)
                break
            elif "cc2 - approximate cc singles and doubles" in ls:
                parser = TMCC2Parser(input_file, *args, **kwargs)
                break
        else:
            raise ValueError("Unrecognized file: parser not found.")

    return parser


__all__ = ["Parser", "QCCSDParser", "QCISParser", "TMCC2Parser"]
