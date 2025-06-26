import re
from pathlib import Path

from xtraee.parser.base import BaseParser
from xtraee.parser.qccsd import QCCSDParser, QCC2Parser
from xtraee.parser.qcis import QCISParser, QCTDDFTParser
from xtraee.parser.qcisd import QCISDParser
from xtraee.parser.tmcc2 import TMCC2Parser


def Parser(input_file: str | Path, *args, **kwargs) -> BaseParser:
    """
    Detects the appropriate parser based on the content of the input file.
    """
    with open(input_file, "r") as f:
        for line in f:
            ls = re.sub(r"\s+", " ", line.strip().lower())
            if "eom-ccsd" in ls:
                parser = QCCSDParser(input_file, *args, **kwargs)
                break
            elif "eom-cc2" in ls:
                parser = QCC2Parser(input_file, *args, **kwargs)
                break
            elif "cis excitation energies" in ls:
                parser = QCISParser(input_file, *args, **kwargs)
                break
            elif "tddft" in ls:
                parser = QCTDDFTParser(input_file, *args, **kwargs)
                break
            elif "doing genuine cisd calculations" in ls or "method cisd" in ls:
                parser = QCISDParser(input_file, *args, **kwargs)
                break
            elif "cc2 - approximate cc singles and doubles" in ls:
                parser = TMCC2Parser(input_file, *args, **kwargs)
                break
        else:
            raise ValueError(f"Could not find suitable parser for {input_file}")

    return parser


__all__ = ["Parser", "QCCSDParser", "QCISParser", "QCISDParser", "TMCC2Parser"]
