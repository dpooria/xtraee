import logging
import re
from enum import Enum
from typing import Optional

from xtraee.irrep import Irrep
from xtraee.parser.base import BaseParser, DatasetType, PathType
from xtraee.trblock import CC2TransitionBlock, TransitionBlock

meta_patterns = {
    "trblock": re.compile(
        r"^\s*\|\s*type:\s*\S+\s+symmetry:\s*(\S+)\s+state:\s*(\d+)\s*\|$"
    ),
    "table": re.compile(
        r"^\s*\|\s*(\w+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*[+-]?\d+\.\d+\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*[+-]?\d+\.\d+\s*\|\s*([+-]?\d+\.\d+)\s*\|\s*([+-]?\d+\.\d+)\s*\|$"
    ),
}

ee_table = (
    "| sym | multi | state |          CC2 excitation energies       |  %t1   |  %t2   |"
)
table_end = "Energy:"


class Block(Enum):
    null = 0
    ee = 1
    tr = 2


class TMCC2Parser(BaseParser):
    name = "CC2TM"

    def __init__(
        self,
        input_file: str | PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.block = Block.null
        self.parser = {
            Block.null: lambda line: None,
            Block.ee: self.process_table,
            Block.tr: self.process_trblocks,
        }
        self._singlets_processed = 0
        self._triplets_processed = 0
        self._current_trblock: Optional[TransitionBlock] = None

    def detect_block(self, line):
        if ee_table in line:
            self.block = Block.ee

    def process_table(self, line):
        if table_end in line:
            self.block = Block.tr
            self.irrep_singlets = [
                irrep
                for irrep in self.irreps_dict.values()
                if irrep.ee_type == "singlet"
            ]
            self.irrep_triplets = [
                irrep
                for irrep in self.irreps_dict.values()
                if irrep.ee_type == "triplet"
            ]
            for irrep in self.irrep_singlets:
                irrep.n_states = len(irrep.trblocks)
                # irrep.sort()
                self.N_singlets += irrep.n_states
            for irrep in self.irrep_triplets:
                irrep.n_states = len(irrep.trblocks)
                # irrep.sort()
                self.N_triplets += irrep.n_states

        elif (m := meta_patterns["table"].match(line)) is not None:
            irrep = m.group(1)
            multi = int(m.group(2))
            if multi == 1:
                multi = "singlet"
            elif multi == 3:
                multi = "triplet"
            else:
                raise ValueError(f"Unknown multiplicity {multi}")
            id_number = int(m.group(3))
            cc2_energy = float(m.group(4))
            t1 = float(m.group(5))
            t2 = float(m.group(6))
            transition_block = CC2TransitionBlock(
                id_number, irrep, multi, cc2_energy, t1, t2
            )
            k = f"{multi}-{irrep}"
            if k not in self.irreps_dict:
                self.irreps_dict[k] = Irrep(irrep, multi, 0)
            self.irreps_dict[k].append(transition_block)

    def process_trblocks(self, line):
        if (m := meta_patterns["trblock"].match(line)) is not None:
            if self._singlets_processed < self.N_singlets:
                multi = "singlet"
                self._singlets_processed += 1
            else:
                multi = "triplet"
                self._triplets_processed += 1
            irrep = m.group(1)
            state = int(m.group(2))
            self._current_trblock = self.irreps_dict[f"{multi}-{irrep}"].trblocks[
                state - 1
            ]
        elif self._current_trblock is not None:
            self._current_trblock.add_data(line)
            # if self.current_trblock.add_data(line):
            #     # self.current_trblock.sort()
            #     self.current_trblock = None
