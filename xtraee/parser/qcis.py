import re
from enum import Enum
from typing import Optional

from xtraee.irrep import Irrep
from xtraee.parser.base import BaseParser, DatasetType, PathType
from xtraee.trblock import CISTransitionBlock, TransitionBlock

start_indicators = {"cisee": "CIS Excitation Energies", "mo": "Orbital Energies (a.u.)"}

trblock_begin_pattern = re.compile(
    r"^Excited\s+state\s+(\d+)\s*:\s*excitation\s+energy\s+\(eV\)\s*=\s*([-+]?\d+.\d+)\s*$"
)


class Block(Enum):
    null = 0
    ee = 1
    mo = 2


class QCISParser(BaseParser):
    name = "CIS"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.block = Block.null
        self.parser = {
            Block.null: lambda line: None,
            Block.ee: self.process_trblock,
            Block.mo: lambda line: None,
        }
        self.vsccsd: DatasetType = {}
        self._inside_eomee = False
        self._current_transition: Optional[TransitionBlock] = None
        self.irrep_singlets = Irrep("A", "singlet", 0, "CIS")
        self.irrep_triplets = Irrep("A", "triplet", 0, "CIS")
        self.homo = 0

    def detect_block(self, line: str) -> None:
        if start_indicators["cisee"] in line:
            self.block = Block.ee
        elif start_indicators["mo"] in line:
            self.process_irreps()
            self.block = Block.mo

    def process_trblock(self, line: str) -> None:
        if (m := trblock_begin_pattern.match(line)) is not None:
            if self._current_transition is not None:
                self._current_transition.sort()
                if self._current_transition.excitation == "singlet":
                    self.irrep_singlets.append(self._current_transition)
                elif self._current_transition.excitation == "triplet":
                    self.irrep_triplets.append(self._current_transition)
                else:
                    raise ValueError(f"Unknown excitation {self._current_transition}")
            self._current_transition = CISTransitionBlock(int(m.group(1)))
            self._current_transition.excitation_energy = float(m.group(2))
        elif self._current_transition is not None:
            self._current_transition.add_data(line)

    def process_irreps(self) -> None:
        self.irrep_singlets.sort()
        self.irrep_triplets.sort()
        self.irrep_singlets.n_states = len(self.irrep_singlets)
        self.irrep_triplets.n_states = len(self.irrep_triplets)
        singlet_trblocks = self.irrep_singlets.trblocks
        triplet_trblocks = self.irrep_triplets.trblocks
        self.homo = 0
        for tr in singlet_trblocks:
            self.homo = max(
                *[int(id[1]) for tr_ in tr.transitions for id in tr_.id_i], self.homo
            )
        for tr in triplet_trblocks:
            self.homo = max(
                *[int(id[1]) for tr_ in tr.transitions for id in tr_.id_i], self.homo
            )
        for tr in singlet_trblocks:
            tr.homo = self.homo
        for tr in triplet_trblocks:
            tr.homo = self.homo
        self.irreps_dict = {
            "singlet": self.irrep_singlets,
            "triplet": self.irrep_triplets,
        }
