from xtraee.irrep import Irrep
from xtraee.lazypattern import LP
from xtraee.parser.base import BaseParser, Block, PathType
from xtraee.trblock import CISTransitionBlock

start_indicators = {"cisee": "CIS Excitation Energies", "mo": "Orbital Energies (a.u.)"}

trblock_begin_pattern = LP(
    r"^Excited\s+state\s+(\d+)\s*:\s*excitation\s+energy\s+\(eV\)\s*=\s*([-+]?\d+.\d+)\s*$"
)


class QCISParser(BaseParser):
    name = "CIS"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.parser.update(
            {Block.ee: self.process_trblock, Block.mo: lambda line: None}
        )
        # not supporting symmetry for CIS
        self.N_singlets = 1
        self.N_triplets = 1
        self.irrep_singlets = Irrep("A", "singlet", 0, "CIS")
        self.irrep_triplets = Irrep("A", "triplet", 0, "CIS")

    def reset(self):
        super().reset()
        self._current_trblock: None | CISTransitionBlock = None

    def detect_block(self, line: str) -> None:
        if start_indicators["cisee"] in line:
            self.block = Block.ee
        elif start_indicators["mo"] in line:
            # move this as detect_block should only detect
            self.process_irreps()
            self.block = Block.mo

    def process_trblock(self, line: str) -> None:
        if (m := trblock_begin_pattern.match(line)) is not None:
            if self._current_trblock is not None:
                self._current_trblock.sort()
                if self._current_trblock.ee_type == "singlet":
                    self.irrep_singlets.append(self._current_trblock)
                elif self._current_trblock.ee_type == "triplet":
                    self.irrep_triplets.append(self._current_trblock)
                else:
                    raise ValueError(f"Unknown excitation {self._current_trblock}")
            self._current_trblock = CISTransitionBlock(int(m.group(1)))
            self._current_trblock.excitation_energy = float(m.group(2))
        elif self._current_trblock is not None:
            self._current_trblock.add_data(line)

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
                *[id_.orb_num for tr_ in tr.transitions for id_ in tr_.id_i], self.homo
            )
        for tr in triplet_trblocks:
            self.homo = max(
                *[id_.orb_num for tr_ in tr.transitions for id_ in tr_.id_i], self.homo
            )
        self.irrep_singlets.scatter_attr("homo", self.homo)
        self.irrep_triplets.scatter_attr("homo", self.homo)
        self.irreps = {
            f"singlet-{self.irrep_singlets.name}": self.irrep_singlets,
            f"triplet-{self.irrep_triplets.name}": self.irrep_triplets,
        }
