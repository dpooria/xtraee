from typing import Optional

from xtraee.irrep import Irrep
from xtraee.lazypattern import LP, VERBOSE
from xtraee.parser.base import BaseParser, Block, PathType
from xtraee.parser.qccsd import meta_patterns as qccsd_meta_patterns
from xtraee.trblock import CISDTransitionBlock

trblock_begin_pattern = LP(
    r"""
    ^Root\s+                      # literal “Root ” at start
    (?P<root>\d+)                 # 1) capture root number
    \s+Conv-d\s+yes\s+Tot\s+Ene=  # skip the fixed text
    \s*[-+]?\d+\.\d+\s+hartree    # skip the total energy value
    \s*\(Ex\s+Ene\s+              # skip “(Ex Ene ”
    (?P<ex>[-+]?\d+\.\d+)         # 2) capture Ex Ene
    \s+eV\),\s*                   # skip “ eV),”
    U0\^2=\s*(?P<U0>[-+]?\d+\.\d+),\s*  # 3) capture U0
    U1\^2=\s*(?P<U1>[-+]?\d+\.\d+),\s*  # 4) capture U1
    U2\^2=\s*(?P<U2>[-+]?\d+\.\d+)      # 5) capture U2
    """,
    VERBOSE,
)

meta_patterns = dict(
    irreps=LP(
        r"^\s*(?P<n_roots>\d+)\s+lowest\s+(?P<multi>LOWSPIN|HIGHSPIN)"
        r"\s+roots\s+of\s+symmetry\s+(?P<irrep>.+?)\s*:\s*$"
    ),
)

start_indicators = {
    Block.null: None,
    Block.input: "$rem",
    Block.irrep: "roots of symmetry",
    Block.mo: "Orbital Energies (a.u.)",
}


class QCISDParser(BaseParser):
    name = "CISD"

    def __init__(
        self,
        input_file: PathType,
        threshold: float = 0.0,
    ):
        super().__init__(input_file, threshold)
        self.parser = {
            Block.null: lambda line: None,
            Block.input: self.process_input,
            Block.irrep: self.process_irreps,
            Block.mo: self.finalize_irreps,
        }
        self._current_transition: Optional[CISDTransitionBlock] = None
        self._current_irrep: Optional[Irrep] = None

    def detect_block(self, line: str) -> None:
        for block, indicator in start_indicators:
            if indicator is not None and indicator in line:
                self.block = block
                break

    def process_input_block(self, line: str) -> None:
        if self.N_singlets == 0:
            if (m := qccsd_meta_patterns["ee_singlets"].match(line)) is not None:
                self.ee_singlets = list(map(int, m.group(1).strip().split(",")))
                self.N_singlets = len(self.ee_singlets)
        if self.N_triplets == 0:
            if (m := qccsd_meta_patterns["ee_triplets"].match(line)) is not None:
                self.ee_triplets = list(map(int, m.group(1).strip().split(",")))
                self.N_triplets = len(self.ee_triplets)

    def process_trblock(self, line: str) -> None:
        if (m := trblock_begin_pattern.match(line)) is not None:
            if self._current_transition is not None:
                self._current_transition.sort()
                if self._current_transition.ee_type == "singlet":
                    self.irrep_singlets.append(self._current_transition)
                elif self._current_transition.ee_type == "triplet":
                    self.irrep_triplets.append(self._current_transition)
                else:
                    raise ValueError(f"Unknown excitation {self._current_transition}")
            self._current_transition = CISDTransitionBlock(int(m.group(1)))
            self._current_transition.excitation_energy = float(m.group(2))
        elif self._current_transition is not None:
            self._current_transition.add_data(line)

    def process_irreps(self, multi: Block) -> None:
        pass

    def finalize_irreps(self) -> None:
        for irrep in self.irreps_dict.values():
            irrep.sort()
            assert irrep.n_states == len(
                irrep
            ), f"Inconsitent number of states in {irrep}"

            for tr in irrep:
                self.homo = max(
                    *[int(id[1]) for tr_ in tr.transitions for id in tr_.id_i],
                    self.homo,
                )
        # broadcast homo
        for irrep in self.irreps_dict.values():
            for tr in irrep:
                tr.homo = self.homo
