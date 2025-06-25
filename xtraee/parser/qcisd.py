from xtraee.irrep import Irrep
from xtraee.lazypattern import LP, VERBOSE
from xtraee.parser.base import BaseParser, Block, PathType
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
        self.parser.update(
            {Block.irrep: self.process_irreps, Block.mo: self.finalize_irreps}
        )

    def reset(self):
        super().reset()
        self._current_trblock: None | CISDTransitionBlock = None
        self._current_irrep: None | Irrep = None
        self._irrep_counter: int = 0

    def detect_block(self, line: str) -> None:
        for block, indicator in start_indicators:
            if indicator and indicator in line:
                self.block = block
                break

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
            self._current_trblock = CISDTransitionBlock(int(m.group(1)))
            self._current_trblock.excitation_energy = float(m.group(2))
        elif self._current_trblock is not None:
            self._current_trblock.add_data(line)

    def process_irreps(self, line: str) -> None:
        if m := meta_patterns["irreps"].match(line):
            multi = "single" if m["multi"] == "LOWSPIN" else "triple"
            irrep = Irrep(m["irrep"], multi, int(m["n_roots"]), "CISD")
            self.irreps_dict[irrep.identifier] = irrep
            self._current_irrep = irrep
        elif m := trblock_begin_pattern.match(line):
            if self._current_irrep is None:
                raise ValueError("trblock before irrep")
            trblock = CISDTransitionBlock(m["root"])

    def finalize_irreps(self) -> None:
        homo = 0
        for irrep in self.irreps_dict.values():
            irrep.sort()
            assert irrep.n_states == len(
                irrep
            ), f"Inconsitent number of states in {irrep}"

            for tr in irrep:
                homo = max(
                    *[int(id_[1]) for tr_ in tr.transitions for id_ in tr_.id_i],
                    homo,
                )
        self.homo = homo
        for irrep in self.irreps_dict.values():
            irrep.scatter_attr("homo", self.homo)
