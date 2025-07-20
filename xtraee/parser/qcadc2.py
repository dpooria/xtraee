from xtraee.irrep import Irrep
from xtraee.lazypattern import LP, VERBOSE
from xtraee.parser.base import BaseParser, Block, PathType
from xtraee.trblock import ADC2TransitionBlock


trblock_begin_pattern = LP(
    r"""
    ^\s*Excited\s+state\s+          # literal
    (?P<root>\d+)\s+                # capture root
    \((?P<multi>triplet|singlet)    # capture multi
    ,\s+(?P<irrep>[\w'"]+)\)        # capture irrep
    \s+\[converged\]$               # literal, maybe we should also consider not converged case
    """,
    VERBOSE,
)

start_indicators = {
    Block.null: None,
    Block.input: "$rem",
    Block.irrep: "Excited State Summary",
    Block.mo: "Orbital Energies (a.u.)",
}


class QCADC2Parser(BaseParser):
    name = "ADC2"

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
        self._current_trblock: None | ADC2TransitionBlock = None
        self._irrep_counter: int = 0

    def detect_block(self, line: str) -> None:
        for block, indicator in start_indicators.items():
            if indicator and indicator in line:
                self.block = block
                break

    def process_irreps(self, line: str) -> None:
        if m := trblock_begin_pattern.match(line):
            key = f"{m['multi']}-{m['irrep']}"
            irrep = self.irreps.get(key)
            if irrep is None:
                irrep = Irrep(m["irrep"], m["multi"], 0, self.name)
                self.irreps[key] = irrep
            trblock = ADC2TransitionBlock(int(m["root"]), irrep.name, irrep.multi)
            irrep.append(trblock)
            self._current_trblock = trblock
        elif self._current_trblock is not None:
            if self._current_trblock.add_data(line):
                self.log.debug(f"Completed {self._current_trblock}")

    def finalize_irreps(self, line: str = "") -> None:
        homo = 0
        for irrep in self.irreps.values():
            irrep.sort()
            irrep.n_states = len(irrep)

            for tr in irrep:
                if tr.transitions:
                    homo = max(
                        *[id_.orb_num for tr_ in tr.transitions for id_ in tr_.id_i],
                        homo,
                    )
        homo += 1
        for irrep in self.irreps.values():
            irrep.scatter_attr("homo", homo)
