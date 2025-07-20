from xtraee.irrep import Irrep
from xtraee.lazypattern import LP, VERBOSE
from xtraee.parser.base import BaseParser, Block, PathType
from xtraee.trblock import CISDTransitionBlock

trblock_begin_pattern = LP(
    r"""
    ^Root\s+                      # literal “Root ” at start
    (?P<root>\d+)                 # 1) capture root number
    \s+Conv-d\s+yes\s+Tot\s+Ene=  # skip the fixed text
    \s*(?P<tote>[-+]?\d+\.\d+)\s+hartree    # skip the total energy value
    \s*\(Ex\s+Ene\s+              # skip “(Ex Ene ”
    (?P<ex>[-+]?\d+\.\d+)         # 2) capture Ex Ene
    \s+eV\),\s*                   # skip “ eV),”
    (U0\^2=\s*(?P<U0>[-+]?\d+\.\d+),\s*)?  # 3) capture U0
    (U1\^2=\s*(?P<U1>[-+]?\d+\.\d+),\s*)?  # 4) capture U1
    (U2\^2=\s*(?P<U2>[-+]?\d+\.\d+))?      # 5) capture U2
    """,
    VERBOSE,
)

meta_patterns = dict(
    irreps=LP(
        r"SOLVE\s+LINEAR\s+RESPONSE\s+EQUATIONS\s+FOR\s+(?P<multi>LOWSPIN|HIGHSPIN)"
        r"\s+SINGLE\s+(AND\s+DOUBLE\s+)?EXCITATIONS\s+OF\s+(?P<irrep>.+)\s+IRREP"
    ),
)

start_indicators = {
    Block.null: None,
    Block.input: "$rem",
    Block.irrep: "SOLVE LINEAR RESPONSE EQUATIONS FOR",
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
        for block, indicator in start_indicators.items():
            if indicator and indicator in line:
                self.block = block
                break

    def process_irreps(self, line: str) -> None:
        if m := meta_patterns["irreps"].match(line):
            multi = "singlet" if m["multi"] == "LOWSPIN" else "triplet"
            irrep = Irrep(m["irrep"], multi, 0, self.name)
            self.irreps[irrep.identifier] = irrep
            self._current_irrep = irrep
        elif m := trblock_begin_pattern.match(line):
            if self._current_irrep is None:
                raise ValueError("trblock before irrep")
            irrep = self._current_irrep
            trblock = CISDTransitionBlock(
                int(m["root"]),
                irrep.name,
                irrep.multi,
                self.name,
                float(m["ex"]),
                float(m["tote"]),
            )
            trblock.meta_data.update(
                {
                    "U0": float(m["U0"] or "nan"),
                    "U1": float(m["U1"] or "nan"),
                    "U2": float(m["U2"] or "nan"),
                }
            )
            irrep.append(trblock)
            self._current_trblock = trblock
        elif self._current_trblock is not None:
            if self._current_trblock.add_data(line):
                self.log.debug(f"Completed {self._current_trblock}")

    def finalize_irreps(self, line: str = "") -> None:
        homo = 0
        for irrep in self.irreps.values():
            irrep.sort()
            # assert irrep.n_states == len(irrep), (
            #     f"Inconsitent number of states in {irrep}"
            # )
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


class QCIS_D_Parser(QCISDParser):
    name = "CIS_D_"
