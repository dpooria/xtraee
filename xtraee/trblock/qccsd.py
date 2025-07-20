from xtraee.lazypattern import LP, VERBOSE
from xtraee.transition import CCSDTransition
from xtraee.utils import Ha, nan
from .trblock import TransitionBlock


class CCSDTransitionBlock(TransitionBlock):
    ee_pattern = LP(
        r"""
        \s*Total\ energy\s*=\s*
        (?P<tot>[-+]?\d+\.\d+)
        .*?
        Excitation\ energy\s*=\s*
        (?P<ee>[-+]?\d*\.?\d+)
        \s*eV\.
        """,
        VERBOSE,
    )
    r2_pattern = LP(
        r"^.*R0\^2\s*=\s*(\d*.\d+)\s*R1\^2\s*=\s*([-+]?\d*\.?\d+)\s*R2\^2\s*=\s*([-+]?\d*\.?\d+).*$"
    )

    occ_frontier_no = "Occupation of frontier NOs:"
    frontier_no_pattern = LP(r"^\s*([-+]?[0-9]+\.[0-9]+)\s+([-+]?[0-9]+\.[0-9]+)\s*$")
    unpaired_no_pattern = LP(
        r"^\s*Number of unpaired electrons:\s*n_u\s*=\s*([-+]?[0-9]+\.[0-9]+),\s*n_u,nl\s*=\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )
    unpaired_no = "Number of unpaired electrons:"
    prno_ind = "NO participation ratio (PR_NO):"
    prno_pattern = LP(
        r"^\s*NO participation ratio \(PR_NO\):\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )

    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        multi: str = "",
        name="CCSD",
        excitation_energy: float = nan,
        total_energy: float = nan,
        oscillator_strength: float = nan,
    ):
        super().__init__(
            id_number,
            irrep,
            multi,
            name,
            excitation_energy,
            total_energy,
            oscillator_strength,
        )
        self.transitions: list[CCSDTransition] = []
        self.tr_cls = CCSDTransition
        self.tr_indicator = "->"
        self.end_trblock = "Summary of significant orbitals:"
        self.wait_frontier_no = False
        self.meta_data.update(
            {
                "R0": nan,
                "R1": nan,
                "R2": nan,
                "gamma": nan,
                "omega": nan,
                "loc": nan,
                "phe": nan,
                "rhre": nan,
                "alphabeta": nan,
                "corr_coef": nan,
                "froniter_no": [],
                "nu": nan,
                "nl": nan,
                "prno": nan,
            }
        )

    def extras(self, line: str) -> None:
        meta = self.meta_data
        if not self.completed:
            if (m := self.r2_pattern.match(line)) is not None:
                meta["R0"] = float(m.group(1))
                meta["R1"] = float(m.group(2))
                meta["R2"] = float(m.group(3))
            elif (m := self.ee_pattern.match(line)) is not None:
                self.excitation_energy = float(m["ee"])
                self.total_energy = float(m["tot"]) * Ha
        elif not self.completed_extras:
            if self.occ_frontier_no in line:
                self.wait_frontier_no = True
            elif self.wait_frontier_no:
                if (m := self.frontier_no_pattern.match(line)) is not None:
                    self.wait_frontier_no = False
                    meta["froniter_no"] = [float(m.group(1)), float(m.group(2))]
                else:
                    raise ValueError(
                        f"Could not match the Occupations of the frontier NO for {self}"
                    )
            elif self.unpaired_no in line:
                if m := self.unpaired_no_pattern.match(line):
                    meta["nu"] = m.group(1)
                    meta["nl"] = m.group(2)
                else:
                    raise ValueError(
                        f"Could not match the number of unpaired electrons for {self}"
                    )
            elif self.prno_ind in line:
                if (m := self.prno_pattern.match(line)) is not None:
                    meta["prno"] = m.group(1)
                else:
                    raise ValueError(
                        f"Could not match the participation number for {self}"
                    )
                self.completed_extras = True
        else:
            raise ValueError(f"Transition {self} is already completed!")

    def compare(
        self, other: TransitionBlock, method: str
    ) -> tuple[float, float, float]:
        # CCSD is usually the reference
        if isinstance(other, CCSDTransitionBlock):
            return super().compare(other, method)
        else:
            return other.compare(self, method)
