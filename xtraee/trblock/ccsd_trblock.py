import re

from xtraee.transition import CCSDTransition

from .trblock import TransitionBlock


class CCSDTransitionBlock(TransitionBlock):
    ee_pattern = re.compile(r"^.*Excitation energy\s*=\s*([-+]?\d*\.?\d+)\s*eV\.\s*$")
    r2_pattern = re.compile(
        r"^.*R0\^2\s*=\s*(\d*.\d+)\s*R1\^2\s*=\s*([-+]?\d*\.?\d+)\s*R2\^2\s*=\s*([-+]?\d*\.?\d+).*$"
    )

    occ_frontier_no = "Occupation of frontier NOs:"
    frontier_no_pattern = re.compile(
        r"^\s*([-+]?[0-9]+\.[0-9]+)\s+([-+]?[0-9]+\.[0-9]+)\s*$"
    )
    unpaired_no_pattern = re.compile(
        r"^\s*Number of unpaired electrons:\s*n_u\s*=\s*([-+]?[0-9]+\.[0-9]+),\s*n_u,nl\s*=\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )
    unpaired_no = "Number of unpaired electrons:"
    prno_ind = "NO participation ratio (PR_NO):"
    prno_pattern = re.compile(
        r"^\s*NO participation ratio \(PR_NO\):\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transitions: list[CCSDTransition] = []
        self.tr_cls = CCSDTransition
        self.tr_indicator = "->"
        self.end_trblock = "Summary of significant orbitals:"
        self.R0 = 0.0
        self.R1 = 0.0
        self.R2 = 0.0
        self.gamma = 0.0
        self.omega = 0.0
        self.loc = 0.0
        self.phe = 0.0
        self.rhre = 0.0
        self.alphabeta = 0.0
        self.corr_coef = 0.0
        self.froniter_no = []
        self.wait_frontier_no = False
        self.nu = 0.0
        self.nl = 0.0
        self.prno = 0.0

    def extras(self, line: str) -> None:
        if not self.completed:
            if (m := self.r2_pattern.match(line)) is not None:
                self.R0 = float(m.group(1))
                self.R1 = float(m.group(2))
                self.R2 = float(m.group(3))
            elif (m := self.ee_pattern.match(line)) is not None:
                self.excitation_energy = float(m.group(1))
        elif not self.completed_extras:
            if self.occ_frontier_no in line:
                self.wait_frontier_no = True
            elif self.wait_frontier_no:
                if (m := self.frontier_no_pattern.match(line)) is not None:
                    self.wait_frontier_no = False
                    self.froniter_no = [float(m.group(1)), float(m.group(2))]
                else:
                    raise ValueError(
                        f"Could not match the Occupations of the frontier NO for {self}"
                    )
            elif self.unpaired_no in line:
                if m := self.unpaired_no_pattern.match(line):
                    self.nu = m.group(1)
                    self.nl = m.group(2)
                else:
                    raise ValueError(
                        f"Could not match the number of unpaired electrons for {self}"
                    )
            elif self.prno_ind in line:
                if (m := self.prno_pattern.match(line)) is not None:
                    self.prno = m.group(1)
                else:
                    raise ValueError(
                        f"Could not match the participation number for {self}"
                    )
                self.completed_extras = True
        else:
            raise ValueError(f"Transition {self} is already completed!")

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"EOMEE transition {self.excitation} {self.id_number}/{self.irrep}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            f"R0^2: {self.R0:.4f} R1^2: {self.R1:.4f} R2^2: {self.R2:.4f}\n"  # noqa
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
            f"Oscillator strength (a.u.): {self.oscillator_strength:.6f},"
            f"omega (Mulliken): {self.omega:.4f}\n"
        )

    def compare(
        self, other: TransitionBlock, method: str
    ) -> tuple[float, float, float]:
        if isinstance(other, CCSDTransitionBlock):
            return super().compare(other, method)
        else:
            return other.compare(self, method)
