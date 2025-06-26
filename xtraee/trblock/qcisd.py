from xtraee.lazypattern import LP
from xtraee.transition import CISDTransition

from .trblock import TransitionBlock


class CISDTransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        ee_type: str = "",
        excitation_energy: float = 0.0,
        oscillator_strength: float = 0.0,
        name="CISD",
    ):
        super().__init__(
            id_number, irrep, ee_type, excitation_energy, oscillator_strength, name=name
        )
        self.transitions: list[CISDTransition] = []
        self.tr_cls = CISDTransition
        self.end_trblock = "\n"
        # -0.6395                 12(   A1) B   ->    0(   A1) B
        self.tr_indicator = LP(
            r"^\s*[-+]?\d+\.\d+"
            r"\s+\d+\s*\(.+\)\s+[AB]"
            r"\s*->\s*"
            r"\d+\s*\(.+\)\s*[AB]\s*$"
        )
        self.R0 = 0.0
        self.R1 = 0.0
        self.R2 = 0.0

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"{self.name} transition {self.id_number}/{self.irrep} {self.ee_type}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            f"U0={self.R0:.4f}, U1={self.R1:.4f}, U2={self.R2:.4f}\n"
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
        )
