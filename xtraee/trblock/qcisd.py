from xtraee.lazypattern import LP
from xtraee.transition import CISDTransition

from .trblock import TransitionBlock


class CISDTransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        excitation: str = "",
        excitation_energy: float = 0.0,
        oscillator_strength: float = 0.0,
    ):
        super().__init__(
            id_number, irrep, excitation, excitation_energy, oscillator_strength
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

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"CISD transition {self.id_number}/{self.irrep} {self.ee_type}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
            # f"Oscillator strength (a.u.): {self.oscillator_strength:.6f}, \n"
        )
