from xtraee.lazypattern import LP
from xtraee.transition import CISTransition

from .trblock import TransitionBlock

oscillator_pattern = LP(r"Strength\s+:\s*([-+]?\d*\.?\d+)\s*")
multplicity_pattern = LP(r"Multiplicity:\s*(Singlet|Triplet)")


class CISTransitionBlock(TransitionBlock):
    def __init__(
        self,
        *args,
        name="CIS",
        **kwargs,
    ):
        super().__init__(*args, **kwargs, name=name)
        self.transitions: list[CISTransition] = []
        self.tr_cls = CISTransition
        self.tr_indicator = "-->"
        self.end_trblock = "\n"

    def extras(self, line: str):
        if (m := oscillator_pattern.match(line)) is not None:
            self.oscillator_strength = float(m.group(1))
        elif (m := multplicity_pattern.match(line)) is not None:
            self.ee_type = m.group(1).lower()

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"{self.name} transition {self.id_number}/{self.irrep} {self.ee_type}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
            f"Oscillator strength (a.u.): {self.oscillator_strength:.6f}, \n"
        )
