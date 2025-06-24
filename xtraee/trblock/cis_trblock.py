import re
from functools import partial

from xtraee.transition import CISTransition

from .ccsd_trblock import CCSDTransitionBlock
from .trblock import TransitionBlock

oscillator_pattern = re.compile(r"Strength\s+:\s*([-+]?\d*\.?\d+)\s*")
multplicity_pattern = re.compile(r"Multiplicity:\s*(Singlet|Triplet)")


class CISTransitionBlock(TransitionBlock):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transitions: list[CISTransition] = []
        self.tr_cls = CISTransition
        self.tr_indicator = "-->"
        self.end_trblock = "\n"

    def extras(self, line: str):
        if (m := oscillator_pattern.match(line)) is not None:
            self.oscillator_strength = float(m.group(1))
        elif (m := multplicity_pattern.match(line)) is not None:
            self.excitation = m.group(1).lower()

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"CIS transition {self.id_number}/{self.irrep} {self.excitation}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
            f"Oscillator strength (a.u.): {self.oscillator_strength:.6f}, \n"
        )
