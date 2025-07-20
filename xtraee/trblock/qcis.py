from xtraee.lazypattern import LP, flp
from xtraee.transition import CISTransition
from xtraee.utils import Ha, nan

from .trblock import TransitionBlock

oscillator_pattern = LP(r"Strength\s+:\s*([-+]?\d*\.?\d+)\s*")
multplicity_pattern = LP(r"Multiplicity:\s*(Singlet|Triplet)")


class CISTransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        multi: str = "",
        name="CIS",
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
        self.transitions: list[CISTransition] = []
        self.tr_cls = CISTransition
        self.tr_indicator = "-->"
        self.end_trblock = "\n"

    def extras(self, line: str):
        if f"Total energy for state  {self.id_number}" in line:
            self.total_energy = float(flp.search(line).group(0)) * Ha
        if (m := oscillator_pattern.match(line)) is not None:
            self.oscillator_strength = float(m.group(1))
        elif (m := multplicity_pattern.match(line)) is not None:
            self.multi = m.group(1).lower()
