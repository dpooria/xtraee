from xtraee.lazypattern import LP
from xtraee.transition import CISDTransition

from .trblock import TransitionBlock, nan


class CISDTransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        multi: str = "",
        name="CISD",
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
        self.meta_data.update({"U0": nan, "U1": nan, "U2": nan})
