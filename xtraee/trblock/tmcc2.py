from xtraee.lazypattern import LP
from xtraee.transition import Transition, TMCC2Transition

from .trblock import TransitionBlock
from xtraee.utils import nan


class TMCC2TransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        multi: str = "",
        name="CC2",
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
        self.transitions: list[Transition] = []
        self.tr_cls = TMCC2Transition
        self.tr_indicator: str | LP = LP(
            r"^\s*\|\s*\d+\s+\w\s+\d+\s*\|\s*\d+\s+\w\s+\d+\s*\|\s*[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s*\|$"
        )
        self.end_trblock = "norm of printed elements:"
        self.meta_data.update({"%t1": nan, "%t2": nan})
