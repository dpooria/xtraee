from xtraee.lazypattern import LP
from xtraee.transition import TMCC2Transition
from xtraee.trblock.trblock import TransitionBlock


class TMCC2TransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        ee_type: str = "",
        excitation_energy: float = 0.0,
        t1: float = 0.0,
        t2: float = 0.0,
    ):
        super().__init__(id_number, irrep, ee_type, excitation_energy)
        self.t1 = t1
        self.t2 = t2
        self.transitions: list[TMCC2Transition] = []
        self.tr_cls = TMCC2Transition
        self.tr_indicator = LP(
            r"^\s*\|\s*\d+\s+\w\s+\d+\s*\|\s*\d+\s+\w\s+\d+\s*\|\s*[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s*\|$"
        )
        self.end_trblock = "norm of printed elements:"

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"CC2 transition {self.id_number}/{self.irrep} {self.ee_type},\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV,\n"
            f"%t1: {self.t1}, %t2: {self.t2}.\n"
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
        )
