from xtraee.lazypattern import LP
from xtraee.transition import Transition, CISDTransition
from xtraee.utils import Ha, nan
from .trblock import TransitionBlock

# for CIS(D)
cis_d_exen_pattern = LP(r"E_ex\s*=\s*([-+]?\d+\.\d+)\s*eV")
cis_d_toten_pattern = LP(r"E_CIS\(D\)\s*=\s*([-+]?\d+\.\d+)\s+hartree")


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
        self.transitions: list[Transition] = []
        self.tr_cls = CISDTransition
        self.end_trblock = "\n"
        # -0.6395                 12(   A1) B   ->    0(   A1) B
        self.tr_indicator: str | LP = LP(
            r"^\s*[-+]?\d+\.\d+"
            r"\s+\d+\s*\(.+\)\s+[AB]"
            r"\s*->\s*"
            r"\d+\s*\(.+\)\s*[AB]\s*$"
        )
        self.meta_data.update({"U0": nan, "U1": nan, "U2": nan})

    def extras(self, line):
        if self.name == "CIS_D_":
            if f"Root {self.id_number} CIS(D) corr=" in line:
                if (m := cis_d_exen_pattern.search(line)) is not None:
                    self.excitation_energy = float(m.group(1))
                if (m := cis_d_toten_pattern.search(line)) is not None:
                    self.total_energy = float(m.group(1)) * Ha
