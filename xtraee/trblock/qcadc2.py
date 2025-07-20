from xtraee.lazypattern import LP, flp
from xtraee.transition import ADC2Transition
from xtraee.utils import Ha, nan

from .trblock import TransitionBlock

tr_indicator = LP(r"(\s*\d+\s+\([\w'\"]+\)\s+[AB]?\s*)+[-+]?\d+\.\d+")


class ADC2TransitionBlock(TransitionBlock):
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        multi: str = "",
        name="ADC2",
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
        self.transitions: list[ADC2Transition] = []
        self.tr_cls = ADC2Transition
        self.tr_indicator = tr_indicator
        self.end_trblock = "\n"
        self.meta_data.update({"V1^2": nan, "V2^2": nan})

    def extras(self, line: str):
        if "Total energy:" in line:
            self.total_energy = float(flp.search(line).group(0)) * Ha
        elif "Excitation energy:" in line:
            self.excitation_energy = float(flp.search(line).group(0))
        elif "V1^2" in line:
            for i, m in enumerate(flp.finditer(line)):
                if i == 0:
                    self.meta_data["V1^2"] = float(m.group(0))
                elif i == 1:
                    self.meta_data["V2^2"] = float(m.group(0))
        elif "Osc. strength:" in line:
            self.oscillator_strength = float(flp.search(line).group(0))
