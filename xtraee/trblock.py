import re
from xtraee.transition import CCSDTransition, Transition, CISTransition
from abc import abstractmethod, ABC
from typing import List, Optional


class TransitionBlock(ABC):
    TrTYPE = Transition

    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        excitation: str = "",
        excitation_energy: float = 0.0,
        oscillator_strength: float = 0.0,
    ):
        self.transitions: List[Transition] = []
        self.id_number = id_number
        self.irrep = irrep
        self.excitation = excitation
        self.completed = False
        self.excitation_energy = excitation_energy
        self.identifier = f"{self.excitation}-{self.id_number}/{self.irrep}"
        self.oscillator_strength = oscillator_strength

    @abstractmethod
    def extras(self, line: str) -> None:
        raise NotImplementedError

    def add_data(self, line: str) -> bool:
        if self.END_TRANSITIONBLOCK in line:
            self.completed = True
            return True
        if self.TR_INDICATOR in line:
            self.transitions.append(self.TrTYPE.from_str(line))
        else:
            self.extras(line)
        return False

    def sort(self) -> None:
        self.transitions.sort(key=lambda t: abs(t.amplitude), reverse=True)

    def __getitem__(self, key) -> Transition:
        return self.transitions[key]


class EOMEETransitionBlock(TransitionBlock):
    EE_PATTERN = re.compile(r"^.*Excitation energy\s*=\s*([-+]?\d*\.?\d+)\s*eV\.\s*$")
    R_PATTERN = re.compile(
        r"^.*R0\^2\s*=\s*(\d*.\d+)\s*R1\^2\s*=\s*([-+]?\d*\.?\d+)\s*R2\^2\s*=\s*([-+]?\d*\.?\d+).*$"
    )
    END_TRANSITIONBLOCK = "Summary of significant orbitals:"
    TR_INDICATOR = "->"
    TrTYPE = CCSDTransition

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.R0 = 0.0
        self.R1 = 0.0
        self.R2 = 0.0
        self.omega = 0.0

    def extras(self, line: str) -> None:
        if (m := self.R_PATTERN.match(line)) is not None:
            self.R0 = float(m.group(1))
            self.R1 = float(m.group(2))
            self.R2 = float(m.group(3))
        elif (m := self.EE_PATTERN.match(line)) is not None:
            self.excitation_energy = float(m.group(1))

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"EOMEE transition {self.excitation} {self.id_number}/{self.irrep}\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV.\n"
            f"R0^2: {self.R0:.4f} R1^2: {self.R1:.4f} R2^2: {self.R2:.4f}\n"  # noqa
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
            f"Oscillator strength (a.u.): {self.oscillator_strength:.6f},"
            f"omega (Mulliken): {self.omega:.4f}\n"
        )


class CISTransitionBlock(TransitionBlock):
    MULTPLICITY_PATTERN = re.compile(r"Multiplicity:\s*(Singlet|Triplet)")
    END_TRANSITIONBLOCK = "\n"
    TR_INDICATOR = "-->"
    OSCILLATOR_PATTERN = re.compile(r"Strength\s+:\s*([-+]?\d*\.?\d+)\s*")
    TrTYPE = CISTransition

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transitions_eomee: List[CCSDTransition] = []

    def extras(self, line: str):
        if (m := self.OSCILLATOR_PATTERN.match(line)) is not None:
            self.oscillator_strength = float(m.group(1))
        elif (m := self.MULTPLICITY_PATTERN.match(line)) is not None:
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

    def generate_eomee(self, homo: int) -> List[CCSDTransition]:
        self.transitions_eomee = list(
            map(
                lambda t: t.to_ccsd(homo),
                self.transitions,
            )
        )
        return self.transitions_eomee

    def is_equal_eomee(self, other: EOMEETransitionBlock, homo: int) -> bool:
        if len(self.transitions_eomee) != len(self.transitions):
            self.generate_eomee(homo)
        for o_tr in other.transitions:
            is_in = False
            for tr in self.transitions_eomee:
                if tr.is_equal(o_tr):
                    is_in = True
                    break
            if not is_in:
                return False
        return True
