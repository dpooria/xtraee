import re
from abc import abstractmethod
from typing import List, Tuple

from xtraee.transition import CCSDTransition, CISTransition, Transition


class TransitionBlock:
    TrTYPE = Transition
    END_TRANSITIONBLOCK = "\n"
    TR_INDICATOR = "->"

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
        self.oscillator_strength = oscillator_strength

    @property
    def identifier(self) -> str:
        return f"{self.excitation}-{self.id_number}/{self.irrep}"

    @abstractmethod
    def extras(self, line: str) -> None:
        raise NotImplementedError

    @property
    def squared_sum(self) -> float:
        return sum(t.amplitude**2 for t in self.transitions)

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

    @abstractmethod
    def compare(self, other) -> bool:
        raise NotImplementedError

    def _compare_acc(
        self, other, transitions: List[Transition]
    ) -> Tuple[float, float, float]:
        amps = []
        other_sq_sum = other.squared_sum
        my_sq_sum = self.squared_sum
        N_tr = min(len(other.transitions), len(self.transitions))
        for o_tr in other.transitions:
            for tr in transitions:
                if tr.is_equal(o_tr):
                    amps.append(
                        (o_tr.amplitude**2 / other_sq_sum, tr.amplitude**2 / my_sq_sum)
                    )
                    break
        N_pos = len(amps)
        if N_pos == 0:
            return (0, 0, 0)
        mae = sum([abs(o_amp - m_amp) for o_amp, m_amp in amps]) / N_pos
        acc = 1.0 - mae
        return N_pos / N_tr * acc, acc, N_pos / N_tr

    def _compare_err(
        self, other, transitions: List[Transition]
    ) -> Tuple[float, float, float]:
        amps = []
        N_tr = min(len(other.transitions), len(self.transitions))
        for o_tr in other.transitions:
            for tr in transitions:
                if tr.is_equal(o_tr):
                    amps.append((o_tr.amplitude**2, tr.amplitude**2))
                    break
        N_pos = len(amps)
        if N_pos == 0:
            return (None, None, 0)
        mae = sum([abs(o_amp - m_amp) for o_amp, m_amp in amps]) / N_pos
        return N_tr / N_pos * mae, mae, N_pos / N_tr

    def _compare_err_tot(
        self, other, transitions: List[Transition]
    ) -> Tuple[float, float, float]:
        amps = []
        N_tr = transitions
        N_pos = 0
        other_sq_sum = other.squared_sum
        my_sq_sum = self.squared_sum
        for tr in transitions:
            matched = False
            for o_tr in other.transitions:
                if tr.is_equal(o_tr):
                    amps.append((tr.amplitude**2 / my_sq_sum, o_tr.amplitude**2 / other_sq_sum))
                    matched = True
                    N_pos += 1
                    break
            if not matched:
                amps.append((tr.amplitude**2 / my_sq_sum, 0))

        mae = sum([abs(m_amp - o_amp) for m_amp, o_amp in amps]) / N_tr
        acc = 1 - mae
        return acc, mae, N_pos / N_tr

    def _compare_error_self_ref(
        self, other, transitions: List[Transition]
    ) -> Tuple[float, float, float]:
        amps = []
        N_tr = len(transitions)
        N_pos = 0
        my_sq_sum = self.squared_sum
        # I am the reference :)
        for tr in transitions:
            matched = False
            for o_tr in other.transitions:
                if tr.is_equal(o_tr):
                    amps.append((tr.amplitude**2, o_tr.amplitude**2))
                    matched = True
                    N_pos += 1
                    break
            if not matched:
                amps.append((tr.amplitude**2, 0))

        mae = sum([m_amp * abs(m_amp - o_amp) / my_sq_sum for m_amp, o_amp in amps])
        acc = 1.0 - mae
        return acc, mae, N_pos / N_tr


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

    def compare(
        self, other: TransitionBlock, method: str
    ) -> Tuple[float, float, float]:
        if isinstance(other, CISTransitionBlock):
            return other.compare(self, method)
        else:
            if method == "acc":
                return self._compare_acc(other, self.transitions)
            elif method == "err_tot":
                return self._compare_err_tot(other, self.transitions)
            elif method == "err":
                return self._compare_err(other, self.transitions)
            elif method == "err_self_ref":
                return self._compare_error_self_ref(other, self.transitions)
            else:
                raise ValueError(f"Method not recognized {method}")


class CISTransitionBlock(TransitionBlock):
    MULTPLICITY_PATTERN = re.compile(r"Multiplicity:\s*(Singlet|Triplet)")
    END_TRANSITIONBLOCK = "\n"
    TR_INDICATOR = "-->"
    OSCILLATOR_PATTERN = re.compile(r"Strength\s+:\s*([-+]?\d*\.?\d+)\s*")
    TrTYPE = CISTransition

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transitions_eomee: List[CCSDTransition] = []
        self.homo = 0

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

    def generate_eomee(self) -> List[CCSDTransition]:
        self.transitions_eomee = [tr.to_ccsd(self.homo) for tr in self.transitions]
        return self.transitions_eomee

    def is_equal_eomee(self, other: EOMEETransitionBlock) -> bool:
        if len(self.transitions_eomee) != len(self.transitions):
            self.generate_eomee()
        for o_tr in other.transitions:
            is_in = False
            for tr in self.transitions_eomee:
                if tr.is_equal(o_tr):
                    is_in = True
                    break
            if not is_in:
                return False

    def compare(
        self, other: TransitionBlock, method: str
    ) -> Tuple[float, float, float]:
        if isinstance(other, EOMEETransitionBlock):
            if len(self.transitions_eomee) != len(self.transitions):
                self.generate_eomee()
            transitions = self.transitions_eomee
            if method == "err_self_ref" or method == "err_tot":
                tr_bkp = self.transitions.copy()
                self.transitions = transitions
                res = other._compare_error_self_ref(self, other.transitions)
                self.transitions = tr_bkp
                return res

        else:
            transitions = self.transitions
        if method == "acc":
            return self._compare_acc(other, transitions)
        elif method == "err":
            return self._compare_err(other, transitions)
        elif method == "err_self_ref":
            return self._compare_error_self_ref(other, transitions)
        elif method == "err_tot":
            return self._compare_err_tot(other, transitions)
        else:
            raise ValueError(f"Method not recognized {method}")
