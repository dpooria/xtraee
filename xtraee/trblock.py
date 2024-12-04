import logging
import re
from abc import abstractmethod
from copy import deepcopy
from typing import List, Tuple

from xtraee.config import debug
from xtraee.transition import CCSDTransition, CISTransition, Transition


def lsq_fit(
    transitions: List[Transition],
    other_transitions: List[Transition],
    amps: List[Tuple[float, float]],
) -> None:
    import numpy as np
    from scipy.optimize import minimize

    residu = 1.0 - sum(
        [tr.probability for tr in transitions]
    )  # the residual transitions
    o_nm_amps = []  # not matched amplitudes
    for o_tr in other_transitions:
        if o_tr not in transitions:
            o_nm_amps.append(abs(o_tr.amplitude))
    o_nm_amp = np.asarray(o_nm_amps)

    # maximize dot product
    def objective(V):
        return -np.dot(o_nm_amp, np.abs(V))

    # the remaining probabilities should be equal to the residual
    def constraint(V):
        return np.dot(V, V) - residu

    con = {"type": "eq", "fun": constraint}
    result = minimize(
        objective,
        np.zeros_like(o_nm_amps),
        constraints=con,
        method="SLSQP",
    )
    if not result.success:
        logging.debug("Unable to retrieve the missing transitions", result)
        return
    for m_amp, o_amps in zip(result.x, o_nm_amps):
        amps.append((abs(m_amp), o_amps))


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

    @property
    def utrs(self) -> List[Transition]:
        ut: List[Transition] = []
        for tr in self.transitions:
            if tr in ut:
                new_tr = deepcopy(ut[ut.index(tr)])
                new_tr.probability += tr.probability
                new_tr.amplitude = new_tr.probability**0.5
                ut[ut.index(tr)] = new_tr
            else:
                ut.append(tr)

        if debug:
            for tr in ut:
                assert tr.probability >= 0.0
                assert tr.probability <= 1.0
                assert abs(tr.amplitude - tr.probability**0.5) < 1e-6
            assert sum([tr.probability for tr in ut]) <= 1.0
        return ut

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

    @abstractmethod
    def compare(self, other, method: str) -> Tuple[float, float, float]:
        raise NotImplementedError

    def _compare_acc(
        self,
        transitions: List[Transition],
        other_transitions: List[Transition],
    ) -> Tuple[float, float, float]:
        probs = []
        N_tr = len(transitions)
        N_pos = 0
        my_sum = sum([tr.probability for tr in transitions])
        # I am the reference :)
        for tr in transitions:
            matched = False
            for o_tr in other_transitions:
                if tr.is_equal(o_tr):
                    probs.append((tr.probability, o_tr.probability))
                    matched = True
                    N_pos += 1
                    break
            if not matched:
                probs.append((tr.probability, 0.0))

        mae = sum([m_prob * abs(m_prob - o_prob) / my_sum for m_prob, o_prob in probs])
        acc = 1.0 - mae
        return acc, mae, N_pos / N_tr

    def _compare_innerprod(
        self,
        transitions: List[Transition],
        other_transitions: List[Transition],
        retreive: bool = False,
    ) -> Tuple[float, float, float]:
        amps = []
        N_tr = len(transitions)
        N_pos = 0
        # I am the reference
        for tr in transitions:
            matched = False
            for o_tr in other_transitions:
                if tr.is_equal(o_tr):
                    amps.append((tr.amplitude, o_tr.amplitude))
                    matched = True
                    N_pos += 1
                    break
            if not matched:
                amps.append((tr.amplitude, 0.0))

        if retreive and len(transitions) < len(other_transitions):
            # construct the missing transitions that maximize the accuracy
            lsq_fit(transitions, other_transitions, amps)

        acc = sum([abs(m_amp) * abs(o_amp) for m_amp, o_amp in amps])
        mae = 1.0 - acc
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
        self.transitions: List[CCSDTransition] = []
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
            if method == "1":
                return self._compare_acc(self.utrs, other.utrs)
            elif method == "2":
                return self._compare_innerprod(self.utrs, other.utrs)
            elif method == "3":
                return self._compare_innerprod(self.utrs, other.utrs, True)
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
        self.transitions: List[CISTransition] = []
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
        return True

    def compare(
        self, other: TransitionBlock, method: str
    ) -> Tuple[float, float, float]:
        is_eomee = isinstance(other, EOMEETransitionBlock)
        if is_eomee:
            if len(self.transitions_eomee) != len(self.transitions):
                self.generate_eomee()
            transitions = self.transitions_eomee
        else:
            transitions = self.transitions
        if method == "1":
            comp = self._compare_acc
        elif method == "2":
            comp = self._compare_innerprod
        elif method == "3":
            comp = lambda tr, o_tr: self._compare_innerprod(tr, o_tr, True)
        else:
            raise ValueError(f"Method not recognized {method}")
        if is_eomee:
            return comp(
                other.utrs,
                transitions,
            )
        else:
            return comp(transitions, other.utrs)
