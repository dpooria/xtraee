import logging
import re
from abc import abstractmethod
from copy import deepcopy
from typing import List, Tuple
from functools import partial

from xtraee.config import debug
from xtraee.transition import CCSDTransition, CISTransition, Transition, CC2Transition


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
        self.completed_extras = False
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
        if isinstance(self.TR_INDICATOR, re.Pattern):
            if self.TR_INDICATOR.search(line):
                self.transitions.append(self.TrTYPE.from_str(line))
        elif self.TR_INDICATOR in line:
            self.transitions.append(self.TrTYPE.from_str(line))
        else:
            self.extras(line)
        return False

    def sort(self) -> None:
        self.transitions.sort(key=lambda t: abs(t.amplitude), reverse=True)

    def __getitem__(self, key) -> Transition:
        return self.transitions[key]

    def _compare_acc(
        self,
        transitions: List[Transition],
        other_transitions: List[Transition],
        check_spin: bool = False,
    ) -> Tuple[float, float, float]:
        probs = []
        N_tr = len(transitions)
        N_pos = 0
        my_sum = sum([tr.probability for tr in transitions])
        # I am the reference :)
        for tr in transitions:
            matched = False
            for o_tr in other_transitions:
                if tr.is_equal(o_tr, check_spin=check_spin):
                    probs.append((tr.probability, o_tr.probability))
                    matched = True
                    N_pos += 1
                    break
            if not matched:
                probs.append((tr.probability, 0.0))

        mae = sum([m_prob * abs(m_prob - o_prob) / my_sum for m_prob, o_prob in probs])
        acc = 1.0 - mae
        return acc, mae, N_pos / N_tr

    def _compare_pearson(
        self,
        transitions: List[Transition],
        other_transitions: List[Transition],
        check_spin: bool = False,
    ) -> Tuple[float, float, float]:
        import numpy as np

        probs = []
        N_tr = len(transitions)
        N_pos = 0
        for tr in transitions:
            for o_tr in other_transitions:
                if tr.is_equal(o_tr, check_spin=check_spin):
                    probs.append((tr.probability, o_tr.probability))
                    N_pos += 1
                    break
            else:
                probs.append((tr.probability, 0.0))
        probs = np.array(probs)
        if np.sum(probs[:, 1]) == 0.0:
            acc = 0.0
        else:
            acc = np.corrcoef(probs.T)[0, 1]
        return acc, 1 - acc**2, N_pos / N_tr

    def _compare_innerprod(
        self,
        transitions: List[Transition],
        other_transitions: List[Transition],
        retreive: bool = False,
        check_spin: bool = False,
    ) -> Tuple[float, float, float]:
        amps = []
        N_tr = len(transitions)
        N_pos = 0
        # I am the reference
        for tr in transitions:
            matched = False
            for o_tr in other_transitions:
                if tr.is_equal(o_tr, check_spin=check_spin):
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

    def compare(
        self, other: "TransitionBlock", method: str, check_spin: bool = False
    ) -> Tuple[float, float, float]:
        if method == "1":
            return self._compare_acc(self.utrs, other.utrs, check_spin=check_spin)
        elif method == "2":
            return self._compare_innerprod(self.utrs, other.utrs, check_spin=check_spin)
        elif method == "3":
            return self._compare_innerprod(
                self.utrs, other.utrs, True, check_spin=check_spin
            )
        elif method == "4":
            return self._compare_pearson(self.utrs, other.utrs, check_spin=check_spin)
        else:
            raise ValueError(f"Method not recognized {method}")


class EOMEETransitionBlock(TransitionBlock):
    EE_PATTERN = re.compile(r"^.*Excitation energy\s*=\s*([-+]?\d*\.?\d+)\s*eV\.\s*$")
    R_PATTERN = re.compile(
        r"^.*R0\^2\s*=\s*(\d*.\d+)\s*R1\^2\s*=\s*([-+]?\d*\.?\d+)\s*R2\^2\s*=\s*([-+]?\d*\.?\d+).*$"
    )

    END_TRANSITIONBLOCK = "Summary of significant orbitals:"
    TR_INDICATOR = "->"
    OCCUPATION_FRONTIER_NO = "Occupation of frontier NOs:"
    FRONTIER_NO_PATTERN = re.compile(
        r"^\s*([-+]?[0-9]+\.[0-9]+)\s+([-+]?[0-9]+\.[0-9]+)\s*$"
    )
    UNPAIRED_NO_PATTERN = re.compile(
        r"^\s*Number of unpaired electrons:\s*n_u\s*=\s*([-+]?[0-9]+\.[0-9]+),\s*n_u,nl\s*=\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )
    UNPAIRED_NO = "Number of unpaired electrons:"
    PARTICIPATION_RATIO_NO = "NO participation ratio (PR_NO):"
    PRNO_PATTERN = re.compile(
        r"^\s*NO participation ratio \(PR_NO\):\s*([-+]?[0-9]+\.[0-9]+)\s*$"
    )

    TrTYPE = CCSDTransition

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transitions: List[CCSDTransition] = []
        self.R0 = 0.0
        self.R1 = 0.0
        self.R2 = 0.0
        self.gamma = 0.0
        self.omega = 0.0
        self.loc = 0.0
        self.phe = 0.0
        self.rhre = 0.0
        self.alphabeta = 0.0
        self.corr_coef = 0.0
        self.froniter_no = []
        self.wait_frontier_no = False
        self.nu = 0.0
        self.nl = 0.0
        self.prno = 0.0

    def extras(self, line: str) -> None:
        if not self.completed:
            if (m := self.R_PATTERN.match(line)) is not None:
                self.R0 = float(m.group(1))
                self.R1 = float(m.group(2))
                self.R2 = float(m.group(3))
            elif (m := self.EE_PATTERN.match(line)) is not None:
                self.excitation_energy = float(m.group(1))
        elif not self.completed_extras:
            if self.OCCUPATION_FRONTIER_NO in line:
                self.wait_frontier_no = True
            elif self.wait_frontier_no:
                if (m := self.FRONTIER_NO_PATTERN.match(line)) is not None:
                    self.wait_frontier_no = False
                    self.froniter_no = [float(m.group(1)), float(m.group(2))]
                else:
                    raise ValueError(
                        f"Could not match the Occupations of the frontier NO for {self}"
                    )
            elif self.UNPAIRED_NO in line:
                if m := self.UNPAIRED_NO_PATTERN.match(line):
                    self.nu = m.group(1)
                    self.nl = m.group(2)
                else:
                    raise ValueError(
                        f"Could not match the number of unpaired electrons for {self}"
                    )
            elif self.PARTICIPATION_RATIO_NO in line:
                if (m := self.PRNO_PATTERN.match(line)) is not None:
                    self.prno = m.group(1)
                else:
                    raise ValueError(
                        f"Could not match the participation number for {self}"
                    )
                self.completed_extras = True
        else:
            raise ValueError(f"Transition {self} is already completed!")

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
        if isinstance(other, EOMEETransitionBlock):
            return super().compare(other, method)
        else:
            return other.compare(self, method)


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
            comp = partial(self._compare_innerprod, retreive=True)
        elif method == "4":
            comp = self._compare_pearson
        else:
            raise ValueError(f"Method not recognized {method}")
        if is_eomee:
            return comp(
                other.utrs,
                transitions,
            )
        else:
            return comp(transitions, other.utrs)


class CC2TransitionBlock(TransitionBlock):
    END_TRANSITIONBLOCK = "norm of printed elements:"
    TR_INDICATOR = re.compile(
        r"^\s*\|\s*\d+\s+\w\s+\d+\s*\|\s*\d+\s+\w\s+\d+\s*\|\s*[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s*\|$"
    )
    TrTYPE = CC2Transition

    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        excitation: str = "",
        excitation_energy: float = 0.0,
        t1: float = 0.0,
        t2: float = 0.0,
    ):
        super().__init__(id_number, irrep, excitation, excitation_energy)
        self.t1 = t1
        self.t2 = t2
        self.transitions: List[CC2Transition] = []
        self.homo = 0

    def extras(self, line: str):
        pass

    def __repr__(self) -> str:
        line = "\n".join(map(str, self.transitions))
        return (
            f"CC2 transition {self.id_number}/{self.irrep} {self.excitation},\n"  # noqa
            f"EE: {self.excitation_energy:.4f} eV,\n"
            f"%t1: {self.t1}, %t2: {self.t2}.\n"
            "Amplitude Transitions between orbitals\n"
            f"{line}\n"
        )

    def compare(
        self,
        other: TransitionBlock,
        method: str,
        check_spin=True,
    ) -> Tuple[float, float, float]:
        assert isinstance(other, CC2TransitionBlock), "Not implemented"
        return super().compare(other, method, check_spin)
