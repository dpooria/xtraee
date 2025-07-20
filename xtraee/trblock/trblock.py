from copy import deepcopy
from functools import partial

from xtraee.config import debug, get_logger
from xtraee.lazypattern import LP
from xtraee.transition import Transition
from xtraee.utils import nan

log = get_logger("trblock", debug)


def lsq_fit(
    transitions: list[Transition],
    other_transitions: list[Transition],
    amps: list[tuple[float, float]],
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
        log.debug("Unable to retrieve the missing transitions", result)
        return
    for m_amp, o_amps in zip(result.x, o_nm_amps):
        amps.append((abs(m_amp), o_amps))


class TransitionBlock:
    def __init__(
        self,
        id_number: int,
        irrep: str = "",
        multi: str = "",
        name="Transition",
        excitation_energy: float = nan,
        total_energy: float = nan,
        oscillator_strength: float = nan,
    ):
        self.name = name
        self.transitions: list[Transition] = []
        self.std_transitions: list[Transition] = []
        self.id_number = id_number
        self.irrep = irrep
        self.multi = multi
        self.completed = False
        self.completed_extras = False
        self.homo = 0
        self.excitation_energy = excitation_energy
        self.total_energy = total_energy
        self.oscillator_strength = oscillator_strength
        self.tr_cls = Transition
        self.tr_indicator = "->"
        # this means the transition never gets terminated as the line is stripped
        self.end_trblock = "\n"
        self.meta_data = {}

    @property
    def identifier(self) -> str:
        return f"{self.multi}-{self.id_number}/{self.irrep}"

    def utrs(self, transitions: list[Transition] | None = None) -> list[Transition]:
        if transitions is None:
            transitions = self.std_transitions
        ut: list[Transition] = []
        for tr in transitions:
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
                assert abs(tr.amplitude) - tr.probability**0.5 < 1e-6
            assert sum([tr.probability for tr in ut]) <= 1.0
        return ut

    def extras(self, line: str) -> None:
        pass

    def add_data(self, line: str) -> bool:
        if self.end_trblock in line:
            self.completed = True
            return True
        if isinstance(self.tr_indicator, LP):
            if self.tr_indicator.search(line):
                self.transitions.append(self.tr_cls.from_str(line))
                return False
        elif self.tr_indicator in line:
            self.transitions.append(self.tr_cls.from_str(line))
            return False
        self.extras(line)
        return False

    def sort(self) -> None:
        self.transitions.sort(key=lambda t: abs(t.amplitude), reverse=True)

    def __repr__(self) -> str:
        transitions = "\n".join(f"  {t}" for t in self.transitions)
        meta = ", ".join([f"{k}: {v}" for k, v in self.meta_data.items()])
        return (
            f"{self.name} transition {self.id_number}/{self.irrep} {self.multi},\n"
            f"Excitation energy: {self.excitation_energy:>9.4f} eV\n"
            f"Total energy: {self.total_energy:>9.6f} eV\n"
            f"Oscillator strenght: {self.oscillator_strength:>9.6f} (a.u.)\n"
            f"{meta}\n"
            "Amplitude Transitions between orbitals\n"
            f"{transitions}\n"
        )

    def __getitem__(self, key) -> Transition:
        return self.transitions[key]

    def _compare_wabs(
        self,
        transitions: list[Transition],
        other_transitions: list[Transition],
        shallow: bool = False,
    ) -> tuple[float, float, float]:
        probs = []
        N_tr = len(transitions)
        N_pos = 0
        my_sum = sum([tr.probability for tr in transitions])
        # I am the reference :)
        for tr in transitions:
            matched = False
            for o_tr in other_transitions:
                if tr.is_equal(o_tr, shallow=shallow):
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
        transitions: list[Transition],
        other_transitions: list[Transition],
        shallow: bool = False,
    ) -> tuple[float, float, float]:
        import numpy as np

        probs = []
        N_tr = len(transitions)
        N_pos = 0
        for tr in transitions:
            for o_tr in other_transitions:
                if tr.is_equal(o_tr, shallow=shallow):
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
        transitions: list[Transition],
        other_transitions: list[Transition],
        retreive: bool = False,
        shallow: bool = False,
    ) -> tuple[float, float, float]:
        amps = []
        N_tr = len(transitions)
        N_pos = 0
        # I am the reference
        for tr in transitions:
            matched = False
            for o_tr in other_transitions:
                if tr.is_equal(o_tr, shallow=shallow):
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
        self, other: "TransitionBlock", method: str, shallow: bool = True
    ) -> tuple[float, float, float]:
        if method == "w-abs":
            comp = self._compare_wabs
        elif method == "inner-prod":
            comp = self._compare_innerprod
        elif method == "innerp-retrieve":
            comp = partial(self._compare_innerprod, retreive=True)
        elif method == "pearson":
            comp = self._compare_pearson
        else:
            raise ValueError(f"Method not recognized {method}")
        if not self.std_ready():
            self.generate_std()
        if not other.std_ready():
            other.generate_std()
        if len(self.std_transitions) == 0 or len(other.std_transitions) == 0:
            return 0.0, 1.0, 0.0
        return comp(self.utrs(), other.utrs(), shallow=shallow)

    def generate_std(self) -> list[Transition]:
        self.std_transitions = [tr.to_std(self.homo) for tr in self.transitions]
        return self.std_transitions

    def std_ready(self) -> bool:
        return len(self.std_transitions) == len(self.transitions)

    def is_equal_std(self, other: "TransitionBlock") -> bool:
        if not self.std_ready():
            self.generate_std()
        if not other.std_ready():
            other.generate_std()

        if len(self.std_transitions) == 0 or len(other.std_transitions) == 0:
            log.warning("The transition block is empty!")
            return False

        for o_tr in other.std_transitions:
            for tr in self.std_transitions:
                if tr.is_equal(o_tr):
                    break
            else:
                return False
        return True
