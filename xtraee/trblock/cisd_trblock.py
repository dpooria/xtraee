import re
from functools import partial

from xtraee.transition import CCSDTransition, CISDTransition

from .ccsd_trblock import CCSDTransitionBlock
from .trblock import TransitionBlock


class CISDTransitionBlock(TransitionBlock):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transitions_eomee: list[CCSDTransition] = []
        self.transitions: list[CISDTransition] = []
        self.homo = 0
        self.tr_cls = CISDTransition
        self.tr_indicator = re.compile(
            r"""^Root\s+\d+\s+Conv-d\s+yes\s+Tot\s+Ene=\s*
        (?P<tot>[+-]?\d+\.\d+)\s+hartree\s+
        \(Ex\s+Ene\s+(?P<ex>[+-]?\d+\.\d+)\s+eV\),\s+
        U0\^2=(?P<u0>[+-]?\d+\.\d+),\s+
        U1\^2=(?P<u1>[+-]?\d+\.\d+),\s+
        U2\^2=(?P<u2>[+-]?\d+\.\d+)\s+\|\|Res\|\|=
        (?P<res>[+-]?\d+\.\d+(?:[Ee][+-]?\d+))$
    """,
            re.VERBOSE,
        )

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

    def is_equal_eomee(self, other: CCSDTransitionBlock) -> bool:
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
    ) -> tuple[float, float, float]:
        is_eomee = isinstance(other, CCSDTransitionBlock)
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
