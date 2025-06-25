from typing import List, Tuple

from xtraee.lazypattern import LP, VERBOSE


class Transition:
    name = "Transition"
    pattern = LP(r".+")

    def __init__(
        self,
        amplitude: float,
        initial: str,
        final: str,
    ):
        self.amplitude = amplitude
        self.probability = amplitude**2
        self.initial = initial.strip()
        self.final = final.strip()
        self.id_i: List[Tuple[str, ...]] = []
        self.id_f: List[Tuple[str, ...]] = []
        self._spin_ind = 2
        self.is_double = False
        if not self.parse():
            raise ValueError(f"Cannot parse {initial} -> {final}")

    @classmethod
    def from_str(cls, line: str):
        raise NotImplementedError

    def parse(
        self,
    ) -> bool:
        if m_l := self.pattern.findall(self.initial):
            self.id_i = m_l
        else:
            return False

        if m_l := self.pattern.findall(self.final):
            self.id_f = m_l
        else:
            return False
        if len(self.id_i) > 1 or len(self.id_f) > 1:
            self.is_double = True
        else:
            self.is_double = False
        return True

    def is_equal(self, other, check_amp=False, check_spin=False) -> bool:
        if not isinstance(other, Transition):
            return NotImplemented
        if self.is_double != other.is_double:
            return False
        if check_spin:
            for id_i, other_id_i in zip(self.id_i, other.id_i):
                if id_i != other_id_i:
                    return False
            for id_f, other_id_f in zip(self.id_f, other.id_f):
                if id_f != other_id_f:
                    return False
        else:
            for id_i, other_id_i, id_f, other_id_f in zip(
                self.id_i, other.id_i, self.id_f, other.id_f
            ):
                if id_i[: self._spin_ind] != other_id_i[: self._spin_ind]:
                    return False
                if id_f[: self._spin_ind] != other_id_f[: self._spin_ind]:
                    return False
        return True

    def to_std(self, homo: int = 0) -> "CCSDTransition":
        return self

    def __eq__(self, other) -> bool:
        return self.is_equal(other)

    def __repr__(self) -> str:
        return f"{self.amplitude:.4f}\t {self.initial} -> {self.final}"


# CCSDTransition is the standard format
class CCSDTransition(Transition):
    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CCSDTransition"
        self.pattern = LP(r"\s*(\d+)\s*\(([^\s]+)\)\s*(\w*)\s*")
        super().__init__(amplitude, initial, final)

    @classmethod
    def from_str(cls, line: str):
        s = line.split("->")
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        lhs = s[0].strip()
        rhs = s[1].strip()
        # match with the float
        if (m := LP(r"([-+]?\d+\.\d+)").match(lhs)) is not None:
            amp_str = m.group(1)
            amplitude = float(amp_str)
            # remove the matched part
            idx = lhs.index(amp_str) + len(amp_str)
            lhs = lhs[idx:].strip()
        else:
            print(f"Cannot match amplitude {line}")
            amplitude = 0.0
        return cls(amplitude, lhs, rhs)


class CISTransition(Transition):
    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CISTransition"
        self.pattern = LP(r"\s*([^\s]+)\s*\(\s*(\d+)\s*\)\s*")
        super().__init__(amplitude, initial, final)

    @classmethod
    def from_str(cls, line: str):
        s = line.split("-->")
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        lhs = s[0].strip()
        rhs = s[1].strip()
        idx = rhs.index("=")
        amplitude = float(rhs[idx + 1 :])
        idx = rhs.find("amplitude")
        rhs = rhs[:idx]
        return cls(amplitude, lhs, rhs)

    def to_std(self, homo: int) -> CCSDTransition:
        # not considering symmetry and there is no double excitation in CIS
        initial = f"{int(self.id_i[0][1])} (A)"
        final = f"{int(self.id_f[0][1]) + homo} (A)"
        return CCSDTransition(self.amplitude, initial, final)


class CISDTransition(Transition):
    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CISDTransition"
        self.pattern = LP(r"\s*(\d+)\s*\(\s*([^\s]+)\s*\)\s*([AB])\s*")
        super().__init__(amplitude, initial, final)

    @classmethod
    def from_str(cls, line: str):
        s = line.split("->")
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        lhs = s[0].strip()
        rhs = s[1].strip()
        amp_str = lhs.split()[0]
        amplitude = float(amp_str)
        lhs = lhs[lhs.find(amp_str) + len(amp_str) :].lstrip()
        return cls(amplitude, lhs, rhs)

    def to_ccsd(self, homo: int, symmetry: bool) -> CCSDTransition:
        initials = []
        finals = []
        # CISD orbitals numbering starts from 0
        for id_i, id_f in zip(self.id_i, self.id_f):
            if symmetry:
                initial = f"{int(id_i[0]) + 1} ({id_i[1]}) {id_i[2]}"
                final = f"{int(id_f[0]) + homo + 1} ({id_f[1]}) {id_f[2]}"
            else:
                initial = f"{int(id_i[0]) + 1} (A) {id_i[2]})"
                final = f"{int(id_f[0]) + homo + 1} (A) {id_f[2]}"
            initials.append(initial)
            finals.append(final)
        return CCSDTransition(self.amplitude, "\t".join(initials), "\t".join(finals))

    def to_std(self, homo):
        return self.to_ccsd(homo, False)


class TMCC2Transition(Transition):
    amppattern = LP(r"\s*([+-]?\d+\.\d+)\s+[+-]?\d+\.\d+\s*")

    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CC2Transition"
        self.pattern = LP(
            r"""
            ^\s*
            (\d+)                # group 1: the first integer
            \s+
            (?=[ab]\s+(\d+))     # LOOKAHEAD: assert “letter + spaces + digits” next,
                                 # and capture that digit as group 2
            ([ab])               # group 3: the letter
            \s+
            \d+                  # match the digit again, but don’t capture it
            \s*$
        """,
            VERBOSE,
        )
        super().__init__(amplitude, initial, final)
        self._spin_ind = 1

    @classmethod
    def from_str(cls, line: str):
        s = line.split("|")
        if len(s) < 4:
            raise ValueError(f"cannot match {line}")
        lhs = s[1].strip()
        rhs = s[2].strip()
        if (m := TMCC2Transition.amppattern.match(s[3].strip())) is not None:
            amplitude = float(m.group(1))
        else:
            raise ValueError(f"could not match the amplitude {line}")
        return cls(amplitude, lhs, rhs)

    def to_std(self, homo: int = 0) -> CCSDTransition:
        # the double transitions is not tested
        initials = []
        finals = []
        for id_i, id_f in zip(self.id_i, self.id_f):
            initial = f"{int(id_i[0])} (A) {id_i[2].upper()}"
            final = f"{int(id_f[0])} (A) {id_f[2].upper()}"
        initials.append(initial)
        finals.append(final)
        return CCSDTransition(self.amplitude, "\t".join(initials), "\t".join(finals))
