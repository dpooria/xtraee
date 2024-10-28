import re
from typing import List, Tuple


class Transition:
    PATTERN = re.compile(r".+")

    def __init__(
        self,
        amplitude: float,
        initial: str,
        final: str,
    ):
        self.amplitude = amplitude
        self.probability = amplitude ** 2
        self.initial = initial.strip()
        self.final = final.strip()
        self.id_i: List[Tuple[str]] = []
        self.id_f: List[Tuple[str]] = []
        self.is_double = False
        if not self.parse():
            raise ValueError(f"Cannot parse {initial} -> {final}")

    @classmethod
    def from_str(cls, line: str):
        raise NotImplementedError

    def parse(
        self,
    ) -> bool:
        if m_l := self.PATTERN.findall(self.initial):
            self.id_i = m_l
        else:
            return False

        if m_l := self.PATTERN.findall(self.final):
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
            for id_i, other_id_i, id_f, other_id_f in zip(self.id_i, other.id_i,
                                                          self.id_f, other.id_f):
                if id_i[:2] != other_id_i[:2]:
                    return False
                if id_f[:2] != other_id_f[:2]:
                    return False
        return True

    def __eq__(self, other) -> bool:
        return self.is_equal(other, check_spin=False)

    def __repr__(self) -> str:
        return f"{self.amplitude:.4f}\t {self.initial} -> {self.final}"


class CCSDTransition(Transition):
    PATTERN = re.compile(r"\s*(\d+)\s*\(([^\s]+)\)\s*(\w*)\s*")
    NAME = "CCSDTransition"

    @classmethod
    def from_str(cls, line: str):
        s = line.split("->")
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        lhs = s[0].strip()
        rhs = s[1].strip()
        # match with the float
        if (m := re.match(r"([-+]?\d+\.\d+)", lhs)) is not None:
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
    PATTERN = re.compile(r"\s*([^\s]+)\s*\(\s*(\d+)\s*\)\s*")
    NAME = "CISTransition"

    @classmethod
    def from_str(cls, line):
        s = line.split("-->")
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        lhs = s[0].strip()
        rhs = s[1].strip()
        idx = rhs.index("=")
        amplitude = float(rhs[idx + 1:])
        idx = rhs.find("amplitude")
        rhs = rhs[:idx]
        return cls(amplitude, lhs, rhs)

    def to_ccsd(self, homo: int) -> CCSDTransition:
        initial = f"{int(self.id_i[0][1])} (A)"
        final = f"{int(self.id_f[0][1]) + homo} (A)"
        return CCSDTransition(
            self.amplitude,
            initial,
            final,
        )
