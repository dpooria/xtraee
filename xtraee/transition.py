from typing import NamedTuple

from xtraee.lazypattern import LP, VERBOSE


class TrID(NamedTuple):
    orb_num: int
    irrep: str
    multi: str

    def __repr__(self):
        return f"{self.orb_num} ({self.irrep}) {self.multi}"


class Transition:
    name = "Transition"
    pattern = LP(
        r"""
    (?P<orb_num>\d+)       # integer
    \s*
    \(\s*
    (?P<irrep>[A-Za-z0-9'"]{1,3})      # allow letters, digits, ' and "
    \s*\)
    \s+
    (?P<multi>[AB])                  # A or B
    """,
        VERBOSE,
    )

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
        self.id_i: list[TrID] = []
        self.id_f: list[TrID] = []
        self._spin_ind = 2
        self.is_double = False
        if not self.parse(self.initial, self.final):
            raise ValueError(f"Cannot parse {initial} -> {final}")

    @classmethod
    def from_str(cls, line: str):
        raise NotImplementedError

    def trid_match(self, tr_str: str):
        trid = []
        for m in self.pattern.finditer(tr_str):
            md = m.groupdict()
            trid.append(
                TrID(int(md["orb_num"]), md.get("irrep", "A"), md.get("multi", "A"))
            )
        return trid

    def parse(self, initial, final) -> bool:
        self.id_i = self.trid_match(initial)
        self.id_f = self.trid_match(final)
        if len(self.id_i) > 1 or len(self.id_f) > 1:
            self.is_double = True
        else:
            self.is_double = False
        return True

    def is_equal(self, other: "Transition", shallow=False) -> bool:
        # currently we don't check for spin matching at all (use trblock.utrs)
        if self.is_double != other.is_double:
            return False

        equal_orbnum = True
        for id_i, id_f, other_id_i, other_id_f in zip(
            self.id_i, self.id_f, other.id_i, other.id_f
        ):
            if id_i.orb_num != other_id_i.orb_num or id_f.orb_num != other_id_f.orb_num:
                equal_orbnum = False
                break
        if (not equal_orbnum) or shallow:
            return equal_orbnum
        full_equal = True
        for id_i, id_f, other_id_i, other_id_f in zip(
            self.id_i, self.id_f, other.id_i, other.id_f
        ):
            if id_i.irrep != other_id_i.irrep or id_f.irrep != other_id_f.irrep:
                full_equal = False
                break
        return full_equal

    def to_std(self, homo: int = 0) -> "CCSDTransition":
        return self

    def __eq__(self, other) -> bool:
        return self.is_equal(other)

    def __repr__(self) -> str:
        return f"{self.amplitude:.4f}\t{self.initial} -> {self.final}"


# CCSDTransition is the standard format
class CCSDTransition(Transition):
    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CCSDTransition"
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
        super().__init__(amplitude, initial, final)

    @staticmethod
    def make_standard(t: str) -> str:
        p_open = t.find("(")
        p_close = t.find(")")
        if p_open == -1 or p_close == -1:
            raise ValueError(f"cannot convert {t} to the standard form")
        return f"{int(t[p_open+1:p_close])} (A) A"

    @classmethod
    def from_str(cls, line: str):
        s = line.split("-->")
        if len(s) != 2:
            raise ValueError(f"cannot match {line}")
        initial = cls.make_standard(s[0].strip())
        rhs = s[1].strip()
        idx = rhs.index("=")
        amplitude = float(rhs[idx + 1 :])
        final = cls.make_standard(rhs)
        return cls(amplitude, initial, final)

    def to_std(self, homo: int) -> CCSDTransition:
        # not considering symmetry and there is no double excitation in CIS
        initial = f"{self.id_i[0].orb_num} (A) A"
        final = f"{self.id_f[0].orb_num + homo} (A) A"
        return CCSDTransition(self.amplitude, initial, final)


class CISDTransition(Transition):
    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CISDTransition"
        # self.pattern = LP(r"\s*(\d+)\s*\(\s*([^\s]+)\s*\)\s*([AB])\s*")
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

    def to_ccsd(self, homo: int, preserve_irreps: bool = False) -> CCSDTransition:
        initials = []
        finals = []
        # CISD orbitals numbering starts from 0
        for id_i, id_f in zip(self.id_i, self.id_f):
            if preserve_irreps:
                initial = f"{id_i.orb_num + 1} ({id_i.irrep}) {id_i.multi}"
                final = f"{id_f.orb_num + homo + 1} ({id_f.irrep}) {id_f.multi}"
            else:
                initial = f"{id_i.orb_num + 1} (A) {id_i.multi})"
                final = f"{id_f.orb_num + homo + 1} (A) {id_f.multi}"
            initials.append(initial)
            finals.append(final)
        return CCSDTransition(self.amplitude, "\t".join(initials), "\t".join(finals))

    def to_std(self, homo):
        return self.to_ccsd(homo, False)


class TMCC2Transition(Transition):
    amppattern = LP(r"\s*([+-]?\d+\.\d+)\s+[+-]?\d+\.\d+\s*")
    nonstd_pattern = LP(
        r"""
        \s*
        (?P<orb_num>\d+)
        \s+
        (?P<irrep>[a-z]{1,3})
        \s+
        \d+
        \s*
    """,
        VERBOSE,
    )

    def __init__(self, amplitude: float, initial: str, final: str):
        self.name = "CC2Transition"
        super().__init__(amplitude, initial, final)
        self._spin_ind = 1

    @classmethod
    def from_str(cls, line: str):
        s = line.split("|")
        if len(s) < 4:
            raise ValueError(f"cannot match {line}")
        lhs = cls.make_standard(s[1].strip())
        rhs = cls.make_standard(s[2].strip())
        if (m := cls.amppattern.match(s[3].strip())) is not None:
            amplitude = float(m.group(1))
        else:
            raise ValueError(f"could not match the amplitude {line}")
        return cls(amplitude, lhs, rhs)

    @classmethod
    def make_standard(cls, t: str) -> str:
        if (m := cls.nonstd_pattern.match(t)) is not None:
            orb_num = m["orb_num"]
            irrep = m["irrep"].upper()  # ?
            # TODO: add spin stuff
            return f"{orb_num} ({irrep}) A"
        else:
            raise ValueError(f"cannot convert {t} to the standard form")
