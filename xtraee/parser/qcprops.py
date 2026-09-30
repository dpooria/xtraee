"""Excited-state properties from a Q-Chem TDDFT/CIS output.

Complements :class:`xtraee.parser.qcis.QCTDDFTParser` (states, energies,
amplitudes) with the post-TDDFT property sections:

* spin-orbit couplings, sublevel resolved and complex, from the
  ``calc_soc`` job:
    ``<S0|H_SO|T_n, m_s>``      -> :attr:`soc_ground`   ``[m_s][n]``
    ``<S_m|H_SO|T_n, m_s>``     -> :attr:`soc_singlet`  ``[m_s][m, n]``
  in cm^-1, as printed.  Row/column k is the k-th state of that multiplicity
  in Q-Chem's labelling (``S1``, ``T1``, ...).
* state and transition dipole moments (atomic units), from the
  ``Electron Dipole Moments ...`` and ``Transition Moments Between ...``
  tables, keyed by Q-Chem's overall state numbers (0 = ground state):
    :attr:`state_dipoles`       ``{number: (x, y, z)}``  (0 included)
    :attr:`transition_dipoles`  ``{(i, j): (x, y, z)}``  with i < j

``multiplicity`` and ``energies`` map each overall state number to
"singlet"/"triplet" and to its excitation energy (eV), and
:meth:`numbers` gives, per multiplicity, the overall numbers in label order,
so ``numbers("triplet")[k]`` is the state Q-Chem calls ``T{k+1}``.

The sign of every transition quantity is an arbitrary wavefunction phase, but
it is *consistent within one output*, which is what products such as
SOC x dipole need.
"""

from pathlib import Path

import numpy as np

from xtraee.lazypattern import LP

PathType = str | Path

_state_header = LP(
    r"^Excited\s+state\s+(\d+)\s*:\s*excitation\s+energy\s+\(eV\)\s*=\s*([-+]?\d+\.\d+)"
)
_multiplicity = LP(r"^Multiplicity:\s*(Singlet|Triplet)")
_soc_ground_header = LP(
    r"^SOC between the singlet ground state and excited triplet states \(ms=(-?\d)\):"
)
_soc_singlet_header = LP(
    r"^SOC between the S(\d+) state and excited triplet states \(ms=(-?\d)\):"
)
_soc_line = LP(
    r"^T(\d+)\(ms=-?\d\)\s+([-+]?\d+\.\d+)\s*([-+])\s*(\d+\.\d+)i\s+cm-1"
)
_dipole_header = LP(r"^Electron Dipole Moments of (Ground State|Singlet Excited State|Triplet Excited State)")
_transition_header = LP(r"^Transition Moments Between (.+)$")
_float = r"([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)"
_dipole_row = LP(rf"^(\d+)\s+{_float}\s+{_float}\s+{_float}\s*$")
_transition_row = LP(rf"^(\d+)\s+(\d+)\s+{_float}\s+{_float}\s+{_float}\s+\S+\s*$")


class QCPropertiesParser:
    """Parse SOCs and dipole tables; see the module docstring."""

    def __init__(self, input_file: PathType):
        self.input_file = input_file
        self.multiplicity: dict[int, str] = {}
        self.energies: dict[int, float] = {}
        self._soc_ground: dict[int, dict[int, complex]] = {}
        self._soc_singlet: dict[int, dict[tuple[int, int], complex]] = {}
        self.state_dipoles: dict[int, tuple[float, float, float]] = {}
        self.transition_dipoles: dict[tuple[int, int], tuple[float, float, float]] = {}

    # ------------------------------------------------------------ parsing
    def run(self) -> "QCPropertiesParser":
        mode, key, state = None, None, None
        with open(self.input_file, "r") as handle:
            for raw in handle:
                line = raw.strip()
                if (m := _state_header.match(line)) is not None:
                    state = int(m.group(1))
                    self.energies[state] = float(m.group(2))
                    mode = None
                    continue
                if state is not None and (m := _multiplicity.match(line)) is not None:
                    self.multiplicity[state] = m.group(1).lower()
                    state = None
                    continue
                if (m := _soc_ground_header.match(line)) is not None:
                    mode, key = "soc_ground", int(m.group(1))
                    self._soc_ground.setdefault(key, {})
                    continue
                if (m := _soc_singlet_header.match(line)) is not None:
                    mode, key = "soc_singlet", (int(m.group(1)), int(m.group(2)))
                    self._soc_singlet.setdefault(key[1], {})
                    continue
                if (m := _dipole_header.match(line)) is not None:
                    mode, key = "dipole", None
                    continue
                if (m := _transition_header.match(line)) is not None:
                    mode, key = "transition", None
                    continue
                if mode in ("soc_ground", "soc_singlet"):
                    if (m := _soc_line.match(line)) is None:
                        mode = None
                        continue
                    imag = float(m.group(4)) * (1 if m.group(3) == "+" else -1)
                    value = complex(float(m.group(2)), imag)
                    triplet = int(m.group(1))
                    if mode == "soc_ground":
                        self._soc_ground[key][triplet] = value
                    else:
                        self._soc_singlet[key[1]][(key[0], triplet)] = value
                elif mode == "dipole":
                    if (m := _dipole_row.match(line)) is not None:
                        self.state_dipoles[int(m.group(1))] = tuple(
                            float(m.group(i)) for i in (2, 3, 4)
                        )
                elif mode == "transition":
                    if (m := _transition_row.match(line)) is not None:
                        i, j = int(m.group(1)), int(m.group(2))
                        vector = tuple(float(m.group(k)) for k in (3, 4, 5))
                        self.transition_dipoles[(min(i, j), max(i, j))] = vector
        return self

    # ------------------------------------------------------------ access
    def numbers(self, multi: str) -> list[int]:
        """Overall state numbers of one multiplicity, in label order (S1, S2 ...)."""
        return sorted(n for n, m in self.multiplicity.items() if m == multi)

    def excitation_energies(self, multi: str) -> np.ndarray:
        """Excitation energies (eV) of one multiplicity, in label order."""
        return np.array([self.energies[n] for n in self.numbers(multi)])

    @property
    def sublevels(self) -> list[int]:
        return sorted(self._soc_ground)

    @property
    def soc_ground(self) -> dict[int, np.ndarray]:
        """``{m_s: array[n_triplets]}`` of <S0|H_SO|T_n, m_s>, cm^-1."""
        out = {}
        for ms, row in self._soc_ground.items():
            size = max(row)
            out[ms] = np.array([row.get(n, np.nan) for n in range(1, size + 1)])
        return out

    @property
    def soc_singlet(self) -> dict[int, np.ndarray]:
        """``{m_s: array[n_singlets, n_triplets]}`` of <S_m|H_SO|T_n, m_s>, cm^-1."""
        out = {}
        for ms, table in self._soc_singlet.items():
            rows = max(m for m, _ in table)
            cols = max(n for _, n in table)
            array = np.full((rows, cols), np.nan, dtype=complex)
            for (m, n), value in table.items():
                array[m - 1, n - 1] = value
            out[ms] = array
        return out

    def dipole(self, i: int, j: int) -> np.ndarray:
        """<i|mu|j> in a.u. between overall state numbers (i == j: state dipole)."""
        if i == j:
            return np.asarray(self.state_dipoles[i], float)
        return np.asarray(self.transition_dipoles[(min(i, j), max(i, j))], float)
