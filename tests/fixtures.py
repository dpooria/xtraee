"""Small synthetic output files, one per supported format.

Each builder reproduces the line layout the corresponding parser reads, so the
block detection, the transition matching and the dataset code can be exercised
without shipping several-hundred-kilobyte calculation logs.  Two real logs are
committed too (``test_qccsd/input2.log`` from Q-Chem, ``test_tmcc2/ricc2.out``
from Turbomole) and are used where fidelity matters.
"""

from pathlib import Path


def write(tmp_path: Path, name: str, lines: list[str]) -> Path:
    path = Path(tmp_path) / name
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------- Q-Chem
def _rem(singlets: int | None, triplets: int | None) -> list[str]:
    lines = ["$molecule", "0 1", "O   0.000000   0.000000   0.000000",
             "H   0.000000   0.700000   0.600000", "$end", "$rem",
             "METHOD EOM-CCSD", "BASIS cc-pVDZ"]
    if singlets is not None:
        lines.append(f"EE_SINGLETS [{singlets}]")
    if triplets is not None:
        lines.append(f"EE_TRIPLETS [{triplets}]")
    lines += ["$end", ""]
    return lines


def _ccsd_state(
    irrep: str,
    number: int,
    energy: float,
    total: float,
    transitions: list[tuple[float, int, int]],
) -> list[str]:
    lines = [
        f" EOMEE transition {number}/{irrep}",
        f" Total energy = {total:.8f} a.u.  Excitation energy = {energy:.4f} eV.",
        " R0^2 = 0.0000  R1^2 = 0.9000  R2^2 = 0.1000  Res^2 = 5.10e-06",
        " Conv-d = yes",
        "",
        " Amplitude    Transitions between orbitals",
    ]
    for amp, initial, final in transitions:
        lines.append(
            f" {amp:8.4f}       {initial:2d} (A) A                  "
            f"->    {final:2d} (A) A                   "
        )
    lines += ["", "              Summary of significant orbitals:", ""]
    return lines


def _ccsd_props(
    irrep: str,
    number: int,
    frontier: tuple[float, float],
    nu: float,
    nl: float,
    prno: float,
) -> list[str]:
    return [
        f" Excited state properties for  EOMEE-CCSD transition {number}/{irrep}",
        " Dipole moment (a.u.): 0.000000 (X 0.000000, Y 0.000000, Z 0.000000)",
        "",
        " NOs",
        "   Occupation of frontier NOs:",
        f"      {frontier[0]:.4f}   {frontier[1]:.4f}",
        "   Number of electrons: 30.000000",
        f"   Number of unpaired electrons: n_u =  {nu:.5f}, n_u,nl =  {nl:.5f}",
        f"   NO participation ratio (PR_NO):  {prno:.6f}",
        "",
    ]


def ccsd_log(
    singlets: list[tuple[int, float, float, list[tuple[float, int, int]]]] | None = None,
    triplets: list[tuple[int, float, float, list[tuple[float, int, int]]]] | None = None,
    props: bool = True,
    trprops: bool = True,
    requested_singlets: int | None = None,
    requested_triplets: int | None = None,
    omit_rem_counts: bool = False,
    irreps: tuple[str, ...] = ("A",),
) -> list[str]:
    """An EOM-CCSD output.

    ``singlets``/``triplets`` are ``(number, excitation_energy_eV,
    total_energy_hartree, [(amplitude, initial_orbital, final_orbital)])`` per
    state.  ``requested_*`` overrides the ``$rem`` counts so a mismatch between
    the request and the printed states can be reproduced.
    """
    if singlets is None:
        singlets = [
            (1, 8.0000, -76.00000000, [(0.9000, 5, 6), (0.2000, 4, 7)]),
            (2, 7.5000, -76.10000000, [(0.6000, 5, 8)]),
        ]
    if triplets is None:
        triplets = [
            (1, 6.0000, -76.20000000, [(0.7000, 5, 6), (0.4000, 4, 8)]),
            (2, 7.0000, -76.30000000, [(0.5000, 4, 7)]),
        ]
    n_singlets = requested_singlets
    if n_singlets is None and singlets is not None:
        n_singlets = len(singlets)
    n_triplets = requested_triplets
    if n_triplets is None and triplets is not None:
        n_triplets = len(triplets)

    lines = _rem(None if omit_rem_counts else n_singlets,
                 None if omit_rem_counts else n_triplets)
    for irrep in irreps:
        if singlets:
            lines.append(f"        Solving for EOMEE-CCSD {irrep} singlet states.")
            lines.append("")
            for number, energy, total, transitions in singlets:
                lines += _ccsd_state(irrep, number, energy, total, transitions)
                if props:
                    lines += _ccsd_props(irrep, number, (0.9432, 1.0423), 2.1, 2.0, 2.46)
        if triplets:
            lines.append(f"        Solving for EOMEE-CCSD {irrep} triplet states.")
            lines.append("")
            for number, energy, total, transitions in triplets:
                lines += _ccsd_state(irrep, number, energy, total, transitions)
                if props:
                    lines += _ccsd_props(irrep, number, (1.0000, 1.0000), 2.0, 2.0, 3.00)
    if trprops:
        lines.append(" Start computing the transition properties")
        for number, *_ in singlets or []:
            lines += [
                f" State B: eomee_ccsd/rhfref/singlets: {number}/A",
                " Oscillator strength (a.u.): 0.7500",
                " ||gamma^AB||*||gamma^BA||: 12.5000",
                " omega = 0.8000",
                " 2<alpha|beta> = 0.9000",
                " LOC = 0.5000",
                " <Phe> = -0.1000",
                " |<r_e - r_h>| [Ang]: 0.0000",
                " Correlation coefficient: 0.0200",
                "",
            ]
    return lines


# ------------------------------------------------------------------ Q-Chem CIS
def _cis_state(
    number: int,
    energy: float,
    total: float,
    multi: str,
    strength: float,
    dipole: tuple[float, float, float],
    transitions: list[tuple[float, int, int]],
) -> list[str]:
    lines = [
        f"Excited state {number}: excitation energy (eV) = {energy:.4f}",
        f"Total energy for state {number:2d}: {total:.8f} au",
        f"Multiplicity: {multi}",
        f"Strength   : {strength:.6f}",
        f"Trans. Mom.: {dipole[0]:.4f} X   {dipole[1]:.4f} Y   {dipole[2]:.4f} Z",
    ]
    for amp, initial, final in transitions:
        lines.append(
            f"     D({initial:5d}) --> V({final:4d}) amplitude = {amp:.4f}"
        )
    return lines


def cis_log(
    states: list[tuple[int, float, float, str, float, tuple[float, float, float],
                        list[tuple[float, int, int]]]] | None = None,
    with_orbitals: bool = True,
    method: str = "CIS",
) -> list[str]:
    """A Q-Chem CIS (or TDDFT) output with singlet and triplet states."""
    if states is None:
        states = [
            (1, 7.5000, -75.90000000, "Singlet", 0.0451, (0.0, 0.1, 0.2),
             [(0.9100, 5, 1), (-0.2000, 4, 2)]),
            (2, 8.2000, -75.80000000, "Triplet", 0.0010, (0.0, 0.0, 0.1),
             [(0.8000, 5, 1), (0.3000, 3, 2)]),
        ]
    lines = [
        "$molecule", "0 1", "O   0.000000   0.000000   0.000000",
        "H   0.000000   0.700000   0.600000", "$end", "$rem", f"METHOD {method}",
        "$end", "",
        f" {method} Excitation Energies",
        "",
    ]
    for number, energy, total, multi, strength, dipole, transitions in states:
        lines += _cis_state(number, energy, total, multi, strength, dipole,
                            transitions)
        lines.append("")
    if with_orbitals:
        lines += [" Orbital Energies (a.u.)", "", " 1  -20.5500", " 5   -0.5000"]
    return lines


# ---------------------------------------------------------------- Q-Chem CISD
def cisd_log(
    roots: list[tuple[int, float, float, list[tuple[float, int, int]]]] | None = None,
    multi: str = "LOWSPIN",
    irrep: str = "A1",
) -> list[str]:
    """A Q-Chem CISD output (``Root N Conv-d yes Tot Ene=...`` blocks)."""
    if roots is None:
        roots = [
            (1, 7.5000, -75.90000000, [(0.9100, 5, 1), (-0.2000, 4, 2)]),
        ]
    lines = [
        "$molecule", "0 1", "O   0.000000   0.000000   0.000000", "$end",
        "",
        " Doing genuine CISD calculations",
        f" SOLVE LINEAR RESPONSE EQUATIONS FOR {multi} SINGLE AND DOUBLE "
        f"EXCITATIONS OF {irrep} IRREP",
        "",
    ]
    for number, energy, total, transitions in roots:
        lines += [
            f" Root {number} Conv-d yes Tot Ene= {total:.6f} hartree "
            f"(Ex Ene  {energy:.4f} eV), U0^2= 0.9000, U1^2= 0.0900, U2^2= 0.0100",
        ]
        for amp, initial, final in transitions:
            lines.append(
                f" {amp:8.4f}                  {initial}(   A1) A   ->     "
                f"{final}(   A1) A"
            )
        lines.append("")
    lines += [" Orbital Energies (a.u.)", ""]
    return lines


# ---------------------------------------------------------------- Q-Chem ADC(2)
def adc2_log(
    states: list[tuple[int, str, str, float, float, float,
                        list[tuple[float, int, int]]]] | None = None,
) -> list[str]:
    """A Q-Chem ADC(2) ``Excited State Summary`` block."""
    if states is None:
        states = [
            (1, "singlet", "a1", 7.5000, -75.90000000, 0.0451,
             [(0.9100, 5, 1)]),
            (2, "triplet", "a1", 8.2000, -75.80000000, 0.0010,
             [(0.8000, 5, 1)]),
        ]
    lines = ["$molecule", "0 1", "O   0.000000   0.000000   0.000000", "$end",
             "", " Running an ADC(2) calculation", " Excited State Summary", ""]
    for number, multi, irrep, energy, total, strength, transitions in states:
        lines += [
            f" Excited state {number} ({multi}, {irrep}) [converged]",
            f"  Excitation energy: {energy:.4f} eV",
            f"  Total energy: {total:.6f} au",
            f"  Osc. strength: {strength:.4f}",
        ]
        for amp, initial, final in transitions:
            lines.append(f"  {initial} ({irrep}) A {final} ({irrep}) A {amp:.4f}")
        lines.append("")
    return lines
