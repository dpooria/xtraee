#!/usr/bin/env python

from __future__ import annotations

from argparse import ArgumentParser
from typing import TYPE_CHECKING
from math import sqrt

import numpy as np
import pandas as pd

from pyscf import gto, scf, ao2mo
from xtraee.parser import Parser, QCISParser

if TYPE_CHECKING:
    from argparse import Namespace
    from typing import Any

    from xtraee.trblock import TransitionBlock

    AtomsT = list[tuple[str, tuple[float, ...]]]
    pyscf_t = Any

Ha = 27.211386245988


def run_hf(atoms: AtomsT, basis: str = 'cc-PVDZ') -> tuple[pyscf_t, pyscf_t]:
    mol = gto.M(
        atom=atoms,
        basis=basis,
        unit="Angstrom",
    )

    mf = scf.RHF(mol).run()

    return mol, mf


def run_cis(mf: pyscf_t, n_singlets: int, n_triplets: int) -> tuple[pyscf_t, pyscf_t]:
    cis_s = mf.TDA()
    cis_s.singlet = True
    cis_s.nstates = n_singlets
    cis_s.kernel()

    cis_t = mf.TDA()
    cis_t.singlet = False
    cis_t.nstates = n_triplets
    cis_t.kernel()

    return cis_s, cis_t


def extract_singlet_triplet(parser: QCISParser) \
        -> tuple[TransitionBlock, TransitionBlock]:

    singlet = parser.data['S1']['singlet']
    scores = parser.scores['A_A']
    t_match_ind = scores.loc[singlet.identifier].argmax()
    t_match = scores.columns[t_match_ind]
    triplet = parser.irreps['triplet-A'].trblocks_dict[t_match]

    return singlet, triplet


def get_iac(trblock: TransitionBlock, homo: int,
            td: pyscf_t, use_pyscf_indices: bool) \
        -> tuple[list[int], list[int], list[float]]:
    i = []
    a = []
    c = []
    if use_pyscf_indices:
        raise NotImplementedError
    else:
        for transition in trblock.transitions:
            i.append(transition.id_i[0].orb_num - 1)
            a.append(transition.id_f[0].orb_num + homo - 1)
            c.append(transition.amplitude)

    i = np.array(i, dtype=int)
    a = np.array(a, dtype=int)
    c = np.array(c, dtype=float)
    return i, a, c


def exchange_integral(mol: pyscf_t, mf: pyscf_t, i: list[int],
                      a: list[int], j: list[int], b: list[int]) \
        -> np.ndarray:
    """(ia|jb)"""

    C = mf.mo_coeff

    N1 = len(i)
    N2 = len(j)
    assert len(a) == N1 and len(b) == N2

    K_iajb = ao2mo.general(
        mol,
        (C[:, i], C[:, a],
         C[:, j], C[:, b]),
        compact=False,
    ).reshape(N1, N1, N2, N2)

    return K_iajb * Ha


def calculate_CIS_exchange(mol: pyscf_t, mf: pyscf_t,
                           c1: list[float], i: list[int],
                           a: list[int], phase1: list[float],
                           c2: list[float], j: list[int],
                           b: list[int], phase2: list[float]) -> float:
    """sum_{pq}{ (phase_p * c_p) (i_p,a_p|j_q,b_q) (phase_q * c_q) }"""

    K_iajb = exchange_integral(mol, mf, i, a, j, b)

    # trace
    K = np.einsum('ppqq->pq', K_iajb)

    c1 = phase1 * np.array(c1)
    c2 = phase2 * np.array(c2)

    return float(c1 @ K @ c2)


def get_phase_amplitudes_energy(amplitudes: list[float], td: pyscf_t, E: float,
                                i: list[int], a: list[int], homo: int) \
        -> tuple[list[float], list[float], float]:
    root = np.argmin(np.abs(E - td.e * Ha))
    energy = td.e[root] * Ha
    td_amps = td.xy[root][0]  # nocc x nvirt
    phase = []
    cis_amplitudes = []
    sqrt2 = sqrt(2)
    for p, c in enumerate(amplitudes):
        td_amp = td_amps[i[p], a[p] - homo]
        cis_amplitudes.append(float(td_amp) * sqrt2)
        if td_amp != 0.0:
            phase.append(np.sign(td_amp) * np.sign(c))
        else:
            phase.append(1.0)
    phase = np.array(phase, dtype=float)
    cis_amplitudes = np.array(cis_amplitudes, dtype=float)
    return phase, cis_amplitudes, energy


def parse_args() -> Namespace:
    parser = ArgumentParser(description="Calculate the CIS exchange integral using Pyscf from a QChem CIS calculation.")
    parser.add_argument('cis_file', type=str)
    parser.add_argument('singlet', type=str, help='e.g. singlet-3/')
    parser.add_argument('triplet', type=str, help='e.g. triplet-2/')
    parser.add_argument('-o', '--out', default=None, type=str, help='name of the output .csv file')
    parser.add_argument('--threshold', type=float, default=1e-5, help='minimum amplitude^2 to consider from the logfile')
    parser.add_argument('--ref-pyscf', action='store_true', help='not implemented yet!')
    return parser.parse_args()


def main() -> int:

    # ----------- Args ------------
    args = parse_args()
    cis_file = args.cis_file
    outfile = args.out or cis_file + '.csv'
    s_label = args.singlet
    t_label = args.triplet
    use_pyscf_indices = args.ref_pyscf
    threshold = args.threshold

    # ----------- Parse inputfile ------------
    parser = Parser(cis_file, threshold=threshold)
    assert isinstance(
        parser, QCISParser), f"This is not a CIS calculation: {cis_file}: {type(parser)}"
    parser.run()
    homo = parser.homo

    # singlet, triplet = extract_singlet_triplet(parser, s_label, t_label)
    singlet = parser.irreps['singlet-A'].trblocks_dict[s_label]
    triplet = parser.irreps['triplet-A'].trblocks_dict[t_label]

    E_S = singlet.excitation_energy
    E_T = triplet.excitation_energy
    n_singlets = parser.irreps['singlet-A'].n_states
    n_triplets = parser.irreps['triplet-A'].n_states

    # ----------- Pyscf run ------------
    mol, mf = run_hf(parser.atoms)
    cis_s, cis_t = run_cis(mf, n_singlets, n_triplets)

    # ----------- extract transition indices ------------
    i, a, c_s = get_iac(singlet, homo, cis_s, use_pyscf_indices)
    j, b, c_t = get_iac(triplet, homo, cis_t, use_pyscf_indices)

    phase_s, cis_amp_s, ECIS_S = get_phase_amplitudes_energy(
        c_s, cis_s, E_S, i, a, homo)
    phase_t, cis_amp_t, ECIS_T = get_phase_amplitudes_energy(
        c_t, cis_t, E_T, j, b, homo)

    print("pyscf amplitudes norm: ", np.linalg.norm(
        cis_amp_s), np.linalg.norm(cis_amp_t))
    print("QChem amplitudes norm: ", np.linalg.norm(
        c_s), np.linalg.norm(c_t))

    # ----------- QChem CIS ------------
    K_S = calculate_CIS_exchange(
        mol, mf, c_s, i, a, phase_s, c_s, i, a, phase_s)

    K_T = calculate_CIS_exchange(
        mol, mf, c_t, j, b, phase_t, c_t, j, b, phase_t)

    K_avg = (K_S + K_T) / 2.0

    # ----------- HL ------------
    homo_ind = homo - 1
    lumo_ind = homo_ind + 1

    ind_hl_s = (i == homo_ind) & (a == lumo_ind)
    if np.any(ind_hl_s):
        c_hl = c_s[ind_hl_s]
        phase_hl = phase_s[ind_hl_s]
        KHL_S = calculate_CIS_exchange(mol, mf,
                                       c_hl, [homo_ind], [lumo_ind], phase_hl,
                                       c_hl, [homo_ind], [lumo_ind], phase_hl)
    else:
        c_hl = 0.0
        KHL_S = 0.0

    print('Singlet HL ampitude^2:', c_hl**2, 'HL transition: ', np.where(ind_hl_s)[0])

    ind_hl_t = (j == homo_ind) & (b == lumo_ind)
    if np.any(ind_hl_t):
        c_hl = c_t[ind_hl_t]
        phase_hl = phase_t[ind_hl_t]
        KHL_T = calculate_CIS_exchange(mol, mf,
                                       c_hl, [homo_ind], [lumo_ind], phase_hl,
                                       c_hl, [homo_ind], [lumo_ind], phase_hl)
    else:
        c_hl = 0.0
        KHL_T = 0.0
    print('Triplet HL ampitude^2:', c_hl**2, 'HL transition: ', np.where(ind_hl_t)[0])

    KHL_avg = (KHL_S + KHL_T) / 2.0

    # ----------- PYSCF CIS ------------
    KCIS_S = calculate_CIS_exchange(
        mol, mf, cis_amp_s, i, a, 1.0, cis_amp_s, i, a, 1.0)

    KCIS_T = calculate_CIS_exchange(
        mol, mf, cis_amp_t, j, b, 1.0, cis_amp_t, j, b, 1.0)

    KCIS_avg = (KCIS_S + KCIS_T) / 2.0

    print("CIS transition amplitudes from pyscf:")
    print("singlets: ", [f"{c:.3f}" for c in cis_amp_s])
    print("triplets: ", [f"{c:.3f}" for c in cis_amp_t])

    # ----------- output ------------
    data = {
        "singlet": singlet.identifier, "triplet": triplet.identifier,
        "E_S": E_S, "E_T": E_T,
        "E_pyscf_S": ECIS_S, "E_pyscf_T": ECIS_T,
        "K_S": K_S, "K_T": K_T,
        "K_pyscfCIS_S": KCIS_S, "K_pyscfCIS_T": KCIS_T,
        "dE_ST": E_S - E_T, "2K_avg": 2 * K_avg,
        "dE_pyscf_ST": ECIS_S - ECIS_T, "2K_pyscf_avg": 2 * KCIS_avg,
        "KHL_S": KHL_S, "KHL_T": KHL_T, "2KHL_avg": 2 * KHL_avg,
    }
    df = pd.DataFrame([data])
    print(df)
    df.to_csv(outfile)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
