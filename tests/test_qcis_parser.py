"""CIS / TDDFT parsing (Q-Chem ``Excited state N: excitation energy`` blocks)."""

import fixtures
import pytest

from xtraee.parser import Parser


def parse(tmp_path, lines, threshold=0.0, name="cis.log"):
    p = Parser(fixtures.write(tmp_path, name, lines), threshold)
    assert p is not None
    p.process_file()
    return p


def test_cis_states_and_amplitudes(tmp_path):
    parser = parse(tmp_path, fixtures.cis_log())
    assert parser.name == "CIS"
    assert {key: len(irrep) for key, irrep in parser.irreps.items()} == {
        "singlet-A": 1,
        "triplet-A": 1,
    }

    state = parser.find_trblock("singlet-1/")
    assert state.excitation_energy == pytest.approx(7.5)
    assert state.oscillator_strength == pytest.approx(0.0451)
    assert state.transition_dipole == pytest.approx((0.0, 0.1, 0.2))
    assert [round(tr.amplitude, 4) for tr in state.transitions] == [0.91, -0.2]


def test_cis_reports_the_total_energy_in_ev(tmp_path):
    parser = parse(tmp_path, fixtures.cis_log())
    state = parser.find_trblock("singlet-1/")
    assert state.total_energy == pytest.approx(-75.9 * 27.211386024367243)


def test_cis_homo_shifts_the_virtual_orbitals_to_ccsd_numbering(tmp_path):
    parser = parse(tmp_path, fixtures.cis_log())
    state = parser.find_trblock("singlet-1/")
    state.generate_std()
    assert str(state.std_transitions[0].id_i[0]) == "5 (A) A"
    assert str(state.std_transitions[0].id_f[0]) == "6 (A) A"
    assert parser.homo == 5


def test_cis_dataset_only_lists_singlets(tmp_path):
    data = parse(tmp_path, fixtures.cis_log()).create_dataset()
    assert list(data) == ["S1"]
    # the triplet state exists, and shares the 5 -> 1 transition
    matches = data["S1"]["matched_triplets"]
    # CIS blocks carry no irrep label, hence the double space
    assert matches["0.9100\t5 (A) A -> 1 (A) A"] == ["2  0.8000 "]
    assert matches["-0.2000\t4 (A) A -> 2 (A) A"] == []


def test_tddft_is_detected_and_parsed(tmp_path):
    parser = parse(tmp_path, fixtures.cis_log(method="TDDFT"), name="tddft.log")
    assert parser.name == "TDDFT"
    assert len(parser.find_trblock("singlet-1/").transitions) == 2
