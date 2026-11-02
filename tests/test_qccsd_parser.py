"""The EOM-CCSD parser, against the committed Q-Chem log and synthetic ones."""

from pathlib import Path

import fixtures
import pytest

from xtraee.parser import QCCSDParser
from xtraee.utils import Ha

REAL_LOG = Path(__file__).parent / "test_qccsd" / "input2.log"


@pytest.fixture(scope="module")
def parser() -> QCCSDParser:
    p = QCCSDParser(REAL_LOG, 0.0)
    p.process_file()
    return p


def test_real_log_parses_every_requested_state(parser):
    assert {key: len(irrep) for key, irrep in parser.irreps.items()} == {
        "singlet-A": 5,
        "triplet-A": 10,
    }
    assert all(irrep.n_states == len(irrep) for irrep in parser.irreps.values())


def test_real_log_energies_and_transitions(parser):
    state = parser.find_trblock("singlet-1/A")
    assert state.excitation_energy == pytest.approx(6.9970)
    # the total energy is printed in hartree and reported in eV
    assert state.total_energy == pytest.approx(-155.25749592 * Ha)
    assert state.oscillator_strength == pytest.approx(0.796701)
    amplitudes = [round(tr.amplitude, 4) for tr in state.transitions]
    assert amplitudes == [0.6732, 0.6732, -0.0408, -0.0408]
    assert str(state.transitions[0]) == "0.6732\t15 (A) A -> 16 (A) A"


def test_real_log_states_are_sorted_by_excitation_energy(parser):
    energies = [block.excitation_energy for block in parser.irreps["singlet-A"]]
    assert energies == sorted(energies)
    energies = [block.excitation_energy for block in parser.irreps["triplet-A"]]
    assert energies == sorted(energies)


def test_real_log_state_properties(parser):
    meta = parser.find_trblock("singlet-1/A").meta_data
    assert meta["R1"] == pytest.approx(0.9426)
    assert meta["R2"] == pytest.approx(0.0574)
    assert meta["gamma"] == pytest.approx(0.887439)
    assert meta["corr_coef"] == pytest.approx(0.015467)
    assert meta["frontier_no"] == pytest.approx([0.9432, 1.0423])
    # these used to be stored as strings, which made the descriptor CSV mixed-type
    assert isinstance(meta["nu"], float)
    assert meta["nu"] == pytest.approx(2.11014)
    assert meta["prno"] == pytest.approx(2.461606)


def test_descriptors_are_extracted_for_every_state(parser):
    rows = parser.gather_descriptors()
    assert len(rows) == 15
    first = rows[0]
    assert first["identifier"] == "singlet-1/A"
    assert first["|r_e-r_h|"] == pytest.approx(0.0)
    assert first["frontier_no_1"] == pytest.approx(0.9432)
    assert first["frontier_no_2"] == pytest.approx(1.0423)
    assert "R0" not in first and "R1" not in first
    assert "frontier_no" not in first


def test_descriptors_survive_states_without_a_property_section(tmp_path):
    """A state printed without the NO analysis used to break ``descriptors``."""
    log = fixtures.write(
        tmp_path,
        "no-props.log",
        fixtures.ccsd_log(props=False, trprops=False),
    )
    p = QCCSDParser(log, 0.0)
    p.process_file()
    rows = p.gather_descriptors()
    assert len(rows) == 4
    assert rows[0]["frontier_no_1"] != rows[0]["frontier_no_1"]  # nan


def test_dataset_of_the_real_log(parser):
    data = parser.create_dataset()
    assert list(data) == ["S1", "S2", "S3", "S4", "S5"]
    matches = data["S1"]["matched_triplets"]
    assert matches["0.6732\t15 (A) A -> 16 (A) A"] == ["1 A 0.6612 "]


def test_compare_all_returns_a_matrix_per_pair(parser):
    scores = parser.compare_all()
    assert list(scores) == ["A_A"]
    frame = scores["A_A"]
    assert frame.shape == (5, 10)
    assert frame.index[0] == "singlet-1/A"
    assert frame.columns[0] == "triplet-1/A"
