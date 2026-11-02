"""CISD and ADC(2) parsing."""

import fixtures
import pytest

from xtraee.parser import Parser, QCADC2Parser, QCISDParser


def test_cisd_root_block(tmp_path):
    path = fixtures.write(tmp_path, "cisd.log", fixtures.cisd_log())
    parser = Parser(path, 0.0)
    assert isinstance(parser, QCISDParser)
    assert parser.name == "CISD"
    parser.process_file()

    (key, irrep), = parser.irreps.items()
    assert key == "CISD-A1.singlet"
    assert irrep.n_states == 1

    state = irrep[0]
    assert state.excitation_energy == pytest.approx(7.5)
    assert state.total_energy == pytest.approx(-75.9 * 27.211386024367243)
    assert [round(tr.amplitude, 4) for tr in state.transitions] == [0.91, -0.2]
    assert state.meta_data["U0"] == pytest.approx(0.9)
    assert state.meta_data["U1"] == pytest.approx(0.09)
    assert state.meta_data["U2"] == pytest.approx(0.01)


def test_cisd_homo_is_one_past_the_highest_occupied(tmp_path):
    path = fixtures.write(tmp_path, "cisd.log", fixtures.cisd_log())
    parser = QCISDParser(path, 0.0)
    parser.process_file()
    state = parser.find_trblock("singlet-1/A1")
    assert parser.homo == 0  # not set on the parser by this path
    assert state.homo == 6


def test_adc2_states(tmp_path):
    path = fixtures.write(tmp_path, "adc2.log", fixtures.adc2_log())
    parser = Parser(path, 0.0)
    assert isinstance(parser, QCADC2Parser)
    assert parser.name == "ADC2"
    parser.process_file()

    assert {key: len(irrep) for key, irrep in parser.irreps.items()} == {
        "singlet-a1": 1,
        "triplet-a1": 1,
    }
    state = parser.find_trblock("singlet-1/a1")
    assert state.excitation_energy == pytest.approx(7.5)
    assert state.oscillator_strength == pytest.approx(0.0451)
    assert [round(tr.amplitude, 4) for tr in state.transitions] == [0.91]
    assert state.transitions[0].parsed
