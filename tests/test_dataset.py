"""Dataset building (the "happy family") and its edge cases."""

import fixtures
import pytest

from xtraee.irrep import Irrep
from xtraee.parser import Parser, QCCSDParser
from xtraee.parser.base import BaseParser


def parse(tmp_path, lines, threshold=0.0, name="input.log"):
    path = fixtures.write(tmp_path, name, lines)
    parser = Parser(path, threshold)
    assert isinstance(parser, BaseParser)
    parser.process_file()
    return parser


def test_every_lowest_singlet_appears_even_without_a_triplet_match(tmp_path):
    """A singlet whose transitions match no triplet used to vanish.

    The entry was only created inside the matching loop, so a state with no
    triplet counterpart (common for Rydberg or charge-transfer states) was
    missing from the output entirely instead of showing an empty match list.
    """
    lines = fixtures.ccsd_log(
        singlets=[
            (1, 8.0, -76.0, [(0.9, 5, 6)]),
            (2, 7.5, -76.1, [(0.6, 9, 12)]),
        ],
        triplets=[(1, 6.0, -76.2, [(0.7, 5, 6)])],
    )
    parser = parse(tmp_path, lines)
    data = parser.create_dataset()

    # the dataset is ordered by excitation energy, not by the log's numbering
    assert list(data) == ["S1", "S2"]
    assert data["S1"]["singlet"].identifier == "singlet-2/A"
    assert data["S1"]["matched_triplets"] == {"0.6000\t9 (A) A -> 12 (A) A": []}
    assert data["S2"]["singlet"].identifier == "singlet-1/A"
    assert data["S2"]["matched_triplets"]["0.9000\t5 (A) A -> 6 (A) A"] == [
        "1 A 0.7000 "
    ]


def test_triplet_only_calculation_yields_an_empty_dataset(tmp_path):
    parser = parse(tmp_path, fixtures.ccsd_log(singlets=[], triplets=None))
    assert "singlet-" not in "".join(parser.irreps)
    assert parser.create_dataset() == {}


def test_dataset_records_the_state_index_not_the_transition_index(tmp_path):
    lines = fixtures.ccsd_log(
        singlets=[
            (1, 9.0, -76.0, [(0.9, 5, 6), (0.3, 4, 7)]),
            (2, 8.0, -76.1, [(0.6, 5, 8)]),
        ],
        triplets=[(1, 6.0, -76.2, [(0.7, 5, 6)])],
    )
    data = parse(tmp_path, lines).create_dataset()
    assert [data[key]["idx"] for key in ("S1", "S2")] == [0, 1]


def test_write_dataset_lists_matches_and_states_the_empty_case(tmp_path):
    parser = parse(
        tmp_path,
        fixtures.ccsd_log(
            singlets=[(1, 8.0, -76.0, [(0.9, 5, 6)]), (2, 7.0, -76.1, [(0.2, 9, 12)])],
            triplets=[(1, 6.0, -76.2, [(0.7, 5, 6)])],
        ),
        # the 0.2 transition of state 2 is filtered out, so it has nothing left
        threshold=0.45,
    )
    out = tmp_path / "happy.txt"
    parser.write_dataset(out)
    text = out.read_text()

    assert text.startswith("---- Happy family ----\n")
    assert text.endswith("--- End of the happy family :) ---\n")
    assert "------------S1------------" in text
    assert "Matched triplets:\n" in text
    assert "Matched triplets: none" in text
    assert "1 A 0.7000" in text


def test_select_lowest_excitations_handles_missing_states():
    assert BaseParser.select_lowest_excitations([]) == []

    short = Irrep("A", "singlet", 5, parent="x")
    short.n_states = 5
    assert BaseParser.select_lowest_excitations([short]) == []


def test_select_lowest_excitations_uses_the_shortest_irrep():
    """The depth is what was parsed, not what the input asked for.

    ``EE_SINGLETS [5]`` with only three converged roots used to raise
    ``IndexError`` (or ``ValueError: min() arg is an empty sequence`` when a
    whole multiplicity was missing) instead of comparing what is available.
    """
    from xtraee.trblock import CCSDTransitionBlock

    def make(name, energies, requested=5):
        irrep = Irrep(name, "singlet", requested, parent="x")
        for i, energy in enumerate(energies, start=1):
            block = CCSDTransitionBlock(i, name, "singlet")
            block.excitation_energy = energy
            irrep.append(block)
        return irrep

    irrep_a = make("A", [8.0, 7.0, 6.0])
    irrep_b = make("B", [7.5])
    lowest = BaseParser.select_lowest_excitations([irrep_a, irrep_b])

    assert len(lowest) == 1
    assert lowest[0].excitation_energy == 7.5


def test_requested_state_count_mismatch_warns_instead_of_raising(tmp_path, logs):
    lines = fixtures.ccsd_log(
        singlets=[(1, 8.0, -76.0, [(0.9, 5, 6)])],
        triplets=[(1, 6.0, -76.2, [(0.7, 5, 6)])],
        requested_singlets=3,
    )
    parser = parse(tmp_path, lines)
    assert isinstance(parser, QCCSDParser)
    assert "3 states were requested but 1 were parsed" in logs.text


def test_missing_rem_counts_fall_back_to_the_parsed_count(tmp_path):
    """``EE_SINGLETS`` absent (or too short) used to raise IndexError."""
    parser = parse(tmp_path, fixtures.ccsd_log(omit_rem_counts=True))
    assert {key: (len(irrep), irrep.n_states) for key, irrep in parser.irreps.items()} == {
        "singlet-A": (2, 2),
        "triplet-A": (2, 2),
    }


def test_several_irreps_with_a_single_requested_count(tmp_path):
    parser = parse(tmp_path, fixtures.ccsd_log(irreps=("A", "B"), requested_singlets=2))
    assert set(parser.irreps) == {
        "singlet-A",
        "triplet-A",
        "singlet-B",
        "triplet-B",
    }
    assert all(len(irrep) == 2 for irrep in parser.irreps.values())
    assert all(irrep.n_states == 2 for irrep in parser.irreps.values())


@pytest.mark.parametrize("threshold, expected", [(0.0, 2), (0.25, 1)])
def test_threshold_filters_ccsd_transitions(tmp_path, threshold, expected):
    """``--threshold`` was only wired into the CIS blocks.

    An EOM-CCSD run parsed every printed amplitude whatever ``--threshold``
    said, so the documented flag silently did nothing for the main method.
    """
    parser = parse(
        tmp_path,
        fixtures.ccsd_log(
            singlets=[(1, 8.0, -76.0, [(0.9, 5, 6), (0.3, 4, 7)])],
            triplets=[(1, 6.0, -76.2, [(0.7, 5, 6)])],
        ),
        threshold=threshold,
    )
    block = parser.find_trblock("singlet-1/A")
    assert len(block.transitions) == expected
    assert block.threshold == threshold
