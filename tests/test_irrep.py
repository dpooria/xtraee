"""Irrrep helpers: sorting, identifiers and the state-by-state matrix."""

import pytest

from xtraee.irrep import Irrep
from xtraee.parser.base import BaseParser
from xtraee.transition import CCSDTransition
from xtraee.trblock import CCSDTransitionBlock


def make_irrep(name="A", multi="singlet", energies=(8.0, 7.0), n_states=None):
    irrep = Irrep(name, multi, len(energies) if n_states is None else n_states,
                  parent="CCSD")
    for i, energy in enumerate(energies, start=1):
        block = CCSDTransitionBlock(i, name, multi)
        block.excitation_energy = energy
        block.transitions = [
            CCSDTransition.from_str(f"0.9000  {i} (A) A -> {i + 10} (A) A")
        ]
        irrep.append(block)
    return irrep


def test_identifier_carries_the_parent():
    assert make_irrep("A1").identifier == "CCSD-A1.singlet"


def test_sort_orders_states_by_excitation_energy():
    irrep = make_irrep(energies=(8.0, 7.0))
    irrep.sort()
    assert [block.identifier for block in irrep.trblocks] == [
        "singlet-2/A",
        "singlet-1/A",
    ]


def test_compare_builds_a_state_matrix():
    left = make_irrep(energies=(8.0, 7.0))
    right = make_irrep(energies=(7.5,))
    frame = left.compare(right, "cosine")
    assert frame.shape == (2, 1)
    assert frame.columns == ["singlet-1/A"]
    assert frame.index[0] == "singlet-1/A"


@pytest.mark.parametrize("empty_side", ["self", "other"])
def test_compare_without_states_returns_an_empty_labelled_frame(empty_side):
    """An irrep that produced no states used to raise ``UnboundLocalError``.

    ``columns`` was only assigned inside the loop over the other irrep, so a
    ``compareall`` on a file whose triplet section was empty (or a singlet-only
    calculation) died with ``UnboundLocalError: cannot access local variable
    'columns'`` instead of writing an empty comparison.
    """
    full = make_irrep(energies=(8.0,))
    empty = Irrep("A", "triplet", 0, parent="CCSD")

    left, right = (empty, full) if empty_side == "self" else (full, empty)
    frame = left.compare(right, "cosine")
    assert frame.empty
    if empty_side == "self":
        assert frame.index.empty and list(frame.columns) == ["singlet-1/A"]
    else:
        assert list(frame.index) == ["singlet-1/A"] and frame.columns.empty


def test_select_lowest_excitations_requires_no_states():
    assert BaseParser.select_lowest_excitations([]) == []
