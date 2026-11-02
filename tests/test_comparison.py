"""The accuracy metrics and the guards around them."""

import pytest

from xtraee.transition import CCSDTransition
from xtraee.trblock import CCSDTransitionBlock


def block(*transitions: str, threshold: float = 0.0) -> CCSDTransitionBlock:
    tb = CCSDTransitionBlock(1, "A", "singlet", threshold=threshold)
    tb.transitions = [CCSDTransition.from_str(line) for line in transitions]
    return tb


def test_cosine_scores_identical_states_exactly_one():
    """``inner-prod`` does not: it sums |a_i||b_i|, i.e. sum(a_i**2) here."""
    state = block(
        "0.9000  5 (A) A -> 6 (A) A",
        "0.2000  4 (A) A -> 7 (A) A",
    )
    accuracy, error, matched = state.compare(state, "cosine")
    assert accuracy == pytest.approx(1.0)
    assert error == pytest.approx(0.0)
    assert matched == 1.0

    inner, _, _ = state.compare(state, "inner-prod")
    assert inner == pytest.approx(0.9**2 + 0.2**2)
    assert inner < 1.0


def test_cosine_ignores_the_absolute_scale_of_the_printed_amplitudes():
    left = block("0.9000  5 (A) A -> 6 (A) A")
    right = block("0.4500  5 (A) A -> 6 (A) A")
    accuracy, _, _ = left.compare(right, "cosine")
    # same direction, different magnitude
    assert accuracy == pytest.approx(1.0)


def test_cosine_of_orthogonal_states_is_zero():
    left = block("1.0000  5 (A) A -> 6 (A) A")
    right = block("1.0000  4 (A) A -> 7 (A) A")
    accuracy, error, matched = left.compare(right, "cosine")
    assert accuracy == pytest.approx(0.0)
    assert error == pytest.approx(1.0)
    assert matched == 0.0


def test_cosine_uses_the_shared_transitions_only():
    left = block(
        "0.9000  5 (A) A -> 6 (A) A",
        "0.4000  4 (A) A -> 7 (A) A",
    )
    right = block("0.9000  5 (A) A -> 6 (A) A")
    accuracy, _, matched = left.compare(right, "cosine")
    assert matched == pytest.approx(0.5)
    # the unmatched 0.4 transition still counts as a zero in the other vector
    expected = 0.81 / (0.9 * (0.81 + 0.16) ** 0.5)
    assert accuracy == pytest.approx(expected)


def test_wabs_relies_on_the_shared_weight():
    left = block("0.9000  5 (A) A -> 6 (A) A", "0.2000  4 (A) A -> 7 (A) A")
    right = block("0.9000  5 (A) A -> 6 (A) A")
    accuracy, error, _ = left.compare(right, "w-abs")
    assert 0.0 < accuracy < 1.0
    assert error == pytest.approx(1.0 - accuracy)


def test_pearson_does_not_return_nan_for_a_single_shared_transition():
    left = block("0.9000  5 (A) A -> 6 (A) A", "0.2000  4 (A) A -> 7 (A) A")
    right = block("0.9000  5 (A) A -> 6 (A) A")
    accuracy, _, _ = left.compare(right, "pearson")
    assert accuracy == accuracy  # not nan


@pytest.mark.parametrize("method", ["inner-prod", "cosine", "w-abs", "pearson"])
def test_empty_blocks_compare_to_nothing(method):
    empty = CCSDTransitionBlock(1, "A", "singlet")
    full = block("0.9000  5 (A) A -> 6 (A) A")
    assert empty.compare(full, method) == (0.0, 1.0, 0.0)
    assert full.compare(empty, method) == (0.0, 1.0, 0.0)


def test_unknown_method_raises():
    left = block("0.9000  5 (A) A -> 6 (A) A")
    with pytest.raises(ValueError, match="Method not recognized"):
        left.compare(left, "made-up")


def test_threshold_is_honoured_when_building_the_block():
    tb = block(
        "0.9000  5 (A) A -> 6 (A) A",
        "0.0400  4 (A) A -> 7 (A) A",
        threshold=0.2,
    )
    kept = [tr for tr in tb.transitions if tr.probability >= tb.threshold]
    assert len(kept) == 1


def test_compare_with_a_cisd_reference_uses_the_reference_as_state_one():
    """CCSD blocks always act as the reference of a cross-method comparison."""
    from xtraee.trblock import CISDTransitionBlock
    from xtraee.transition import CISDTransition

    ccsd = block("0.9000  15 (A) A -> 16 (A) A")
    cisd = CISDTransitionBlock(1, "A", "singlet", "CISD", threshold=0.0)
    cisd.homo = 14
    cisd.transitions = [
        CISDTransition.from_str("0.9000   14(   A) A   ->    1(   A) A")
    ]
    forward = ccsd.compare(cisd, "cosine")
    backward = cisd.compare(ccsd, "cosine")
    assert forward == backward
