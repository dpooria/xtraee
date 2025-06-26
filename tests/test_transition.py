import re

import pytest
from xtraee.transition import CCSDTransition, CISTransition, CISDTransition


@pytest.mark.parametrize(
    "transition_str, amplitude, transition_str_clean",
    [
        (
            "  0.6228       13 (A) A    ->    16 (A) A                   ",
            0.6228,
            "0.6228 13 (A) A -> 16 (A) A",
        ),
    ],
)
def test_ccsdtransition_fromstr(transition_str, amplitude, transition_str_clean):
    transition = CCSDTransition.from_str(transition_str)
    assert transition.amplitude == pytest.approx(amplitude)
    transition_repr = re.sub(r"\s+", " ", transition.__repr__()).strip()
    print(transition_repr, "\t==\t", transition_str_clean)
    assert transition_repr == transition_str_clean


@pytest.mark.parametrize(
    "transition_str, amplitude, transition_str_clean",
    [
        (
            "     D(   14) --> V(    2) amplitude = -0.1002",
            -0.1002,
            "-0.1002 14 (A) A -> 2 (A) A",
        ),
    ],
)
def test_cistransition_fromstr(transition_str, amplitude, transition_str_clean):
    transition = CISTransition.from_str(transition_str)
    assert transition.amplitude == pytest.approx(amplitude)
    transition_repr = re.sub(r"\s+", " ", transition.__repr__()).strip()
    print(transition_repr, "\t==\t", transition_str_clean)
    assert transition_repr == transition_str_clean


def test_ccsdtransition_double():
    transition = CCSDTransition.from_str(
        "  0.4088       1 (B2u) A     1 (B2u) B   ->    1 (B3g) A     1 (B3g) B"
    )
    assert transition.id_f == [(1, "B3g", "A"), (1, "B3g", "B")]
    assert transition.id_i == [(1, "B2u", "A"), (1, "B2u", "B")]
    assert transition.is_double


def test_cisd_transition():
    transition = CISDTransition.from_str(
        "  0.6397                 12(   A1) B   ->    0(   A1) B "
    )
    assert transition.id_i == [(12, "A1", "B")]
    assert transition.id_f == [(0, "A1", "B")]
    assert transition.amplitude == pytest.approx(0.6397)
