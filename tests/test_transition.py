import re

import pytest
from xtraee.transition import CCSDTransition, CISTransition


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
    print(transition.id_i)
    print(transition.id_f)
    print(len(transition.id_i))
    print(len(transition.id_f))
    assert transition.id_f == [(1, "B3g", "A"), (1, "B3g", "B")]
    assert transition.id_i == [(1, "B2u", "A"), (1, "B2u", "B")]
    assert transition.is_double


# def test_cistransition_double():
#     transition = CISTransition.from_str(
#         "  D(   15) H (   2  )--> V(    1)  V(16      ) amplitude =  0.9698        "
#     )
#     print()
#     print(transition.id_i)
#     print(transition.id_f)
#     assert transition.amplitude == pytest.approx(0.9698)
#     assert transition.id_i == [(15, "A", "A"), ("H", "2")]
#     assert transition.id_f == [("V", "1"), ("G", "16")]
#     assert transition.is_double
