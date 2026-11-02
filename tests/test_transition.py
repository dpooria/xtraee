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


def test_unparsed_transition_is_not_equal_to_anything():
    """A line without recognizable orbitals used to match *every* transition.

    ``is_equal`` walked the (empty) orbital lists with ``zip``, so it reported
    equality for any two unparsed transitions and the accuracy of a state was
    computed against transitions that were never matched.
    """
    broken = CCSDTransition.from_str("0.5000        ->    16 (A) A")
    other = CCSDTransition.from_str("0.5000        ->    16 (A) A")
    good = CCSDTransition.from_str("0.5000       15 (A) A    ->    16 (A) A")

    assert broken.parsed is False
    assert broken.is_equal(other) is False
    assert broken.is_equal(good) is False
    assert good.is_equal(good) is True


def test_double_and_single_never_match():
    single = CCSDTransition.from_str("0.5  15 (A) A -> 16 (A) A")
    double = CCSDTransition.from_str(
        "0.5  15 (A) A  17 (A) B -> 16 (A) A  18 (A) B"
    )
    assert single.is_double is False
    assert double.is_double is True
    assert single.is_equal(double) is False


def test_cistransition_to_std_shifts_empty_orbitals_by_homo():
    transition = CISTransition.from_str(
        "     D(    5) --> V(    1) amplitude = 0.9100"
    )
    standard = transition.to_std(homo=5)
    assert str(standard.id_i[0]) == "5 (A) A"
    assert str(standard.id_f[0]) == "6 (A) A"


def test_cistransition_to_std_does_not_crash_without_orbitals():
    """``to_std`` used to index ``id_i[0]`` and raise IndexError."""
    transition = CISTransition(0.5, "no orbital here", "nor here")
    assert transition.parsed is False
    standard = transition.to_std(homo=5)
    assert standard.parsed is False


def test_cisd_to_ccsd_string_has_no_stray_parenthesis():
    transition = CISDTransition.from_str(
        "  0.6397                 12(   A1) B   ->    0(   A1) B "
    )
    standard = transition.to_ccsd(homo=6)
    assert standard.initial.count("(") == standard.initial.count(")")
    assert str(standard.id_i[0]) == "13 (A) B"
    assert str(standard.id_f[0]) == "7 (A) B"

