"""
Tests for ame2020.isomers (NUBASE2020-based isomer/excited-state lookup).

Run against the bundled placeholder nubase.mas20 (A=1..14 sample). Like
test_ame2020.py, these are written so they keep passing once the full
NUBASE2020 file is swapped in.
"""

import math

import pytest

from ame2020 import (
    Isomer,
    IsomerNotFoundError,
    get_excitation_energy,
    get_isomer_mass,
    list_states,
)


def test_ground_state_default():
    gs = Isomer("Li", 10)
    assert gs.level_index == 0
    assert gs.state_label == ""
    assert gs.is_isomer is False
    assert math.isclose(gs.mass_excess, 33053.0, rel_tol=1e-9)
    assert math.isnan(gs.excitation_energy)


def test_ground_state_explicit_gs_string():
    gs1 = Isomer("Li", 10)
    gs2 = Isomer("Li", 10, state="gs")
    gs3 = Isomer("Li", 10, state="ground")
    gs4 = Isomer("Li", 10, state=0)
    assert gs1.mass_excess == gs2.mass_excess == gs3.mass_excess == gs4.mass_excess


def test_first_isomer_m():
    m = Isomer("Li", 10, state="m")
    assert m.level_index == 1
    assert m.state_label == "m"
    assert m.is_isomer is True
    assert math.isclose(m.mass_excess, 33250.0, rel_tol=1e-9)
    assert math.isclose(m.excitation_energy, 200.0, rel_tol=1e-9)
    assert math.isclose(m.excitation_energy_unc, 40.0, rel_tol=1e-9)
    assert m.half_life == "3.7"
    assert m.half_life_unit == "zs"


def test_second_isomer_n():
    n = Isomer("Li", 10, state="n")
    assert n.level_index == 2
    assert n.state_label == "n"
    assert n.is_isomer is True
    assert math.isclose(n.excitation_energy, 480.0, rel_tol=1e-9)


def test_isomer_by_integer_level_index():
    by_label = Isomer("Li", 10, state="m")
    by_index = Isomer("Li", 10, state=1)
    assert by_label.mass_excess == by_index.mass_excess


def test_symbol_property():
    gs = Isomer("Li", 10)
    m = Isomer("Li", 10, state="m")
    assert gs.symbol == "10Li"
    assert m.symbol == "10Li-m"


def test_mass_property_consistency():
    # mass[u] = A + mass_excess[keV] / U_TO_KEV
    m = Isomer("Be", 12, state="m")
    expected = 12 + 27328.8 / 931494.10242
    assert math.isclose(m.mass, expected, rel_tol=1e-9)


def test_levels_excluded_from_is_isomer():
    # 6Li has a level_index=8 'i' row (a level, not a true isomer)
    level = Isomer("Li", 6, state="i")
    assert level.is_isomer is False
    assert level.level_index == 8


def test_case_insensitive_element_and_state():
    a = Isomer("li", 10, state="M")
    b = Isomer("LI", 10, state="m")
    assert a.mass_excess == b.mass_excess
    assert a.symbol == b.symbol


def test_unknown_state_raises_with_helpful_message():
    with pytest.raises(IsomerNotFoundError) as excinfo:
        Isomer("Li", 10, state="q")
    assert "Available states" in str(excinfo.value)


def test_unknown_nuclide_raises():
    with pytest.raises(IsomerNotFoundError):
        Isomer("Xx", 999)


def test_unknown_nuclide_is_keyerror_subclass():
    with pytest.raises(KeyError):
        Isomer("Xx", 999)


def test_list_states_returns_all_rows():
    df = list_states("Li", 10)
    assert len(df) == 3
    assert set(df["state_label"]) == {"", "m", "n"}


def test_list_states_unknown_nuclide_raises():
    with pytest.raises(IsomerNotFoundError):
        list_states("Xx", 999)


def test_functional_get_isomer_mass_default_state_m():
    mass = get_isomer_mass("Li", 10)  # default state='m'
    assert mass == Isomer("Li", 10, state="m").mass


def test_functional_get_isomer_mass_with_flag():
    mass, est = get_isomer_mass("Li", 10, state="n", with_flag=True)
    assert isinstance(est, bool)
    assert mass == Isomer("Li", 10, state="n").mass


def test_functional_get_excitation_energy():
    exc = get_excitation_energy("Li", 10, state="m")
    assert math.isclose(exc, 200.0, rel_tol=1e-9)


def test_repr_contains_symbol():
    m = Isomer("Li", 10, state="m")
    r = repr(m)
    assert "10Li-m" in r


def test_single_state_nuclide_has_no_isomer():
    # 1H in the sample only has a ground state -- requesting 'm' should fail clearly
    with pytest.raises(IsomerNotFoundError):
        Isomer("H", 1, state="m")
