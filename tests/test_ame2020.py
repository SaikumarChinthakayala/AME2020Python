"""
Tests for ame2020. Run with: pytest

These tests run against the bundled data file. Since the placeholder data
currently only contains n, H, He, Li, Be (A=1..6), tests are written against
that range. If/when the full mass.mas20 is swapped in, these should still
pass unchanged (they don't assert anything about isotopes outside that
range).
"""

import math

import pytest

from ame2020 import (
    Isotope,
    IsotopeNotFoundError,
    get_beta_decay_energy,
    get_binding_energy,
    get_binding_per_A,
    get_mass,
    get_mass_excess,
)


def test_basic_lookup_normal_value():
    # 2H (deuteron): mass_excess = 13135.722895 keV, not estimated
    iso = Isotope("H", 2)
    assert iso.A == 2
    assert iso.Z == 1
    assert iso.N == 1
    assert math.isclose(iso.mass_excess, 13135.722895, rel_tol=1e-9)
    assert iso.mass_excess_estimated is False
    assert math.isclose(iso.mass_excess_unc, 0.000015, rel_tol=1e-6)


def test_atomic_mass_concatenation():
    # 2H atomic mass = 2.014101777844 u (int part "2" + frac "014101.777844")
    iso = Isotope("H", 2)
    assert math.isclose(iso.mass, 2.014101777844, rel_tol=1e-12)
    assert math.isclose(iso.mass_unc, 0.000015 / 1e6, rel_tol=1e-6)


def test_neutron_lookup():
    n = Isotope("n", 1)
    assert n.Z == 0
    assert n.N == 1
    assert math.isclose(n.mass, 1.00866491590, rel_tol=1e-9)


def test_case_insensitivity_with_known_element():
    lower = Isotope("he", 4)
    upper = Isotope("HE", 4)
    mixed = Isotope("He", 4)
    assert lower.mass == upper.mass == mixed.mass


def test_estimated_value_flagged():
    # 5Li mass excess: "11678.887" with unc "50.000" -- NOT marked with '#' itself,
    # but 5Be (A=5) IS estimated. Use 5Be o='x' which has '#' markers.
    be5 = Isotope("Be", 5)
    assert be5.mass_excess_estimated is True
    assert math.isclose(be5.mass_excess, 37139.0, rel_tol=1e-9)


def test_estimated_binding_energy():
    li3 = Isotope("Li", 3)
    assert li3.binding_per_A_estimated is True
    assert math.isclose(li3.binding_per_A, -2267.0, rel_tol=1e-9)


def test_not_calculable_beta_energy_is_nan():
    # 1H beta decay energy is '*' -> NaN
    h1 = Isotope("H", 1)
    assert math.isnan(h1.beta_energy)


def test_negative_beta_energy():
    he3 = Isotope("He", 3)
    assert he3.beta_energy < 0
    assert math.isclose(he3.beta_energy, -13736.0, rel_tol=1e-6)
    assert he3.beta_energy_estimated is True


def test_binding_energy_total():
    li6 = Isotope("Li", 6)
    expected_total = li6.binding_per_A * 6
    assert math.isclose(li6.binding_energy(), expected_total, rel_tol=1e-12)


def test_isotope_not_found_raises():
    with pytest.raises(IsotopeNotFoundError):
        Isotope("Xx", 999)


def test_isotope_not_found_is_keyerror_subclass():
    with pytest.raises(KeyError):
        Isotope("Xx", 999)


def test_symbol_property():
    iso = Isotope("Li", 6)
    assert iso.symbol == "6Li"


def test_is_estimated_aggregate_true():
    be5 = Isotope("Be", 5)
    assert be5.is_estimated is True


def test_is_estimated_aggregate_false():
    h2 = Isotope("H", 2)
    assert h2.is_estimated is False


def test_functional_get_mass_matches_object():
    assert get_mass("Li", 6) == Isotope("Li", 6).mass


def test_functional_get_mass_with_flag():
    mass, est = get_mass("Be", 5, with_flag=True)
    assert est is True
    assert math.isclose(mass, Isotope("Be", 5).mass, rel_tol=1e-12)


def test_functional_get_mass_excess():
    assert get_mass_excess("He", 6) == Isotope("He", 6).mass_excess


def test_functional_get_binding_per_A():
    assert get_binding_per_A("He", 6) == Isotope("He", 6).binding_per_A


def test_functional_get_binding_energy_total():
    iso = Isotope("He", 6)
    assert math.isclose(get_binding_energy("He", 6), iso.binding_per_A * iso.A, rel_tol=1e-12)


def test_functional_get_beta_decay_energy():
    assert get_beta_decay_energy("He", 6) == Isotope("He", 6).beta_energy


def test_repr_contains_symbol_and_mass():
    iso = Isotope("Li", 6)
    r = repr(iso)
    assert "6Li" in r
    assert "estimated" in r


def test_equality_and_hash():
    a = Isotope("Li", 6)
    b = Isotope("Li", 6)
    c = Isotope("He", 6)
    assert a == b
    assert a != c
    assert hash(a) == hash(b)
