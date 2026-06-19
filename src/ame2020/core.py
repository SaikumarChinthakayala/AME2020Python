"""
Core lookup API for AME2020 isotope data.

Two equivalent ways to use this:

    >>> from ame2020 import Isotope
    >>> pd102 = Isotope("Pd", 102)
    >>> pd102.mass            # atomic mass in u (default)
    >>> pd102.mass_excess     # keV
    >>> pd102.binding_per_A   # keV

    >>> from ame2020 import get_mass, get_mass_excess, get_binding_per_A
    >>> get_mass("Pd", 102)
    >>> get_mass_excess("Pd", 102)
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from .parser import load_default_table

# physical constant: 1 u in keV/c^2 (CODATA 2018), used if you want to
# convert mass <-> mass-excess-style quantities yourself.
U_TO_KEV = 931494.10242


class IsotopeNotFoundError(KeyError):
    """Raised when the requested (element, A) combination is not in the table."""


def _normalize_element(element: str) -> str:
    element = element.strip()
    # AME uses 'n' (lowercase) for the free neutron; everything else is
    # Title-cased, e.g. 'pd' -> 'Pd', 'HE' -> 'He'.
    if element.lower() == "n":
        return "n"
    return element[:1].upper() + element[1:].lower()


class Isotope:
    """A single isotope's AME2020 data, looked up by element symbol and mass number.

    Parameters
    ----------
    element : str
        Element symbol, e.g. "Pd", "Rh". Case-insensitive.
    A : int
        Mass number.
    table : pandas.DataFrame, optional
        Use a custom parsed table instead of the bundled default
        (mainly useful for testing).

    Examples
    --------
    >>> iso = Isotope("Pd", 102)
    >>> iso.mass
    101.905609...
    >>> iso.mass_excess, iso.mass_excess_unc, iso.mass_excess_estimated
    (-87.917, 0.0023, False)
    """

    __slots__ = (
        "element", "A", "Z", "N", "origin",
        "mass_excess", "mass_excess_unc", "mass_excess_estimated",
        "binding_per_A", "binding_per_A_unc", "binding_per_A_estimated",
        "beta_type", "beta_energy", "beta_energy_unc", "beta_energy_estimated",
        "mass", "mass_unc", "mass_estimated",
    )

    def __init__(self, element: str, A: int, table: Optional[pd.DataFrame] = None):
        df = table if table is not None else load_default_table()
        element_norm = _normalize_element(element)
        match = df[(df["element"] == element_norm) & (df["A"] == A)]

        if match.empty:
            raise IsotopeNotFoundError(
                f"No AME2020 entry found for element={element!r}, A={A}. "
                "Check the element symbol (e.g. 'Pd' not 'pd') and mass number."
            )
        if len(match) > 1:
            # Shouldn't happen in practice (element+A is unique), but guard anyway.
            match = match.iloc[[0]]

        row = match.iloc[0]

        self.element = row["element"]
        self.A = int(row["A"])
        self.Z = int(row["Z"])
        self.N = int(row["N"])
        self.origin = row["origin"]

        self.mass_excess = row["mass_excess_keV"]
        self.mass_excess_unc = row["mass_excess_unc_keV"]
        self.mass_excess_estimated = bool(row["mass_excess_estimated"])

        self.binding_per_A = row["binding_per_A_keV"]
        self.binding_per_A_unc = row["binding_per_A_unc_keV"]
        self.binding_per_A_estimated = bool(row["binding_per_A_estimated"])

        self.beta_type = row["beta_type"]
        self.beta_energy = row["beta_energy_keV"]
        self.beta_energy_unc = row["beta_energy_unc_keV"]
        self.beta_energy_estimated = bool(row["beta_energy_estimated"])

        self.mass = row["atomic_mass_u"]
        self.mass_unc = row["atomic_mass_unc_u"]
        self.mass_estimated = bool(row["atomic_mass_estimated"])

    @property
    def symbol(self) -> str:
        """Isotope label like '102Pd'."""
        return f"{self.A}{self.element}"

    @property
    def is_estimated(self) -> bool:
        """True if any of the core observables (mass excess, binding energy,
        atomic mass) is flagged as estimated rather than experimental."""
        return self.mass_excess_estimated or self.binding_per_A_estimated or self.mass_estimated

    def binding_energy(self) -> float:
        """Total binding energy in keV (binding_per_A * A)."""
        return self.binding_per_A * self.A

    def __repr__(self) -> str:
        return f"Isotope({self.symbol!r}, mass={self.mass:.6f} u, estimated={self.is_estimated})"

    def __eq__(self, other) -> bool:
        if not isinstance(other, Isotope):
            return NotImplemented
        return self.element == other.element and self.A == other.A

    def __hash__(self) -> int:
        return hash((self.element, self.A))


# ---------------------------------------------------------------------------
# Functional wrappers (thin convenience layer over Isotope)
# ---------------------------------------------------------------------------


def get_mass(element: str, A: int, with_flag: bool = False):
    """Atomic mass in u. If with_flag=True, returns (mass, is_estimated)."""
    iso = Isotope(element, A)
    return (iso.mass, iso.mass_estimated) if with_flag else iso.mass


def get_mass_excess(element: str, A: int, with_flag: bool = False):
    """Mass excess in keV. If with_flag=True, returns (value, is_estimated)."""
    iso = Isotope(element, A)
    return (iso.mass_excess, iso.mass_excess_estimated) if with_flag else iso.mass_excess


def get_binding_per_A(element: str, A: int, with_flag: bool = False):
    """Binding energy per nucleon in keV. If with_flag=True, returns (value, is_estimated)."""
    iso = Isotope(element, A)
    return (
        (iso.binding_per_A, iso.binding_per_A_estimated)
        if with_flag
        else iso.binding_per_A
    )


def get_binding_energy(element: str, A: int, with_flag: bool = False):
    """Total binding energy in keV (binding_per_A * A)."""
    iso = Isotope(element, A)
    total = iso.binding_energy()
    return (total, iso.binding_per_A_estimated) if with_flag else total


def get_beta_decay_energy(element: str, A: int, with_flag: bool = False):
    """Beta-decay energy in keV (sign/type given by .beta_type, e.g. 'B-')."""
    iso = Isotope(element, A)
    return (iso.beta_energy, iso.beta_energy_estimated) if with_flag else iso.beta_energy
