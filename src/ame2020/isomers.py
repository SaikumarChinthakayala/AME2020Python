"""
Isomer-aware lookup, combining:
  - NUBASE2020 (nubase.mas20) for excitation energies and isomer labels
  - AME2020 (mass.mas20), via ame2020.core, for ground-state masses (when
    you want an absolute mass rather than just an excitation energy)

NUBASE's own "Mass Excess" column on an isomer row is already the isomer's
own absolute mass excess (not a delta) -- so for the headline numbers
(mass excess, mass in u) we read it directly from NUBASE. The "Exc" column
(excitation energy above the ground state) is kept alongside as it's
independently useful and is sometimes known more precisely than the
isomer's absolute mass excess.
"""

from __future__ import annotations

import math
from typing import Optional

import pandas as pd

from .core import _normalize_element
from .nubase_parser import load_default_nubase_table

U_TO_KEV = 931494.10242  # CODATA 2018, keV per u; same constant as ame2020.core


class IsomerNotFoundError(KeyError):
    """Raised when the requested (element, A, state) combination is not in NUBASE."""


# Map convenient user-facing keys to NUBASE's level_index / state_label values.
_STATE_ALIASES = {
    "gs": 0,
    "ground": 0,
    "m": "m",
    "n": "n",
    "m1": "m",
    "m2": "n",
    1: 1,
    2: 2,
}


class Isomer:
    """A single nuclear state (ground state, isomer, level, or IAS) from NUBASE2020.

    Parameters
    ----------
    element : str
        Element symbol, e.g. "Pd", "Rh". Case-insensitive.
    A : int
        Mass number.
    state : int, str, or None
        Which state to retrieve:
          - None or "gs"/"ground"/0  -> ground state
          - "m"                       -> first (lowest) isomer
          - "n"                       -> second isomer
          - 1 or 2                    -> isomer by NUBASE level_index directly
          - any other single character (e.g. "p", "i") -> matched against
            NUBASE's raw state_label, for levels/IAS rows
    table : pandas.DataFrame, optional
        Use a custom parsed NUBASE table instead of the bundled default.

    Examples
    --------
    >>> gs = Isomer("Li", 10)            # ground state
    >>> m  = Isomer("Li", 10, state="m") # first isomer (10Li-m)
    >>> m.excitation_energy              # keV above the ground state
    >>> m.mass_excess                    # the ISOMER's own absolute mass excess, keV
    >>> m.half_life, m.half_life_unit
    """

    __slots__ = (
        "element", "A", "Z", "level_index", "state_label",
        "mass_excess", "mass_excess_unc", "mass_excess_estimated",
        "excitation_energy", "excitation_energy_unc", "excitation_energy_estimated",
        "half_life", "half_life_unit", "Jpi", "is_isomer",
    )

    def __init__(self, element: str, A: int, state=None, table: Optional[pd.DataFrame] = None):
        df = table if table is not None else load_default_nubase_table()
        element_norm = _normalize_element(element)

        candidates = df[(df["element"] == element_norm) & (df["A"] == A)]
        if candidates.empty:
            raise IsomerNotFoundError(
                f"No NUBASE2020 entries found for element={element!r}, A={A}."
            )

        row = self._select_state(candidates, state, element, A)

        self.element = row["element"]
        self.A = int(row["A"])
        self.Z = int(row["Z"])
        self.level_index = int(row["level_index"])
        self.state_label = row["state_label"]

        self.mass_excess = row["mass_excess_keV"]
        self.mass_excess_unc = row["mass_excess_unc_keV"]
        self.mass_excess_estimated = bool(row["mass_excess_estimated"])

        self.excitation_energy = row["exc_energy_keV"]
        self.excitation_energy_unc = row["exc_energy_unc_keV"]
        self.excitation_energy_estimated = bool(row["exc_energy_estimated"])

        self.half_life = row["half_life"]
        self.half_life_unit = row["half_life_unit"]
        self.Jpi = row["Jpi"]
        self.is_isomer = bool(row["is_isomer"])

    @staticmethod
    def _select_state(candidates: pd.DataFrame, state, element: str, A: int) -> pd.Series:
        if state is None:
            state = "gs"

        key = state.lower() if isinstance(state, str) else state
        resolved = _STATE_ALIASES.get(key, key)

        if isinstance(resolved, int):
            match = candidates[candidates["level_index"] == resolved]
        else:
            # match against the raw single-character state label (m, n, p, q, i, j, ...)
            match = candidates[candidates["state_label"].str.lower() == str(resolved).lower()]

        if match.empty:
            available = sorted(
                candidates["state_label"].replace("", "gs").unique().tolist()
            )
            raise IsomerNotFoundError(
                f"No state '{state}' found for {element}-{A}. "
                f"Available states: {available}"
            )
        return match.iloc[0]

    @property
    def symbol(self) -> str:
        suffix = f"-{self.state_label}" if self.state_label else ""
        return f"{self.A}{self.element}{suffix}"

    @property
    def has_known_mass(self) -> bool:
        """False when NUBASE has no mass-excess value at all for this state
        (seen in practice for some isomers where only half-life/spin-parity
        are established but the mass excess is blank in the file, often
        alongside an excitation energy marked 'non-exist'). Check this
        before trusting `.mass` / `.mass_excess` -- they'll be NaN otherwise,
        but a NaN alone doesn't tell you whether that's because the state
        truly has no evaluated mass, vs. some other parsing gap.
        """
        return not math.isnan(self.mass_excess)

    @property
    def mass(self) -> float:
        """Absolute mass in u, derived from this state's own mass excess.

        mass[u] = A + mass_excess[keV] / U_TO_KEV
        (matches the standard mass-excess definition: ME = (M - A*u) * c^2)
        """
        return self.A + self.mass_excess / U_TO_KEV

    def __repr__(self) -> str:
        return (
            f"Isomer({self.symbol!r}, mass_excess={self.mass_excess:.3f} keV, "
            f"exc={self.excitation_energy} keV, is_isomer={self.is_isomer})"
        )


def list_states(element: str, A: int, table: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """List all NUBASE2020 states (ground state, isomers, levels, IAS) for a
    given element/A, as a small DataFrame for quick inspection.
    """
    df = table if table is not None else load_default_nubase_table()
    element_norm = _normalize_element(element)
    candidates = df[(df["element"] == element_norm) & (df["A"] == A)]
    if candidates.empty:
        raise IsomerNotFoundError(f"No NUBASE2020 entries found for element={element!r}, A={A}.")
    cols = [
        "level_index", "state_label", "mass_excess_keV", "exc_energy_keV",
        "half_life", "half_life_unit", "Jpi", "is_isomer",
    ]
    return candidates[cols].reset_index(drop=True)


def get_isomer_mass(element: str, A: int, state="m", with_flag: bool = False):
    """Convenience function: absolute mass (u) of a specific isomer."""
    iso = Isomer(element, A, state=state)
    return (iso.mass, iso.mass_excess_estimated) if with_flag else iso.mass


def get_excitation_energy(element: str, A: int, state="m", with_flag: bool = False):
    """Convenience function: excitation energy (keV) above the ground state."""
    iso = Isomer(element, A, state=state)
    return (
        (iso.excitation_energy, iso.excitation_energy_estimated)
        if with_flag
        else iso.excitation_energy
    )
