"""
Fixed-width parser for the AME2020 atomic mass table (mass.mas20).

Format reference (from the file header), Fortran format string:
    a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6

    cc NZ  N  Z  A    el  o     mass  unc binding unc      B  beta  unc    atomic_mass   unc

Quirks handled:
    '#' in place of / after a value -> estimated (non-experimental) value.
      The numeric value is still present and parseable; we strip the '#'
      and flag the value as estimated.
    '*' in place of a value -> not calculable -> NaN.
    The "ATOMIC MASS" column is split into an integer part (i3) and a
      fractional part in micro-u (f13.6); these must be concatenated and
      divided by 1e6 to get a normal mass-excess-like fractional value in u.
    Beta-decay columns can be blank (stable / no beta decay defined).
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import pandas as pd

# 0-indexed [start, end) column slices, derived from the Fortran format spec.
_COLSPECS = {
    "N_minus_Z": (1, 4),
    "N": (4, 9),
    "Z": (9, 14),
    "A": (14, 19),
    "element": (20, 23),
    "origin": (23, 27),
    "mass_excess": (28, 42),
    "mass_excess_unc": (42, 54),
    "binding_per_A": (54, 67),
    "binding_per_A_unc": (68, 78),
    "beta_type": (79, 81),
    "beta_energy": (81, 94),
    "beta_energy_unc": (94, 105),
    "atomic_mass_int": (106, 109),
    "atomic_mass_frac": (110, 123),
    "atomic_mass_unc": (123, 135),
}

# Minimum line length we expect for a genuine data row (avoids picking up
# short header/footer lines).
_MIN_DATA_LINE_LEN = 100


def _slice(line: str, key: str) -> str:
    start, end = _COLSPECS[key]
    # Lines can be shorter than the nominal width if trailing whitespace was
    # stripped by an editor; pad defensively.
    if len(line) < end:
        line = line.ljust(end)
    return line[start:end]


def _is_data_line(line: str) -> bool:
    """Heuristic matching the real layout: A field is numeric, element field
    is alphabetic (allowing the 1-char neutron 'n' as well as 2-3 char
    element symbols)."""
    if len(line) < _MIN_DATA_LINE_LEN:
        return False
    a_field = _slice(line, "A").strip()
    el_field = _slice(line, "element").strip()
    return a_field.isdigit() and el_field.isalpha()


def _parse_value_unc(raw_value: str, raw_unc: str) -> tuple[float, float, bool]:
    """Parse a (value, uncertainty) pair, handling '*' (not calculable) and
    '#' (estimated) markers.

    Returns (value, uncertainty, is_estimated). value/uncertainty are NaN
    if not calculable.
    """
    raw_value = raw_value.strip()
    raw_unc = raw_unc.strip()

    if raw_value == "*" or raw_value == "":
        return float("nan"), float("nan"), False

    is_estimated = "#" in raw_value or "#" in raw_unc
    value = float(raw_value.replace("#", ""))

    if raw_unc == "" or raw_unc == "*":
        unc = float("nan")
    else:
        unc = float(raw_unc.replace("#", ""))

    return value, unc, is_estimated


def _parse_atomic_mass(int_part: str, frac_part: str, unc_part: str) -> tuple[float, float, bool]:
    """The atomic mass is split across an integer-u column (i3) and a
    micro-u fractional column (f13.6). Concatenate as strings (matching the
    Fortran fixed-point convention) then convert to atomic mass units.
    """
    int_part = int_part.strip()
    frac_part = frac_part.strip()
    unc_part = unc_part.strip()

    if frac_part == "" or frac_part == "*":
        return float("nan"), float("nan"), False

    is_estimated = "#" in frac_part or "#" in unc_part
    frac_clean = frac_part.replace("#", "")
    unc_clean = unc_part.replace("#", "") if unc_part not in ("", "*") else "nan"

    # int_part + frac_clean (micro-u, 6 decimal-equivalent digits) concatenated
    # as in the original Fortran fixed-width fields, e.g. "1" + "008664.91590"
    # -> mass = 1 + 008664.91590 / 1e6
    micro_u_value = float(frac_clean)
    int_value = int(int_part) if int_part else 0
    atomic_mass_u = int_value + micro_u_value / 1e6

    unc_u = float(unc_clean) / 1e6 if unc_clean != "nan" else float("nan")

    return atomic_mass_u, unc_u, is_estimated


@dataclass
class _ParsedRow:
    N: int
    Z: int
    A: int
    element: str
    origin: str
    mass_excess_keV: float
    mass_excess_unc_keV: float
    mass_excess_estimated: bool
    binding_per_A_keV: float
    binding_per_A_unc_keV: float
    binding_per_A_estimated: bool
    beta_type: str
    beta_energy_keV: float
    beta_energy_unc_keV: float
    beta_energy_estimated: bool
    atomic_mass_u: float
    atomic_mass_unc_u: float
    atomic_mass_estimated: bool


def _parse_line(line: str) -> _ParsedRow:
    N = int(_slice(line, "N").strip())
    Z = int(_slice(line, "Z").strip())
    A = int(_slice(line, "A").strip())
    element = _slice(line, "element").strip()
    origin = _slice(line, "origin").strip()

    me, me_unc, me_est = _parse_value_unc(
        _slice(line, "mass_excess"), _slice(line, "mass_excess_unc")
    )
    be, be_unc, be_est = _parse_value_unc(
        _slice(line, "binding_per_A"), _slice(line, "binding_per_A_unc")
    )
    beta_type = _slice(line, "beta_type").strip()
    beta_e, beta_e_unc, beta_e_est = _parse_value_unc(
        _slice(line, "beta_energy"), _slice(line, "beta_energy_unc")
    )
    am, am_unc, am_est = _parse_atomic_mass(
        _slice(line, "atomic_mass_int"),
        _slice(line, "atomic_mass_frac"),
        _slice(line, "atomic_mass_unc"),
    )

    return _ParsedRow(
        N=N,
        Z=Z,
        A=A,
        element=element,
        origin=origin,
        mass_excess_keV=me,
        mass_excess_unc_keV=me_unc,
        mass_excess_estimated=me_est,
        binding_per_A_keV=be,
        binding_per_A_unc_keV=be_unc,
        binding_per_A_estimated=be_est,
        beta_type=beta_type,
        beta_energy_keV=beta_e,
        beta_energy_unc_keV=beta_e_unc,
        beta_energy_estimated=beta_e_est,
        atomic_mass_u=am,
        atomic_mass_unc_u=am_unc,
        atomic_mass_estimated=am_est,
    )


def parse_mass_file(path: str | Path) -> pd.DataFrame:
    """Parse an AME2020-format mass table file into a DataFrame.

    Parameters
    ----------
    path : str or Path
        Path to a mass.mas20-format file.

    Returns
    -------
    pandas.DataFrame
        One row per isotope, indexed by nothing in particular (use
        `set_index(['element', 'A'])` or rely on the lookup helpers in
        `ame2020.core` instead of querying this DataFrame directly).
    """
    path = Path(path)
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if _is_data_line(line):
                rows.append(_parse_line(line))

    if not rows:
        raise ValueError(
            f"No data rows parsed from {path}. Is this a valid AME2020 "
            "mass.mas20-format file?"
        )

    df = pd.DataFrame([vars(r) for r in rows])
    # Normalize element symbol casing (AME uses e.g. 'He', 'n' for neutron).
    df["element"] = df["element"].astype(str)
    return df


def default_data_path() -> Path:
    """Path to the mass.mas20 file bundled with this package."""
    return resources.files("ame2020.data").joinpath("mass.mas20")


_CACHED_DF: pd.DataFrame | None = None


def load_default_table(force_reload: bool = False) -> pd.DataFrame:
    """Load (and cache) the bundled AME2020 mass table."""
    global _CACHED_DF
    if _CACHED_DF is None or force_reload:
        _CACHED_DF = parse_mass_file(default_data_path())
    return _CACHED_DF
