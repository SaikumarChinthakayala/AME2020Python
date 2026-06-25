"""
Fixed-width parser for the NUBASE2020 table (nubase.mas20).

Column spec (1-indexed inclusive, as published in the file header), converted
to 0-indexed Python slices below and verified against real sample lines:

      1: 3   AAA           a3       Mass Number (AAA)
      5: 8   ZZZi          a4       Atomic Number (ZZZ); i=0 (gs); i=1,2 (isomers);
                                     i=3,4 (levels); i=5 (resonance); i=8,9 (IAS)
     12:16   A El          a5       A Element
     17:17   s             a1       s=m,n (isomers); s=p,q (levels); s=r (resonance);
                                     s=i,j (IAS)
     19:31   Mass #     f13.6       Mass Excess in keV (# from systematics)
     32:42   dMass #    f11.6       Mass Excess uncertainty in keV
     43:54   Exc #      f12.6       Isomer/level Excitation Energy in keV, RELATIVE
                                     TO THE GROUND STATE (i.e. NOT an absolute mass).
     55:65   dE #       f11.6       Excitation Energy uncertainty in keV
     66:67   Orig          a2       Origin of Excitation Energy
     68:68   Isom.Unc      a1       '*' = gs/isomer ordering uncertain
     69:69   Isom.Inv      a1       '&' = gs/isomer order reversed vs ENSDF
     70:78   T #         f9.4       Half-life; 'stbl'=stable; 'p-unst'=particle unstable
     79:80   unit T        a2       Half-life unit
     82:88   dT            a7       Half-life uncertainty
     89:102  Jpi */#/T=    a14      Spin and parity
    103:104  Ensdf year    a2       ENSDF update year
    115:118  Discovery     a4       Year of discovery
    120:209  BR            a90      Decay modes / intensities / isotopic abundance

IMPORTANT: the "Mass Excess" column on an isomer row is its OWN mass excess
(ground-state mass excess + excitation energy), not a delta. The "Exc" column
is the delta relative to the ground state. In practice the Mass-Excess column
on isomer rows already equals (gs mass excess + Exc) to within rounding, so we
trust it directly rather than re-deriving it, but we keep both available.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Optional

import pandas as pd

_COLSPECS = {
    "AAA": (0, 3),
    "ZZZi": (4, 8),
    "A_El": (11, 16),
    "s": (16, 17),
    "mass_excess": (18, 31),
    "dmass": (31, 42),
    "exc_energy": (42, 54),
    "dexc": (54, 65),
    "orig": (65, 67),
    "isom_unc": (67, 68),
    "isom_inv": (68, 69),
    "half_life": (69, 78),
    "unit_T": (78, 80),
    "dT": (81, 88),
    "Jpi": (88, 102),
    "ensdf_yr": (102, 104),
    "discovery": (114, 118),
    "BR": (119, 209),
}

_MIN_DATA_LINE_LEN = 40  # short lines (e.g. no half-life data) still exceed this


def _slice(line: str, key: str) -> str:
    start, end = _COLSPECS[key]
    if len(line) < end:
        line = line.ljust(end)
    return line[start:end]


def _is_data_line(line: str) -> bool:
    if len(line) < _MIN_DATA_LINE_LEN or line.startswith("#"):
        return False
    aaa = _slice(line, "AAA").strip()
    zzzi = _slice(line, "ZZZi").strip()
    return aaa.isdigit() and zzzi.isdigit()


def _parse_value_unc(raw_value: str, raw_unc: str) -> tuple:
    """Same convention as the AME parser: '#' = estimated, blank = NaN.
    NUBASE has no '*' marker for not-calculable in these numeric fields
    (blank is used instead), but we handle '*' defensively anyway.
    """
    raw_value = raw_value.strip()
    raw_unc = raw_unc.strip()

    if raw_value in ("", "*"):
        return float("nan"), float("nan"), False

    is_estimated = "#" in raw_value or "#" in raw_unc
    value = float(raw_value.replace("#", ""))

    if raw_unc in ("", "*"):
        unc = float("nan")
    else:
        unc = float(raw_unc.replace("#", ""))

    return value, unc, is_estimated


@dataclass
class _ParsedNubaseRow:
    A: int
    Z: int
    level_index: int          # the trailing digit of ZZZi: 0=gs, 1/2=isomers, 3/4=levels, 5=resonance, 8/9=IAS
    element: str
    state_label: str          # the 's' character: '', 'm', 'n', 'p', 'q', 'r', 'i', 'j', 'x'...
    mass_excess_keV: float
    mass_excess_unc_keV: float
    mass_excess_estimated: bool
    exc_energy_keV: float          # excitation above ground state; NaN for the gs row itself
    exc_energy_unc_keV: float
    exc_energy_estimated: bool
    half_life: str             # raw string, e.g. '609.8', 'stbl', 'p-unst'
    half_life_unit: str
    Jpi: str
    is_isomer: bool            # True only for level_index in {1, 2} (the 'm'/'n' true isomers)


def _parse_line(line: str) -> _ParsedNubaseRow:
    aaa = _slice(line, "AAA").strip()
    zzzi = _slice(line, "ZZZi").strip().zfill(4)

    A = int(aaa)
    Z = int(zzzi[:3])
    level_index = int(zzzi[3])

    a_el = _slice(line, "A_El").strip()
    # A_El looks like "102Pd" -- strip the leading mass number to get the element symbol
    element = a_el[len(aaa):].strip() if a_el.startswith(aaa) else "".join(
        ch for ch in a_el if ch.isalpha()
    )

    state_label = _slice(line, "s").strip()

    me, me_unc, me_est = _parse_value_unc(
        _slice(line, "mass_excess"), _slice(line, "dmass")
    )
    exc, exc_unc, exc_est = _parse_value_unc(
        _slice(line, "exc_energy"), _slice(line, "dexc")
    )

    half_life = _slice(line, "half_life").strip()
    half_life_unit = _slice(line, "unit_T").strip()
    jpi = _slice(line, "Jpi").strip()

    is_isomer = level_index in (1, 2)

    return _ParsedNubaseRow(
        A=A,
        Z=Z,
        level_index=level_index,
        element=element,
        state_label=state_label,
        mass_excess_keV=me,
        mass_excess_unc_keV=me_unc,
        mass_excess_estimated=me_est,
        exc_energy_keV=exc,
        exc_energy_unc_keV=exc_unc,
        exc_energy_estimated=exc_est,
        half_life=half_life,
        half_life_unit=half_life_unit,
        Jpi=jpi,
        is_isomer=is_isomer,
    )


def parse_nubase_file(path) -> pd.DataFrame:
    """Parse a NUBASE2020-format file into a DataFrame, one row per nuclide
    state (ground state, isomers, levels, IAS — everything in the file).

    Use the `level_index` / `is_isomer` columns to filter to what you need;
    see `ame2020.isomers` for a higher-level lookup API.
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
            f"No data rows parsed from {path}. Is this a valid NUBASE2020 "
            "(nubase.mas20) format file?"
        )

    df = pd.DataFrame([vars(r) for r in rows])
    df["element"] = df["element"].astype(str)
    return df


def default_data_path() -> Path:
    """Path to the bundled nubase.mas20 file, if present."""
    return resources.files("ame2020.data").joinpath("nubase.mas20")


_CACHED_NUBASE_DF: Optional[pd.DataFrame] = None


def load_default_nubase_table(force_reload: bool = False) -> pd.DataFrame:
    """Load (and cache) the bundled NUBASE2020 table."""
    global _CACHED_NUBASE_DF
    if _CACHED_NUBASE_DF is None or force_reload:
        _CACHED_NUBASE_DF = parse_nubase_file(default_data_path())
    return _CACHED_NUBASE_DF
