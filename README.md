# ame2020
!!!!!!!!! Package made using the help of Claude !!!!!!!!!

Lookup atomic masses, mass excesses, and binding energies from the
[AME2020 atomic mass evaluation](https://www-nds.iaea.org/amdc/) (`mass.mas20`),
without re-parsing the fixed-width file yourself every time.

## ⚠️ Data file placeholder

`src/ame2020/data/mass.mas20` currently contains **only the 16 lightest
isotopes** (n, H, He, Li, Be up to A=6) as a placeholder used to validate the
parser. **Replace this file with the full AME2020 `mass.mas20`** (available
from the AMDC: https://www-nds.iaea.org/amdc/) before relying on this for
Pd/Rh data. The parser itself handles the full file's format already
(estimated values, non-calculable fields, negative N-Z, etc.) — only the
bundled data needs swapping in.

To replace it:
```bash
cp /path/to/your/mass.mas20 src/ame2020/data/mass.mas20
```
Then reinstall (`pip install -e .`) or just restart your Python session —
the table is cached in memory per-process via `load_default_table()`.

## Install

```bash
pip install -e .
```
(`-e` for an editable/development install; drop it for a normal install.)

## Usage

### Object-oriented

```python
from ame2020 import Isotope

iso = Isotope("Pd", 102)
iso.mass                  # atomic mass, u
iso.mass_excess           # keV
iso.mass_excess_unc       # keV
iso.mass_excess_estimated # bool: True if AME flagged this as estimated (the '#' marker)
iso.binding_per_A         # keV
iso.binding_energy()      # keV, = binding_per_A * A
iso.beta_type             # 'B-', 'B+', or '' 
iso.beta_energy           # keV
iso.symbol                # '102Pd'
iso.is_estimated          # True if ANY of mass/mass_excess/binding_per_A is estimated
```

### Functional

```python
from ame2020 import get_mass, get_mass_excess, get_binding_per_A, get_binding_energy

mass = get_mass("Rh", 103)                          # u
mass, is_est = get_mass("Rh", 103, with_flag=True)  # (u, bool)

me = get_mass_excess("Pd", 102)                     # keV
be_per_A = get_binding_per_A("Pd", 102)              # keV
be_total = get_binding_energy("Pd", 102)             # keV
```

### Estimated values

AME2020 marks some values as estimated (non-experimental) with a `#`. These
are parsed normally (the numeric value is kept) but flagged:

```python
iso = Isotope("Li", 5)
iso.mass_excess              # 11678.887 (keV) — still usable
iso.mass_excess_estimated    # True
```

### Not-calculable values

Where AME marks a field as not calculable (`*`), the corresponding value is
`NaN` (e.g. some beta-decay energies for stable isotopes).

### Custom data file

If you don't want to overwrite the bundled file, point the parser at your
own copy directly:

```python
from ame2020.parser import parse_mass_file
from ame2020.core import Isotope

table = parse_mass_file("/path/to/mass.mas20")
iso = Isotope("Pd", 102, table=table)
```

## Notes

- Element symbols are case-insensitive (`"pd"`, `"Pd"`, `"PD"` all work);
  the free neutron is `"n"`.
- Looking up an isotope not present in the table raises
  `ame2020.IsotopeNotFoundError` (a `KeyError` subclass).
- The underlying parsed table is cached per-process; call
  `load_default_table(force_reload=True)` if you've swapped the data file
  mid-session.
