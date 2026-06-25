# ame2020

Lookup atomic masses, mass excesses, and binding energies from the
[AME2020 atomic mass evaluation](https://www-nds.iaea.org/amdc/) (`mass.mas20`),
without re-parsing the fixed-width file yourself every time.

## ⚠️ Data file placeholders

Both bundled data files are currently **small samples used to validate the
parsers**, not the full evaluations:

- `src/ame2020/data/mass.mas20` (AME2020) — only n, H, He, Li, Be up to A=6.
- `src/ame2020/data/nubase.mas20` (NUBASE2020) — only A=1–14.

**Replace both with the full files from the AMDC**
(https://www-nds.iaea.org/amdc/) before relying on this for Pd/Rh data.
Both parsers already handle the full files' format quirks (estimated
values marked `#`, non-calculable fields marked `*`, isomer/level rows,
negative N-Z, etc.) — only the bundled data needs swapping in.

To replace it:
```bash
cp /path/to/your/mass.mas20 src/ame2020/data/mass.mas20
```
Then reinstall (`pip install -e .`) or just restart your Python session —
the table is cached in memory per-process via `load_default_table()`.

## Install

```bash
pip install .
```
(`-e` for an editable/development install; drop it for a normal install.)
A normal install (pip install .) copies your package into site-packages at install time. 

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

## Isomer masses (NUBASE2020)

`mass.mas20` (AME2020) only has ground states. Isomer/excited-state data
comes from a separate file, `nubase.mas20` (NUBASE2020), which this package
also parses — bundled at `src/ame2020/data/nubase.mas20`. **Same placeholder
caveat applies**: it currently only covers A=1–14. Replace it with the full
file from the AMDC before relying on it for Pd/Rh isomers.

```python
from ame2020 import Isomer, list_states

# see what states exist for a nuclide
list_states("Li", 10)
#    level_index state_label  mass_excess_keV  exc_energy_keV half_life ...
# 0            0                      33053.0             NaN        2.0 zs
# 1            1           m          33250.0           200.0        3.7 zs
# 2            2           n          33530.0           480.0       1.35 zs

gs = Isomer("Li", 10)              # ground state (default)
m  = Isomer("Li", 10, state="m")   # first/lowest isomer
n  = Isomer("Li", 10, state="n")   # second isomer

m.mass                 # absolute mass, u (derived from m.mass_excess)
m.mass_excess          # keV, the ISOMER's own mass excess (not a delta)
m.excitation_energy    # keV, above the ground state
m.half_life, m.half_life_unit
m.Jpi                  # spin-parity string, e.g. "1+"
m.is_isomer            # True for state in {m, n}; False for gs/levels/IAS
```

Functional wrappers:
```python
from ame2020 import get_isomer_mass, get_excitation_energy

get_isomer_mass("Li", 10, state="m")            # u
get_excitation_energy("Li", 10, state="m")      # keV above ground state
```

NUBASE also lists higher-lying "levels" and isospin-analogue states
(`state_label` of `p`, `q`, `i`, `j`, ...) for some nuclides — these are
real entries in the file but generally aren't what's meant by "isomer" in
the usual sense; `Isomer.is_isomer` is `False` for these (only `m`/`n`
states are flagged `True`). You can still retrieve them explicitly via
`state="i"` etc. if needed.

Requesting a state that doesn't exist for a given nuclide raises
`IsomerNotFoundError` (also a `KeyError` subclass) and lists which states
*are* available.

## Notes

- Element symbols are case-insensitive (`"pd"`, `"Pd"`, `"PD"` all work);
  the free neutron is `"n"`.
- Looking up an isotope not present in the table raises
  `ame2020.IsotopeNotFoundError` (a `KeyError` subclass).
- The underlying parsed table is cached per-process; call
  `load_default_table(force_reload=True)` if you've swapped the data file
  mid-session.
