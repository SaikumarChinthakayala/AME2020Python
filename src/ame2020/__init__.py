"""
ame2020: lookup tools for the AME2020 atomic mass evaluation.

    >>> from ame2020 import Isotope
    >>> iso = Isotope("Pd", 102)
    >>> iso.mass
    101.905609...

    >>> from ame2020 import get_mass
    >>> get_mass("Rh", 103)
    102.905504...
"""

from .core import (
    Isotope,
    IsotopeNotFoundError,
    get_beta_decay_energy,
    get_binding_energy,
    get_binding_per_A,
    get_mass,
    get_mass_excess,
)
from .isomers import (
    Isomer,
    IsomerNotFoundError,
    get_excitation_energy,
    get_isomer_mass,
    list_states,
)
from .parser import load_default_table, parse_mass_file
from .nubase_parser import load_default_nubase_table, parse_nubase_file

__all__ = [
    "Isotope",
    "IsotopeNotFoundError",
    "get_mass",
    "get_mass_excess",
    "get_binding_per_A",
    "get_binding_energy",
    "get_beta_decay_energy",
    "load_default_table",
    "parse_mass_file",
    "Isomer",
    "IsomerNotFoundError",
    "get_isomer_mass",
    "get_excitation_energy",
    "list_states",
    "load_default_nubase_table",
    "parse_nubase_file",
]

__version__ = "0.1.0"
