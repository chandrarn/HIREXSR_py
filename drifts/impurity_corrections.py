"""Compatibility re-exports for impurity correction models.

Loop-voltage correction models live in impurity_corrections_loop_voltage.py.
Ti-gradient correction models live in impurity_corrections_ti_gradient.py.
"""

from .impurity_corrections_loop_voltage import (
    estimate_main_ion_self_collision_time_s,
    impurity_main_ion_delta_f_loop_voltage_first_term_hz,
    impurity_main_ion_delta_f_loop_voltage_hz,
    impurity_main_ion_delta_v_loop_voltage,
    impurity_main_ion_delta_v_loop_voltage_first_term,
)
from .impurity_corrections_ti_gradient import (
    impurity_main_ion_delta_f_ti_gradient_hz,
    impurity_main_ion_delta_v_ti_gradient,
    impurity_main_ion_regime_parameter,
    ion_poloidal_larmor_radius_from_amu_m,
    ion_poloidal_larmor_radius_m,
)

__all__ = [
    "estimate_main_ion_self_collision_time_s",
    "impurity_main_ion_delta_f_loop_voltage_first_term_hz",
    "impurity_main_ion_delta_f_loop_voltage_hz",
    "impurity_main_ion_delta_f_ti_gradient_hz",
    "impurity_main_ion_delta_v_loop_voltage",
    "impurity_main_ion_delta_v_loop_voltage_first_term",
    "impurity_main_ion_delta_v_ti_gradient",
    "impurity_main_ion_regime_parameter",
    "ion_poloidal_larmor_radius_from_amu_m",
    "ion_poloidal_larmor_radius_m",
]
