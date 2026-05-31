"""Ti-gradient impurity/main-ion toroidal correction models."""

from __future__ import annotations

import warnings

import numpy as np

E_CHARGE = 1.602176634e-19  # [C]
TWO_PI = 2.0 * np.pi
AMU_TO_KG = 1.66053906660e-27  # [kg]


def _velocity_to_toroidal_frequency_hz(
    v_tor_m_s: np.ndarray, r_major_m: np.ndarray
) -> np.ndarray:
    """Convert toroidal velocity [m/s] to toroidal rotation frequency [Hz]."""
    r_safe = np.maximum(np.asarray(r_major_m, dtype=float), 1e-6)
    return np.asarray(v_tor_m_s, dtype=float) / (TWO_PI * r_safe)


def ion_poloidal_larmor_radius_m(
    ti_j: np.ndarray,
    mi_kg: float,
    b_pol_T: np.ndarray,
) -> np.ndarray:
    """Compute the ion poloidal Larmor radius [m]."""
    ti_joules = np.maximum(np.asarray(ti_j, dtype=float), 0.0)
    mi_safe = max(float(mi_kg), 1e-30)
    b_pol = np.maximum(np.abs(np.asarray(b_pol_T, dtype=float)), 1e-12)
    return np.sqrt(2.0 * mi_safe * ti_joules) / (E_CHARGE * b_pol)


def ion_poloidal_larmor_radius_from_amu_m(
    ti_j: np.ndarray,
    mi_amu: float,
    b_pol_T: np.ndarray,
) -> np.ndarray:
    """Convenience wrapper for [m] using ion mass in atomic mass units."""
    return ion_poloidal_larmor_radius_m(
        ti_j=ti_j,
        mi_kg=float(mi_amu) * AMU_TO_KG,
        b_pol_T=b_pol_T,
    )


def impurity_main_ion_regime_parameter(
    n_imp_m3: np.ndarray,
    z_imp: float,
    m_imp_amu: float,
    t_main_eV: np.ndarray,
    n_main_m3: np.ndarray,
    z_main: float,
    m_main_amu: float,
    t_imp_eV: np.ndarray,
) -> np.ndarray:
    """Compute the Testa et al. ordering parameter for Eq. 1 validity."""
    n_imp = np.maximum(np.asarray(n_imp_m3, dtype=float), 1e-30)
    n_main = np.maximum(np.asarray(n_main_m3, dtype=float), 1e-30)
    t_main = np.maximum(np.asarray(t_main_eV, dtype=float), 1e-12)
    t_imp = np.maximum(np.asarray(t_imp_eV, dtype=float), 1e-12)

    z_imp_safe = max(float(z_imp), 1e-12)
    z_main_safe = max(float(z_main), 1e-12)
    m_imp_safe = max(float(m_imp_amu), 1e-12)
    m_main_safe = max(float(m_main_amu), 1e-12)

    density_ratio = (n_imp * z_imp_safe * z_imp_safe) / (
        n_main * z_main_safe * z_main_safe
    )
    mass_ratio = np.sqrt(m_imp_safe / m_main_safe)
    temp_ratio = (t_main / t_imp) ** 1.5
    return density_ratio * mass_ratio * temp_ratio


def _impurity_main_ion_ti_gradient_correction_factor(
    inv_lt_main_1_m: np.ndarray,
    inv_lp_main_1_m: np.ndarray,
    inv_lp_imp_1_m: np.ndarray,
    b_phi_T: np.ndarray,
    b_sq_flux_avg_T2: np.ndarray,
    k2: float = 1.0,
    z_main: float = 1.0,
    z_imp: float = 16.0,
    t_main_eV: np.ndarray | None = None,
    t_imp_eV: np.ndarray | None = None,
) -> np.ndarray:
    """Return Eq. (1) multiplicative correction factor for model 2."""
    inv_lt = np.asarray(inv_lt_main_1_m, dtype=float)
    inv_lp_main = np.asarray(inv_lp_main_1_m, dtype=float)
    inv_lp_imp = np.asarray(inv_lp_imp_1_m, dtype=float)
    b_phi = np.asarray(b_phi_T, dtype=float)
    b2_avg = np.maximum(np.asarray(b_sq_flux_avg_T2, dtype=float), 1e-12)

    t_main = (
        np.asarray(t_main_eV, dtype=float)
        if t_main_eV is not None
        else np.ones_like(inv_lt, dtype=float)
    )
    t_imp = (
        np.asarray(t_imp_eV, dtype=float)
        if t_imp_eV is not None
        else np.asarray(t_main, dtype=float)
    )

    inv_lt_safe = np.where(np.abs(inv_lt) > 1e-12, inv_lt, np.nan)
    ltd_over_lpd = inv_lp_main / inv_lt_safe
    ltd_over_lpc = inv_lp_imp / inv_lt_safe

    b_ratio = (b_phi * b_phi) / b2_avg
    species_term = 1.0 - (
        (float(z_main) * t_imp * ltd_over_lpc)
        / np.maximum(float(z_imp) * np.asarray(t_main, dtype=float), 1e-12)
    )

    core = float(k2) * (1.0 - 3.0 * b_ratio) - ltd_over_lpd * species_term
    return core * (1.0 - b_ratio)


def _impurity_main_ion_ti_gradient_legacy_factor(
    r_minor_m: np.ndarray,
    r_major_m: np.ndarray,
    b_pol_T: np.ndarray,
    species_factor: float = 2.25,
) -> np.ndarray:
    """Legacy model-2 multiplicative factor for Ti-gradient correction."""
    r_minor = np.maximum(np.asarray(r_minor_m, dtype=float), 1e-6)
    r_major = np.maximum(np.asarray(r_major_m, dtype=float), 1e-6)
    b_pol = np.asarray(b_pol_T, dtype=float)

    geom = np.sqrt(np.clip(r_minor / r_major, 0.0, None))
    denom = E_CHARGE * np.where(np.abs(b_pol) > 1e-8, b_pol, np.nan)
    return float(species_factor) * geom / denom


def _run_eq1_regime_checks(
    ti_eV: np.ndarray,
    ti_j: np.ndarray,
    ni_m3: np.ndarray,
    b_pol_T: np.ndarray,
    inv_lt_main_1_m: np.ndarray,
    inv_lp_main_1_m: np.ndarray,
    z_main: float,
    z_imp: float,
    t_imp_eV: np.ndarray,
    impurity_dilution_for_checks: float,
    m_main_amu_for_checks: float,
    m_imp_amu_for_checks: float,
) -> None:
    """Run Eq. (1) validity checks and emit RuntimeWarnings on violations."""
    ti = np.asarray(ti_eV, dtype=float)
    ni = np.asarray(ni_m3, dtype=float)
    b_pol = np.asarray(b_pol_T, dtype=float)
    inv_lt = np.asarray(inv_lt_main_1_m, dtype=float)
    inv_lp_main = np.asarray(inv_lp_main_1_m, dtype=float)

    interior = np.zeros_like(ti, dtype=bool)
    if ti.shape[0] > 2:
        interior[1:-1, ...] = True
    else:
        interior[...] = True

    n_imp_est = np.maximum(float(impurity_dilution_for_checks), 0.0) * ni
    ordering_param = impurity_main_ion_regime_parameter(
        n_imp_m3=n_imp_est,
        z_imp=z_imp,
        m_imp_amu=m_imp_amu_for_checks,
        t_main_eV=ti,
        n_main_m3=ni,
        z_main=z_main,
        m_main_amu=m_main_amu_for_checks,
        t_imp_eV=t_imp_eV,
    )
    ordering_mask = interior & np.isfinite(ordering_param)
    if np.any(ordering_mask) and np.any(ordering_param[ordering_mask] >= 1.0):
        worst = float(np.nanmax(ordering_param[ordering_mask]))
        warnings.warn(
            "Eq. (1) ordering check failed at some interior points: "
            f"max ordering parameter={worst:.3g} (expected < 1).",
            RuntimeWarning,
            stacklevel=3,
        )

    rho_i_theta = ion_poloidal_larmor_radius_from_amu_m(
        ti_j=ti_j,
        mi_amu=m_main_amu_for_checks,
        b_pol_T=b_pol,
    )
    lt_i = 1.0 / np.maximum(np.abs(inv_lt), 1e-12)
    lp_i = 1.0 / np.maximum(np.abs(inv_lp_main), 1e-12)
    l_min = np.minimum(lt_i, lp_i)
    scale_mask = interior & np.isfinite(rho_i_theta) & np.isfinite(l_min)
    if np.any(scale_mask) and np.any(rho_i_theta[scale_mask] > l_min[scale_mask]):
        worst_ratio = float(np.nanmax(rho_i_theta[scale_mask] / l_min[scale_mask]))
        warnings.warn(
            "Eq. (1) scale-separation check failed at some interior points: "
            f"max rho_i_theta/L={worst_ratio:.3g} (expected <= 1).",
            RuntimeWarning,
            stacklevel=3,
        )


def impurity_main_ion_delta_v_ti_gradient(
    ti_eV: np.ndarray,
    r_minor_m: np.ndarray,
    r_major_m: np.ndarray,
    b_pol_T: np.ndarray,
    ni_m3: np.ndarray | None = None,
    b_phi_T: np.ndarray | None = None,
    b_sq_flux_avg_T2: np.ndarray | None = None,
    z_main: float = 1.0,
    z_imp: float = 16.0,
    t_imp_eV: np.ndarray | None = None,
    k2: float = 1.0,
    use_eq1: bool = True,
    impurity_dilution_for_checks: float = 1e-4,
    m_main_amu_for_checks: float = 2.0,
    m_imp_amu_for_checks: float = 40.0,
) -> np.ndarray:
    """Model 2: impurity-main ion toroidal velocity difference from Ti gradient."""
    ti = np.asarray(ti_eV, dtype=float)
    r_minor = np.maximum(np.asarray(r_minor_m, dtype=float), 1e-6)
    b_pol = np.asarray(b_pol_T, dtype=float)

    inv_lt = np.gradient((np.maximum(ti, 1e-8)), r_minor, edge_order=2)
    ti_j = np.maximum(ti, 1e-8) * E_CHARGE
    base = (
        ti_j
        / (2.0 * max(float(z_main), 1e-8) * E_CHARGE)
        * (inv_lt / np.where(np.abs(b_pol) > 1e-8, b_pol, np.nan))
    )

    if use_eq1 and ni_m3 is not None and b_phi_T is not None:
        ni = np.asarray(ni_m3, dtype=float)
        p_main = np.maximum(ni * np.maximum(ti, 1e-8), 1e-16)
        inv_lp_main = np.gradient((p_main), r_minor, edge_order=2)

        if t_imp_eV is not None:
            t_imp_arr = np.asarray(t_imp_eV, dtype=float)
        else:
            t_imp_arr = ti

        p_imp_proxy = np.maximum(
            np.maximum(ni, 1e-16) * np.maximum(t_imp_arr, 1e-8), 1e-16
        )
        inv_lp_imp = np.gradient((p_imp_proxy), r_minor, edge_order=2)

        if b_sq_flux_avg_T2 is None:
            b_sq_flux_avg_T2 = np.asarray(b_phi_T, dtype=float) ** 2 + b_pol**2
        _run_eq1_regime_checks(
            ti_eV=ti,
            ti_j=ti_j,
            ni_m3=ni,
            b_pol_T=b_pol,
            inv_lt_main_1_m=inv_lt,
            inv_lp_main_1_m=inv_lp_main,
            z_main=z_main,
            z_imp=z_imp,
            t_imp_eV=t_imp_arr,
            impurity_dilution_for_checks=impurity_dilution_for_checks,
            m_main_amu_for_checks=m_main_amu_for_checks,
            m_imp_amu_for_checks=m_imp_amu_for_checks,
        )

        correction = _impurity_main_ion_ti_gradient_correction_factor(
            inv_lt_main_1_m=inv_lt,
            inv_lp_main_1_m=inv_lp_main,
            inv_lp_imp_1_m=inv_lp_imp,
            b_phi_T=np.asarray(b_phi_T, dtype=float),
            b_sq_flux_avg_T2=np.asarray(b_sq_flux_avg_T2, dtype=float),
            k2=k2,
            z_main=z_main,
            z_imp=z_imp,
            t_main_eV=ti,
            t_imp_eV=t_imp_arr,
        )
        return base * correction

    dti_dr_j_m = np.gradient(ti * E_CHARGE, r_minor, edge_order=2)
    return (
        _impurity_main_ion_ti_gradient_legacy_factor(
            r_minor_m=r_minor_m,
            r_major_m=r_major_m,
            b_pol_T=b_pol_T,
            species_factor=2.25,
        )
        * dti_dr_j_m
    )


def impurity_main_ion_delta_f_ti_gradient_hz(
    ti_eV: np.ndarray,
    r_minor_m: np.ndarray,
    r_major_m: np.ndarray,
    b_pol_T: np.ndarray,
    ni_m3: np.ndarray | None = None,
    b_phi_T: np.ndarray | None = None,
    b_sq_flux_avg_T2: np.ndarray | None = None,
    z_main: float = 1.0,
    z_imp: float = 16.0,
    t_imp_eV: np.ndarray | None = None,
    k2: float = 0.6,
    use_eq1: bool = False,
) -> np.ndarray:
    """Model 2 correction converted to toroidal rotation frequency [Hz]."""
    delta_v = impurity_main_ion_delta_v_ti_gradient(
        ti_eV=ti_eV,
        r_minor_m=r_minor_m,
        r_major_m=r_major_m,
        b_pol_T=b_pol_T,
        ni_m3=ni_m3,
        b_phi_T=b_phi_T,
        b_sq_flux_avg_T2=b_sq_flux_avg_T2,
        z_main=z_main,
        z_imp=z_imp,
        t_imp_eV=t_imp_eV,
        k2=k2,
        use_eq1=use_eq1,
    )
    return _velocity_to_toroidal_frequency_hz(delta_v, r_major_m)
