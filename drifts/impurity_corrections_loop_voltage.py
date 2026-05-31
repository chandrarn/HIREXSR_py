"""Loop-voltage impurity/main-ion toroidal correction models."""

from __future__ import annotations

import numpy as np

E_CHARGE = 1.602176634e-19  # [C]
TWO_PI = 2.0 * np.pi
AMU_TO_KG = 1.66053906660e-27  # [kg]


def load_vres_from_zeff_neo(
    shot: int,
    dt: float = 0.1,
    trange: tuple[float, float] = (0.5, 1.5),
    verbose: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Load v_res(t) from zeff_neo, where v_res = P_oh / I_p.

    Returns
    -------
    t_vres_s : ndarray
        Time grid used by zeff_neo for v_res.
    vres_v : ndarray
        Effective resistive voltage [V].
    """
    from ..zeff_neo_python import zeff_neo

    out = zeff_neo(
        shot=shot,
        dt=dt,
        trange=trange,
        plot=False,
        verbose=verbose,
        strict_diagnostics=False,
        return_vres=True,
    )
    if len(out) < 4:
        raise RuntimeError("zeff_neo did not return vres outputs as expected.")
    _zeff, _times, vres_v, t_vres_s = out[:4]
    return np.asarray(t_vres_s, dtype=float), np.asarray(vres_v, dtype=float)


def _velocity_to_toroidal_frequency_hz(
    v_tor_m_s: np.ndarray, r_major_m: np.ndarray
) -> np.ndarray:
    """Convert toroidal velocity [m/s] to toroidal rotation frequency [Hz]."""
    r_safe = np.maximum(np.asarray(r_major_m, dtype=float), 1e-6)
    return np.asarray(v_tor_m_s, dtype=float) / (TWO_PI * r_safe)


def _parallel_electric_field_from_loop_voltage_v_m(
    loop_voltage_V: float | np.ndarray, r_major_m: np.ndarray
) -> np.ndarray:
    """Approximate parallel electric field from loop voltage [V/m].

    Uses the same geometric relation as the ephi term in zeff_neo:
      E_parallel ~= E_phi ~= V_loop / (2 pi R).
    """
    r_safe = np.maximum(np.asarray(r_major_m, dtype=float), 1e-6)
    return np.asarray(loop_voltage_V, dtype=float) / (TWO_PI * r_safe)


def _effective_resistive_voltage_v(
    loop_voltage_V: float | None = None,
    ohmic_power_W: float | np.ndarray | None = None,
    plasma_current_A: float | np.ndarray | None = None,
) -> np.ndarray:
    """Return effective resistive voltage, preferring P_oh / I_p when available.

    Mirrors zeff_neo construction where
      vres = poh / ip.
    If both ohmic_power_W and plasma_current_A are provided, vres is computed
    from these inputs. Any non-finite or near-zero-current points fall back to
    loop_voltage_V when provided.
    """
    if ohmic_power_W is not None and plasma_current_A is not None:
        poh = np.asarray(ohmic_power_W, dtype=float)
        ip = np.asarray(plasma_current_A, dtype=float)
        vres = np.divide(
            poh,
            ip,
            out=np.full_like(poh, np.nan, dtype=float),
            where=np.abs(ip) > 0.0,
        )
        if loop_voltage_V is not None:
            fallback = np.asarray(loop_voltage_V, dtype=float)
            vres = np.where(np.isfinite(vres), vres, fallback)
        return vres

    if loop_voltage_V is None:
        raise ValueError(
            "Provide loop_voltage_V, or provide both ohmic_power_W and plasma_current_A."
        )

    return np.asarray(loop_voltage_V, dtype=float)


def _coulomb_logarithm_ii(
    ti_eV: np.ndarray,
    ni_m3: np.ndarray,
) -> np.ndarray:
    """Return a robust 2023 NRL-style ion-ion Coulomb logarithm estimate."""
    ti = np.maximum(np.asarray(ti_eV, dtype=float), 1e-3)
    ni = np.maximum(np.asarray(ni_m3, dtype=float), 1e10)
    ni_cm3 = np.maximum(ni * 1e-6, 1e4)
    ln_lambda = 23.0 - np.log(np.sqrt(2.0) * np.sqrt(ni_cm3) / (ti**1.5))
    return np.clip(ln_lambda, 2.0, 40.0)


def estimate_main_ion_self_collision_time_s(
    ti_eV: np.ndarray,
    ni_m3: np.ndarray,
    z_main: float = 1.0,
    m_main_amu: float = 2.0,
    coulomb_log: np.ndarray | None = None,
) -> np.ndarray:
    """Estimate main-ion self-collision time tau_ii [s].

    Uses the NRL-style ion-ion collision frequency scaling:
      nu_ii [1/s] = 4.80e-8 * Z^4 * n_i[cm^-3] * lnLambda / (sqrt(A_i) * T_i[eV]^(3/2))
    and returns tau_ii = 1 / nu_ii.
    """
    ti = np.maximum(np.asarray(ti_eV, dtype=float), 1e-3)
    ni_cm3 = np.maximum(np.asarray(ni_m3, dtype=float) * 1e-6, 1e4)
    z_main_safe = max(float(z_main), 1e-8)
    m_main_safe = max(float(m_main_amu), 1e-8)
    ln_lambda = (
        _coulomb_logarithm_ii(ti, ni_m3)
        if coulomb_log is None
        else np.maximum(np.asarray(coulomb_log, dtype=float), 1.0)
    )

    nu_ii = (
        4.80e-8
        * (z_main_safe**4)
        * ni_cm3
        * ln_lambda
        / (np.sqrt(m_main_safe) * np.maximum(ti, 1e-8) ** 1.5)
    )
    return 1.0 / np.maximum(nu_ii, 1e-30)


def _impurity_main_ion_loop_voltage_correction_factor(
    z_imp: float,
    mu_imp_amu: float,
    dilution_factor: float = 1.0,
) -> float:
    """Return legacy model-1 prefactor for loop-voltage correction."""
    mu_safe = max(float(mu_imp_amu), 1e-8)
    return -4.19e3 * float(dilution_factor) * (float(z_imp) / np.sqrt(mu_safe))


def impurity_main_ion_delta_v_loop_voltage(
    ti_eV: np.ndarray,
    ni_m3: np.ndarray,
    r_major_m: np.ndarray,
    loop_voltage_V: float,
    z_imp: float,
    mu_imp_amu: float,
    dilution_factor: float = 1.0,
) -> np.ndarray:
    """Legacy approximation: impurity-main ion delta v from loop voltage.

    The legacy Hutchinson-side ``dilution_factor`` is intended to be order unity.
    If you start from a small impurity fraction f, pass ``1 - f`` here.
    """
    """ Source: I. Hutchinson 2000, Eqn 1"""
    ti_keV = np.maximum(np.asarray(ti_eV, dtype=float) * 1e-3, 0.0)
    n20 = np.maximum(np.asarray(ni_m3, dtype=float) * 1e-20, 1e-8)
    r_major_safe = np.maximum(np.asarray(r_major_m, dtype=float), 1e-6)

    coeff = _impurity_main_ion_loop_voltage_correction_factor(
        z_imp=z_imp,
        mu_imp_amu=mu_imp_amu,
        dilution_factor=dilution_factor,
    )
    return coeff * (float(loop_voltage_V) / r_major_safe) * ((ti_keV**1.5) / n20)


def impurity_main_ion_delta_f_loop_voltage_hz(
    ti_eV: np.ndarray,
    ni_m3: np.ndarray,
    r_major_m: np.ndarray,
    loop_voltage_V: float,
    z_imp: float,
    mu_imp_amu: float,
    dilution_factor: float = 1.0,
) -> np.ndarray:
    """Legacy approximation converted to toroidal rotation frequency [Hz]."""
    delta_v = impurity_main_ion_delta_v_loop_voltage(
        ti_eV=ti_eV,
        ni_m3=ni_m3,
        r_major_m=r_major_m,
        loop_voltage_V=loop_voltage_V,
        z_imp=z_imp,
        mu_imp_amu=mu_imp_amu,
        dilution_factor=dilution_factor,
    )
    return _velocity_to_toroidal_frequency_hz(delta_v, r_major_m)


def impurity_main_ion_delta_v_loop_voltage_first_term(
    ni_m3: np.ndarray,
    r_major_m: np.ndarray,
    loop_voltage_V: float | None = None,
    ohmic_power_W: float | np.ndarray | None = None,
    plasma_current_A: float | np.ndarray | None = None,
    z_imp: float = 18.0,
    impurity_dilution_factor: float = 1e-4,
    alpha: np.ndarray | float | None = None,
    tau_ii_s: np.ndarray | None = None,
    ti_eV_for_tau: np.ndarray | None = None,
    m_main_amu: float = 2.0,
    m_imp_amu: float = 39.948,
    z_main: float = 1.0,
) -> np.ndarray:
    """First-term-only Y. Kim 1991 Eq.56-Eq.57 closed-form estimate of v_II,I - v_II,i.

        Implements the user-specified simplified expression:
      v_II,I - v_II,i = -(tau_ii Z_i e E_|| / m_i)
        * ((Z_I - Z_i) / Z_I)
        * ((sqrt(2) + 13 alpha / 4) / ((1 + alpha) * (sqrt(2) + alpha)))
        * ((n_i m_i - n_I m_I) / (n_i m_i + n_I m_I))

        with n_I = impurity_fraction * n_i by default, and
    alpha = n_I Z_I^2 / (n_i Z_i^2) unless alpha is explicitly provided.

        For comparison against the Hutchinson-style legacy form, the corresponding
        order-unity factor is ``1 - impurity_fraction``.

        The electric-field drive uses vres from zeff_neo-style inputs when present:
            vres = P_oh / I_p,
        then E_|| ~= vres / (2 pi R).
    """

    ni = np.maximum(np.asarray(ni_m3, dtype=float), 1e-30)
    r_major = np.asarray(r_major_m, dtype=float)
    vres_v = _effective_resistive_voltage_v(
        loop_voltage_V=loop_voltage_V,
        ohmic_power_W=ohmic_power_W,
        plasma_current_A=plasma_current_A,
    )
    e_par = _parallel_electric_field_from_loop_voltage_v_m(vres_v, r_major)
    z_main_safe = max(float(z_main), 1e-8)
    z_imp_safe = max(float(z_imp), 1e-8)
    dilution = max(float(impurity_dilution_factor), 0.0)
    n_imp = dilution * ni

    if tau_ii_s is None:
        if ti_eV_for_tau is None:
            raise ValueError(
                "Provide tau_ii_s directly, or provide ti_eV_for_tau to estimate tau_ii."
            )
        tau_ii = estimate_main_ion_self_collision_time_s(
            ti_eV=np.asarray(ti_eV_for_tau, dtype=float),
            ni_m3=ni,
            z_main=z_main,
            m_main_amu=m_main_amu,
        )
    else:
        tau_ii = np.maximum(np.asarray(tau_ii_s, dtype=float), 0.0)

    m_main_kg = max(float(m_main_amu), 1e-8) * AMU_TO_KG
    m_imp_kg = max(float(m_imp_amu), 1e-8) * AMU_TO_KG
    if alpha is None:
        alpha_use = (n_imp * z_imp_safe * z_imp_safe) / (ni * z_main_safe * z_main_safe)
    else:
        alpha_use = np.maximum(np.asarray(alpha, dtype=float), 1e-12)

    sqrt2 = np.sqrt(2.0)
    charge_factor = (z_imp_safe - z_main_safe) / z_imp_safe
    collisional_factor = (sqrt2 + 13.0 * alpha_use / 4.0) / (
        (1.0 + alpha_use) * (sqrt2 + alpha_use)
    )
    density_moment_factor = (ni * m_main_kg - n_imp * m_imp_kg) / np.maximum(
        ni * m_main_kg + n_imp * m_imp_kg, 1e-30
    )

    return (
        -tau_ii
        * z_main_safe
        * E_CHARGE
        * e_par
        / m_main_kg
        * charge_factor
        * collisional_factor
        * density_moment_factor
    )


def impurity_main_ion_delta_f_loop_voltage_first_term_hz(
    ni_m3: np.ndarray,
    r_major_m: np.ndarray,
    loop_voltage_V: float | None = None,
    ohmic_power_W: float | np.ndarray | None = None,
    plasma_current_A: float | np.ndarray | None = None,
    z_imp: float = 18.0,
    impurity_dilution_factor: float = 1e-4,
    alpha: np.ndarray | float | None = None,
    tau_ii_s: np.ndarray | None = None,
    ti_eV_for_tau: np.ndarray | None = None,
    m_main_amu: float = 2.0,
    m_imp_amu: float = 39.948,
    z_main: float = 1.0,
) -> np.ndarray:
    """First-term-only Eq.56-Eq.57 correction as toroidal frequency [Hz]."""
    delta_v = impurity_main_ion_delta_v_loop_voltage_first_term(
        ni_m3=ni_m3,
        r_major_m=r_major_m,
        loop_voltage_V=loop_voltage_V,
        ohmic_power_W=ohmic_power_W,
        plasma_current_A=plasma_current_A,
        z_imp=z_imp,
        impurity_dilution_factor=impurity_dilution_factor,
        alpha=alpha,
        tau_ii_s=tau_ii_s,
        ti_eV_for_tau=ti_eV_for_tau,
        m_main_amu=m_main_amu,
        m_imp_amu=m_imp_amu,
        z_main=z_main,
    )
    return _velocity_to_toroidal_frequency_hz(delta_v, r_major_m)
