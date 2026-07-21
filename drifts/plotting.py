"""Plotting helpers for diamagnetic drift and rotation comparisons."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .compute import (
    E_CHARGE,
    DiamagneticResult,
    EquilibriumData,
    ProfileData,
)
from .impurity_corrections_loop_voltage import (
    impurity_main_ion_delta_f_loop_voltage_first_term_hz,
)
from .impurity_corrections_ti_gradient import impurity_main_ion_delta_f_ti_gradient_hz


def plot_diamagnetic_vs_q_times(
    result: DiamagneticResult,
    profiles: ProfileData,
    equilibrium: EquilibriumData,
    shot: int,
    selected_times_s: list[float],
    doSave: str = "",
    f_lims_omega: list[float] = [-50, 50],
    f_lims_f: list[float] = [-15, 15],
    f_lims_f_tor: list[float] = [-10, 10],
    line: int = 2,
    tht: int = 0,
    max_err_omega: float = 10.0,
    impurity_corr1_params: dict[str, float] | None = None,
    ti_gradient_core_edge_guard_points: int = 1,
    show_plot: bool = True,
) -> None:
    """Production summary plot for diamagnetic drift and toroidal rotation."""
    if len(selected_times_s) == 0:
        raise ValueError("selected_times_s must contain at least one time.")

    t_res = result.time_s
    t_diag = np.asarray(profiles.time_diag, dtype=float).squeeze()
    t_eq = np.asarray(equilibrium.time_s, dtype=float).squeeze()

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(20, 5),
        num=f"Diamagnetic Drift Summary: shot {shot}",
    )
    ax_p = axes[0]
    ax_f = axes[1]
    ax_f_tor = axes[2]
    ax_omg = axes[3]
    ax_q_omg = ax_omg.twinx()

    cmap = plt.get_cmap("tab10")

    has_omega = (
        profiles.omega_tor_rad_s is not None
        and profiles.omega_tor_err_rad_s is not None
        and profiles.psi_hx is not None
    )
    tor_y_min = np.inf
    tor_y_max = -np.inf

    for i_sel, t_sel in enumerate(selected_times_s):
        it_res = int(np.argmin(np.abs(t_res - t_sel)))
        it_diag = int(np.argmin(np.abs(t_diag - t_sel)))
        color = cmap(i_sel % 10)
        time_label = f"t={t_res[it_res]:.2f} s"

        psi_col = result.psi_n[:, it_res]
        pe_col = result.ne_m3[:, it_res] * result.te_eV[:, it_res] * E_CHARGE
        pi_col = result.ni_m3[:, it_res] * result.ti_eV[:, it_res] * E_CHARGE
        order_psi = np.argsort(psi_col)
        psi_s = psi_col[order_psi]
        ax_p.plot(
            psi_s,
            pe_col[order_psi] / 101325,
            color=color,
            lw=2.0,
            ls="-",
            label=rf"$p_e$ ({time_label})",
        )
        ax_p.plot(
            psi_s,
            pi_col[order_psi] / 101325,
            color=color,
            lw=2.0,
            ls="--",
            label=rf"$p_i$ ({time_label})",
        )

        q_col = result.q[:, it_res]
        order_q = np.argsort(q_col)
        q_s = q_col[order_q]
        valid_q = np.isfinite(q_s)
        fe_khz = result.f_star_e_Hz[:, it_res][order_q] / 1e3
        fi_khz = result.f_star_i_Hz[:, it_res][order_q] / 1e3
        ax_f.plot(
            q_s,
            fe_khz,
            color=color,
            lw=2.0,
            ls="-",
            label=rf"$f_{{*e}}$ ({time_label})",
        )
        ax_f.plot(
            q_s,
            fi_khz,
            color=color,
            lw=2.0,
            ls="--",
            label=rf"$f_{{*i}}$ ({time_label})",
        )

        fe_tor_khz = result.f_star_e_tor_Hz[:, it_res][order_q] / 1e3
        fi_tor_khz = result.f_star_i_tor_Hz[:, it_res][order_q] / 1e3
        mask_fe_tor = valid_q & np.isfinite(fe_tor_khz)
        mask_fi_tor = valid_q & np.isfinite(fi_tor_khz)
        ax_f_tor.plot(
            q_s[mask_fe_tor],
            fe_tor_khz[mask_fe_tor],
            color=color,
            lw=2.0,
            ls="-",
            label=rf"$f_{{*e}}^{{\mathrm{{tor}}}}$ ({time_label})",
        )
        ax_f_tor.plot(
            q_s[mask_fi_tor],
            fi_tor_khz[mask_fi_tor],
            color=color,
            lw=2.0,
            ls="--",
            label=rf"$f_{{*i}}^{{\mathrm{{tor}}}}$ ({time_label})",
        )
        if np.any(mask_fe_tor):
            tor_y_min = min(tor_y_min, float(np.nanmin(fe_tor_khz[mask_fe_tor])))
            tor_y_max = max(tor_y_max, float(np.nanmax(fe_tor_khz[mask_fe_tor])))
        if np.any(mask_fi_tor):
            tor_y_min = min(tor_y_min, float(np.nanmin(fi_tor_khz[mask_fi_tor])))
            tor_y_max = max(tor_y_max, float(np.nanmax(fi_tor_khz[mask_fi_tor])))

        r_minor_col = np.asarray(equilibrium.r_minor_q_m[:, it_res], dtype=float)
        r_major_col = np.asarray(equilibrium.r_major_q_m[:, it_res], dtype=float)
        b_pol_col = np.asarray(equilibrium.b_pol_q_T[:, it_res], dtype=float)
        b_phi_col = np.interp(
            r_major_col,
            np.asarray(equilibrium.r_major_full_m, dtype=float),
            np.asarray(equilibrium.b_t_t[:, it_res], dtype=float),
        )
        b_sq_avg_col = b_phi_col * b_phi_col + b_pol_col * b_pol_col
        ti_col = np.asarray(result.ti_eV[:, it_res], dtype=float)
        delta_f2_khz = (
            impurity_main_ion_delta_f_ti_gradient_hz(
                ti_eV=ti_col,
                r_minor_m=r_minor_col,
                r_major_m=r_major_col,
                b_pol_T=b_pol_col,
                ni_m3=np.asarray(result.ni_m3[:, it_res], dtype=float),
                b_phi_T=b_phi_col,
                b_sq_flux_avg_T2=b_sq_avg_col,
                z_main=1.0,
                z_imp=16.0,
                t_imp_eV=ti_col,
                k2=1.0,
                use_eq1=False,
            )[order_q]
            / 1e3
        )
        mask_df2 = valid_q & np.isfinite(delta_f2_khz)
        if ti_gradient_core_edge_guard_points > 0:
            guard_idx_df2 = np.flatnonzero(mask_df2)
            if guard_idx_df2.size > ti_gradient_core_edge_guard_points:
                mask_df2[guard_idx_df2[:ti_gradient_core_edge_guard_points]] = False
        ax_f_tor.plot(
            q_s[mask_df2],
            delta_f2_khz[mask_df2],
            color=color,
            lw=1.8,
            ls=":",
            alpha=0.9,
            zorder=4,
            label=rf"$\Delta f_{{I-i}}^{{(\nabla T_i)}}$ ({time_label})",
        )
        if np.any(mask_df2):
            tor_y_min = min(tor_y_min, float(np.nanmin(delta_f2_khz[mask_df2])))
            tor_y_max = max(tor_y_max, float(np.nanmax(delta_f2_khz[mask_df2])))

        has_vres_traces = impurity_corr1_params is not None and all(
            k in impurity_corr1_params for k in ["corr_time_s", "vres_v_t"]
        )
        has_loop_voltage = (
            impurity_corr1_params is not None
            and "loop_voltage_V" in impurity_corr1_params
        )
        if (
            impurity_corr1_params is not None
            and "z_imp" in impurity_corr1_params
            and (has_vres_traces or has_loop_voltage)
        ):
            t_now = float(t_res[it_res])
            loop_voltage_fallback = float(
                impurity_corr1_params.get("loop_voltage_V", 1.0)
            )
            core_edge_guard_points = max(
                0,
                int(impurity_corr1_params.get("core_edge_guard_points", 2)),
            )
            impurity_fraction = float(
                impurity_corr1_params.get(
                    "impurity_fraction",
                    impurity_corr1_params.get("dilution_factor", 1e-4),
                )
            )

            if has_vres_traces:
                corr_t = np.asarray(impurity_corr1_params["corr_time_s"], dtype=float)
                vres_t = np.asarray(impurity_corr1_params["vres_v_t"], dtype=float)
                vres_now = float(np.interp(t_now, corr_t, vres_t))
            else:
                vres_now = loop_voltage_fallback

            delta_f1_khz = (
                impurity_main_ion_delta_f_loop_voltage_first_term_hz(
                    ni_m3=np.asarray(result.ni_m3[:, it_res], dtype=float),
                    r_major_m=r_major_col,
                    loop_voltage_V=vres_now,
                    z_imp=float(impurity_corr1_params["z_imp"]),
                    impurity_dilution_factor=impurity_fraction,
                    ti_eV_for_tau=ti_col,
                    m_main_amu=2.0,
                    z_main=1.0,
                )[order_q]
                / 1e3
            )
            mask_df1 = valid_q & np.isfinite(delta_f1_khz)
            if core_edge_guard_points > 0:
                guard_idx = np.flatnonzero(mask_df1)
                if guard_idx.size > core_edge_guard_points:
                    mask_df1[guard_idx[:core_edge_guard_points]] = False
            ax_f_tor.plot(
                q_s[mask_df1],
                delta_f1_khz[mask_df1],
                color=color,
                lw=2.8,
                ls="-.",
                marker="o",
                markevery=max(1, max(1, int(np.count_nonzero(mask_df1))) // 12),
                ms=3.5,
                alpha=1.0,
                zorder=6,
                label=rf"$\Delta f_{{I-i}}^{{(V_l,\,\mathrm{{Eq56-57\ first\ term}})}}$ ({time_label})",
            )
            if np.any(mask_df1):
                tor_y_min = min(tor_y_min, float(np.nanmin(delta_f1_khz[mask_df1])))
                tor_y_max = max(tor_y_max, float(np.nanmax(delta_f1_khz[mask_df1])))

        if has_omega:
            psi_hx_arr = profiles.psi_hx
            omg_arr = profiles.omega_tor_rad_s
            omg_err_arr = profiles.omega_tor_err_rad_s
            if (
                psi_hx_arr is not None
                and omg_arr is not None
                and omg_err_arr is not None
            ):
                psi_hx_col = np.asarray(psi_hx_arr[:-1, it_diag], dtype=float)
                omg_col = np.asarray(omg_arr[:, it_diag], dtype=float)
                omg_err_col = np.asarray(omg_err_arr[:, it_diag], dtype=float)
                order_hx = np.argsort(psi_hx_col)
                psi_hx_s = psi_hx_col[order_hx]
                omg_s = omg_col[order_hx]
                omg_err_s = omg_err_col[order_hx]
                valid_mask = omg_err_s <= max_err_omega
                ax_omg.errorbar(
                    psi_hx_s[valid_mask],
                    omg_s[valid_mask],
                    yerr=omg_err_s[valid_mask],
                    color=color,
                    lw=1.8,
                    elinewidth=0.9,
                    capsize=3,
                    label=rf"$\omega_{{\phi}}$ ({time_label})",
                )

        it_eq = int(np.argmin(np.abs(t_eq - t_sel)))
        q95_val = float(np.asarray(equilibrium.q95, dtype=float).squeeze().flat[it_eq])
        psi_eq_col = np.asarray(equilibrium.psi_n_q[:, it_eq], dtype=float)
        q_eq_col = np.asarray(equilibrium.q_profile[:, it_eq], dtype=float)
        order_eq = np.argsort(psi_eq_col)
        psi_eq_s = psi_eq_col[order_eq]
        q_eq_s = q_eq_col[order_eq]
        mask_q95 = q_eq_s <= q95_val
        ax_q_omg.plot(
            psi_eq_s[mask_q95],
            q_eq_s[mask_q95],
            color=color,
            lw=1.4,
            ls=":",
            label=rf"$q$ ({time_label})",
        )

    if f_lims_f:
        ax_f.set_ylim(f_lims_f)
    if f_lims_f_tor:
        y0, y1 = float(f_lims_f_tor[0]), float(f_lims_f_tor[1])
        if np.isfinite(tor_y_min) and np.isfinite(tor_y_max):
            if tor_y_min < y0 or tor_y_max > y1:
                span = max(tor_y_max - tor_y_min, 1.0)
                pad = 0.05 * span
                y0 = min(y0, tor_y_min - pad)
                y1 = max(y1, tor_y_max + pad)
        ax_f_tor.set_ylim([y0, y1])
    if f_lims_omega:
        ax_omg.set_ylim(f_lims_omega)

    ax_p.set_title(r"Pressure Profiles")
    ax_p.set_xlabel(r"$\psi_N$")
    ax_p.set_ylabel(r"$p$ [atm]")
    ax_p.axhline(0.0, color="k", ls=":", lw=0.8)
    ax_p.grid(alpha=0.3)
    ax_p.legend(fontsize=9, ncol=1)

    ax_f.set_title(r"Diamagnetic Drift Frequencies")
    ax_f.set_xlabel(r"$q(\psi_N,\,t)$")
    ax_f.set_ylabel(r"$f_*$ [kHz]")
    ax_f.axhline(0.0, color="k", ls="--", lw=1)
    ax_f.grid(alpha=0.3)
    ax_f.legend(fontsize=9, ncol=1)

    ax_f_tor.set_title(r"Toroidal Diamagnetic Drift Frequencies ($B_p$)")
    ax_f_tor.set_xlabel(r"$q(\psi_N,\,t)$")
    ax_f_tor.set_ylabel(r"$f_*^{\mathrm{tor}}$ [kHz]")
    ax_f_tor.axhline(0.0, color="k", ls="--", lw=1)
    ax_f_tor.grid(alpha=0.3)
    ax_f_tor.legend(fontsize=9, ncol=1)

    if has_omega:
        ax_omg.set_title(r"Toroidal Rotation (HIREXSR) \& $q(\psi_N)$")
        ax_omg.set_xlabel(r"$\psi_N$")
        ax_omg.set_ylabel(r"$\omega_\phi$ [kHz]")
        ax_q_omg.set_ylabel(r"$q$")
        ax_omg.axhline(0.0, color="k", ls=":", lw=0.8)
        ax_omg.grid(alpha=0.3)
        h1, l1 = ax_omg.get_legend_handles_labels()
        h2, l2 = ax_q_omg.get_legend_handles_labels()
        ax_omg.legend(h1 + h2, l1 + l2, fontsize=9, ncol=1)
    else:
        ax_omg.set_visible(False)

    fig.suptitle(f"Shot {shot}: diamagnetic drift and toroidal rotation")
    fig.tight_layout()

    if doSave:
        save_path = doSave + f"diamagnetic_drift_shot_{shot}_tht_{tht}_line_{line}.pdf"
        fig.savefig(save_path, transparent=True)
        print(f"Saved figure to {save_path}")

    if show_plot:
        plt.show(block=True)
    else:
        plt.close(fig)


def plot_ti_grid_diagnostic(
    profiles: ProfileData,
    equilibrium: EquilibriumData,
    result: DiamagneticResult,
    selected_times_s: list[float],
    shot: int | None = None,
    doSave: str = "",
    show_plot: bool = True,
) -> None:
    """Diagnostic plot comparing the native HIREX-SR Ti grid versus the
    Ti profile after mapping onto the Thomson diagnostic and q-profile grids.

    Three panels per selected time:
      Left  – Ti(R_major): HIREX source, mapped-to-Thomson, and mapped-to-q.
      Centre – dTi/dr(R_major): derivative on each grid.
      Right  – dTi/dr(q): same derivatives plotted against q, which is what
               the impurity correction formula sees.
    """
    from scipy.interpolate import PchipInterpolator
    from .compute import _interp_columns_to_grid

    t_res = np.asarray(result.time_s, dtype=float)
    t_diag = np.asarray(profiles.time_diag, dtype=float).squeeze()
    cmap = plt.get_cmap("tab10")

    fig, axes = plt.subplots(
        len(selected_times_s),
        3,
        figsize=(15, 4 * len(selected_times_s)),
        num=(
            f"Ti Grid Diagnostic: shot {shot}"
            if shot is not None
            else "Ti Grid Diagnostic"
        ),
        squeeze=False,
    )

    for row, t_sel in enumerate(selected_times_s):
        it_res = int(np.argmin(np.abs(t_res - t_sel)))
        it_diag = int(np.argmin(np.abs(t_diag - t_sel)))
        t_actual = float(t_res[it_res])
        color = cmap(row % 10)

        ax_ti = axes[row, 0]
        ax_dti = axes[row, 1]
        ax_dti_q = axes[row, 2]

        # --- Source HIREX grid (native sparse grid) ---
        rho_src = np.asarray(profiles.rho_hx[:, it_diag], dtype=float)
        ti_src_raw = np.asarray(profiles.ti_eV[:, it_diag], dtype=float)
        valid_hx = np.isfinite(rho_src) & np.isfinite(ti_src_raw) & (ti_src_raw > 0)
        rho_s = rho_src[valid_hx]
        ti_s = ti_src_raw[valid_hx]
        ord_s = np.argsort(rho_s)
        rho_s = rho_s[ord_s]
        ti_s = ti_s[ord_s]
        keep = np.concatenate(([True], np.diff(rho_s) > 1e-8))
        rho_s = rho_s[keep]
        ti_s = ti_s[keep]

        dti_s = np.gradient(ti_s, rho_s, edge_order=2)

        # --- Mapped-to-Thomson-diagnostic grid (dense, ~17→ Thomson n-channels) ---
        r_diag_col = np.asarray(profiles.r_major_diag_m[:, it_diag], dtype=float)
        ti_mapped_lin = _interp_columns_to_grid(
            profiles.ti_eV[:, it_diag : it_diag + 1],
            profiles.rho_hx[:, it_diag : it_diag + 1],
            r_diag_col[:, np.newaxis],
            method="linear",
        )[:, 0]
        ti_mapped_pch = _interp_columns_to_grid(
            profiles.ti_eV[:, it_diag : it_diag + 1],
            profiles.rho_hx[:, it_diag : it_diag + 1],
            r_diag_col[:, np.newaxis],
            method="pchip",
        )[:, 0]

        r_diag_sorted_idx = np.argsort(r_diag_col)
        r_diag_s = r_diag_col[r_diag_sorted_idx]
        ti_lin_s = ti_mapped_lin[r_diag_sorted_idx]
        ti_pch_s = ti_mapped_pch[r_diag_sorted_idx]
        dti_lin_diag = np.gradient(ti_lin_s, r_diag_s, edge_order=2)
        dti_pch_diag = np.gradient(ti_pch_s, r_diag_s, edge_order=2)

        # --- Mapped-to-q grid ---
        r_q_col = np.asarray(equilibrium.r_major_q_m[:, it_res], dtype=float)
        ti_q_col = np.asarray(result.ti_eV[:, it_res], dtype=float)
        r_minor_col = np.asarray(equilibrium.r_minor_q_m[:, it_res], dtype=float)
        q_col = np.asarray(result.q[:, it_res], dtype=float)

        ord_rq = np.argsort(r_q_col)
        r_q_s = r_q_col[ord_rq]
        ti_q_s = ti_q_col[ord_rq]
        r_minor_q_s = r_minor_col[ord_rq]
        q_s_rq = q_col[ord_rq]
        dti_q = np.gradient(ti_q_s, r_minor_q_s, edge_order=2)

        # --- Plot: Ti vs R_major ---
        ax_ti.plot(rho_s, ti_s * 1e-3, "o", color=color, ms=6, label="HIREX source")
        if rho_s.size >= 3:
            r_dense = np.linspace(rho_s[0], rho_s[-1], 300)
            ax_ti.plot(
                r_dense,
                PchipInterpolator(rho_s, ti_s)(r_dense) * 1e-3,
                color=color,
                lw=1.4,
                ls="--",
                alpha=0.6,
                label="PCHIP fit (dense)",
            )
        ax_ti.plot(
            r_diag_s,
            ti_lin_s * 1e-3,
            "^",
            color="steelblue",
            ms=4,
            alpha=0.7,
            label="mapped→Thomson (linear)",
        )
        ax_ti.plot(
            r_diag_s,
            ti_pch_s * 1e-3,
            "s",
            color="tomato",
            ms=4,
            alpha=0.7,
            label="mapped→Thomson (pchip)",
        )
        ax_ti.plot(
            r_q_s,
            ti_q_s * 1e-3,
            "-",
            color="gray",
            lw=1.2,
            alpha=0.5,
            label="mapped→q-grid",
        )
        ax_ti.set_title(
            rf"$T_i$ vs $R_\mathrm{{maj}}$ — $t={t_actual:.2f}\,\mathrm{{s}}$"
        )
        ax_ti.set_xlabel(r"$R_\mathrm{major}$ [m]")
        ax_ti.set_ylabel(r"$T_i$ [keV]")
        ax_ti.legend(fontsize=8)
        ax_ti.grid(alpha=0.3)

        # --- Plot: dTi/dr vs R_major ---
        ax_dti.plot(
            rho_s, dti_s * 1e-3, "o", color=color, ms=6, label="HIREX source (FD)"
        )
        ax_dti.plot(
            r_diag_s,
            dti_lin_diag * 1e-3,
            "^",
            color="steelblue",
            ms=4,
            alpha=0.7,
            label=r"$d T_i/dr$ Thomson (linear)",
        )
        ax_dti.plot(
            r_diag_s,
            dti_pch_diag * 1e-3,
            "s",
            color="tomato",
            ms=4,
            alpha=0.7,
            label=r"$d T_i/dr$ Thomson (pchip)",
        )
        ax_dti.plot(
            r_q_s[np.isfinite(q_s_rq)],
            dti_q[np.isfinite(q_s_rq)] * 1e-3,
            "-",
            color="gray",
            lw=1.2,
            alpha=0.5,
            label=r"$d T_i/dr$ q-grid",
        )
        ax_dti.set_title(
            rf"$dT_i/dr$ vs $R_\mathrm{{maj}}$ — $t={t_actual:.2f}\,\mathrm{{s}}$"
        )
        ax_dti.set_xlabel(r"$R_\mathrm{major}$ [m]")
        ax_dti.set_ylabel(r"$dT_i/dr$ [keV/m]")
        ax_dti.legend(fontsize=8)
        ax_dti.grid(alpha=0.3)
        ax_dti.axhline(0.0, color="k", ls=":", lw=0.8)

        # --- Plot: dTi/dr vs q (what correction formula sees) ---
        ord_q = np.argsort(q_col)
        q_sorted = q_col[ord_q]
        ti_q_ord = ti_q_col[ord_q]
        r_minor_ord = r_minor_col[ord_q]
        valid_q = np.isfinite(q_sorted)
        dti_vs_q = np.gradient(ti_q_ord[valid_q], r_minor_ord[valid_q], edge_order=2)
        ax_dti_q.plot(
            q_sorted[valid_q],
            dti_vs_q * 1e-3,
            "-",
            color=color,
            lw=1.8,
            label="q-grid (pchip Ti mapped)",
        )
        ax_dti_q.set_title(rf"$dT_i/dr$ vs $q$ — $t={t_actual:.2f}\,\mathrm{{s}}$")
        ax_dti_q.set_xlabel(r"$q(\psi_N,\,t)$")
        ax_dti_q.set_ylabel(r"$dT_i/dr$ [keV/m]")
        ax_dti_q.legend(fontsize=8)
        ax_dti_q.grid(alpha=0.3)
        ax_dti_q.axhline(0.0, color="k", ls=":", lw=0.8)

    fig.tight_layout()

    if doSave:
        save_path = doSave + (
            f"ti_grid_diagnostic_shot_{shot}.pdf"
            if shot is not None
            else "ti_grid_diagnostic.pdf"
        )
        fig.savefig(save_path, transparent=True)
        print(f"Saved Ti diagnostic to {save_path}")

    if show_plot:
        plt.show(block=True)
    else:
        plt.close(fig)
