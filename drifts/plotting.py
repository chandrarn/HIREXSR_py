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
from .impurity_corrections import (
    impurity_main_ion_delta_f_loop_voltage_hz,
    impurity_main_ion_delta_f_ti_gradient_hz,
)


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
        ax_f_tor.plot(
            q_s,
            fe_tor_khz,
            color=color,
            lw=2.0,
            ls="-",
            label=rf"$f_{{*e}}^{{\mathrm{{tor}}}}$ ({time_label})",
        )
        ax_f_tor.plot(
            q_s,
            fi_tor_khz,
            color=color,
            lw=2.0,
            ls="--",
            label=rf"$f_{{*i}}^{{\mathrm{{tor}}}}$ ({time_label})",
        )

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
                use_eq1=True,
            )[order_q]
            / 1e3
        )
        ax_f_tor.plot(
            q_s,
            delta_f2_khz,
            color=color,
            lw=1.8,
            ls=":",
            label=rf"$\Delta f_{{I-i}}^{{(\nabla T_i)}}$ ({time_label})",
        )

        if impurity_corr1_params is not None and all(
            k in impurity_corr1_params
            for k in ["loop_voltage_V", "z_imp", "mu_imp_amu"]
        ):
            delta_f1_khz = (
                impurity_main_ion_delta_f_loop_voltage_hz(
                    ti_eV=ti_col,
                    ni_m3=np.asarray(result.ni_m3[:, it_res], dtype=float),
                    r_major_m=r_major_col,
                    loop_voltage_V=float(impurity_corr1_params["loop_voltage_V"]),
                    z_imp=float(impurity_corr1_params["z_imp"]),
                    mu_imp_amu=float(impurity_corr1_params["mu_imp_amu"]),
                    dilution_factor=float(
                        impurity_corr1_params.get("dilution_factor", 1.0)
                    ),
                )[order_q]
                / 1e3
            )
            ax_f_tor.plot(
                q_s,
                delta_f1_khz,
                color=color,
                lw=1.8,
                ls="-.",
                label=rf"$\Delta f_{{I-i}}^{{(V_l)}}$ ({time_label})",
            )

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
        ax_f_tor.set_ylim(f_lims_f_tor)
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
