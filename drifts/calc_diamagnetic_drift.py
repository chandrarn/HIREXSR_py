#!/usr/bin/env python3
"""CLI entrypoint for diamagnetic drift computation and plotting."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Support direct script execution as well as `python -m HIREXSR_py.drifts...`.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute time-dependent electron/ion diamagnetic frequencies and plot selected time slices versus q."
    )
    parser.add_argument("--shot", type=int, default=1120906030, help="Shot number")
    parser.add_argument("--tree", type=str, default="efit20", help="EFIT tree name")
    parser.add_argument(
        "--zi", type=float, default=1.0, help="Ion charge state Z_i for q_i = Z_i e"
    )
    parser.add_argument(
        "--plot-times",
        nargs="+",
        type=float,
        default=[0.6, 1, 1.3],
        help="Selected times [s] to plot (nearest available times are used)",
    )
    parser.add_argument(
        "--demo",
        default=False,
        action="store_true",
        help="Use fully synthetic time-dependent demo data (0.0:0.1:1.0 s)",
    )
    parser.add_argument(
        "--f_lims_omega",
        nargs=2,
        type=float,
        default=[-50, 50],
        help="Y-axis limits for omega_* plot (kHz)",
    )
    parser.add_argument(
        "--doSave", default="", type=str, help="Save plots to file path"
    )
    parser.add_argument(
        "--tht", type=int, default=0, help="THT branch number for HIREXSR MHD+ Tree"
    )
    parser.add_argument(
        "--line",
        type=int,
        default=2,
        help="Diagnostic line name for HIREXSR MHD+ Tree (e.g. 'core', 'edge')",
    )
    parser.add_argument(
        "--max-omega-err",
        type=float,
        default=10.0,
        help="Maximum allowed error in omega_* [kHz]",
    )
    parser.add_argument(
        "--corr1-loop-voltage-v",
        type=float,
        default=1.0,
        help="Loop voltage V_l [V] for impurity-main-ion correction model 1",
    )
    parser.add_argument(
        "--corr1-z-imp",
        type=float,
        default=16,
        help="Impurity charge state Z_i for correction model 1",
    )
    parser.add_argument(
        "--corr1-mu-imp-amu",
        type=float,
        default=18,
        help="Impurity atomic mass mu [amu] for correction model 1",
    )
    parser.add_argument(
        "--corr1-dilution",
        type=float,
        default=1.0,
        help="Impurity dilution factor f for correction model 1",
    )
    parser.add_argument(
        "--headless",
        default=False,
        action="store_true",
        help="Force non-interactive matplotlib backend (Agg)",
    )
    parser.add_argument(
        "--no-show",
        default=False,
        action="store_true",
        help="Do not open interactive plot windows (still saves figures)",
    )
    return parser.parse_args()


def _configure_matplotlib_backend(force_headless: bool = False) -> bool:
    """Configure matplotlib backend before importing pyplot-based modules.

    Returns True if non-interactive Agg backend is selected.
    """
    auto_headless = (
        os.environ.get("DISPLAY") is None and os.environ.get("WAYLAND_DISPLAY") is None
    )
    use_agg = force_headless or auto_headless
    if use_agg:
        os.environ["MPLBACKEND"] = "Agg"
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.use("Agg", force=True)
        plt.switch_backend("Agg")
    return use_agg


def main() -> None:
    args = parse_args()

    using_headless_backend = _configure_matplotlib_backend(force_headless=args.headless)

    from HIREXSR_py.drifts.compute import (  # type: ignore[reportMissingImports]
        _demo_equilibrium,
        _demo_profiles,
        compute_diamagnetic_drift_frequencies,
        load_equilibrium_for_shot,
        load_profiles_for_shot,
    )
    from HIREXSR_py.drifts.plotting import plot_diamagnetic_vs_q_times  # type: ignore[reportMissingImports]

    corr1_params = None
    if (
        args.corr1_loop_voltage_v is not None
        and args.corr1_z_imp is not None
        and args.corr1_mu_imp_amu is not None
    ):
        corr1_params = {
            "loop_voltage_V": float(args.corr1_loop_voltage_v),
            "z_imp": float(args.corr1_z_imp),
            "mu_imp_amu": float(args.corr1_mu_imp_amu),
            "dilution_factor": float(args.corr1_dilution),
        }

    if args.demo:
        profiles = _demo_profiles()
        equilibrium = _demo_equilibrium(profiles.time_diag)
    else:
        profiles = load_profiles_for_shot(args.shot, line=args.line, tht=args.tht)
        equilibrium = load_equilibrium_for_shot(args.shot, tree=args.tree)

    result = compute_diamagnetic_drift_frequencies(
        profiles=profiles,
        equilibrium=equilibrium,
        ion_charge_state=args.zi,
        selected_times_s=args.plot_times,
        shot=args.shot,
        do_diagnostic_plot=not (args.no_show or using_headless_backend),
    )

    plot_diamagnetic_vs_q_times(
        result,
        profiles,
        equilibrium,
        args.shot,
        args.plot_times,
        doSave=args.doSave,
        f_lims_omega=args.f_lims_omega,
        tht=args.tht,
        line=args.line,
        max_err_omega=args.max_omega_err,
        impurity_corr1_params=corr1_params,
        show_plot=not (args.no_show or using_headless_backend),
    )


if __name__ == "__main__":
    main()
    print("Done")
