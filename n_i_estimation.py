from zeff_neo_python import zeff_neo
from get_Te_ne import YAG
import numpy as np
import matplotlib.pyplot as plt


def estimate_n_i(
    shotno: int,
    dropChansMain: list = [3],
    dropChansEdge: list = [0, 1, 2, 3],
    doPlot: bool = False,
    doSave: str = "",
):
    # Get YAG data
    yag = YAG(shotno)

    # Trim core channels
    keep_core = np.delete(np.arange(yag.Te.shape[0]), dropChansMain)
    te_core = yag.Te[keep_core, :]
    ne_core = yag.Ne[keep_core, :]
    r_core = yag.R_Map[keep_core, :]

    # check if edge Thomson data is active
    has_edge = all(
        hasattr(yag, attr) for attr in ["Te_Edge", "Ne_Edge", "R_Map_Edge", "time_Edge"]
    )

    if has_edge:
        # t_edge = np.asarray(yag.time_Edge, dtype=float).squeeze()
        # nt_edge = t_edge.size

        te_edge = np.asarray(yag.Te_Edge, dtype=float)
        ne_edge = np.asarray(yag.Ne_Edge, dtype=float)
        r_edge = np.asarray(yag.R_Map_Edge, dtype=float)

        keep_edge = np.delete(np.arange(te_edge.shape[0]), dropChansEdge)
        te_edge = te_edge[keep_edge, :]
        ne_edge = ne_edge[keep_edge, :]
        r_edge = r_edge[keep_edge, :]

        te_all = np.concatenate((te_core, te_edge), axis=0)
        ne_all = np.concatenate((ne_core, ne_edge), axis=0)
        r_all = np.concatenate((r_core, r_edge), axis=0)
    else:
        te_all = te_core
        ne_all = ne_core
        r_all = r_core

    # Get Zeff and Neo
    zeff, zeff_times = zeff_neo(
        shotno,
        dt=np.mean(np.diff(yag.time)),
        trange=[yag.time[0], yag.time[-1]],
        verbose=False,
    )

    # Ensure that the timebases of the YAG data and Zeff data are compatible
    zeff_interpolated = np.interp(yag.time, zeff_times, zeff)

    # Estimate ni using ne and Zeff

    ni_estimate = ne_all / zeff_interpolated  # Simplified estimation

    if doPlot:
        make_plots(
            shotno,
            yag.time,
            te_all,
            ne_all,
            r_all,
            zeff_interpolated,
            ni_estimate,
            doSave,
        )

    return ni_estimate, te_all, ne_all, r_all, zeff_interpolated, yag.time


#######################################################
def make_plots(
    shotno, yag_time, te_all, ne_all, r_all, zeff_interpolated, ni_estimate, doSave
):
    fig, ax = plt.subplots(4, 1, sharex=True, layout="constrained")

    time_2d = np.broadcast_to(yag_time, te_all.shape)
    te_masked = te_all
    te_masked[r_all <= 0] = (
        np.nan
    )  # Mask out invalid measurements (noted as negative R values)
    ne_masked = ne_all
    ne_masked[r_all <= 0] = (
        np.nan
    )  # Mask out invalid measurements (noted as negative R values)
    ni_estimate_masked = ni_estimate
    ni_estimate_masked[r_all <= 0] = np.nan  # Mask out
    ax[0].contourf(time_2d, r_all, te_masked, levels=100, cmap="viridis", zorder=-5)

    ax[1].contourf(
        time_2d, r_all, ne_masked * 1e-20, levels=100, cmap="viridis", zorder=-5
    )

    ax[2].contourf(
        time_2d,
        r_all,
        ni_estimate_masked * 1e-20,
        levels=100,
        cmap="viridis",
        zorder=-5,
    )

    ax[3].plot(
        yag_time, zeff_interpolated, label="shotno={}".format(shotno), color="orange"
    )

    for i in range(3):
        ax[i].set_ylabel("R [m]")
        ax[i].set_rasterization_zorder(-1)
        ax[i].set_ylim(np.nanmin(r_all[r_all > 0]), np.nanmax(r_all[r_all > 0]))
        cbar = plt.colorbar(ax[i].collections[0], ax=ax[i])
        if i == 0:
            cbar.set_label(r"T$_\mathrm{e}$ [keV]")
        elif i == 1:
            cbar.set_label(r"n$_\mathrm{e}$ [$10^{20}\,\mathrm{m}^{-3}$]")
        elif i == 2:
            cbar.set_label(r"n$_\mathrm{i}$ [$10^{20}\,\mathrm{m}^{-3}$]")

    ax[3].set_xlabel("Time [s]")
    ax[3].set_ylabel("$Z_{eff}$")
    ax[3].legend()
    ax[3].grid()
    ax[3].set_xlim([0, yag_time[-1]])

    if doSave:
        plt.savefig(doSave + "ni_estimate_{}.pdf".format(shotno), transparent=True)
        print(f"Plot saved as {doSave}ni_estimate_{shotno}.pdf")

    plt.show()


##################################3
if __name__ == "__main__":
    shotno = 1140221012
    ni_estimate, te_all, ne_all, r_all, zeff_interpolated, yag_time = estimate_n_i(
        shotno=shotno, doPlot=True
    )
    print(
        f"Estimated ion density for shot {shotno} at time {yag_time[0]:.2f} s: {ni_estimate[:, 0]} m^-3"
    )
