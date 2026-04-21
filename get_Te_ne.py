import numpy as np
import matplotlib.pyplot as plt
import mdsthin as mds


###############################################################################
def openTree(shotno, treeName="CMOD"):
    # Connect to data tree
    conn = mds.Connection("alcdata")
    conn.openTree(treeName, shotno)
    return conn


###############################################################################
def currentShot(conn):
    return conn.get('current_shot("cmod")').data()


###############################################################################
class YAG:
    # Nd:YAG Thomson scattering laser Te, ne profiles
    def __init__(self, shotno: int, debug: bool = False):
        if debug:
            print("Loading Thomson (YAG) Signal")

        conn = openTree(shotno)
        self.shotno = shotno if shotno != 0 else currentShot(conn)

        # only the old system is availible early on
        old = shotno < 1020000000
        edge = shotno > 1000000000  # edge TS only avilible here and newer

        ne_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.GLOBAL.PROFILE:NE_RZ_T"
            if old
            else r"\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:NE_RZ"
        )
        ne_err_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.GLOBAL.PROFILE:NE_ERR_ZT"
            if old
            else r"\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:NE_ERR"
        )
        te_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.GLOBAL.PROFILE:TE_RZ_T"
            if old
            else r"\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:TE_RZ"
        )
        te_err_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.GLOBAL.PROFILE:TE_ERR_ZT"
            if old
            else r"\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:TE_ERR"
        )
        r_mapped_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.GLOBAL.PROFILE:R_MID_T"
            if old
            else r"\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:R_MID_T"
        )
        z_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.GLOBAL.PROFILE:Z_SORTED"
            if old
            else r"\ELECTRONS::TOP.YAG_NEW.RESULTS.PROFILES:Z_SORTED"
        )
        r_node = (
            r"\ELECTRONS::TOP.YAG.RESULTS.PARAM:R"
            if old
            else r"\ELECTRONS::TOP.YAG.RESULTS.PARAM:R"
        )

        self.Ne = np.array(conn.get(ne_node).data())
        self.Ne_Err = np.array(conn.get(ne_err_node).data())
        self.Te = np.array(conn.get(te_node).data())
        self.Te_Err = np.array(conn.get(te_err_node).data())
        self.R_Map = np.array(conn.get(r_mapped_node).data())
        self.Z = np.array(conn.get(z_node).data())
        self.R = np.array(conn.get(r_node).data())
        self.time = np.array(conn.get("dim_of(" + ne_node + ")"))

        # If edge exists, pull it
        if edge:
            self.Ne_Edge = np.array(
                conn.get(r"\ELECTRONS::TOP.YAG_EDGETS.RESULTS:NE").data()
            )
            self.Ne_Err_Edge = np.array(
                conn.get(r"\ELECTRONS::TOP.YAG_EDGETS.RESULTS:NE:ERROR").data()
            )
            self.Te_Edge = (
                np.array(conn.get(r"\ELECTRONS::TOP.YAG_EDGETS.RESULTS:TE").data())
                * 1e-3
            )
            self.Te_Err_Edge = (
                np.array(
                    conn.get(r"\ELECTRONS::TOP.YAG_EDGETS.RESULTS:TE:ERROR").data()
                )
                * 1e-3
            )
            self.R_Map_Edge = np.array(
                conn.get(r"\ELECTRONS::TOP.YAG_EDGETS.RESULTS:RMID").data()
            )
            self.Z_Edge = np.array(
                conn.get(r"\ELECTRONS::TOP.YAG_EDGETS.DATA:FIBER_Z").data()
            )
            self.R_Edge = np.array(
                conn.get(r"\ELECTRONS::TOP.YAG.RESULTS.PARAM:R").data()
            )
            self.time_Edge = np.array(
                conn.get(r"dim_of(\ELECTRONS::TOP.YAG_EDGETS.RESULTS:NE)").data()
            )

        conn.closeAllTrees()

    def return_Profile(
        self,
        timePoint: float,
        dropChansMain: list = [3],
        dropChansEdge: list = [0, 1, 2, 3],
    ):
        # Select desired R points
        r_inds = np.delete(np.arange(len(self.R_Map)), dropChansMain)
        r_inds_Edge = np.delete(np.arange(len(self.R_Map_Edge)), dropChansEdge)

        tInd = np.argmin((self.time - timePoint) ** 2)
        tInd_Edge = np.argmin((self.time_Edge - timePoint) ** 2)

        # Check if timepoint is valid
        if self.R_Map[0, tInd] < 0:
            raise SyntaxError("Thomson Data Invalid at %2.2fs" % timePoint)

        Te = np.concatenate(
            (self.Te[r_inds, tInd], self.Te_Edge[r_inds_Edge, tInd_Edge])
        )
        Ne = np.concatenate(
            (self.Ne[r_inds, tInd], self.Ne_Edge[r_inds_Edge, tInd_Edge])
        )
        R = np.concatenate(
            (self.R_Map[r_inds, tInd], self.R_Map_Edge[r_inds_Edge, tInd_Edge])
        )

        r_sort = np.argsort(R)

        return Te[r_sort], Ne[r_sort], R[r_sort]

    def makePlot(
        self,
        time=1,
        dropChansMain=[3],
        dropChansEdge=[0, 1, 2, 3],
        doSave="",
        eqdsk=None,
    ):

        plt.close("Thomson_%d_%1.1f" % (self.shotno, time))
        fig, ax = plt.subplots(
            1,
            1,
            num="Thomson_%d_%1.1f" % (self.shotno, time),
            tight_layout=True,
            figsize=(4.5, 2.0),
        )

        tInd = np.argmin((self.time - time) ** 2)
        r_inds = np.delete(np.arange(len(self.R_Map)), dropChansMain)

        r_plot = (
            self.R_Map[r_inds, tInd]
            if eqdsk is None
            else eqdsk.rz2psinorm(
                self.R_Map[r_inds, tInd], 0, time, sqrt=True, make_grid=True
            )[0]
        )

        ax.errorbar(
            r_plot, self.Te[r_inds, tInd], fmt="*", yerr=self.Te_Err[r_inds, tInd]
        )
        ax1 = ax.twinx()
        ax1.errorbar(
            r_plot,
            self.Ne[r_inds, tInd] * 1e-20,
            yerr=self.Ne_Err[r_inds, tInd] * 1e-20,
            c=plt.get_cmap("tab10")(1),
            fmt="*",
        )

        tInd = np.argmin((self.time_Edge - time) ** 2)
        r_inds_Edge = np.delete(np.arange(len(self.R_Map_Edge)), dropChansEdge)

        r_plot_edge = (
            self.R_Map_Edge[r_inds_Edge, tInd]
            if eqdsk is None
            else eqdsk.rz2psinorm(
                self.R_Map_Edge[r_inds_Edge, tInd], 0, time, sqrt=True, make_grid=True
            )[0]
        )

        ax.errorbar(
            r_plot_edge,
            self.Te_Edge[r_inds_Edge, tInd],
            fmt="*",
            yerr=self.Te_Err_Edge[r_inds_Edge, tInd],
            c=plt.get_cmap("tab10")(0),
            alpha=0.7,
        )
        ax1.errorbar(
            r_plot_edge,
            self.Ne_Edge[r_inds_Edge, tInd] * 1e-20,
            yerr=self.Ne_Err_Edge[r_inds_Edge, tInd] * 1e-20,
            c=plt.get_cmap("tab10")(1),
            fmt="^",
            alpha=0.7,
        )

        p1 = ax.plot(0.8, 1, "k*", label="TS Core", ms=3)
        p2 = ax.plot(0.8, 1, "k^", label="TS Edge", ms=1)
        p3 = ax.plot(
            0.8, 1, "*", label=r"$\mathrm{T_e}$", ms=3, c=plt.get_cmap("tab10")(0)
        )
        p4 = ax.plot(
            0.8, 1, "*", label=r"$\mathrm{n_e}$", ms=3, c=plt.get_cmap("tab10")(1)
        )

        ax.legend(
            fontsize=8,
            title="%d: %1.1f s" % (self.shotno, time),
            title_fontsize=9,
            loc="lower left",
            ncol=2,
            columnspacing=0.5,
        )
        p1[0].remove()
        p2[0].remove()
        p3[0].remove()
        p4[0].remove()
        yl1 = ax.get_ylim()
        ax.set_ylim([0, yl1[1]])
        yl2 = ax1.get_ylim()
        ax1.set_ylim([0, yl2[1]])
        ax.grid()
        ax.set_xlabel("R [m]" if eqdsk is None else r"$\sqrt{\psi_n}$")
        ax.set_ylabel(r"T$_\mathrm{e}$ [keV]")
        ax1.set_ylabel(r"n$_\mathrm{e}$ [$10^{20}\,\mathrm{m}^{-3}$]")

        if doSave:
            fig.savefig(
                doSave + fig.canvas.manager.get_window_title() + ".pdf",
                transparent=True,
            )
        # safe_show()

        return ax, ax1


###############################################################################
