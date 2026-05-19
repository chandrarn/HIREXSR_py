"""C-Mod profile helpers for HIREXSR drift calculations.

This module intentionally duplicates the minimal YAG/openTree/currentShot logic
from Synthetic_Mirnov C-Mod utilities so the drifts submodule does not depend on
get_Cmod_Data.py.
"""

from __future__ import annotations

import numpy as np
import mdsthin as mds


def openTree(shotno: int, treeName: str = "CMOD"):
    """Connect to a CMOD MDS tree and return an open connection."""
    conn = mds.Connection("alcdata")
    conn.openTree(treeName, shotno)
    return conn


def currentShot(conn) -> int:
    """Return the active CMOD shot on a live MDS connection."""
    return int(conn.get('current_shot("cmod")').data())


class YAG:
    """Nd:YAG Thomson scattering Te/ne profiles for C-Mod."""

    def __init__(self, shotno: int, debug: bool = False):
        if debug:
            print("Loading Thomson (YAG) Signal")

        conn = openTree(shotno)
        self.shotno = shotno if shotno != 0 else currentShot(conn)

        old = shotno < 1020000000
        edge = shotno > 1105000000

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
        r_node = r"\ELECTRONS::TOP.YAG.RESULTS.PARAM:R"

        self.Ne = np.array(conn.get(ne_node).data())
        self.Ne_Err = np.array(conn.get(ne_err_node).data())
        self.Te = np.array(conn.get(te_node).data())
        self.Te_Err = np.array(conn.get(te_err_node).data())
        self.R_Map = np.array(conn.get(r_mapped_node).data())
        self.Z = np.array(conn.get(z_node).data())
        self.R = np.array(conn.get(r_node).data())
        self.time = np.array(conn.get("dim_of(" + ne_node + ")"))

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
