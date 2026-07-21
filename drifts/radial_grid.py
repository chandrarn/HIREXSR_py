"""Radial-grid utilities for per-time profile operations."""

from __future__ import annotations

import numpy as np


def sort_columns_by_x(
    x_st: np.ndarray,
    *y_st: np.ndarray,
) -> tuple[np.ndarray, ...]:
    """Sort each time column by x and apply the same order to companion arrays.

    Parameters
    ----------
    x_st : ndarray, shape [nspace, nt]
        Radial coordinate per time slice.
    *y_st : ndarrays, each shape [nspace, nt]
        Profile arrays defined on x_st.

    Returns
    -------
    tuple
        (x_sorted, y1_sorted, y2_sorted, ...)
    """
    x = np.asarray(x_st, dtype=float)
    if x.ndim != 2:
        raise ValueError(f"Expected x_st to be 2D, got shape {x.shape}.")

    ys = [np.asarray(y, dtype=float) for y in y_st]
    for y in ys:
        if y.shape != x.shape:
            raise ValueError(
                f"Companion array shape mismatch: expected {x.shape}, got {y.shape}."
            )

    x_sorted = np.empty_like(x, dtype=float)
    y_sorted = [np.empty_like(y, dtype=float) for y in ys]

    for it in range(x.shape[1]):
        order = np.argsort(x[:, it])
        x_sorted[:, it] = x[:, it][order]
        for j, y in enumerate(ys):
            y_sorted[j][:, it] = y[:, it][order]

    return (x_sorted, *y_sorted)


def merge_nearby_channels_by_mean_radius(
    x_st: np.ndarray,
    *y_st: np.ndarray,
    min_dr_m: float = 1e-3,
) -> tuple[np.ndarray, ...]:
    """Merge channels with nearly identical mean radius.

    This is intended for overlapping core/edge Thomson channels that map to
    almost the same major radius but may have different values.
    """
    x = np.asarray(x_st, dtype=float)
    if x.ndim != 2:
        raise ValueError(f"Expected x_st to be 2D, got shape {x.shape}.")

    ys = [np.asarray(y, dtype=float) for y in y_st]
    for y in ys:
        if y.shape != x.shape:
            raise ValueError(
                f"Companion array shape mismatch: expected {x.shape}, got {y.shape}."
            )

    if min_dr_m <= 0.0 or x.shape[0] <= 1:
        return (x, *ys)

    r_mean = np.nanmean(x, axis=1)
    order = np.argsort(r_mean)
    xs = x[order, :]
    ys_s = [y[order, :] for y in ys]
    rm = r_mean[order]

    groups: list[list[int]] = [[0]]
    for i in range(1, xs.shape[0]):
        if (
            np.isfinite(rm[i])
            and np.isfinite(rm[i - 1])
            and (rm[i] - rm[i - 1]) < float(min_dr_m)
        ):
            groups[-1].append(i)
        else:
            groups.append([i])

    x_out = np.empty((len(groups), xs.shape[1]), dtype=float)
    y_out = [np.empty((len(groups), xs.shape[1]), dtype=float) for _ in ys_s]

    for ig, g in enumerate(groups):
        x_out[ig, :] = np.nanmean(xs[g, :], axis=0)
        for j, y in enumerate(ys_s):
            y_out[j][ig, :] = np.nanmean(y[g, :], axis=0)

    return (x_out, *y_out)
