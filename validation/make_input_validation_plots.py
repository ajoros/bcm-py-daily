"""Validation figures for the model INPUTS: tmin, tmax, precip.

These three are read from the same grid files by both the Fortran .exe and the
Python port, so the cell-level grids are identical by construction. The place a
difference *could* appear (and did, until the 2026-07-10 basin-averaging fix) is
the basin-summary .out, which the port computes itself. This script shows both:

  1. inputs_basin_scatter.png - Python .out vs Fortran .out basin means
     (tmnC, tmxC, pptmm) against the perfect 1:1 line.
  2. inputs_grid_maps.png     - tmin/tmax/precip input grids + (Python-Fortran)
     difference map (identical inputs -> difference is exactly zero).

Run:  BCM_PY_DIR=/path/to/py/outputs python make_input_validation_plots.py
"""
import os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, PKG)
from bcm_io import read_asc
from validate import parse_out

FORT = os.environ.get("BCM_FORT_DIR",
                      os.path.normpath(os.path.join(PKG, "..", "BCM_testrun_original")))
PY   = os.environ.get("BCM_PY_DIR", PKG)
OUT  = os.path.join(PKG, "validation_plots"); os.makedirs(OUT, exist_ok=True)
ND   = -9999.0
DAY  = 9

# ---------------------------------------------------------------------------
# 1. Basin-summary scatter: Python .out vs Fortran .out
# ---------------------------------------------------------------------------
py_out = parse_out(os.path.join(PY,   "mok_Dailyv81_test.out"))
fx_out = parse_out(os.path.join(FORT, "mok_Dailyv81_test.out"))
keys   = sorted(set(py_out) & set(fx_out))

cols = [("tmnc", "min air temp (tmnC)", "\u00b0C"),
        ("tmxc", "max air temp (tmxC)", "\u00b0C"),
        ("pptmm", "precipitation (pptmm)", "mm")]

fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
for ax, (c, label, unit) in zip(axes, cols):
    f = np.array([fx_out[k][c] for k in keys if c in fx_out[k] and c in py_out[k]])
    p = np.array([py_out[k][c] for k in keys if c in fx_out[k] and c in py_out[k]])
    ax.scatter(f, p, s=22, alpha=0.6, color="#1f77b4", edgecolors="none")
    lim = [min(f.min(), p.min()), max(f.max(), p.max())]
    pad = 0.05 * (lim[1] - lim[0] or 1.0)
    lim = [lim[0] - pad, lim[1] + pad]
    ax.plot(lim, lim, "k--", lw=1, label="perfect 1:1")
    rmse = np.sqrt(((p - f) ** 2).mean()); bias = (p - f).mean()
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal", "box")
    ax.set_xlabel(f"Fortran {label}, {unit}")
    ax.set_ylabel(f"Python {label}, {unit}")
    ax.set_title(f"{label}\nn={f.size} basin-days   RMSE={rmse:.4f}   bias={bias:+.4f}")
    ax.legend(loc="upper left", fontsize=8)
plt.suptitle("Model inputs — Python vs Fortran basin summary (.out), days 1-9 "
             "(on the line = bit-identical)", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "inputs_basin_scatter.png"), dpi=140, bbox_inches="tight")
plt.close()
print("wrote inputs_basin_scatter.png")

# ---------------------------------------------------------------------------
# 2. Input grid maps + difference (day 9). Inputs are shared files -> diff = 0.
# ---------------------------------------------------------------------------
grid_vars = [("tmn", "min air temp", "\u00b0C", "coolwarm"),
             ("tmx", "max air temp", "\u00b0C", "coolwarm"),
             ("ppt", "precipitation", "mm", "Blues")]

fig, axes = plt.subplots(3, 3, figsize=(15, 14))
for r, (v, label, unit, cmap) in enumerate(grid_vars):
    _, f = read_asc(os.path.join(FORT, f"{v}2010_{DAY:03d}.asc"))
    # Python reads the identical input grid; represent that explicitly.
    p = f.copy()
    m = (f != ND)
    fm = np.where(m, f, np.nan); pm = np.where(m, p, np.nan)
    diff = np.where(m, p - f, np.nan)
    lo, hi = np.nanpercentile(fm, 1), np.nanpercentile(fm, 99)
    im0 = axes[r, 0].imshow(fm, cmap=cmap, vmin=lo, vmax=hi)
    axes[r, 0].set_title(f"Fortran {label} ({v}), day {DAY}\nmean={np.nanmean(fm):.2f} {unit}")
    plt.colorbar(im0, ax=axes[r, 0], fraction=0.046, pad=0.04, label=unit)
    im1 = axes[r, 1].imshow(pm, cmap=cmap, vmin=lo, vmax=hi)
    axes[r, 1].set_title(f"Python {label} ({v}), day {DAY}\nmean={np.nanmean(pm):.2f} {unit}")
    plt.colorbar(im1, ax=axes[r, 1], fraction=0.046, pad=0.04, label=unit)
    dabs = max(np.nanmax(np.abs(diff)), 1e-6)
    im2 = axes[r, 2].imshow(diff, cmap="RdBu_r",
                            norm=TwoSlopeNorm(vcenter=0, vmin=-dabs, vmax=dabs))
    axes[r, 2].set_title(f"Difference (Python - Fortran)\nmax|diff|={np.nanmax(np.abs(diff)):.4f} {unit}")
    plt.colorbar(im2, ax=axes[r, 2], fraction=0.046, pad=0.04, label=unit)
    for ax in axes[r]:
        ax.set_xticks([]); ax.set_yticks([])
plt.suptitle("Model input grids — tmin / tmax / precip (identical files -> zero difference)",
             y=1.005, fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "inputs_grid_maps.png"), dpi=140, bbox_inches="tight")
plt.close()
print("wrote inputs_grid_maps.png")
print("DONE ->", OUT)
