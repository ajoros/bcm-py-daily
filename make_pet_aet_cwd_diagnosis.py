"""PET / AET / CWD diagnosis figures for Michelle.

Findings this figure documents (2026-07-14, Mokelumne WY2010 days 1-9):
  * PET, AET: Python matches the exe closely (small per-cell scatter from the
    solar-radiation model; see difference maps).
  * CWD: the compiled .exe writes CWD *identical to PET* on every cell/day,
    i.e. the AET subtraction is missing in the exe. The .f90 source (line 1983)
    and the Python port both compute the textbook CWD = PET - AET.

Run:  BCM_PY_DIR=/path/to/py/outputs python make_pet_aet_cwd_diagnosis.py
"""
import os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from bcm_io import read_asc
FORT = os.environ.get("BCM_FORT_DIR",
                      os.path.normpath(os.path.join(HERE, "..", "BCM_testrun_original")))
PY   = os.environ.get("BCM_PY_DIR", HERE)
OUT  = os.path.join(HERE, "validation_plots"); os.makedirs(OUT, exist_ok=True)
ND   = -9999.0
DAY  = 9

def load(folder, var, day):
    _, a = read_asc(os.path.join(folder, f"{var}2010_{day:03d}.asc"))
    return a.astype(float)

def r2(x, y):
    ss_res = np.sum((y - x) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

# ---------------------------------------------------------------------------
# 1. CWD diagnosis: is the exe's CWD = PET, or = PET - AET?
# ---------------------------------------------------------------------------
pe, ae, ce = load(FORT,"pet",DAY), load(FORT,"aet",DAY), load(FORT,"cwd",DAY)
pp, ap, cp = load(PY,"pet",DAY),   load(PY,"aet",DAY),   load(PY,"cwd",DAY)
m = (pe!=ND)&(ae!=ND)&(ce!=ND)&(pp!=ND)&(ap!=ND)&(cp!=ND)

fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
panels = [
    ("exe CWD  vs  exe PET",            ce[m], pe[m], "exe PET (mm)",       "exe CWD (mm)"),
    ("exe CWD  vs  exe (PET - AET)",    ce[m], (pe-ae)[m], "exe PET - AET (mm)", "exe CWD (mm)"),
    ("Python CWD  vs  Python (PET-AET)",cp[m], (pp-ap)[m], "Python PET - AET (mm)", "Python CWD (mm)"),
]
for ax, (title, y, x, xl, yl) in zip(axes, panels):
    ax.scatter(x, y, s=3, alpha=0.25, color="#1f77b4", edgecolors="none")
    lim = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(lim, lim, "k--", lw=1, label="1:1")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal","box")
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(f"{title}\nR\u00b2 = {r2(x, y):.4f}   RMSE = {np.sqrt(((y-x)**2).mean()):.4f} mm")
    ax.legend(loc="upper left", fontsize=8)
plt.suptitle("CWD definition check (day 9): the compiled .exe writes CWD \u2261 PET; "
             "the .f90 source and Python port compute CWD = PET \u2212 AET", y=1.03, fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "cwd_definition_diagnosis.png"), dpi=140, bbox_inches="tight")
plt.close(); print("wrote cwd_definition_diagnosis.png")

# ---------------------------------------------------------------------------
# 2. PET & AET difference maps (Fortran / Python / diff), day 9
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15, 9.5))
for r, (v, label) in enumerate([("pet","potential ET"), ("aet","actual ET")]):
    f = load(FORT, v, DAY); p = load(PY, v, DAY)
    mm = (f!=ND)&(p!=ND)
    fm = np.where(mm, f, np.nan); pm = np.where(mm, p, np.nan); dm = np.where(mm, p-f, np.nan)
    vmax = np.nanpercentile(fm, 99)
    im0 = axes[r,0].imshow(fm, cmap="viridis", vmin=0, vmax=vmax)
    axes[r,0].set_title(f"Fortran {label} ({v})\nmean={np.nanmean(fm):.3f} mm")
    plt.colorbar(im0, ax=axes[r,0], fraction=0.046, pad=0.04, label="mm")
    im1 = axes[r,1].imshow(pm, cmap="viridis", vmin=0, vmax=vmax)
    axes[r,1].set_title(f"Python {label} ({v})\nmean={np.nanmean(pm):.3f} mm")
    plt.colorbar(im1, ax=axes[r,1], fraction=0.046, pad=0.04, label="mm")
    dabs = max(np.nanpercentile(np.abs(dm[mm]), 99), 1e-3)
    im2 = axes[r,2].imshow(dm, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-dabs, vmax=dabs))
    axes[r,2].set_title(f"Difference (Python - Fortran)\nRMSE={np.sqrt(np.nanmean(dm[mm]**2)):.3f}  bias={np.nanmean(dm[mm]):+.3f} mm")
    plt.colorbar(im2, ax=axes[r,2], fraction=0.046, pad=0.04, label="mm")
    for ax in axes[r]: ax.set_xticks([]); ax.set_yticks([])
plt.suptitle(f"PET and AET: Fortran vs Python, day {DAY} "
             "(differences are the solar-radiation-model residual, RMSE < 0.07 mm)", y=1.0, fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "pet_aet_difference_maps.png"), dpi=140, bbox_inches="tight")
plt.close(); print("wrote pet_aet_difference_maps.png")
print("DONE ->", OUT)
