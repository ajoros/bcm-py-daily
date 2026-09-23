"""Generate Python-vs-Fortran validation figures for the 9-day Mokelumne test.

Outputs (PNG) into ./validation_plots:
  1. rmse_by_variable.png   - avg cell-level RMSE per variable (log scale)
  2. str_maps.png           - Fortran vs Python soil-storage maps + difference
  3. scatter_str_aet.png    - cell-by-cell Python vs Fortran (1:1 line, R2)
"""
import os, sys, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, PKG)
from bcm_io import read_asc
# Fortran (.exe) reference grids; Python grids default to the package folder but can be
# pointed at a fresh run dir via BCM_PY_DIR (keeps this cross-platform).
FORT = os.environ.get("BCM_FORT_DIR",
                      os.path.normpath(os.path.join(PKG, "..", "BCM_testrun_original")))
PY   = os.environ.get("BCM_PY_DIR", PKG)
OUT  = os.path.join(PKG, "validation_plots"); os.makedirs(OUT, exist_ok=True)
ND   = -9999.0

def load(folder, var, day):
    fp = os.path.join(folder, f"{var}2010_{day:03d}.asc")
    if not os.path.exists(fp): return None
    _, a = read_asc(fp)
    return a.astype(float)

def valid_mask(*arrs):
    m = np.ones_like(arrs[0], dtype=bool)
    for a in arrs:
        m &= (a != ND)
    return m

# ----------------------------------------------------------------------------
# 1. Error metrics per variable — Michelle's prioritized set, incl. temp/precip
#    Outputs -> cell-level grids; temp/precip inputs -> basin-summary .out
#    (identical input files, so the meaningful validated metric is the .out).
# ----------------------------------------------------------------------------
from validate import parse_out
TOL = 0.01

def _metrics(f, p):
    d = p - f
    ss_tot = np.sum((f - f.mean())**2)
    return {"rmse": float(np.sqrt((d**2).mean())),
            "bias": float(d.mean()),
            "maxabs": float(np.abs(d).max()),
            "pct": float((np.abs(d) <= TOL).mean()*100.0),
            "r2": float(1 - np.sum(d**2)/ss_tot) if ss_tot > 0 else 1.0,
            "n": int(f.size)}

def cell_metrics(v):
    fs, ps = [], []
    for d in range(1, 10):
        f = load(FORT, v, d); p = load(PY, v, d)
        if f is None or p is None: return None
        m = valid_mask(f, p)
        if m.sum() == 0: continue
        fs.append(f[m]); ps.append(p[m])
    return _metrics(np.concatenate(fs), np.concatenate(ps)) if fs else None

def out_metrics(col):
    py = parse_out(os.path.join(PY, "mok_Dailyv81_test.out"))
    fx = parse_out(os.path.join(FORT, "mok_Dailyv81_test.out"))
    keys = sorted(set(py) & set(fx))
    f = np.array([fx[k][col] for k in keys if col in fx[k] and col in py[k]])
    p = np.array([py[k][col] for k in keys if col in fx[k] and col in py[k]])
    return _metrics(f, p) if f.size else None

# (var, label, unit, source) — Michelle's equally-important list + temp/precip.
# mlt/sbl/exc are priority too but the .exe wrote no reference grids for them
# (map flags off in the .ctl), so they can't be graded here yet — noted below.
SPEC = [("tmn", "min air temp", "\u00b0C", "out:tmnc"),
        ("tmx", "max air temp", "\u00b0C", "out:tmxc"),
        ("ppt", "precipitation", "mm", "out:pptmm"),
        ("snw", "snowfall",      "mm", "grid"),
        ("pck", "snowpack",      "mm", "grid"),
        ("rch", "recharge",      "mm", "grid"),
        ("run", "runoff",        "mm", "grid"),
        ("aet", "actual ET",     "mm", "grid"),
        ("pet", "PET",           "mm", "grid"),
        ("cwd", "clim. water deficit", "mm", "grid"),
        ("str", "soil storage",  "mm", "grid")]
NO_REF = ["mlt (snowmelt)", "sbl (sublimation)", "exc (excess = PPT-PET)"]

met, unit_of, lab_of = {}, {}, {}
for v, lab, unit, src in SPEC:
    m = out_metrics(src.split(":")[1]) if src.startswith("out:") else cell_metrics(v)
    if m is None: continue
    met[v] = m; unit_of[v] = unit; lab_of[v] = lab

print("\nError metrics vs Fortran (Michelle's priority variables)")
print(f"{'variable':>26} {'unit':>4} {'RMSE':>10} {'bias':>10} {'maxabs':>10} {'%match':>8} {'R2':>8}")
for v in met:
    m = met[v]
    print(f"{lab_of[v]+' ('+v+')':>26} {unit_of[v]:>4} {m['rmse']:>10.4f} "
          f"{m['bias']:>+10.4f} {m['maxabs']:>10.4f} {m['pct']:>7.1f}% {m['r2']:>8.4f}")
print("  no reference grid yet (exe map flags off): " + ", ".join(NO_REF))

order = sorted(met, key=lambda k: met[k]["rmse"])
vals   = [max(met[k]["rmse"], 1e-4) for k in order]
names  = [f"{lab_of[k]} ({k})  [{unit_of[k]}]" for k in order]
colors = ["#2ca02c" if met[k]["rmse"] < 1 else ("#ff9900" if met[k]["rmse"] < 20 else "#d62728") for k in order]

fig, ax = plt.subplots(figsize=(10.5, 6.2))
bars = ax.barh(names, vals, color=colors)
ax.set_xscale("log")
ax.set_xlabel("RMSE vs Fortran, log scale  (mm, except temp in \u00b0C)")
ax.set_title("Python vs Fortran — Michelle's priority variables, 9-day Mokelumne test\n"
             "(smaller = closer; green <1, orange <20; temp/precip from basin-summary .out)")
for b, k in zip(bars, order):
    m = met[k]
    txt = ("0.000" if m["rmse"] < 0.001 else f"{m['rmse']:.3f}") + f"  (bias {m['bias']:+.3f})"
    xtext = max(b.get_width()*1.25, 1.1e-3)   # keep zero-bars' labels off the y-axis
    ax.text(xtext, b.get_y()+b.get_height()/2, txt, va="center", fontsize=8)
ax.set_xlim(left=8e-4, right=max(vals)*8)
ax.text(0.99, 0.02, "not yet graded (no exe reference grid): " + ", ".join(NO_REF),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="#666")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "rmse_by_variable.png"), dpi=140)
plt.close()
print("wrote rmse_by_variable.png")

# ----------------------------------------------------------------------------
# 1b. Same metrics as a standalone table image (for slides/reports)
# ----------------------------------------------------------------------------
TABLE_ORDER = ["tmn", "tmx", "ppt", "snw", "rch", "run", "aet", "pet", "pck", "cwd", "str"]
headers = ["Variable", "Unit", "RMSE", "Bias", "MaxAbs", "%Match", "R\u00b2"]
rows = []
for k in TABLE_ORDER:
    if k not in met: continue
    m = met[k]
    rows.append([f"{lab_of[k]} ({k})", unit_of[k], f"{m['rmse']:.4f}",
                 f"{m['bias']:+.4f}", f"{m['maxabs']:.3f}", f"{m['pct']:.1f}%",
                 f"{m['r2']:.4f}"])

fig, ax = plt.subplots(figsize=(11, 0.5 + 0.42*(len(rows)+1)))
ax.axis("off")
tbl = ax.table(cellText=rows, colLabels=headers, cellLoc="center", loc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(11); tbl.scale(1, 1.5)
col_w = [0.26, 0.09, 0.13, 0.13, 0.13, 0.13, 0.13]  # wide 1st col so names fit
for (r, c), cell in tbl.get_celld().items():
    cell.set_width(col_w[c])
for c in range(len(headers)):                      # header styling
    cell = tbl[0, c]; cell.set_facecolor("#33475b")
    cell.set_text_props(color="w", fontweight="bold")
for r, k in enumerate(TABLE_ORDER, start=1):        # row color by RMSE band
    if k not in met: continue
    band = "#e7f4e4" if met[k]["rmse"] < 1 else ("#fde6cf" if met[k]["rmse"] < 20 else "#f9d6d5")
    for c in range(len(headers)):
        tbl[r, c].set_facecolor(band)
        if c == 0: tbl[r, c].set_text_props(ha="left")
ax.set_title("Python vs Fortran error metrics \u2014 Michelle's priority variables\n"
             "9-day Mokelumne test (temp/precip from basin-summary .out)",
             fontsize=12, pad=14)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "metrics_table.png"), dpi=160, bbox_inches="tight")
plt.close()
print("wrote metrics_table.png")

# ----------------------------------------------------------------------------
# 2. Maps: Fortran, Python, difference (day 9) for the flagship metrics
#    CWD (climatic water deficit) and STR (soil storage) — Michelle's key outputs.
# ----------------------------------------------------------------------------
day = 9
MAP_VARS = [("cwd", "clim. water deficit", "YlOrRd"),
            ("str", "soil storage", "YlGnBu")]
for v, unit, cmap in MAP_VARS:
    f = load(FORT, v, day); p = load(PY, v, day)
    m = valid_mask(f, p)
    fm = np.where(m, f, np.nan); pm = np.where(m, p, np.nan); dm = np.where(m, p-f, np.nan)
    vmax = np.nanpercentile(np.concatenate([fm[m], pm[m]]), 99)
    dabs = max(np.nanpercentile(np.abs(dm[m]), 99), 1e-2)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    im0 = axes[0].imshow(fm, cmap=cmap, vmin=0, vmax=vmax)
    axes[0].set_title(f"Fortran {unit} ({v}), day {day}\nmean = {np.nanmean(fm):.2f} mm")
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04, label="mm")
    im1 = axes[1].imshow(pm, cmap=cmap, vmin=0, vmax=vmax)
    axes[1].set_title(f"Python {unit} ({v}), day {day}\nmean = {np.nanmean(pm):.2f} mm")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label="mm")
    im2 = axes[2].imshow(dm, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-dabs, vmax=dabs))
    axes[2].set_title(f"Difference (Python − Fortran)\nRMSE = {np.sqrt(np.nanmean(dm[m]**2)):.3f} mm  "
                      f"bias = {np.nanmean(dm[m]):+.3f} mm")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04, label="mm")
    for ax in axes: ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"{v}_maps.png"), dpi=140, bbox_inches="tight")
    plt.close()
    print(f"wrote {v}_maps.png")

# ----------------------------------------------------------------------------
# 3. Scatter Python vs Fortran (1:1 line) for the flagship metrics, day 9:
#    PET (driver) -> AET -> CWD, plus STR.
# ----------------------------------------------------------------------------
def r2(x, y):
    ss_res = np.sum((y-x)**2); ss_tot = np.sum((y-np.mean(y))**2)
    return 1 - ss_res/ss_tot if ss_tot > 0 else 1.0

scatter_vars = [("pet", "potential ET"), ("aet", "actual ET"),
                ("cwd", "clim. water deficit"), ("str", "soil storage")]
fig, axes = plt.subplots(1, 4, figsize=(20, 5.4))
for ax, (v, unit) in zip(axes, scatter_vars):
    f = load(FORT, v, day); p = load(PY, v, day); m = valid_mask(f, p)
    fx = f[m]; py = p[m]
    ax.scatter(fx, py, s=3, alpha=0.25, color="#1f77b4", edgecolors="none")
    lim = [min(fx.min(), py.min()), max(fx.max(), py.max())]
    ax.plot(lim, lim, "k--", lw=1, label="perfect 1:1")
    ax.set_xlabel(f"Fortran {unit} ({v}), mm")
    ax.set_ylabel(f"Python {unit} ({v}), mm")
    ax.set_title(f"{unit} ({v}), day {day}\nR\u00b2 = {r2(fx, py):.4f}   RMSE = {np.sqrt(((py-fx)**2).mean()):.3f} mm")
    ax.legend(loc="upper left", fontsize=8); ax.set_aspect("equal", "box")
plt.suptitle("Python vs Fortran on Michelle's flagship metrics (PET → AET → CWD, and soil storage)", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "scatter_flagship.png"), dpi=140, bbox_inches="tight")
plt.close()
print("wrote scatter_flagship.png")

# ----------------------------------------------------------------------------
# 4. Scatter Python vs Fortran for ALL 9 days (rows = aet, cwd; cols = day)
# ----------------------------------------------------------------------------
panel_vars = [("aet", "actual ET"), ("cwd", "clim. water deficit")]
days = list(range(1, 10))

# Per-variable common axis limits across all days (for fair comparison)
lims = {}
for v, _ in panel_vars:
    lo, hi = np.inf, -np.inf
    for d in days:
        f = load(FORT, v, d); p = load(PY, v, d)
        if f is None or p is None: continue
        m = valid_mask(f, p)
        if m.sum() == 0: continue
        lo = min(lo, f[m].min(), p[m].min()); hi = max(hi, f[m].max(), p[m].max())
    pad = 0.03 * (hi - lo)
    lims[v] = (lo - pad, hi + pad)

fig, axes = plt.subplots(len(panel_vars), len(days), figsize=(26, 6.2),
                         sharex="row", sharey="row")
for r, (v, unit) in enumerate(panel_vars):
    for c, d in enumerate(days):
        ax = axes[r, c]
        f = load(FORT, v, d); p = load(PY, v, d); m = valid_mask(f, p)
        fx = f[m]; py = p[m]
        ax.scatter(fx, py, s=2, alpha=0.20, color="#1f77b4", edgecolors="none")
        ax.plot(lims[v], lims[v], "k--", lw=0.8)
        ax.set_xlim(lims[v]); ax.set_ylim(lims[v])
        ax.set_aspect("equal", "box")
        ax.set_title(f"day {d}\nR\u00b2={r2(fx, py):.3f}  RMSE={np.sqrt(((py-fx)**2).mean()):.2f}",
                     fontsize=8)
        ax.tick_params(labelsize=6)
        if c == 0:
            ax.set_ylabel(f"Python\n{unit} ({v}), mm", fontsize=9)
        if r == len(panel_vars) - 1:
            ax.set_xlabel("Fortran, mm", fontsize=8)
plt.suptitle("Python vs Fortran, every grid cell, all 9 days (dashed = perfect 1:1)",
             y=1.01, fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "scatter_flagship_days1-9.png"), dpi=140, bbox_inches="tight")
plt.close()
print("wrote scatter_flagship_days1-9.png")
print("DONE ->", OUT)
