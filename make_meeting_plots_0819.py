"""Meeting plots for 2026-08-19 — same Fortran-vs-Python table as the scorecard."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bcm_io import read_asc
from validate import parse_out

FORT = os.path.normpath(os.path.join(HERE, "..", "BCM_testrun_original"))
PY = os.environ.get("BCM_PY_DIR", os.path.join(os.environ.get("TEMP", ""), "bcm_snowmaps_0819"))
OUT = os.path.normpath(os.path.join(HERE, "..", "08192026meetingplots"))
os.makedirs(OUT, exist_ok=True)
ND = -9999.0
TOL = 0.01
DAY = 9

# status: exact / residual / exefloor / design / lowpri
ROWS = [
    ("tmn", "min air temp (tmn)", "out:tmnc", "exact"),
    ("tmx", "max air temp (tmx)", "out:tmxc", "exact"),
    ("ppt", "precipitation (ppt)", "out:pptmm", "exact"),
    ("snw", "snowfall (snw)", "grid", "exact"),
    ("rch", "recharge (rch)", "grid", "exact"),
    ("run", "runoff (run)", "grid", "exact"),
    ("exc", "excess (exc)", "grid", "residual"),
    ("pet", "PET (pet)", "grid", "residual"),
    ("sbl", "sublimation (sbl)", "grid", "residual"),
    ("aet", "actual ET (aet)", "grid", "residual"),
    ("pck", "snowpack (pck)", "grid", "residual"),
    ("mlt", "snowmelt (mlt)", "grid", "residual"),
    ("str", "soil storage (str)", "grid", "exefloor"),
    ("cwd", "clim. water deficit (cwd)", "grid", "design"),
    ("smd", "soil deficit (smd)", "grid", "lowpri"),
    ("smr", "soil moisture (smr)", "grid", "lowpri"),
    ("rvs", "rain/snow (rvs)", "grid", "lowpri"),
]
STATUS_COLOR = {
    "exact": "#2e7d32",
    "residual": "#1565c0",
    "exefloor": "#e65100",
    "design": "#6a1b9a",
    "lowpri": "#9e9e9e",
}


def load(folder, var, day):
    fp = os.path.join(folder, f"{var}2010_{day:03d}.asc")
    if not os.path.exists(fp):
        return None
    return read_asc(fp)[1].astype(float)


def metrics_grid(v):
    fs, ps = [], []
    for d in range(1, 10):
        f, p = load(FORT, v, d), load(PY, v, d)
        if f is None or p is None:
            return None
        m = (f != ND) & (p != ND)
        if m.any():
            fs.append(f[m]); ps.append(p[m])
    if not fs:
        return None
    f, p = np.concatenate(fs), np.concatenate(ps)
    d = p - f
    sst = np.sum((f - f.mean()) ** 2)
    return {
        "rmse": float(np.sqrt((d ** 2).mean())),
        "bias": float(d.mean()),
        "maxabs": float(np.abs(d).max()),
        "pct": float((np.abs(d) <= TOL).mean() * 100),
        "r2": float(1 - np.sum(d ** 2) / sst) if sst > 0 else 1.0,
    }


def metrics_out(col):
    py = parse_out(os.path.join(PY, "mok_Dailyv81_test.out"))
    fx = parse_out(os.path.join(FORT, "mok_Dailyv81_test.out"))
    keys = sorted(set(py) & set(fx))
    f = np.array([fx[k][col] for k in keys if col in fx[k] and col in py[k]])
    p = np.array([py[k][col] for k in keys if col in fx[k] and col in py[k]])
    if not f.size:
        return None
    d = p - f
    sst = np.sum((f - f.mean()) ** 2)
    return {
        "rmse": float(np.sqrt((d ** 2).mean())),
        "bias": float(d.mean()),
        "maxabs": float(np.abs(d).max()),
        "pct": float((np.abs(d) <= TOL).mean() * 100),
        "r2": float(1 - np.sum(d ** 2) / sst) if sst > 0 else 1.0,
    }


met = {}
for key, lab, src, status in ROWS:
    m = metrics_out(src.split(":")[1]) if src.startswith("out:") else metrics_grid(key)
    if m is None:
        print("SKIP", key)
        continue
    m["lab"], m["status"] = lab, status
    met[key] = m

print(f"py={PY}\nfort={FORT}\nout={OUT}")
print(f"{'var':6} {'rmse':>8} {'status'}")
for k, m in met.items():
    print(f"{k:6} {m['rmse']:8.4f} {m['status']}")

# --- 1. RMSE scoreboard ---
order = [k for k, *_ in ROWS if k in met]
fig, ax = plt.subplots(figsize=(11.2, 7.4))
vals = [max(met[k]["rmse"], 1e-4) for k in order]
names = [met[k]["lab"] for k in order]
colors = [STATUS_COLOR[met[k]["status"]] for k in order]
bars = ax.barh(names[::-1], vals[::-1], color=colors[::-1])
ax.set_xscale("log")
ax.set_xlabel("RMSE vs shipped Fortran exe  (mm; temp in °C)")
ax.set_title("Python vs Fortran RMSE — SNOW=0 (maps), full-precision storage\n"
             "green = 0   blue = residual   orange = exe still floors   "
             "purple = exe CWD bug   gray = not scored")
for b, k in zip(bars, order[::-1]):
    m = met[k]
    txt = "0.000" if m["rmse"] < 0.001 else f"{m['rmse']:.3f}"
    ax.text(max(b.get_width() * 1.2, 1.2e-3), b.get_y() + b.get_height() / 2,
            txt, va="center", fontsize=8)
ax.set_xlim(8e-4, max(vals) * 12)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "01_rmse_scoreboard.png"), dpi=150)
plt.close()

# --- 2. Table image ---
headers = ["Variable", "Michelle", "RMSE", "Note"]
note = {
    "tmn": "same input files", "tmx": "same input files", "ppt": "same input files",
    "snw": "100% of cells", "rch": "RMSE = 0", "run": "RMSE = 0",
    "exc": "max(ppt−AET, 0); follows AET",
    "pet": "compiler FP residual (~0.005)",
    "sbl": "uses PET; more snow → a bit more sbl",
    "aet": "PET residual × Kv",
    "pck": "thin-pack timing; tracks melt",
    "mlt": "same leftover as pack (maxabs 2.38)",
    "str": "exe floors ~1 mm/day; Python full precision",
    "cwd": "exe writes CWD≡PET; Python keeps PET−AET",
    "smd": "tracks storage (fc−str)",
    "smr": "tracks storage (por−str)",
    "rvs": "cold-day split only",
}
michelle = {k: ("not worried" if met[k]["status"] == "lowpri"
                else ("inputs" if k in ("tmn", "tmx", "ppt") else "priority"))
            for k in met}
rows = []
row_status = []
for k in order:
    m = met[k]
    rows.append([m["lab"], michelle[k],
                 "0.000" if m["rmse"] < 0.001 else f"{m['rmse']:.3f}", note[k]])
    row_status.append(m["status"])

fig, ax = plt.subplots(figsize=(13.2, 0.48 * (len(rows) + 2)))
ax.axis("off")
tbl = ax.table(cellText=rows, colLabels=headers, loc="center", cellLoc="left")
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.45)
widths = [0.28, 0.12, 0.10, 0.50]
for (r, c), cell in tbl.get_celld().items():
    cell.set_width(widths[c])
    if r == 0:
        cell.set_facecolor("#263238")
        cell.set_text_props(color="w", fontweight="bold", ha="center")
    else:
        cell.set_facecolor(STATUS_COLOR[row_status[r - 1]] + "22")
ax.set_title("Fortran vs Python  —  19 Aug 2026  —  SNOW=0 (maps), full-precision storage",
             fontsize=12, pad=10)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "02_scorecard_table.png"), dpi=160, bbox_inches="tight")
plt.close()

# --- 3. 1:1 scatters, day 9, priority + CWD ---
scatter = [
    ("snw", "snowfall"), ("rch", "recharge"), ("run", "runoff"),
    ("pet", "PET"), ("aet", "AET"), ("exc", "excess"),
    ("sbl", "sublimation"), ("pck", "snowpack"), ("mlt", "snowmelt"),
    ("str", "soil storage"), ("cwd", "CWD (exe bug)"),
]
fig, axes = plt.subplots(3, 4, figsize=(16, 12))
axes = axes.ravel()
for i, (v, lab) in enumerate(scatter):
    ax = axes[i]
    f, p = load(FORT, v, DAY), load(PY, v, DAY)
    if f is None or p is None:
        ax.set_title(f"{lab}  (no map)")
        ax.axis("off")
        continue
    msk = (f != ND) & (p != ND)
    fx, py = f[msk], p[msk]
    ax.scatter(fx, py, s=3, alpha=0.22, c="#1565c0", edgecolors="none")
    lo = min(fx.min(), py.min()); hi = max(fx.max(), py.max())
    pad = 0.04 * (hi - lo if hi > lo else 1)
    lim = [lo - pad, hi + pad]
    ax.plot(lim, lim, "k--", lw=1)
    rmse = float(np.sqrt(((py - fx) ** 2).mean()))
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal", "box")
    ax.set_xlabel(f"Fortran {v} (mm)", fontsize=8)
    ax.set_ylabel(f"Python {v} (mm)", fontsize=8)
    ax.set_title(f"{lab}   RMSE={rmse:.3f}", fontsize=10)
    ax.tick_params(labelsize=7)
axes[-1].axis("off")
axes[-1].text(0.05, 0.7, "Dashed line: Python = Fortran\nDay 9, every valid cell\n\n"
              "CWD sits off the line because\nthe exe writes CWD = PET.\nPython writes PET − AET.",
              fontsize=10, va="top", transform=axes[-1].transAxes)
fig.suptitle("Python vs Fortran  —  day 9  —  SNOW=0 (maps)", fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "03_scatter_1to1_day9.png"), dpi=140, bbox_inches="tight")
plt.close()


def three_maps(v, title, cmap, fname):
    f, p = load(FORT, v, DAY), load(PY, v, DAY)
    msk = (f != ND) & (p != ND)
    fm = np.where(msk, f, np.nan)
    pm = np.where(msk, p, np.nan)
    dm = np.where(msk, p - f, np.nan)
    vmax = np.nanpercentile(np.concatenate([fm[msk], pm[msk]]), 99)
    dabs = max(np.nanpercentile(np.abs(dm[msk]), 99), 1e-2)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.0))
    im0 = axes[0].imshow(fm, cmap=cmap, vmin=0, vmax=vmax)
    axes[0].set_title(f"Fortran {title}\nmean {np.nanmean(fm):.2f} mm")
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04, label="mm")
    im1 = axes[1].imshow(pm, cmap=cmap, vmin=0, vmax=vmax)
    axes[1].set_title(f"Python {title}\nmean {np.nanmean(pm):.2f} mm")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label="mm")
    im2 = axes[2].imshow(dm, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-dabs, vmax=dabs))
    axes[2].set_title(f"Python − Fortran\nRMSE {np.sqrt(np.nanmean(dm[msk]**2)):.3f} mm")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04, label="mm")
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{title}  —  day 9  —  SNOW=0 (maps)", fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, fname), dpi=140, bbox_inches="tight")
    plt.close()


three_maps("str", "soil storage (str)", "YlGnBu", "04_str_maps.png")
three_maps("cwd", "clim. water deficit (cwd)  —  exe writes CWD=PET", "YlOrRd", "05_cwd_maps.png")
three_maps("pet", "PET (pet)", "YlOrRd", "06_pet_maps.png")
three_maps("aet", "actual ET (aet)", "YlGn", "07_aet_maps.png")
three_maps("mlt", "snowmelt (mlt)", "Blues", "08_mlt_maps.png")
three_maps("exc", "excess (exc)", "Purples", "09_exc_maps.png")

# --- CWD diagnosis: exe CWD vs exe PET should be 1:1 ---
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
f_cwd, f_pet = load(FORT, "cwd", DAY), load(FORT, "pet", DAY)
p_cwd, p_pet, p_aet = load(PY, "cwd", DAY), load(PY, "pet", DAY), load(PY, "aet", DAY)
m1 = (f_cwd != ND) & (f_pet != ND)
axes[0].scatter(f_pet[m1], f_cwd[m1], s=3, alpha=0.2, c="#6a1b9a", edgecolors="none")
lim = [0, max(f_pet[m1].max(), f_cwd[m1].max())]
axes[0].plot(lim, lim, "k--", lw=1)
axes[0].set_xlabel("Fortran PET (mm)"); axes[0].set_ylabel("Fortran CWD (mm)")
axes[0].set_title(f"Fortran: CWD vs PET, day {DAY}\nR² = 1.0 — exe writes CWD ≡ PET")
axes[0].set_aspect("equal", "box")
m2 = (p_cwd != ND) & (p_pet != ND) & (p_aet != ND)
expect = np.maximum(p_pet[m2] - p_aet[m2], 0)
axes[1].scatter(expect, p_cwd[m2], s=3, alpha=0.2, c="#1565c0", edgecolors="none")
hi = max(expect.max(), p_cwd[m2].max())
axes[1].plot([0, hi], [0, hi], "k--", lw=1)
axes[1].set_xlabel("Python max(PET−AET, 0) (mm)")
axes[1].set_ylabel("Python CWD (mm)")
axes[1].set_title(f"Python: CWD vs PET−AET, day {DAY}\nmatches the correct definition")
axes[1].set_aspect("equal", "box")
fig.suptitle("CWD issue is the exe map, not the port  —  19 Aug 2026", fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "10_cwd_definition.png"), dpi=140, bbox_inches="tight")
plt.close()

# --- str RMSE grows because the exe still floors ---
days, str_rmse, str_bias = [], [], []
for d in range(1, 10):
    f, p = load(FORT, "str", d), load(PY, "str", d)
    msk = (f != ND) & (p != ND)
    dd = p[msk] - f[msk]
    days.append(d)
    str_rmse.append(float(np.sqrt((dd ** 2).mean())))
    str_bias.append(float(dd.mean()))
fig, ax = plt.subplots(figsize=(8.4, 4.6))
ax.plot(days, str_rmse, "o-", color="#e65100", lw=2, label="RMSE")
ax.plot(days, str_bias, "s--", color="#1565c0", lw=1.6, label="bias (Python − Fortran)")
ax.set_xlabel("Day of water year 2010")
ax.set_ylabel("mm")
ax.set_title("Soil storage: exe still floors, Python does not\n"
             f"gap grows ~1 mm/day  —  9-day RMSE = {met['str']['rmse']:.2f} mm")
ax.set_xticks(days)
ax.legend(frameon=False)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUT, "11_str_floor_drift.png"), dpi=140)
plt.close()

print("DONE ->", OUT)
print("files:", ", ".join(sorted(os.listdir(OUT))))
