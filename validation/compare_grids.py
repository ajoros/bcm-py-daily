"""
compare_grids.py
Compares Python (BCM_testrun_python_v1) vs Fortran (BCM_testrun_original)
output .asc grids cell-by-cell and prints a summary table.
"""
import os, sys, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, PKG)
from bcm_io import read_asc

# Python maps live in the package folder. Fortran maps sit next to that package.
PY   = PKG
FORT = os.path.normpath(os.path.join(PKG, "..", "BCM_testrun_original"))

# prefixes that exist in both runs
VARS = ['aet','cwd','hdi','mwt','pck','pet','rad','rch','run','rvs','smd','smr','snw','str','ati']
DAYS = range(1, 10)
NODATA = -9999.0

print("=" * 90)
print(f"{'Variable':<6}  {'Day':<4}  {'Fort_mean':>10}  {'Py_mean':>10}  "
      f"{'RMSE':>10}  {'MaxAbsErr':>10}  {'% cells matching':>16}  Status")
print("=" * 90)

summary = {}   # var → list of RMSE

for var in VARS:
    rmse_list = []
    for day in DAYS:
        fname = f"{var}2010_{day:03d}.asc"
        fp = os.path.join(FORT, fname)
        pp = os.path.join(PY,   fname)
        if not (os.path.exists(fp) and os.path.exists(pp)):
            continue

        _, fa = read_asc(fp)
        _, pa = read_asc(pp)

        if fa.shape != pa.shape:
            print(f"{var:<6}  {day:<4}  SHAPE MISMATCH {fa.shape} vs {pa.shape}")
            continue

        valid = (fa != NODATA) & (pa != NODATA)
        if valid.sum() == 0:
            continue

        fv = fa[valid]; pv = pa[valid]
        diff = pv - fv
        rmse = np.sqrt((diff**2).mean())
        mae  = np.abs(diff).max()
        tol  = 0.01  # 0.01 mm tolerance for "matching"
        pct  = (np.abs(diff) < tol).mean() * 100
        fmn  = fv.mean()
        pmn  = pv.mean()
        rmse_list.append(rmse)

        status = "OK" if rmse < 1.0 else ("WARN" if rmse < 10.0 else "DIFF")
        print(f"{var:<6}  {day:<4}  {fmn:>10.3f}  {pmn:>10.3f}  "
              f"{rmse:>10.3f}  {mae:>10.3f}  {pct:>15.1f}%  {status}")

    if rmse_list:
        summary[var] = np.mean(rmse_list)

print()
print("=" * 90)
print("PER-VARIABLE AVERAGE RMSE ACROSS ALL 9 DAYS:")
print("-" * 50)
for var, rmse in sorted(summary.items(), key=lambda x: -x[1]):
    bar = "█" * min(int(rmse * 5), 40)
    print(f"  {var:<5}  avg RMSE = {rmse:8.4f} mm  {bar}")
print()
print("NOTE: RMSE computed over valid (non-nodata) cells only.")
print("      Differences >10mm are flagged DIFF, 1-10mm as WARN, <1mm as OK.")
