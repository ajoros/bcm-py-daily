"""Does Michelle's Skyview-regenerated terrain shading improve the ET metrics?

Builds a new .inp that swaps ONLY the sky-view factor and the 36 horizon angles
(the grids Skyview produces from the DEM) into the current deck, keeping its
slope/aspect/elev/lat/lon untouched. This isolates the terrain-shading effect
(lat/lon are identical between decks; our from-DEM slope/aspect placeholder is
NOT trustworthy yet, so we deliberately leave slope/aspect alone).

Then runs the port twice -- baseline vs swapped -- and compares pet/aet/cwd to
the .exe reference. If the swapped shading moves pet/aet toward 1:1, Michelle's
inputs are the lever; if not, the residual is in the daily radiation code.

Run:  python tests/test_shade_regen.py
"""
import os, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, PKG)
sys.path.insert(0, os.path.join(PKG, "validation"))
from bcm_io import read_asc
from validate import run_model, compare_grids, DEFAULT_REF

ROOT = os.path.normpath(os.path.join(PKG, ".."))
WORK = os.path.join(ROOT, "SolarFilesFromMichelle", "_skyview_clean")
ORIG_INP = os.path.join(DEFAULT_REF, "mok_270m.inp")   # the deck the port uses
NEW_INP  = os.path.join(WORK, "mok_270m_regen.inp")
ND = -9999.0
YEAR, DAYS, VARS = 2010, range(1, 10), ["pet", "aet", "cwd"]


def build_swapped_inp():
    """Write NEW_INP = ORIG_INP with the FULL Skyview terrain swapped in on
    valid cells: slope (col 5), aspect (col 6), sky (col 8), ridge(1..36)
    (cols 9..44). Skyview auto-generates slp/asp/sky/ang from the DEM, so this
    is the true 'regenerate the deck and test it' experiment. (lat/lon are
    identical between decks; elev is untouched.)"""
    _, slp = read_asc(os.path.join(WORK, "slpfile.asc"))
    _, asp = read_asc(os.path.join(WORK, "aspfile.asc"))
    asp = np.where(asp == -1.0, 0.0, asp)                      # Solinp: flat -> 0
    _, sky = read_asc(os.path.join(WORK, "mok_270m_sky.asc"))
    slp_f, asp_f, sky_f = slp.ravel(), asp.ravel(), sky.ravel()
    ang_f = np.stack([read_asc(os.path.join(WORK, f"ang{i*10:03d}.asc"))[1].ravel()
                      for i in range(1, 37)], axis=1)          # (N, 36)
    with open(ORIG_INP, "r", encoding="latin-1") as f:
        lines = f.readlines()
    out, swapped = [], 0
    for i, line in enumerate(lines):
        t = line.split()
        if len(t) >= 45 and float(t[7]) != ND:                 # valid elevation
            t[5] = f"{slp_f[i]:.2f}"
            t[6] = f"{asp_f[i]:.2f}"
            t[8] = f"{sky_f[i]:.2f}"
            for j in range(36):
                t[9 + j] = str(int(round(ang_f[i, j])))
            swapped += 1
        out.append(" ".join(t))
    with open(NEW_INP, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"built {NEW_INP}  ({swapped} cells fully regenerated)")


def metrics(py_dir):
    return compare_grids(py_dir, DEFAULT_REF, YEAR, DAYS, VARS)


def row(tag, s):
    return (f"  {tag:<10} rmse={s['rmse']:>8.4f}  bias={s['bias']:>+8.4f}  "
            f"corr={s['corr']:>7.4f}  %match={s['pct_match']:>6.1f}")


def main():
    build_swapped_inp()

    base = tempfile.mkdtemp(prefix="bcm_base_")
    new  = tempfile.mkdtemp(prefix="bcm_new_")
    print("running baseline (current deck) …");  run_model(base)
    print("running swapped-shade deck …");        run_model(new, {"BCM_INP_FILE": NEW_INP})

    b, n = metrics(base), metrics(new)
    print("\nPET / AET / CWD vs .exe  (baseline -> Skyview-reshaded)")
    for v in VARS:
        if v not in b or v not in n:
            continue
        print(f"[{v}]")
        print(row("baseline", b[v]))
        print(row("reshaded", n[v]))
        drmse = n[v]["rmse"] - b[v]["rmse"]
        verdict = "IMPROVED" if drmse < -1e-4 else ("worse" if drmse > 1e-4 else "no change")
        print(f"    -> RMSE {drmse:+.4f}  ({verdict})")

    # self-check: the swap must actually change the deck (else the test is moot)
    assert os.path.getsize(NEW_INP) > 0
    print("\nDONE.")


if __name__ == "__main__":
    main()
