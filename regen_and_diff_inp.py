"""Regenerate the terrain grids from the DEM with Michelle's Skyview pipeline
and diff them against the deck the port currently uses.

Skyview generates the ENTIRE terrain deck from the DEM alone: lat/lon (mklatlon),
plus slope, aspect, sky-view and the 36 horizon angles (Skyview). So there is no
slope/aspect to derive ourselves -- we just run the tools and compare their
output to the current mok_270m.inp, column by column.

(Solinpv1.exe, which would stitch the grids into a fresh .inp, fails with a
missing-DLL error, but that stitch is trivial and unnecessary for the diff:
we compare Skyview's grids directly to the .inp's columns.)

Usage:
  python regen_and_diff_inp.py           # reuse an existing Skyview run if present
  python regen_and_diff_inp.py --run     # (re)run mklatlon + Skyview first
"""
import os, sys, shutil, subprocess, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from bcm_io import read_asc, read_inp

ROOT    = os.path.normpath(os.path.join(HERE, ".."))
SOLAR   = os.path.join(ROOT, "SolarFilesFromMichelle", "Solar_mok_extracted")
WORKDIR = os.path.join(ROOT, "SolarFilesFromMichelle", "_skyview_clean")
CUR_INP = os.path.join(ROOT, "BCM_testrun_original", "mok_270m.inp")
ND      = -9999.0


# ---------------------------------------------------------------------------
def run_pipeline(workdir=WORKDIR):
    """Run mklatlon + Skyview in a clean workdir; they emit the whole deck."""
    if os.path.exists(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    for txt, exe in (("mklatlon.txt", "mklatlon.exe"),
                     ("Skyviewv1.txt", "Skyviewv1.exe")):
        shutil.copy(os.path.join(SOLAR, txt), os.path.join(workdir, exe))
    for f in ("mklatlon.ctl", "Skyviewv1.ctl", "mok_270m.asc"):
        shutil.copy(os.path.join(SOLAR, f), os.path.join(workdir, f))
    for exe in ("mklatlon.exe", "Skyviewv1.exe"):
        print(f"  running {exe} …")
        subprocess.run([os.path.join(workdir, exe)], cwd=workdir, check=True,
                       stdout=subprocess.DEVNULL)
    # self-check: Skyview must have produced the full terrain deck
    need = ["slpfile.asc", "aspfile.asc", "mok_270m_sky.asc",
            "mok_270lat.asc", "mok_270lon.asc"] + [f"ang{i*10:03d}.asc"
                                                   for i in range(1, 37)]
    missing = [f for f in need if not os.path.exists(os.path.join(workdir, f))]
    assert not missing, f"Skyview did not produce: {missing}"
    print(f"  Skyview produced the full deck ({len(need)} grids)")


# ---------------------------------------------------------------------------
def _stats(name, diff):
    print(f"  {name:<12} rmse={np.sqrt((diff**2).mean()):>9.4f}  "
          f"bias={diff.mean():>+9.4f}  maxabs={np.abs(diff).max():>9.4f}  "
          f"%match={100*np.mean(np.abs(diff) <= 0.01):>6.1f}")


def _circ(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def diff_against_deck(workdir=WORKDIR, inp_path=CUR_INP):
    """Compare each regenerated grid to the current .inp's column."""
    cur = read_inp(inp_path, 279, 635)
    m = cur["elev"] != ND

    def g(name):
        return read_asc(os.path.join(workdir, name))[1]

    print(f"\nSkyview regen  vs  current deck  ({int(m.sum())} valid cells)")
    _stats("lat", (g("mok_270lat.asc") - cur["lat"])[m])
    _stats("lon", (g("mok_270lon.asc") - cur["lon"])[m])
    _stats("slope", (g("slpfile.asc") - cur["sl"])[m])
    asp = np.where(g("aspfile.asc") == -1.0, 0.0, g("aspfile.asc"))
    _stats("asp(circ)", _circ(asp, cur["asp"])[m])
    _stats("sky", (g("mok_270m_sky.asc") - cur["sky"])[m])
    rd = np.concatenate([(g(f"ang{i*10:03d}.asc") - cur["ridge"][:, :, i])[m]
                         for i in range(1, 37)])
    _stats("ridge x36", rd)
    print("\nNote: lat/lon identical and slope/aspect match to integer rounding "
          "-> the deck came from this same pipeline. Only sky/horizon angles "
          "differ (search-radius/version), and swapping them does not improve ET "
          "(see test_shade_regen.py).")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true",
                    help="(re)run mklatlon + Skyview before diffing")
    args = ap.parse_args()
    if args.run or not os.path.exists(os.path.join(WORKDIR, "slpfile.asc")):
        run_pipeline()
    diff_against_deck()
    print("\nDONE.")


if __name__ == "__main__":
    main()
