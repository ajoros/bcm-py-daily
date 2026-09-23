"""Compare mlt/sbl/exc Python vs Fortran (grids + basin .out)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, PKG)
from validate import compare_out, compare_grids

PY = os.environ.get("BCM_PY_DIR", os.path.join(os.environ.get("TEMP", "/tmp"), "bcm_mlt_sbl_exc"))
REF = os.environ.get("BCM_FORT_DIR", os.path.normpath(os.path.join(PKG, "..", "BCM_testrun_original")))


def main():
    print(f"py  = {PY}")
    print(f"ref = {REF}")
    for prefix in ("mlt", "sbl", "exc"):
        py_n = sum(os.path.exists(os.path.join(PY, f"{prefix}2010_{d:03d}.asc")) for d in range(1, 10))
        rf_n = sum(os.path.exists(os.path.join(REF, f"{prefix}2010_{d:03d}.asc")) for d in range(1, 10))
        print(f"  {prefix} grids: python {py_n}/9  fortran {rf_n}/9")

    py_out = os.path.join(PY, "mok_Dailyv81_test.out")
    rf_out = os.path.join(REF, "mok_Dailyv81_test.out")
    if os.path.exists(py_out) and os.path.exists(rf_out):
        r = compare_out(py_out, rf_out, cols=["mltmm", "sblmm", "excmm", "snwmm", "pckmm", "aetmm"])
        print(f"\nBasin .out ({r['keys']} shared basin-days)")
        print(f"{'col':8} {'exe_mean':>10} {'bias':>10} {'rmse':>10} {'maxabs':>10}")
        for c, s in r["cols"].items():
            print(f"{c:8} {s['exe_mean']:10.4f} {s['bias']:+10.4f} {s['rmse']:10.4f} {s['mae']:10.4f}")

    stats = compare_grids(PY, REF, 2010, list(range(1, 10)),
                          variables=["mlt", "sbl", "exc", "snw", "pck", "aet"],
                          tol=0.01)
    if stats:
        print("\nGrid cell-by-cell (days 1-9)")
        print(f"{'var':5} {'n':>8} {'exe_mean':>10} {'py_mean':>10} {'bias':>10} {'rmse':>10} {'maxabs':>10} {'%match':>8}")
        for v, s in stats.items():
            if not s.get("n"):
                print(f"{v:5} (no overlap)")
                continue
            print(f"{v:5} {s['n']:8} {s['exe_mean']:10.4f} {s['py_mean']:10.4f} "
                  f"{s['bias']:+10.4f} {s['rmse']:10.4f} {s['mae']:10.3f} {s['pct_match']:7.1f}%")
    else:
        print("\nNo grid overlap yet (Fortran maps missing until the exe is re-run).")


if __name__ == "__main__":
    main()
