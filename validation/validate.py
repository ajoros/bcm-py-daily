"""
validate.py — BCM daily Python-vs-Fortran(.exe) validation harness.

Quantifies how closely the Python port reproduces the compiled .exe reference
outputs, cell-by-cell and at the basin-summary level, and checks the
water-balance closure invariant. Designed so `run_scenarios.py` can import
`compare_grids()` / `compare_out()` to sweep the calibration open-items
(see Calibration_OpenItems_Memo_2026-07-09.md).

Usage:
    python validate.py                 # run the model fresh, then compare
    python validate.py --no-run        # compare existing outputs in --py-dir
    python validate.py --py-dir DIR --ref-dir DIR --days 1-9

Reference (.exe truth) defaults to ../BCM_testrun_original, which ships the
aet/cwd/pet/rch/run/str/... 2010_00N.asc grids and mok_Dailyv81_test.out.
"""
from __future__ import annotations
import os
import sys
import argparse
import subprocess
import tempfile
import datetime as _dt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, PKG)
from bcm_io import read_asc  # noqa: E402

NODATA = -9999.0
MODEL  = os.path.join(PKG, "BCM_Dailyv81_python.py")
DEFAULT_REF = os.path.normpath(os.path.join(PKG, "..", "BCM_testrun_original"))

# Variables worth grading. The calibration open-items live in pet/str/rch/run.
GRID_VARS = ["pet", "aet", "cwd", "rch", "run", "str", "snw", "rvs",
             "smd", "smr", "mlt", "sbl", "exc"]
KEY_VARS  = ["pet", "str", "rch", "run", "aet", "mlt", "sbl", "exc"]


# ---------------------------------------------------------------------------
# Grid comparison
# ---------------------------------------------------------------------------
def _stats(py: np.ndarray, ex: np.ndarray, tol: float) -> dict | None:
    """Per-grid stats over cells valid (non-nodata) in BOTH grids."""
    if py.shape != ex.shape:
        return {"error": f"shape {py.shape} vs {ex.shape}"}
    both = (py != NODATA) & (ex != NODATA)
    n = int(both.sum())
    # nodata agreement: fraction of cells whose nodata status matches
    nd_match = float(((py == NODATA) == (ex == NODATA)).mean()) * 100.0
    if n == 0:
        return {"n": 0, "nodata_match": nd_match}
    p = py[both]; e = ex[both]
    d = p - e
    r = (np.corrcoef(p, e)[0, 1] if n > 1 and p.std() > 0 and e.std() > 0
         else float("nan"))
    return {
        "n": n,
        "exe_mean": float(e.mean()),
        "py_mean":  float(p.mean()),
        "bias":     float(d.mean()),          # py − exe
        "rmse":     float(np.sqrt((d**2).mean())),
        "mae":      float(np.abs(d).max()),
        "pct_match": float((np.abs(d) < tol).mean() * 100.0),
        "corr":     float(r),
        "nodata_match": nd_match,
    }


def compare_grids(py_dir: str, ref_dir: str, year: int, days,
                  variables=None, tol: float = 0.01) -> dict:
    """Aggregate per-variable stats across the requested days.
    Returns {var: {aggregated stats, 'days': [...per-day...]}}."""
    variables = variables or GRID_VARS
    out = {}
    for var in variables:
        per_day = []
        for d in days:
            fn = f"{var}{year}_{d:03d}.asc"
            pp = os.path.join(py_dir, fn)
            rp = os.path.join(ref_dir, fn)
            if not (os.path.exists(pp) and os.path.exists(rp)):
                continue
            _, pa = read_asc(pp)
            _, ea = read_asc(rp)
            s = _stats(pa, ea, tol)
            if s:
                s["day"] = d
                per_day.append(s)
        if not per_day:
            continue
        good = [s for s in per_day if s.get("n", 0) > 0]
        if not good:
            out[var] = {"n": 0, "days": per_day}
            continue
        # n-weighted aggregate
        ntot = sum(s["n"] for s in good)
        agg = {
            "n": ntot,
            "exe_mean": sum(s["exe_mean"]*s["n"] for s in good)/ntot,
            "py_mean":  sum(s["py_mean"] *s["n"] for s in good)/ntot,
            "bias":     sum(s["bias"]    *s["n"] for s in good)/ntot,
            "rmse":     float(np.sqrt(sum((s["rmse"]**2)*s["n"] for s in good)/ntot)),
            "mae":      max(s["mae"] for s in good),
            "pct_match": sum(s["pct_match"]*s["n"] for s in good)/ntot,
            "corr":     float(np.nanmean([s["corr"] for s in good])),
            "nodata_match": float(np.mean([s["nodata_match"] for s in per_day])),
            "days":     per_day,
        }
        out[var] = agg
    return out


# ---------------------------------------------------------------------------
# Basin summary (.out) comparison
# ---------------------------------------------------------------------------
def parse_out(path: str) -> dict:
    """Parse a BCM basin-summary .out into {(year,doy,basin): {col: val}}.
    Aligns by column header names so Python/exe layouts can differ."""
    with open(path, "r", encoding="latin-1") as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    if not lines:
        return {}
    hdr = lines[0].split()
    # locate the year/doy/basin key columns (case-insensitive, tolerant)
    low = [h.lower() for h in hdr]
    def idx(*names):
        for nm in names:
            if nm in low:
                return low.index(nm)
        return None
    iy, idoy, ib = idx("year"), idx("doy"), idx("basin")
    rows = {}
    for ln in lines[1:]:
        parts = ln.split()
        if len(parts) < len(hdr):
            continue
        try:
            key = (int(float(parts[iy])), int(float(parts[idoy])),
                   int(float(parts[ib])))
        except (ValueError, TypeError):
            continue
        rec = {}
        for j, name in enumerate(hdr):
            try:
                rec[name.lower()] = float(parts[j])
            except ValueError:
                pass
        rows[key] = rec
    return rows


def compare_out(py_out: str, ref_out: str, cols=None) -> dict:
    """Compare shared numeric columns of two .out files over shared keys."""
    py = parse_out(py_out)
    ex = parse_out(ref_out)
    keys = sorted(set(py) & set(ex))
    if not keys:
        return {"keys": 0}
    all_cols = set(py[keys[0]]) & set(ex[keys[0]])
    cols = [c for c in (cols or all_cols) if c in all_cols]
    res = {"keys": len(keys), "cols": {}}
    for c in cols:
        d = np.array([py[k][c] - ex[k][c] for k in keys
                      if c in py[k] and c in ex[k]])
        e = np.array([ex[k][c] for k in keys if c in ex[k]])
        if d.size:
            res["cols"][c] = {
                "bias": float(d.mean()), "rmse": float(np.sqrt((d**2).mean())),
                "mae": float(np.abs(d).max()), "exe_mean": float(e.mean()),
            }
    return res


def water_balance_closure(py_out: str) -> dict:
    """Check ppt − snw + mlt − aet − rch − run − Δstr ≈ 0 per basin/day from
    the Python .out (needs consecutive days for Δstr, so day 1 is skipped)."""
    rows = parse_out(py_out)
    by_basin = {}
    for (yn, dn, b), rec in rows.items():
        by_basin.setdefault(b, []).append((yn, dn, rec))
    resid = []
    for b, seq in by_basin.items():
        seq.sort(key=lambda t: (t[0], t[1]))
        for i in range(1, len(seq)):
            _, _, r  = seq[i]
            _, _, rp = seq[i-1]
            need = ("pptmm", "snwmm", "mltmm", "aetmm", "rchmm", "runmm", "strmm")
            if not all(k in r for k in need) or "strmm" not in rp:
                continue
            dstr = r["strmm"] - rp["strmm"]
            res = (r["pptmm"] - r["snwmm"] + r["mltmm"] - r["aetmm"]
                   - r["rchmm"] - r["runmm"] - dstr)
            resid.append(res)
    if not resid:
        return {"n": 0}
    a = np.array(resid)
    return {"n": int(a.size), "max_abs": float(np.abs(a).max()),
            "rms": float(np.sqrt((a**2).mean())), "mean": float(a.mean())}


# ---------------------------------------------------------------------------
# Model runner
# ---------------------------------------------------------------------------
def run_model(out_dir: str, env_extra: dict | None = None) -> None:
    os.makedirs(out_dir, exist_ok=True)
    env = os.environ.copy()
    env["BCM_OUT_DIR"] = out_dir
    env["PYTHONIOENCODING"] = "utf-8"
    if env_extra:
        env.update({k: str(v) for k, v in env_extra.items()})
    subprocess.run([sys.executable, MODEL], cwd=PKG, env=env,
                   check=True, stdout=subprocess.DEVNULL)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _fmt_grid_table(stats: dict) -> str:
    hdr = (f"| {'var':<5} | {'n':>7} | {'exe_mean':>9} | {'py_mean':>9} | "
           f"{'bias':>8} | {'rmse':>8} | {'maxabs':>8} | {'%match':>7} | "
           f"{'corr':>6} | {'nd%':>5} |")
    sep = "|" + "|".join(["-"*7, "-"*9, "-"*11, "-"*11, "-"*10, "-"*10,
                          "-"*10, "-"*9, "-"*8, "-"*7]) + "|"
    out = [hdr, sep]
    for var, s in stats.items():
        if s.get("n", 0) == 0:
            out.append(f"| {var:<5} | {'0':>7} | (no valid overlap) |")
            continue
        out.append(f"| {var:<5} | {s['n']:>7} | {s['exe_mean']:>9.3f} | "
                   f"{s['py_mean']:>9.3f} | {s['bias']:>8.3f} | {s['rmse']:>8.3f} | "
                   f"{s['mae']:>8.2f} | {s['pct_match']:>6.1f}% | "
                   f"{s['corr']:>6.3f} | {s['nodata_match']:>4.0f}% |")
    return "\n".join(out)


def write_report(report_dir: str, grid_stats: dict, out_cmp: dict,
                 wb: dict, meta: dict) -> str:
    os.makedirs(report_dir, exist_ok=True)
    md = os.path.join(report_dir, "validation_report.md")
    lines = [
        "# BCM Daily — Validation Report",
        "",
        f"- Generated: {meta['when']}",
        f"- Python outputs: `{meta['py_dir']}`",
        f"- Reference (.exe): `{meta['ref_dir']}`",
        f"- Year {meta['year']}, days {meta['days'][0]}–{meta['days'][-1]}, "
        f"match tolerance {meta['tol']} mm",
        "",
        "## Grid comparison (Python − exe, valid cells in both)",
        "*bias = mean(py−exe); rmse n-weighted across days; nd% = nodata-mask "
        "agreement.*",
        "",
        _fmt_grid_table(grid_stats),
        "",
        "## Basin summary (.out) comparison",
    ]
    if out_cmp.get("keys"):
        lines.append(f"Shared basin×day rows: {out_cmp['keys']}")
        lines.append("")
        lines.append(f"| {'col':<8} | {'exe_mean':>9} | {'bias':>9} | {'rmse':>9} | {'maxabs':>9} |")
        lines.append("|" + "|".join(["-"*10]*5) + "|")
        for c, s in sorted(out_cmp["cols"].items()):
            lines.append(f"| {c:<8} | {s['exe_mean']:>9.3f} | {s['bias']:>9.3f} | "
                         f"{s['rmse']:>9.3f} | {s['mae']:>9.3f} |")
    else:
        lines.append("_No shared rows (Python .out vs exe .out) found._")
    lines += [
        "",
        "## Water-balance closure (Python .out, days 2+)",
        (f"- residuals: n={wb['n']}, max|res|={wb['max_abs']:.4g} mm, "
         f"rms={wb['rms']:.4g} mm, mean={wb['mean']:.4g} mm"
         if wb.get("n") else "- not computable (missing columns)"),
        "",
        "## Key-variable headline (RMSE vs exe)",
    ]
    for v in KEY_VARS:
        if v in grid_stats and grid_stats[v].get("n"):
            s = grid_stats[v]
            lines.append(f"- **{v}**: rmse={s['rmse']:.3f} mm, bias={s['bias']:+.3f} mm, "
                         f"corr={s['corr']:.3f}")
    md_text = "\n".join(lines) + "\n"
    with open(md, "w", encoding="utf-8") as f:
        f.write(md_text)
    # CSV
    import csv
    csv_path = os.path.join(report_dir, "grid_stats.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["var", "n", "exe_mean", "py_mean", "bias", "rmse",
                    "maxabs", "pct_match", "corr", "nodata_match"])
        for var, s in grid_stats.items():
            if s.get("n", 0) == 0:
                continue
            w.writerow([var, s["n"], f"{s['exe_mean']:.5f}", f"{s['py_mean']:.5f}",
                        f"{s['bias']:.5f}", f"{s['rmse']:.5f}", f"{s['mae']:.5f}",
                        f"{s['pct_match']:.2f}", f"{s['corr']:.5f}",
                        f"{s['nodata_match']:.2f}"])
    return md


def _parse_days(spec: str):
    if "-" in spec:
        a, b = spec.split("-"); return list(range(int(a), int(b)+1))
    return [int(x) for x in spec.split(",")]


def main():
    ap = argparse.ArgumentParser(description="BCM Python-vs-exe validation")
    ap.add_argument("--py-dir", default=None,
                    help="dir with Python outputs (default: run fresh to temp)")
    ap.add_argument("--ref-dir", default=DEFAULT_REF)
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--days", default="1-9")
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--no-run", action="store_true",
                    help="don't run the model; compare existing --py-dir")
    ap.add_argument("--report-dir", default=None)
    args = ap.parse_args()

    days = _parse_days(args.days)
    py_dir = args.py_dir
    if not args.no_run:
        py_dir = py_dir or tempfile.mkdtemp(prefix="bcm_val_")
        print(f"Running model → {py_dir}")
        run_model(py_dir)
    if not py_dir:
        sys.exit("--py-dir required with --no-run")

    print(f"Comparing\n  py  = {py_dir}\n  ref = {args.ref_dir}")
    grid_stats = compare_grids(py_dir, args.ref_dir, args.year, days, tol=args.tol)
    py_out  = os.path.join(py_dir, "mok_Dailyv81_test.out")
    ref_out = os.path.join(args.ref_dir, "mok_Dailyv81_test.out")
    out_cmp = (compare_out(py_out, ref_out)
               if os.path.exists(py_out) and os.path.exists(ref_out)
               else {"keys": 0})
    wb = water_balance_closure(py_out) if os.path.exists(py_out) else {"n": 0}

    report_dir = args.report_dir or os.path.join(
        PKG, "validation_reports",
        _dt.datetime.now().strftime("%Y%m%d_%H%M%S"))
    meta = {"when": _dt.datetime.now().isoformat(timespec="seconds"),
            "py_dir": py_dir, "ref_dir": args.ref_dir, "year": args.year,
            "days": days, "tol": args.tol}
    md = write_report(report_dir, grid_stats, out_cmp, wb, meta)

    # console summary
    print("\n" + _fmt_grid_table(grid_stats))
    if wb.get("n"):
        print(f"\nWater-balance closure: max|res|={wb['max_abs']:.4g} mm "
              f"rms={wb['rms']:.4g} mm (n={wb['n']})")
    print(f"\nReport: {md}")


if __name__ == "__main__":
    main()
