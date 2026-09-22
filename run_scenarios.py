"""
run_scenarios.py — sweep the calibration open-items and quantify each one's
effect on the Python-vs-exe error.

For each scenario it runs the model (into its own temp dir, via env overrides)
then compares KEY variables to the .exe reference. Prints a table and writes a
markdown summary. See Calibration_OpenItems_Memo_2026-07-09.md.

Knobs (env vars, all default-preserving):
    BCM_PT_ALPHA        Priestley-Taylor alpha (default 0.95)
    BCM_SOIL_AUG        'full' (default, soild+root) | 'none' (raw soild)
    BCM_RCHRUN_DISABLE  '1' to force the RCH/RUN scaler off

Usage:
    python run_scenarios.py                 # default sweep
    python run_scenarios.py --days 1-9
"""
from __future__ import annotations
import os
import sys
import argparse
import tempfile
import datetime as _dt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import validate as V  # noqa: E402

# (label, env-overrides) — one baseline + one change per open-item.
SCENARIOS = [
    ("baseline (a=0.95, soil=full, rchrun=on)", {}),
    ("PET alpha=0.85",   {"BCM_PT_ALPHA": "0.85"}),
    ("PET alpha=1.05",   {"BCM_PT_ALPHA": "1.05"}),
    ("PET alpha=1.26",   {"BCM_PT_ALPHA": "1.26"}),
    ("soil_aug=none",    {"BCM_SOIL_AUG": "none"}),
    ("rchrun disabled",  {"BCM_RCHRUN_DISABLE": "1"}),
]
KEY = V.KEY_VARS + ["cwd"]


def _parse_days(spec: str):
    if "-" in spec:
        a, b = spec.split("-"); return list(range(int(a), int(b)+1))
    return [int(x) for x in spec.split(",")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--days", default="1-9")
    ap.add_argument("--ref-dir", default=V.DEFAULT_REF)
    ap.add_argument("--year", type=int, default=2010)
    args = ap.parse_args()
    days = _parse_days(args.days)

    results = []   # (label, {var: stats})
    for label, env in SCENARIOS:
        d = tempfile.mkdtemp(prefix="bcm_scn_")
        print(f"[run] {label}  → {d}")
        V.run_model(d, env_extra=env)
        gs = V.compare_grids(d, args.ref_dir, args.year, days, variables=KEY)
        wb = V.water_balance_closure(os.path.join(d, "mok_Dailyv81_test.out"))
        results.append((label, gs, wb))

    # ── table: RMSE (and bias) per key var per scenario ──────────────────────
    def cell(gs, var, field):
        s = gs.get(var)
        return f"{s[field]:+.3f}" if s and s.get("n") else "   -  "

    lines = []
    lines.append(f"Sweep over year {args.year}, days {days[0]}–{days[-1]}, "
                 f"ref={args.ref_dir}")
    lines.append("")
    head = f"| {'scenario':<42} | " + " | ".join(f"{v:>8}" for v in KEY) + " | wb_rms |"
    lines.append(head)
    lines.append("|" + "-"*44 + "|" + "|".join(["-"*10]*len(KEY)) + "|" + "-"*8 + "|")
    lines.append("| " + " "*42 + " | " + " | ".join(" RMSE(mm)" for _ in KEY) + " |        |")
    for label, gs, wb in results:
        row = f"| {label:<42} | " + " | ".join(cell(gs, v, "rmse") for v in KEY)
        row += f" | {wb.get('rms', float('nan')):>6.3f} |" if wb.get("n") else " |   -    |"
        lines.append(row)
    lines.append("")
    lines.append("Bias (py − exe), same layout:")
    lines.append(head)
    lines.append("|" + "-"*44 + "|" + "|".join(["-"*10]*len(KEY)) + "|" + "-"*8 + "|")
    for label, gs, wb in results:
        row = f"| {label:<42} | " + " | ".join(cell(gs, v, "bias") for v in KEY) + " |        |"
        lines.append(row)

    table = "\n".join(lines)
    print("\n" + table + "\n")

    report_dir = os.path.join(HERE, "validation_reports",
                              "scenarios_" + _dt.datetime.now().strftime("%Y%m%d_%H%M%S"))
    os.makedirs(report_dir, exist_ok=True)
    md = os.path.join(report_dir, "scenario_sweep.md")
    with open(md, "w", encoding="utf-8") as f:
        f.write("# BCM Daily — Calibration scenario sweep\n\n")
        f.write(f"Generated {_dt.datetime.now().isoformat(timespec='seconds')}\n\n")
        f.write("Lower RMSE = closer to the exe. Use to quantify the open-items "
                "(PET alpha, soil augmentation, RCH/RUN scaler).\n\n")
        f.write("```\n" + table + "\n```\n")
    print(f"Report: {md}")


if __name__ == "__main__":
    main()
