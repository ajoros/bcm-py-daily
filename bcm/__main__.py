"""python -m bcm run|check|import-ctl"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bcm_config import check_run, import_ctl, load_run  # noqa: E402


def cmd_run(args):
    problems = check_run(args.config)
    fatal = [p for p in problems if "terrain.inp" not in p]
    if fatal:
        print("Not ready:")
        for p in fatal:
            print(f"  {p}")
        sys.exit(1)
    os.environ["BCM_CONFIG"] = str(Path(args.config).resolve())
    from BCM_Dailyv81_python import main
    main()


def cmd_check(args):
    problems = check_run(args.config)
    cfg, indir, outdir = load_run(args.config)
    print(f"run     {Path(args.config).resolve()}")
    print(f"domain  {indir}")
    print(f"output  {outdir}")
    print(f"window  {cfg['yn1']}-{cfg['dn1']:03d} → {cfg['yn2']}-{cfg['dn2']:03d}")
    if problems:
        print("problems:")
        for p in problems:
            print(f"  {p}")
        sys.exit(1)
    print("ok")


def cmd_import(args):
    out = import_ctl(args.ctl, args.out, domain=args.domain)
    print(f"wrote {out}")
    print("tables: geology.csv vegetation.csv next to that file")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m bcm",
        description="BCM Daily in Python. No Fortran control file required.",
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="run a YAML config")
    p_run.add_argument("config")
    p_run.set_defaults(func=cmd_run)

    p_chk = sub.add_parser("check", help="verify inputs exist")
    p_chk.add_argument("config")
    p_chk.set_defaults(func=cmd_check)

    p_imp = sub.add_parser("import-ctl", help="write YAML from a Fortran .ctl")
    p_imp.add_argument("ctl")
    p_imp.add_argument("-o", "--out", default="run.yaml")
    p_imp.add_argument("--domain", default=None)
    p_imp.set_defaults(func=cmd_import)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
