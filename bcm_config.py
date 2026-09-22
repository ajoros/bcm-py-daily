"""Python run config. Replaces BCM_Dailyv81.ctl for new runs.

A run.yaml sets dates, flags, and paths. Geology / veg tables ship in
bcm_tables/. Existing Fortran control files still import.
"""
from __future__ import annotations

import csv
import os
from datetime import date, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLES = HERE / "bcm_tables"

# yaml key -> cfg key used by BCM_Dailyv81_python.py
LAYERS = {
    "dem": "demfile",
    "soil_depth": "soildfile",
    "wilting_point": "wpfile",
    "field_capacity": "fcfile",
    "porosity": "porfile",
    "geology": "geolfile",
    "vegetation": "vegfile",
    "bc_a": "bc_a_file",
    "bc_b": "bc_b_file",
    "bc_c": "bc_c_file",
    "pt_alpha": "pt_alpha_file",
    "snow_accum": "snowaccumfile",
    "mfmax": "mfmaxfile",
    "mfmin": "mfminfile",
    "basins": "areafile",
    "aridity": "aridityfile",
    "terrain": "inpfile",
}

DEFAULT_FILES = {k: (f"{k}.inp" if k == "terrain" else f"{k}.asc") for k in LAYERS}

PRINT_FLAGS = (
    "aet", "cld", "cwd", "evp", "exc", "mlt", "pet", "rad",
    "rch", "run", "sbl", "snw", "smd", "smr", "rvs",
)
DEFAULT_PRINT = ["pet", "aet", "cwd", "rch", "run", "snw", "mlt", "sbl", "pck", "str"]
# pck/str are state maps (always written). listed so the yaml reads naturally.


# ---------------------------------------------------------------------------
# tiny YAML (our files only — scalars, [lists], one-level maps)
# ---------------------------------------------------------------------------
def _parse_val(s: str):
    s = s.strip()
    if s in ("true", "True", "yes"):
        return True
    if s in ("false", "False", "no"):
        return False
    if s.startswith("[") and s.endswith("]"):
        inner = s[1:-1].strip()
        if not inner:
            return []
        return [_parse_val(x) for x in inner.split(",")]
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
        return s[1:-1]
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        return s


def load_simple_yaml(path) -> dict:
    """Read the small YAML subset this project writes. Not a general parser."""
    root: dict = {}
    stack = [(0, root)]
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip())
        line = raw.strip()
        while len(stack) > 1 and indent < stack[-1][0]:
            stack.pop()
        cur = stack[-1][1]
        if ":" not in line:
            raise ValueError(f"{path}: cannot parse: {line}")
        key, _, val = line.partition(":")
        key, val = key.strip(), val.strip()
        if " #" in val:
            val = val.split(" #", 1)[0].rstrip()
        if val == "":
            cur[key] = {}
            stack.append((indent + 1, cur[key]))
        else:
            cur[key] = _parse_val(val)
    return root


def dump_simple_yaml(data: dict, path) -> None:
    lines = []

    def emit(d, pad=""):
        for k, v in d.items():
            if isinstance(v, dict):
                lines.append(f"{pad}{k}:")
                emit(v, pad + "  ")
            elif isinstance(v, list):
                inner = ", ".join(str(x) for x in v)
                lines.append(f"{pad}{k}: [{inner}]")
            elif isinstance(v, bool):
                lines.append(f"{pad}{k}: {'true' if v else 'false'}")
            else:
                lines.append(f"{pad}{k}: {v}")

    emit(data)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
def _table_dir(yaml_dir: Path, domain: Path) -> Path:
    for d in (yaml_dir, domain, TABLES):
        if (d / "geology.csv").is_file() and (d / "vegetation.csv").is_file():
            return d
    raise FileNotFoundError(
        f"need geology.csv and vegetation.csv in {yaml_dir}, {domain}, or {TABLES}"
    )


def load_tables(table_dir: Path) -> tuple[dict, dict]:
    rockks = {}
    with open(table_dir / "geology.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rockks[int(row["id"])] = float(row["ks_mm_day"])
    veg = {}
    months = ["jan", "feb", "mar", "apr", "may", "jun",
              "jul", "aug", "sep", "oct", "nov", "dec"]
    with open(table_dir / "vegetation.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            veg[int(row["id"])] = {
                "init_lai": float(row["init_lai"]),
                "up_limit": float(row["up_limit"]),
                "dn_limit": float(row["dn_limit"]),
                "up_rate": float(row["up_rate"]),
                "dn_rate": float(row["dn_rate"]),
                "root_depth": float(row["root_depth"]),
                "kfactor": [float(row[f"k_{m}"]) for m in months],
            }
    return rockks, veg


def _as_date(val) -> date:
    if isinstance(val, date) and not isinstance(val, datetime):
        return val
    return date.fromisoformat(str(val))


def _doy(d: date) -> int:
    return d.timetuple().tm_yday


def _defaults() -> dict:
    """cfg keys the driver reads, with no-CTL defaults."""
    cfg = {
        "header": "bcm-py-daily",
        "outfile": "basin.out",
        "areatable": "",
        "hstep": 1,
        "std": -120.0,
        "prior": 0,
        "istor": 0.9,
        "soilksfile": "",
        "cwd_dry_file": "",
        "rockksfile": "",
        "aridityfile": "",
        "maskfile": "",
        "hourpetfile": "",
        "climate_dir_file": "",
        "climate_dir_flag": 0,
        "rockks_flag": 0,
        "drydown_flag": 0,
        "drydown_a": 0.65,
        "drydown_b": -10.0,
        "drydown_c": 0.05,
        "alphaprime_flag": 0,
        "rchrun_flag": 0,
        "rchrun_limit": 0.3,
        "urban_flag": 0,
        "urban_veg": 53,
        "urban_soild": 0.1,
        "urban_ks": 0.5,
        "rosflag": 0,
        "rvs_percent": 0.0,
        "snow_flag": 0,
        "snowaccum_t": 1.5,
        "maxmf": 1.8,
        "minmf": 0.1,
        "maf": 7.0,
        "tipm": 0.1,
        "nmf": 0.21,
        "solar_flag": 1,
        "hourly_pet_flag": 0,
        "rainfrac_flag": 0,
        "subl_flag": 0,
        "subl_par": 1.0,
        "mask_flag": 0,
        "ingest_flag": 0,
        "num_ingest": 0,
        "ingest_pairs": [],
        "rad_ave_files": [],
        "ppt_ave_files": [],
    }
    for name in PRINT_FLAGS:
        cfg[f"{name}_flag"] = 0
    return cfg


def _climate_folder(raw, yaml_path: Path, domain: Path) -> str:
    """off = ppt/tmn/tmx live in the domain folder. A path is the other folder."""
    if raw in (None, False, "off", "false", "False", 0, "0", ""):
        return str(domain)
    folder = Path(str(raw))
    if not folder.is_absolute():
        folder = (yaml_path.parent / folder).resolve()
    return str(folder)


def _apply_flags(doc: dict, cfg: dict) -> None:
    """Switches the Python run actually uses. Unknown keys are rejected."""
    flags = dict(doc.get("flags") or {})
    if "rch_run" not in flags and "rch_run_scaler" in doc:
        flags["rch_run"] = doc["rch_run_scaler"]
    known = {
        "rch_run", "solar", "snow", "snow_accum_temp",
        "mfmax", "mfmin", "maf", "tipm", "nmf",
        "drydown", "drydown_a", "drydown_b", "drydown_c",
        "pt_modified", "urban", "urban_veg", "urban_soil_depth", "urban_bedrock_k",
        "rain_fraction", "sublimation", "sublimation_constant", "ingest",
        "ingest_files",
    }
    unknown = set(flags) - known
    if unknown:
        raise ValueError(f"unknown flag {sorted(unknown)[0]}. want one of: {', '.join(sorted(known))}")
    if "rch_run" in flags:
        v = flags["rch_run"]
        if v in (False, "off", "false", "False", 0, "0"):
            cfg["rchrun_flag"] = 0
        else:
            cfg["rchrun_flag"] = 1
            cfg["rchrun_limit"] = float(v)
    if "solar" in flags:
        cfg["solar_flag"] = 1 if flags["solar"] in (True, "on", "yes", "true", "True", 1, "1") else 0
    if "snow" in flags:
        mode = str(flags["snow"]).lower()
        if mode == "maps":
            cfg["snow_flag"] = 0
        elif mode == "scalar":
            cfg["snow_flag"] = 1
        else:
            raise ValueError(f"snow must be maps or scalar, got {flags['snow']}")
    for yaml_key, cfg_key in (
        ("snow_accum_temp", "snowaccum_t"),
        ("mfmax", "maxmf"),
        ("mfmin", "minmf"),
        ("maf", "maf"),
        ("tipm", "tipm"),
        ("nmf", "nmf"),
    ):
        if yaml_key in flags:
            cfg[cfg_key] = float(flags[yaml_key])
    if "drydown" in flags:
        cfg["drydown_flag"] = 1 if flags["drydown"] in (
            True, "on", "yes", "true", "True", 1, "1") else 0
    for yaml_key, cfg_key in (
        ("drydown_a", "drydown_a"),
        ("drydown_b", "drydown_b"),
        ("drydown_c", "drydown_c"),
    ):
        if yaml_key in flags:
            cfg[cfg_key] = float(flags[yaml_key])
    if "pt_modified" in flags:
        cfg["alphaprime_flag"] = 1 if _flag_on(flags["pt_modified"]) else 0
    if "urban" in flags:
        cfg["urban_flag"] = 1 if _flag_on(flags["urban"]) else 0
    if "urban_veg" in flags:
        cfg["urban_veg"] = int(flags["urban_veg"])
    if "urban_soil_depth" in flags:
        cfg["urban_soild"] = float(flags["urban_soil_depth"])
    if "urban_bedrock_k" in flags:
        cfg["urban_ks"] = float(flags["urban_bedrock_k"])
    if "rain_fraction" in flags:
        cfg["rainfrac_flag"] = 1 if _flag_on(flags["rain_fraction"]) else 0
    if "sublimation" in flags:
        cfg["subl_flag"] = 1 if _flag_on(flags["sublimation"]) else 0
    if "sublimation_constant" in flags:
        cfg["subl_par"] = float(flags["sublimation_constant"])
    if "ingest" in flags:
        mode = str(flags["ingest"]).lower()
        try:
            cfg["ingest_flag"] = {"off": 0, "swe": 1, "str": 2, "lai": 3}[mode]
        except KeyError:
            raise ValueError(f"ingest must be off, swe, str, or lai, got {flags['ingest']}")
    if "ingest_files" in flags:
        pairs = flags["ingest_files"]
        if not isinstance(pairs, dict):
            raise ValueError("ingest_files must be file-to-replace: replacement-file")
        cfg["ingest_pairs"] = [(str(old), str(new)) for old, new in pairs.items()]
        cfg["num_ingest"] = len(cfg["ingest_pairs"])


def _flag_on(v) -> bool:
    return v in (True, "on", "yes", "true", "True", 1, "1")


def load_run(path) -> tuple[dict, str, str]:
    """Return (cfg, indir, outdir) from a run.yaml (or .yml)."""
    path = Path(path).resolve()
    doc = load_simple_yaml(path)
    start = _as_date(doc["start"])
    end = _as_date(doc["end"])
    if end < start:
        raise ValueError(f"end {end} is before start {start}")

    domain = Path(doc.get("domain", path.parent / "domain"))
    if not domain.is_absolute():
        domain = (path.parent / domain).resolve()
    outdir = Path(doc.get("output", path.parent / "out"))
    if not outdir.is_absolute():
        outdir = (path.parent / outdir).resolve()
    cfg = _defaults()
    cfg["climate_dir"] = _climate_folder(doc.get("climate", "off"), path, domain)
    cfg["yn1"], cfg["dn1"] = start.year, _doy(start)
    cfg["yn2"], cfg["dn2"] = end.year, _doy(end)
    cfg["prior"] = 1 if doc.get("antecedent") else 0
    if "initial_soil" in doc:
        cfg["istor"] = float(doc["initial_soil"])
    if "timezone" in doc:
        cfg["std"] = float(doc["timezone"])
    if "outfile" in doc:
        cfg["outfile"] = str(doc["outfile"])
    _apply_flags(doc, cfg)

    names = {k: DEFAULT_FILES[k] for k in LAYERS}
    for k, v in (doc.get("layers") or {}).items():
        if k not in LAYERS:
            raise ValueError(f"unknown layer '{k}'. want one of: {', '.join(LAYERS)}")
        names[k] = str(v)
    for yaml_key, cfg_key in LAYERS.items():
        cfg[cfg_key] = names[yaml_key]

    for name in PRINT_FLAGS:
        cfg[f"{name}_flag"] = 0
    for item in doc.get("print", DEFAULT_PRINT):
        key = str(item).strip().lower()
        if key in ("pck", "str", "ati", "hdi", "mwt", "lai"):
            continue
        if key not in PRINT_FLAGS:
            raise ValueError(f"unknown print flag '{item}'")
        cfg[f"{key}_flag"] = 1

    bbox = doc.get("bbox") or {}
    if bbox:
        cfg["bbox"] = {k: float(bbox[k]) for k in ("west", "east", "south", "north")}

    tdir = _table_dir(path.parent, domain)
    rockks, veg = load_tables(tdir)
    cfg["rockks_table"] = rockks
    cfg["veg_params"] = veg
    cfg["n_geol"] = len(rockks)
    cfg["n_veg"] = len(veg)
    cfg["_yaml"] = str(path)
    return cfg, str(domain), str(outdir)


def required_layer_files(cfg: dict) -> list[str]:
    skip = {"terrain"}
    if not cfg.get("drydown_flag"):
        skip.add("aridity")
    return [cfg[LAYERS[k]] for k in LAYERS if k not in skip]


def climate_days(cfg: dict):
    from BCM_Dailyv81_python import calendar_days
    return calendar_days(cfg)


def check_run(path) -> list[str]:
    """Return human-readable problems. Empty list = ready to run."""
    cfg, indir, outdir = load_run(path)
    problems = []
    root = Path(indir)
    if not root.is_dir():
        return [f"domain folder not found: {root}"]
    for name in required_layer_files(cfg):
        fp = root / name
        if not fp.is_file():
            problems.append(f"missing {fp}")
    from bcm_io import day_filename
    croot = Path(cfg.get("climate_dir") or indir)
    if not croot.is_dir():
        problems.append(f"climate folder not found: {croot}")
    for yn, dn, _leap in climate_days(cfg):
        for prefix in ("ppt", "tmn", "tmx"):
            fp = croot / day_filename(prefix, yn, dn)
            if not fp.is_file():
                problems.append(f"missing {fp} (climate for {yn} DOY {dn:03d})")
    if cfg.get("ingest_flag"):
        pairs = cfg.get("ingest_pairs") or []
        if not pairs:
            problems.append("ingest is swe, str, or lai and ingest_files is empty")
        for _replaced, replacement in pairs:
            fp = root / replacement
            if not fp.is_file():
                problems.append(f"missing {fp} (ingest replacement)")
    inp = root / cfg["inpfile"]
    if not inp.is_file() and not cfg.get("bbox"):
        problems.append(
            f"missing {inp.name} and no bbox: in the yaml "
            "(terrain.inp is built from the DEM when bbox is set)"
        )
    return problems


def import_ctl(ctl_path, out_yaml, domain=None) -> Path:
    """Write a run.yaml from a Fortran BCM_Dailyv81.ctl."""
    from bcm_io import parse_ctl

    ctl_path = Path(ctl_path).resolve()
    out_yaml = Path(out_yaml).resolve()
    cfg = parse_ctl(str(ctl_path))
    dom = Path(domain).resolve() if domain else ctl_path.parent
    try:
        domain_s = str(dom.relative_to(out_yaml.parent))
    except ValueError:
        domain_s = str(dom)

    start = date(int(cfg["yn1"]), 1, 1)
    from datetime import timedelta
    start = start + timedelta(days=int(cfg["dn1"]) - 1)
    end = date(int(cfg["yn2"]), 1, 1) + timedelta(days=int(cfg["dn2"]) - 1)

    layers = {}
    for yaml_key, cfg_key in LAYERS.items():
        layers[yaml_key] = cfg.get(cfg_key) or DEFAULT_FILES[yaml_key]

    printed = [n for n in PRINT_FLAGS if cfg.get(f"{n}_flag")]
    doc = {
        "start": start.isoformat(),
        "end": end.isoformat(),
        "antecedent": bool(cfg.get("prior")),
        "initial_soil": cfg.get("istor", 0.9),
        "domain": domain_s,
        "climate": cfg["climate_dir_path"] if cfg.get("climate_dir_flag") else domain_s,
        "output": "./out",
        "timezone": cfg.get("std", -120),
        "outfile": cfg.get("outfile", "basin.out"),
        "print": printed,
        "layers": layers,
        "flags": {
            "drydown": "on" if cfg.get("drydown_flag") else "off",
            "drydown_a": cfg.get("drydown_a", 0.65),
            "drydown_b": cfg.get("drydown_b", -10.0),
            "drydown_c": cfg.get("drydown_c", 0.05),
            "pt_modified": "on" if cfg.get("alphaprime_flag") else "off",
            "urban": "on" if cfg.get("urban_flag") else "off",
            "urban_veg": int(cfg.get("urban_veg", 53)),
            "urban_soil_depth": cfg.get("urban_soild", 0.1),
            "urban_bedrock_k": cfg.get("urban_ks", 0.5),
            "rain_fraction": "on" if cfg.get("rainfrac_flag") else "off",
            "sublimation": "on" if cfg.get("subl_flag") else "off",
            "sublimation_constant": cfg.get("subl_par", 1.0),
            "ingest": {0: "off", 1: "swe", 2: "str", 3: "lai"}.get(int(cfg.get("ingest_flag") or 0), "off"),
            "ingest_files": {old: new for old, new in (cfg.get("ingest_pairs") or [])},
        },
    }
    if cfg.get("rchrun_flag"):
        doc["rch_run_scaler"] = cfg.get("rchrun_limit", 0.3)

    out_yaml.parent.mkdir(parents=True, exist_ok=True)
    dump_simple_yaml(doc, out_yaml)

    # tables from THIS ctl, next to the yaml (domain may differ from shipped)
    gpath = out_yaml.parent / "geology.csv"
    vpath = out_yaml.parent / "vegetation.csv"
    if not gpath.exists():
        with open(gpath, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["id", "ks_mm_day", "name"])
            for gid, ks in cfg["rockks_table"].items():
                w.writerow([gid, ks, ""])
    if not vpath.exists():
        months = ["jan", "feb", "mar", "apr", "may", "jun",
                  "jul", "aug", "sep", "oct", "nov", "dec"]
        with open(vpath, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["id", "init_lai", "up_limit", "dn_limit", "up_rate",
                        "dn_rate", "root_depth", *[f"k_{m}" for m in months], "name"])
            for vid, v in cfg["veg_params"].items():
                w.writerow([vid, v["init_lai"], v["up_limit"], v["dn_limit"],
                            v["up_rate"], v["dn_rate"], v["root_depth"],
                            *v["kfactor"], ""])
    return out_yaml


def _check():
    rockks, veg = load_tables(TABLES)
    assert len(rockks) == 54, len(rockks)
    assert len(veg) == 62, len(veg)
    assert rockks[20] == 15.0
    assert veg[51]["root_depth"] == 0.5
    toy = HERE / "examples" / "toy" / "run.yaml"
    if toy.is_file():
        cfg, indir, _out = load_run(toy)
        assert cfg["yn1"] == 2010 and cfg["dn1"] == 1
        assert cfg["pet_flag"] == 1
        assert Path(indir).name == "domain"
    print("bcm_config check ok")


if __name__ == "__main__":
    _check()
