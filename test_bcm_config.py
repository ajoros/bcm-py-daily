"""Fails if the no-CTL config path breaks."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from bcm_config import load_run, load_simple_yaml, load_tables, TABLES

rockks, veg = load_tables(TABLES)
assert len(rockks) == 54 and len(veg) == 62
assert rockks[20] == 15.0
assert abs(veg[51]["kfactor"][0] - 0.206) < 1e-6

from tempfile import TemporaryDirectory

sample = """
start: 2010-01-01
end: 2010-01-03
antecedent: false
timezone: -120   # pacific meridian
domain: ./domain
output: ./out
print: [pet, aet]
"""
with TemporaryDirectory() as td:
    yp = Path(td) / "run.yaml"
    yp.write_text(sample, encoding="utf-8")
    doc = load_simple_yaml(yp)
    assert doc["start"] == "2010-01-01"
    assert "pet" in doc["print"]
    cfg, indir, outdir = load_run(yp)
assert cfg["yn1"] == 2010 and cfg["dn1"] == 1 and cfg["dn2"] == 3
assert cfg["std"] == -120.0
assert cfg["prior"] == 0
assert cfg["rchrun_flag"] == 0
assert cfg["pet_flag"] == 1 and cfg["aet_flag"] == 1
assert cfg["demfile"] == "dem.asc"
assert Path(indir).name == "domain"
assert Path(outdir).name == "out"

shipped, shipped_dir, _ = load_run(HERE / "run.yaml")[:3]
assert shipped["rchrun_flag"] == 0 and shipped["solar_flag"] == 1
assert shipped["snow_flag"] == 0
assert shipped["maf"] == 7.0 and shipped["tipm"] == 0.1 and shipped["nmf"] == 0.21
assert shipped["maxmf"] == 1.8 and shipped["minmf"] == 0.1
assert shipped["snowaccum_t"] == 1.5
assert shipped["drydown_flag"] == 0
assert shipped["drydown_a"] == 0.65
assert shipped["drydown_b"] == -10.0
assert shipped["drydown_c"] == 0.05
assert shipped["aridityfile"] == "aridity.asc"
assert shipped["alphaprime_flag"] == 0 and shipped["urban_flag"] == 0
assert shipped["urban_veg"] == 53 and shipped["urban_soild"] == 0.1 and shipped["urban_ks"] == 0.5
assert shipped["rainfrac_flag"] == 0 and shipped["subl_flag"] == 0 and shipped["subl_par"] == 1.0
assert shipped["ingest_flag"] == 0
assert shipped["ingest_pairs"] == [
    ("pck2019mar.asc", "aso2019mar.asc"),
    ("pck2019apr.asc", "aso2019apr.asc"),
]
assert Path(shipped["climate_dir"]) == Path(shipped_dir)

from bcm_config import check_run
assert not any("aridity" in p for p in check_run(HERE / "run.yaml"))
assert not any("ingest replacement" in p for p in check_run(HERE / "run.yaml"))

with TemporaryDirectory() as td:
    yp = Path(td) / "run.yaml"
    yp.write_text(sample.replace("domain: ./domain", "domain: ./domain\nclimate: ./climate"), encoding="utf-8")
    cfg, indir, _ = load_run(yp)
    assert Path(cfg["climate_dir"]).name == "climate"
    assert Path(cfg["climate_dir"]) != Path(indir)

# Same driver knobs as parse_ctl when the Fortran CTL is on disk (local, not GitHub).
ctl = HERE.parent / "BCM_testrun_original" / "BCM_Dailyv81.ctl"
if ctl.is_file():
    from tempfile import TemporaryDirectory
    from bcm_config import import_ctl
    from bcm_io import parse_ctl
    with TemporaryDirectory() as td:
        yp = Path(td) / "run.yaml"
        import_ctl(ctl, yp, domain=ctl.parent)
        ycfg, _, _ = load_run(yp)
    ccfg = parse_ctl(str(ctl))
    keys = [
        "dn1", "yn1", "dn2", "yn2", "hstep", "std", "prior", "istor",
        "demfile", "soildfile", "wpfile", "fcfile", "porfile",
        "geolfile", "vegfile", "bc_a_file", "bc_b_file", "bc_c_file",
        "pt_alpha_file", "snowaccumfile", "mfmaxfile", "mfminfile",
        "areafile", "inpfile", "aridityfile", "rchrun_flag", "rchrun_limit",
        "maf", "tipm", "nmf", "solar_flag", "outfile",
        "drydown_flag", "drydown_a", "drydown_b", "drydown_c",
        "alphaprime_flag", "urban_flag", "urban_veg", "urban_soild", "urban_ks",
        "rainfrac_flag", "subl_flag", "subl_par", "ingest_flag", "ingest_pairs",
        "aet_flag", "cwd_flag", "exc_flag", "mlt_flag", "pet_flag",
        "rad_flag", "rch_flag", "run_flag", "sbl_flag", "snw_flag",
    ]
    for k in keys:
        assert ccfg[k] == ycfg[k], (k, ccfg[k], ycfg[k])
    assert ccfg["rockks_table"] == ycfg["rockks_table"]
    assert ccfg["veg_params"] == ycfg["veg_params"]

# v86 meridian: a cell off -120° must not match a run whose STD equals that longitude.
import numpy as np
from BCM_Dailyv81_python import compute_solar_pet, soil_drydown
import math

arid = np.array([[0.2, -9999.0]])
got = soil_drydown(arid, 0.65, -10.0, 0.05, -9999.0)
assert abs(float(got[0, 0]) - (0.65 * math.exp(-10.0 * 0.2) + 0.05)) < 1e-12
assert float(got[0, 1]) == 0.0

def _pet(std, lon):
    z = np.zeros((1, 1))
    one = np.ones((1, 1))
    terrain = {
        "lat": one * 38.4, "lon": one * lon, "sl": z, "asp": z,
        "elev": one * 500.0, "sky": one, "ridge": np.zeros((1, 1, 37)),
    }
    _, _, pet, _ = compute_solar_pet(
        15, 1, {"std": std}, terrain,
        one * 12.0, one * 2.0, one * 0.7, one * 0.01, one * 2.4,
        one, -9999.0, one, one * 0.1, one * 0.2, one * 0.3,
    )
    return float(pet[0, 0])

assert _pet(-120.0, -120.5) != _pet(-120.5, -120.5)

print("test_bcm_config ok")
