"""
Regression test for water / no-geology masking in the BCM daily water balance.

Guards the "Andy TODO": water cells must NOT emit spurious runoff.

The Fortran (BCM_Dailyv81.f90) masks water two independent ways:
  (1) bedrock Ks == 0  → geology "Water" (id 58) / "No Geology" (id 54):
      soild forced to 0 (L1929-1931) and, since maxrch = geolks = 0,
      run/rch/aet forced to 0 (L2067-2072).
  (2) WHR vegetation "Water" (id 57): every soil/snow term zeroed (L2073-2085).

Run directly:   python test_water_masking.py
Or with pytest: pytest test_water_masking.py
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from BCM_Dailyv81_python import water_balance  # noqa: E402

NODATA = -9999.0


def _base_inputs(ncell: int):
    """Land cell with plenty of incoming water so it WOULD generate runoff,
    recharge, and AET unless masked. Returns a dict of (1, ncell) arrays."""
    z = lambda v: np.full((1, ncell), float(v))  # noqa: E731
    return dict(
        ppt=z(500.0), melt=z(0.0), snow=z(0.0), petday=z(5.0),
        prior_str=z(200.0),
        soild_eff=z(1.0),               # 1 m soil
        wpmm=z(100.0), fcmm=z(300.0), pormm=z(400.0),  # wp/fc/por * 1 m * 1000
        geolks=z(10.0),                 # finite bedrock Ks
        kfac=z(0.5),
    )


def _call(inp, **kw):
    return water_balance(
        inp["ppt"], inp["melt"], inp["snow"], inp["petday"], inp["prior_str"],
        inp["soild_eff"], inp["wpmm"], inp["fcmm"], inp["pormm"], inp["geolks"],
        inp["kfac"], NODATA, rchrun_flag=0, rchrun_limit=0.3, **kw,
    )


def test_geology_water_is_fully_masked():
    """geolks==0 with the Fortran's soild=0 → run/rch/aet/str/smd/smr all 0."""
    inp = _base_inputs(1)
    inp["geolks"] = np.zeros((1, 1))     # water / no-geology
    inp["soild_eff"] = np.zeros((1, 1))  # main() forces soild=0 when geolks==0
    inp["wpmm"] = np.zeros((1, 1))
    inp["fcmm"] = np.zeros((1, 1))
    inp["pormm"] = np.zeros((1, 1))
    zero_ks = np.array([[True]])
    wb = _call(inp, zero_ks=zero_ks)
    for key in ("run2day", "rch3day", "vegaet", "storday", "smdday", "smrday"):
        assert np.allclose(wb[key], 0.0), f"{key} should be 0 over water, got {wb[key]}"


def test_geology_water_without_mask_WOULD_spill_runoff():
    """Demonstrates the bug the mask fixes: valid soil over Ks==0 water,
    with NO masking, produces spurious runoff (rejected recharge → runoff)."""
    inp = _base_inputs(1)
    inp["geolks"] = np.zeros((1, 1))     # Ks == 0
    # soil_eff kept > 0 (the un-guarded case) and no zero_ks mask passed
    wb = _call(inp, zero_ks=None)
    assert wb["run2day"][0, 0] > 0.0, (
        "Without the water mask, Ks==0 cells spill rejected recharge into "
        "runoff — this is the spurious runoff the fix must eliminate."
    )
    assert np.isclose(wb["rch3day"][0, 0], 0.0), "recharge is capped at Ks=0"


def test_vegetation_water_is_fully_masked():
    """WHR veg id 57 (Water) → all soil terms 0 even with finite geology."""
    inp = _base_inputs(1)
    water_veg = np.array([[True]])
    wb = _call(inp, water_veg=water_veg)
    for key in ("run2day", "rch3day", "vegaet", "storday",
                "cwd1day", "smdday", "smrday"):
        assert np.allclose(wb[key], 0.0), f"{key} should be 0 over veg-water, got {wb[key]}"


def test_normal_land_cell_unchanged():
    """A normal land cell still produces bounded recharge, some runoff, and AET."""
    inp = _base_inputs(1)
    zero_ks = np.array([[False]])
    water_veg = np.array([[False]])
    wb = _call(inp, zero_ks=zero_ks, water_veg=water_veg)
    assert wb["rch3day"][0, 0] <= inp["geolks"][0, 0] + 1e-9, "recharge must be capped at bedrock Ks"
    assert wb["run2day"][0, 0] > 0.0, "excess water should run off on a saturated land cell"
    assert wb["vegaet"][0, 0] > 0.0, "land cell should transpire"


def test_mixed_grid_only_water_cells_masked():
    """4-cell grid: land / geology-water / veg-water / land — only water masked."""
    inp = _base_inputs(4)
    inp["geolks"] = np.array([[10.0, 0.0, 10.0, 10.0]])
    inp["soild_eff"] = np.array([[1.0, 0.0, 1.0, 1.0]])  # soild=0 where Ks=0
    inp["wpmm"] = np.array([[100.0, 0.0, 100.0, 100.0]])
    inp["fcmm"] = np.array([[300.0, 0.0, 300.0, 300.0]])
    inp["pormm"] = np.array([[400.0, 0.0, 400.0, 400.0]])
    zero_ks = np.array([[False, True, False, False]])
    water_veg = np.array([[False, False, True, False]])
    wb = _call(inp, zero_ks=zero_ks, water_veg=water_veg)
    # water cells (index 1, 2) fully zeroed
    for idx in (1, 2):
        for key in ("run2day", "rch3day", "vegaet", "storday", "smdday", "smrday"):
            assert np.isclose(wb[key][0, idx], 0.0), f"cell {idx} {key} not masked"
    # land cells (index 0, 3) active
    for idx in (0, 3):
        assert wb["run2day"][0, idx] > 0.0, f"land cell {idx} should have runoff"


def _run_all():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"  FAIL  {t.__name__}: {e}")
    print("-" * 60)
    print(f"{len(tests) - failures}/{len(tests)} passed")
    return failures


if __name__ == "__main__":
    sys.exit(1 if _run_all() else 0)
