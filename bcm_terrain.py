"""Build a BCM Daily terrain .inp from a DEM (no Fortran).

Needs the DEM ASCII plus the geographic / projected box so each cell
gets lat/lon. Slope, aspect, sky-view, and 36 horizon angles are
computed in Python (Flint Skyview method).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bcm_io import read_asc, write_inp


def latlon_from_box(nrows, ncols, west, east, south, north, left, right, top, bottom):
    """Lat at the north edge of the cell, lon at the east edge."""
    xs = left + (np.arange(ncols) + 1) * ((right - left) / ncols)
    ys = top - np.arange(nrows) * ((top - bottom) / nrows)
    lon = west + (xs - left) / (right - left) * (east - west)
    lat = south + (ys - bottom)[:, None] / (top - bottom) * (north - south)
    lon = np.broadcast_to(lon, (nrows, ncols)).copy()
    return lat, lon


def xy_from_dem(nrows, ncols, xll, yll, cellsize):
    cols = np.arange(ncols) + 1
    rows = np.arange(nrows, 0, -1)
    east = cols * cellsize + xll + cellsize / 2.0
    north = rows * cellsize + yll + cellsize / 2.0
    return np.broadcast_to(east, (nrows, ncols)).copy(), np.repeat(
        north[:, None], ncols, axis=1
    )


def slope_aspect_from_dem(elev, cellsize, nodata=-9999.0):
    """Horn 8-neighbor slope (deg) and aspect (deg). Flat -> 0."""
    z = elev.astype(np.float64)
    p = np.pad(z, 1, mode="edge")
    center = p[1:-1, 1:-1]

    def nbr(a):
        return np.where(a == nodata, center, a)

    nw, n, ne = nbr(p[:-2, :-2]), nbr(p[:-2, 1:-1]), nbr(p[:-2, 2:])
    w, e = nbr(p[1:-1, :-2]), nbr(p[1:-1, 2:])
    sw, s, se = nbr(p[2:, :-2]), nbr(p[2:, 1:-1]), nbr(p[2:, 2:])
    dzdx = ((ne + 2 * e + se) - (nw + 2 * w + sw)) / (8.0 * cellsize)
    dzdy = ((nw + 2 * n + ne) - (sw + 2 * s + se)) / (8.0 * cellsize)
    sl = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))
    asp = (np.degrees(np.arctan2(dzdx, dzdy)) + 180.0) % 360.0
    asp = np.where(sl < 1e-6, 0.0, asp)
    ocean = elev == nodata
    return np.where(ocean, nodata, sl), np.where(ocean, nodata, asp)


# Flint Skyview lattice (x east, y north). Numpy row increases south.
_SKY_XANG = np.array(
    [1, 1, 2, 1, 3, 2, 3, 1, 3, 2, 3, 1, 2, 1, 1, 0,
     -1, -1, -2, -1, -3, -2, -3, -1, -3, -2, -3, -1, -2, -1, -1, 0],
    dtype=np.intp,
)
_SKY_YANG = np.array(
    [3, 2, 3, 1, 2, 1, 1, 0, -1, -1, -2, -1, -3, -2, -3, -1,
     -3, -2, -3, -1, -2, -1, -1, 0, 1, 1, 2, 1, 3, 2, 3, 1],
    dtype=np.intp,
)


def _interp_32_to_36(a32):
    """32 lattice horizons → 36 ten-degree bins."""
    a = [None] + [a32[:, i] for i in range(32)]
    out = np.empty((a32.shape[0], 36), dtype=np.float64)
    out[:, 0] = (10.0 / 18.43) * (a[1] - a[32]) + a[32]
    out[:, 1] = (20 - 18.43) / (26.57 - 18.43) * (a[2] - a[1]) + a[1]
    out[:, 2] = (30 - 26.57) / (33.69 - 26.57) * (a[3] - a[2]) + a[2]
    out[:, 3] = (40 - 33.69) / (45 - 33.69) * (a[4] - a[3]) + a[3]
    out[:, 4] = (50 - 45) / (56.31 - 45) * (a[5] - a[4]) + a[4]
    out[:, 5] = (60 - 56.31) / (63.43 - 56.31) * (a[6] - a[5]) + a[5]
    out[:, 6] = (70 - 63.43) / (71.57 - 63.43) * (a[7] - a[6]) + a[6]
    out[:, 7] = (80 - 71.57) / (90 - 71.57) * (a[8] - a[7]) + a[7]
    out[:, 8] = a[8]
    out[:, 9] = (80 - 71.57) / (90 - 71.57) * (a[8] - a[9]) + a[9]
    out[:, 10] = (70 - 63.43) / (71.57 - 63.43) * (a[9] - a[10]) + a[10]
    out[:, 11] = (60 - 56.31) / (63.43 - 56.31) * (a[10] - a[11]) + a[11]
    out[:, 12] = (50 - 45) / (56.31 - 45) * (a[11] - a[12]) + a[12]
    out[:, 13] = (40 - 33.69) / (45 - 33.69) * (a[12] - a[13]) + a[13]
    out[:, 14] = (30 - 26.57) / (33.69 - 26.57) * (a[13] - a[14]) + a[14]
    out[:, 15] = (20 - 18.43) / (26.57 - 18.43) * (a[14] - a[15]) + a[15]
    out[:, 16] = (10.0 / 18.43) * (a[15] - a[16]) + a[16]
    out[:, 17] = a[16]
    out[:, 18] = (10.0 / 18.43) * (a[17] - a[16]) + a[16]
    out[:, 19] = (20 - 18.43) / (26.57 - 18.43) * (a[18] - a[17]) + a[17]
    out[:, 20] = (30 - 26.57) / (33.69 - 26.57) * (a[19] - a[18]) + a[18]
    out[:, 21] = (40 - 33.69) / (45 - 33.69) * (a[20] - a[19]) + a[19]
    out[:, 22] = (50 - 45) / (56.31 - 45) * (a[21] - a[20]) + a[20]
    out[:, 23] = (60 - 56.31) / (63.43 - 56.31) * (a[22] - a[21]) + a[21]
    out[:, 24] = (70 - 63.43) / (71.57 - 63.43) * (a[23] - a[22]) + a[22]
    out[:, 25] = (80 - 71.57) / (90 - 71.57) * (a[24] - a[23]) + a[23]
    out[:, 26] = a[24]
    out[:, 27] = (80 - 71.57) / (90 - 71.57) * (a[24] - a[25]) + a[25]
    out[:, 28] = (70 - 63.43) / (71.57 - 63.43) * (a[25] - a[26]) + a[26]
    out[:, 29] = (60 - 56.31) / (63.43 - 56.31) * (a[26] - a[27]) + a[27]
    out[:, 30] = (50 - 45) / (56.31 - 45) * (a[27] - a[28]) + a[28]
    out[:, 31] = (40 - 33.69) / (45 - 33.69) * (a[28] - a[29]) + a[29]
    out[:, 32] = (30 - 26.57) / (33.69 - 26.57) * (a[29] - a[30]) + a[30]
    out[:, 33] = (20 - 18.43) / (26.57 - 18.43) * (a[30] - a[31]) + a[31]
    out[:, 34] = (10.0 / 18.43) * (a[31] - a[32]) + a[32]
    out[:, 35] = a[32]
    return out


def horizon_sky_from_dem(elev, cellsize, radius=75000.0, nodata=-9999.0,
                         sl=None, asp=None):
    """32 lattice rays, 36 interpolated bins, slope-adjusted sky-view."""
    nrows, ncols = elev.shape
    land = elev != nodata
    ir, ic = np.nonzero(land)
    z0 = elev[ir, ic].astype(np.float64)
    nland = int(z0.size)
    nstep = max(1, int(radius / cellsize))
    a32 = np.zeros((nland, 32), dtype=np.float64)
    print(f"  horizons: {nland} cells, {nstep} steps x 32 lattice rays ...")
    z = elev.astype(np.float64)
    for j in range(32):
        dx, dy = int(_SKY_XANG[j]), int(_SKY_YANG[j])
        maxelev = z0.copy()
        dist = np.zeros(nland, dtype=np.float64)
        hmax = np.zeros(nland, dtype=np.float64)
        step_len = float(np.hypot(dx, dy) * cellsize)
        for k in range(1, nstep + 1):
            ri = ir - k * dy
            ci = ic + k * dx
            ok = (ri >= 0) & (ci >= 0) & (ri < nrows) & (ci < ncols)
            if not ok.any():
                break
            zz = np.where(
                ok,
                z[np.clip(ri, 0, nrows - 1), np.clip(ci, 0, ncols - 1)],
                nodata,
            )
            higher = ok & (maxelev < zz)
            maxelev = np.where(higher, zz, maxelev)
            dist = np.where(higher, k * step_len, dist)
            dist_use = np.where(dist == 0.0, 1.0, dist)
            agl = np.degrees(np.arctan((maxelev - z0) / dist_use))
            hmax = np.maximum(hmax, agl)
        a32[:, j] = hmax
    ang2 = np.floor(np.maximum(_interp_32_to_36(a32), 0.0) + 0.5)
    if sl is None or asp is None:
        sl, asp = slope_aspect_from_dem(elev, cellsize, nodata)
    sl_l = sl[ir, ic]
    asp_l = np.where(asp[ir, ic] == -1.0, 0.0, asp[ir, ic])
    az = (np.arange(1, 37) * 10.0)[None, :]
    ridge = 90.0 - ang2
    slr = np.radians(sl_l)[:, None]
    rdr = np.radians(ridge)
    costheta = (
        np.cos(slr) * np.cos(rdr)
        + np.sin(slr) * np.sin(rdr) * np.cos(np.radians(az - asp_l[:, None]))
    )
    viewf = np.minimum(np.degrees(np.arccos(np.clip(costheta, -1.0, 1.0))), 90.0)
    sky = np.full(elev.shape, nodata, dtype=np.float64)
    sky[land] = np.clip(viewf.mean(axis=1) / 90.0, 0.0, 1.0)
    ridge_out = np.zeros((nrows, ncols, 37), dtype=np.float64)
    ridge_out[ir, ic, 1:37] = ang2
    return sky, ridge_out


def parse_box_file(path):
    """Eight numbers: west east north south left right top bottom."""
    nums = []
    with open(path, encoding="latin-1") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            tok = parts[0].strip().strip("'")
            try:
                nums.append(float(tok))
            except ValueError:
                continue
    if len(nums) < 8:
        raise ValueError(f"{path}: need 8 numbers (west east north south left right top bottom)")
    west, east, north, south, left, right, top, bottom = nums[:8]
    return dict(west=west, east=east, south=south, north=north,
                left=left, right=right, top=top, bottom=bottom)


def build_from_dem(dem_path, out_path, box, radius=75000.0):
    """Write a terrain .inp from a DEM ASCII and a lat/lon / projected box."""
    hdr, elev = read_asc(dem_path)
    nrows, ncols = int(hdr["nrows"]), int(hdr["ncols"])
    cs, nd = hdr["cellsize"], hdr["nodata_value"]
    lat, lon = latlon_from_box(nrows, ncols, **box)
    sl, asp = slope_aspect_from_dem(elev, cs, nd)
    land = elev != nd
    lat = np.where(land, lat, nd)
    lon = np.where(land, lon, nd)
    sky, ridge = horizon_sky_from_dem(elev, cs, radius, nd, sl=sl, asp=asp)
    sl = np.where(land & (sl != nd), np.rint(sl), nd)
    asp = np.where(land & (asp != nd), np.rint(asp), nd)
    east, north = xy_from_dem(nrows, ncols, hdr["xllcorner"], hdr["yllcorner"], cs)
    site = np.arange(1, nrows * ncols + 1, dtype=np.float64).reshape(nrows, ncols)
    write_inp(out_path, {
        "site": site, "east": east, "north": north,
        "lat": lat, "lon": lon, "sl": sl, "asp": asp,
        "elev": elev, "sky": sky, "ridge": ridge,
    })
    print(f"  wrote {out_path}")
    return out_path


def _self_check():
    """Tiny hill: sky in (0, 1], peak more open than a lower neighbor."""
    elev = np.array(
        [[100, 100, 100, 100, 100],
         [100, 130, 150, 130, 100],
         [100, 150, 220, 150, 100],
         [100, 130, 150, 130, 100],
         [100, 100, 100, 100, 100]],
        dtype=np.float64,
    )
    sl, asp = slope_aspect_from_dem(elev, 270.0)
    sky, ridge = horizon_sky_from_dem(elev, 270.0, radius=810.0, sl=sl, asp=asp)
    assert np.all((sky > 0.0) & (sky <= 1.0))
    assert sky[2, 2] >= sky[1, 0]
    assert ridge.shape == (5, 5, 37)
    lat, lon = latlon_from_box(
        5, 5, west=-121.0, east=-120.0, south=38.0, north=39.0,
        left=0.0, right=1350.0, top=1350.0, bottom=0.0,
    )
    assert lat[0, 0] > lat[-1, 0]
    assert lon[0, -1] > lon[0, 0]
    print("bcm_terrain check ok")


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Build a BCM Daily terrain .inp from a DEM ASCII.")
    p.add_argument("dem", nargs="?", help="DEM as ESRI ASCII (.asc)")
    p.add_argument("-o", "--out", help="Output .inp path")
    p.add_argument("--west", type=float, help="West longitude")
    p.add_argument("--east", type=float, help="East longitude")
    p.add_argument("--south", type=float, help="South latitude")
    p.add_argument("--north", type=float, help="North latitude")
    p.add_argument("--left", type=float, help="Projected west (same units as DEM)")
    p.add_argument("--right", type=float, help="Projected east")
    p.add_argument("--top", type=float, help="Projected north")
    p.add_argument("--bottom", type=float, help="Projected south")
    p.add_argument("--box-file", help="Text file with 8 numbers: west east north south left right top bottom")
    p.add_argument("--radius", type=float, default=75000.0,
                   help="Horizon search radius in DEM map units (default 75000)")
    p.add_argument("--check", action="store_true", help="Run the built-in self-check and exit")
    args = p.parse_args(argv)
    if args.check:
        _self_check()
        return 0
    if not args.dem or not args.out:
        p.error("dem and -o are required (or use --check)")
    if args.box_file:
        box = parse_box_file(args.box_file)
    else:
        keys = ("west", "east", "south", "north", "left", "right", "top", "bottom")
        missing = [k for k in keys if getattr(args, k) is None]
        if missing:
            p.error("need --box-file or " + " ".join(f"--{k}" for k in missing))
        box = {k: getattr(args, k) for k in keys}
    build_from_dem(args.dem, args.out, box, radius=args.radius)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
