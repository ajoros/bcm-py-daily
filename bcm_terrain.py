"""Python terrain-.inp builder (mklatlon + Skyview + Solinp).

From a DEM + projection box:
  lat/lon, slope/aspect, sky-view, 36 horizon angles, then one .inp.
"""
from __future__ import annotations
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bcm_io import read_asc, write_inp


def latlon_from_box(nrows, ncols, west, east, south, north, left, right, top, bottom):
    """mklatlon layout: lat at the north edge of the cell, lon at the east edge."""
    xs = left + (np.arange(ncols) + 1) * ((right - left) / ncols)
    ys = top - np.arange(nrows) * ((top - bottom) / nrows)
    lon = west + (xs - left) / (right - left) * (east - west)
    lat = south + (ys - bottom)[:, None] / (top - bottom) * (north - south)
    lon = np.broadcast_to(lon, (nrows, ncols)).copy()
    return lat, lon


def xy_solinp(nrows, ncols, xll, yll, cellsize):
    """East/north as Solinp writes them (one cellsize east/north of cell center)."""
    cols = np.arange(ncols) + 1
    rows = np.arange(nrows, 0, -1)
    east = cols * cellsize + xll + cellsize / 2.0
    north = rows * cellsize + yll + cellsize / 2.0
    return np.broadcast_to(east, (nrows, ncols)).copy(), np.repeat(north[:, None], ncols, axis=1)


def slope_aspect_from_dem(elev, cellsize, nodata=-9999.0):
    """Horn 8-neighbor slope (deg) and aspect (deg).

    Aspect is Skyview's convention: atan2(dzdx, dzdy)+180, flat -> 0.
    Nodata neighbors are replaced with the cell's own elevation so every
    land cell still gets a number (we do not write -9999 on land).
    """
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
    sl = np.where(ocean, nodata, sl)
    asp = np.where(ocean, nodata, asp)
    return sl, asp


# Skyviewv1.f90 lattice steps (xang, yang). Their row index increases north;
# numpy row increases south, so a step is drow=-yang, dcol=+xang.
_SKY_XANG = np.array([
    1, 1, 2, 1, 3, 2, 3, 1, 3, 2, 3, 1, 2, 1, 1, 0,
    -1, -1, -2, -1, -3, -2, -3, -1, -3, -2, -3, -1, -2, -1, -1, 0,
], dtype=np.intp)
_SKY_YANG = np.array([
    3, 2, 3, 1, 2, 1, 1, 0, -1, -1, -2, -1, -3, -2, -3, -1,
    -3, -2, -3, -1, -2, -1, -1, 0, 1, 1, 2, 1, 3, 2, 3, 1,
], dtype=np.intp)


def _interp_32_to_36(a32):
    """32 lattice horizons → 36 ten-degree bins. Copied from Skyviewv1.f90."""
    a = [None] + [a32[:, i] for i in range(32)]  # 1-based like Fortran ANGLE
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
    """Horizons + sky-view as in Alan Flint Skyviewv1.f90.

    32 integer lattice rays, interpolate to 36 ten-degree bins, sky-view
    is the slope-adjusted mean (not 1-sin(H)).
    """
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
            zz = np.where(ok, z[np.clip(ri, 0, nrows - 1),
                                np.clip(ci, 0, ncols - 1)], nodata)
            higher = ok & (maxelev < zz)
            maxelev = np.where(higher, zz, maxelev)
            dist = np.where(higher, k * step_len, dist)
            dist_use = np.where(dist == 0.0, 1.0, dist)
            agl = np.degrees(np.arctan((maxelev - z0) / dist_use))
            hmax = np.maximum(hmax, agl)
        a32[:, j] = hmax
    a32 = np.maximum(a32, 0.0)
    ang = _interp_32_to_36(a32)
    ang2 = np.floor(ang + 0.5)
    if sl is None or asp is None:
        sl, asp = slope_aspect_from_dem(elev, cellsize, nodata)
    sl_l = sl[ir, ic]
    asp_l = np.where(asp[ir, ic] == -1.0, 0.0, asp[ir, ic])
    az = (np.arange(1, 37) * 10.0)[None, :]
    ridge = 90.0 - ang2
    slr = np.radians(sl_l)[:, None]
    rdr = np.radians(ridge)
    costheta = (np.cos(slr) * np.cos(rdr)
                + np.sin(slr) * np.sin(rdr)
                * np.cos(np.radians(az - asp_l[:, None])))
    viewf = np.degrees(np.arccos(np.clip(costheta, -1.0, 1.0)))
    viewf = np.minimum(viewf, 90.0)
    sky = np.full(elev.shape, nodata, dtype=np.float64)
    sky[land] = np.clip(viewf.mean(axis=1) / 90.0, 0.0, 1.0)
    ridge_out = np.zeros((nrows, ncols, 37), dtype=np.float64)
    ridge_out[ir, ic, 1:37] = ang2
    return sky, ridge_out


def parse_mklatlon_ctl(path):
    """Read west/east/south/north and the projected box from mklatlon.ctl."""
    nums = []
    with open(path, encoding="latin-1") as f:
        next(f)  # title
        for line in f:
            parts = line.split()
            if not parts:
                continue
            tok = parts[0].strip().strip("'")
            try:
                nums.append(float(tok))
            except ValueError:
                continue
    west, east, north, south, left, right, top, bottom = nums[:8]
    return dict(west=west, east=east, south=south, north=north,
                left=left, right=right, top=top, bottom=bottom)


def build_from_box(dem_path, out_path, west, east, south, north,
                   radius=None, donor_inp=None):
    """Write a .inp from a DEM + geographic box. No mklatlon.ctl."""
    hdr, elev = read_asc(dem_path)
    nrows, ncols = int(hdr["nrows"]), int(hdr["ncols"])
    cs, nd = hdr["cellsize"], hdr["nodata_value"]
    left = hdr["xllcorner"]
    bottom = hdr["yllcorner"]
    right = left + ncols * cs
    top = bottom + nrows * cs
    if radius is None:
        span = max(right - left, top - bottom)
        radius = max(8 * cs, min(75000.0, 2 * span))
    return _write_inp_from_dem(
        hdr, elev, out_path,
        west, east, south, north, left, right, top, bottom,
        radius, donor_inp,
    )


def build_from_dem(dem_path, ctl_path, out_path, donor_inp=None, radius=75000.0):
    """Write a .inp from the DEM (and mklatlon.ctl box).

    If donor_inp is given, sky/ridge are copied from it (debug). Otherwise
    they are computed from the DEM. Slope/aspect are whole degrees (Solinp).
    """
    hdr, elev = read_asc(dem_path)
    nrows, ncols = int(hdr["nrows"]), int(hdr["ncols"])
    cs, nd = hdr["cellsize"], hdr["nodata_value"]
    box = parse_mklatlon_ctl(ctl_path)
    return _write_inp_from_dem(
        hdr, elev, out_path,
        box["west"], box["east"], box["south"], box["north"],
        box["left"], box["right"], box["top"], box["bottom"],
        radius, donor_inp,
    )


def _write_inp_from_dem(hdr, elev, out_path, west, east, south, north,
                        left, right, top, bottom, radius, donor_inp):
    nrows, ncols = int(hdr["nrows"]), int(hdr["ncols"])
    cs, nd = hdr["cellsize"], hdr["nodata_value"]
    lat, lon = latlon_from_box(
        nrows, ncols, west, east, south, north, left, right, top, bottom,
    )
    sl, asp = slope_aspect_from_dem(elev, cs, nd)
    land = elev != nd
    lat = np.where(land, lat, nd)
    lon = np.where(land, lon, nd)
    if donor_inp is None:
        sky, ridge = horizon_sky_from_dem(elev, cs, radius, nd, sl=sl, asp=asp)
    sl = np.where(land & (sl != nd), np.rint(sl), nd)
    asp = np.where(land & (asp != nd), np.rint(asp), nd)
    east, north = xy_solinp(nrows, ncols, hdr["xllcorner"], hdr["yllcorner"], cs)
    site = np.arange(1, nrows * ncols + 1, dtype=np.float64).reshape(nrows, ncols)
    if donor_inp is not None:
        sky, ridge = donor_inp["sky"], donor_inp["ridge"]
    write_inp(out_path, {
        "site": site, "east": east, "north": north,
        "lat": lat, "lon": lon, "sl": sl, "asp": asp,
        "elev": elev, "sky": sky, "ridge": ridge,
    })
    return out_path


def stitch_grids(workdir, out_path, dem_name="mok_270m.asc",
                 lat_name="mok_270lat.asc", lon_name="mok_270lon.asc",
                 slp_name="slpfile.asc", asp_name="aspfile.asc",
                 sky_name="mok_270m_sky.asc"):
    """Python Solinp: ASCII grids → one .inp."""
    hdr, elev = read_asc(os.path.join(workdir, dem_name))
    nrows, ncols = int(hdr["nrows"]), int(hdr["ncols"])
    cs, nd = hdr["cellsize"], hdr["nodata_value"]
    _, lat = read_asc(os.path.join(workdir, lat_name))
    _, lon = read_asc(os.path.join(workdir, lon_name))
    _, sl = read_asc(os.path.join(workdir, slp_name))
    _, asp = read_asc(os.path.join(workdir, asp_name))
    asp = np.where(asp == -1.0, 0.0, asp)
    _, sky = read_asc(os.path.join(workdir, sky_name))
    ridge = np.zeros((nrows, ncols, 37), dtype=np.float64)
    for k in range(1, 37):
        _, ridge[:, :, k] = read_asc(os.path.join(workdir, f"ang{k * 10:03d}.asc"))
    east, north = xy_solinp(nrows, ncols, hdr["xllcorner"], hdr["yllcorner"], cs)
    site = np.arange(1, nrows * ncols + 1, dtype=np.float64).reshape(nrows, ncols)
    write_inp(out_path, {
        "site": site, "east": east, "north": north,
        "lat": lat, "lon": lon, "sl": sl, "asp": asp,
        "elev": elev, "sky": sky, "ridge": ridge,
    })
    return out_path


def _check():
    """Fails if lat/lon/slope/aspect drift off the Skyview / mklatlon grids."""
    root = os.path.normpath(os.path.join(HERE, ".."))
    work = os.path.join(root, "SolarFilesFromMichelle", "_skyview_clean")
    hdr, dem = read_asc(os.path.join(work, "mok_270m.asc"))
    _, sl_sv = read_asc(os.path.join(work, "slpfile.asc"))
    _, asp_sv = read_asc(os.path.join(work, "aspfile.asc"))
    asp_sv = np.where(asp_sv == -1.0, 0.0, asp_sv)
    _, lat_sv = read_asc(os.path.join(work, "mok_270lat.asc"))
    _, lon_sv = read_asc(os.path.join(work, "mok_270lon.asc"))
    nd = hdr["nodata_value"]
    nrows, ncols = dem.shape
    lat, lon = latlon_from_box(
        nrows, ncols,
        -121.6721530, -119.7013026, 38.0328759, 38.6993700,
        -145265.836400, 26184.163600, 77206.665800, 1876.665800,
    )
    sl, asp = slope_aspect_from_dem(dem, hdr["cellsize"], nd)
    _, sky_sv = read_asc(os.path.join(work, "mok_270m_sky.asc"))
    ang_sv = np.stack([
        read_asc(os.path.join(work, f"ang{i * 10:03d}.asc"))[1] for i in range(1, 37)
    ], axis=2)
    sky, ridge = horizon_sky_from_dem(
        dem, hdr["cellsize"], 75000.0, nd, sl=sl, asp=asp)
    land = (dem != nd) & (lat_sv != nd) & (sky_sv != nd)
    finite = land & (sl != nd) & (sl_sv != nd)
    lat_err = np.max(np.abs(lat - lat_sv)[land])
    lon_err = np.max(np.abs(lon - lon_sv)[land])
    sl_rmse = float(np.sqrt(np.mean((sl - sl_sv)[finite] ** 2)))
    dasp = ((asp - asp_sv + 180.0) % 360.0 - 180.0)[finite]
    asp_rmse = float(np.sqrt(np.mean(dasp ** 2)))
    sky_rmse = float(np.sqrt(np.mean((sky - sky_sv)[land] ** 2)))
    rd = np.concatenate([(ridge[:, :, k] - ang_sv[:, :, k - 1])[land] for k in range(1, 37)])
    ang_rmse = float(np.sqrt(np.mean(rd ** 2)))
    print(f"lat maxabs={lat_err:.6f}  lon maxabs={lon_err:.6f}")
    print(f"slope rmse={sl_rmse:.4f}  aspect rmse={asp_rmse:.4f}")
    print(f"sky rmse={sky_rmse:.4f}  horizon rmse={ang_rmse:.4f}")
    assert lat_err < 1e-4, lat_err
    assert lon_err < 1e-4, lon_err
    assert sl_rmse < 0.1, sl_rmse
    assert asp_rmse < 0.5, asp_rmse
    assert sky_rmse < 0.002, sky_rmse
    assert ang_rmse < 0.01, ang_rmse
    print("bcm_terrain check ok")


if __name__ == "__main__":
    _check()
