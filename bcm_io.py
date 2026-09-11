"""
bcm_io.py — BCM Daily v8.1 input/output helpers

Provides:
  parse_ctl(path)        → dict of configuration values (mirrors BCM_Dailyv81.ctl)
  read_asc(path)         → (header_dict, 2-D numpy array)
  write_asc(path, hdr, arr) → write a 2-D numpy array as an ESRI ASCII grid
  read_inp(path, nrows, ncols) → structured dict of per-cell terrain data
"""

from __future__ import annotations
import os
import numpy as np
from typing import Tuple


# ---------------------------------------------------------------------------
# Helper: read one token line from CTL (first 30 characters, strip whitespace)
# ---------------------------------------------------------------------------
def _token(line: str) -> str:
    """Return the first meaningful token from a CTL line (30-char field)."""
    return line[:30].strip()


def _scalar(line: str):
    """Return the first numeric token from a CTL line."""
    return line.split()[0]


# ---------------------------------------------------------------------------
# CTL parser helpers
# ---------------------------------------------------------------------------
def _getnums(line: str, n: int = None) -> list:
    """
    Extract numeric tokens from a CTL line, stopping at the first
    non-numeric token (mirrors Fortran list-directed I/O stopping at
    non-numeric data on the same record).
    Returns up to n values; pads with None if fewer are found.
    """
    result = []
    for tok in line.split():
        try:
            result.append(float(tok))
            if n is not None and len(result) == n:
                break
        except ValueError:
            break
    return result


def _n1(line: str, default=0.0):
    """First numeric token from a CTL line, or default."""
    v = _getnums(line, 1)
    return v[0] if v else default


# ---------------------------------------------------------------------------
# CTL parser
# ---------------------------------------------------------------------------
def parse_ctl(path: str) -> dict:
    """
    Read BCM_Dailyv81.ctl and return a configuration dict.

    The CTL file is read line-by-line, extracting numeric tokens and
    30-character filename fields.  Multi-value switch lines use _getnums()
    which stops at the first non-numeric token, matching the behaviour of
    Fortran list-directed I/O on these lines.

    Two "phantom" reads in the Fortran source (rosflag/rvspercent and
    hourlypetflag) do not correspond to dedicated CTL lines in this file;
    they are set to defaults (0 / 0.0) without consuming a line.
    """
    cfg = {}
    with open(path, "r", encoding="latin-1") as f:
        lines = [l.rstrip("\n") for l in f]

    i = 0
    def nxt():
        nonlocal i
        while i < len(lines):
            l = lines[i]; i += 1
            return l
        return ""

    cfg["header"]        = nxt()     # line 1 – description
    nxt()                            # line 2 – section comment
    cfg["outfile"]       = _token(nxt())
    cfg["areafile"]      = _token(nxt())
    cfg["areatable"]     = _token(nxt())
    cfg["dn1"]           = int(_n1(nxt()))
    cfg["yn1"]           = int(_n1(nxt()))
    cfg["dn2"]           = int(_n1(nxt()))
    cfg["yn2"]           = int(_n1(nxt()))
    cfg["hstep"]         = int(_n1(nxt()))
    cfg["std"]           = _n1(nxt())
    cfg["prior"]         = int(_n1(nxt()))
    cfg["istor"]         = _n1(nxt())
    nxt()                            # section comment
    cfg["demfile"]       = _token(nxt())
    cfg["soildfile"]     = _token(nxt())
    cfg["wpfile"]        = _token(nxt())
    cfg["fcfile"]        = _token(nxt())
    cfg["porfile"]       = _token(nxt())
    cfg["soilksfile"]    = _token(nxt())
    cfg["geolfile"]      = _token(nxt())
    cfg["vegfile"]       = _token(nxt())
    cfg["bc_a_file"]     = _token(nxt())
    cfg["bc_b_file"]     = _token(nxt())
    cfg["bc_c_file"]     = _token(nxt())
    cfg["pt_alpha_file"] = _token(nxt())
    cfg["cwd_dry_file"]  = _token(nxt())
    cfg["inpfile"]       = _token(nxt())
    nxt()                            # section comment (optional layer files)
    cfg["rockksfile"]       = _token(nxt())
    cfg["snowaccumfile"]    = _token(nxt())
    cfg["mfmaxfile"]        = _token(nxt())
    cfg["mfminfile"]        = _token(nxt())
    cfg["aridityfile"]      = _token(nxt())
    cfg["maskfile"]         = _token(nxt())
    cfg["hourpetfile"]      = _token(nxt())
    # Fortran: READ header2 (eats climate-dir path), READ ClimateDirFile
    nxt()                            # eats D:\Mok_BCM_Heggli\... (Fortran header2)
    cfg["climate_dir_file"] = _token(nxt())  # !----- BELOW use an on-off SWITCH
    # ------- SWITCH lines (one per Fortran READ, numeric tokens extracted) -------
    # NOTE: Fortran eats the first switch line (ClimateDirflag) as another
    # header2 before the switch block, so ClimateDirflag comes from line 39
    # and subsequent flags are in comment order below.
    cfg["climate_dir_flag"] = int(_n1(nxt()))   # line 39: Read pet from another dir
    cfg["rockks_flag"]      = int(_n1(nxt()))   # line 40: GEOL Bedrock
    dd = _getnums(nxt(), 4)                     # line 41: SOIL DRYDOWN (0 0.65 -10 0.05)
    cfg["drydown_flag"] = int(dd[0]) if dd else 0
    cfg["drydown_a"]    = dd[1] if len(dd) > 1 else 0.65
    cfg["drydown_b"]    = dd[2] if len(dd) > 2 else -10.0
    cfg["drydown_c"]    = dd[3] if len(dd) > 3 else 0.05
    cfg["alphaprime_flag"]  = int(_n1(nxt()))   # line 42: PT
    rr = _getnums(nxt(), 2)                     # line 43: RCH/RUN (1 0.300)
    cfg["rchrun_flag"]  = int(rr[0]) if rr else 0
    cfg["rchrun_limit"] = rr[1] if len(rr) > 1 else 0.3
    uu = _getnums(nxt(), 4)                     # line 44: URBAN (0 53 0.1 0.5)
    cfg["urban_flag"]   = int(uu[0]) if uu else 0
    cfg["urban_veg"]    = int(uu[1]) if len(uu) > 1 else 53
    cfg["urban_soild"]  = uu[2] if len(uu) > 2 else 0.1
    cfg["urban_ks"]     = uu[3] if len(uu) > 3 else 0.5
    # rosflag/rvspercent: no dedicated CTL line — set defaults
    cfg["rosflag"]      = 0
    cfg["rvs_percent"]  = 0.0
    cfg["snow_flag"]     = int(_n1(nxt()))       # line 45: SNOW
    cfg["snowaccum_t"]   = _n1(nxt())            # line 46: snowaccum temperature
    cfg["maxmf"]         = _n1(nxt())            # line 47: max melt factor
    cfg["minmf"]         = _n1(nxt())            # line 48: min melt factor
    mf_vals = _getnums(nxt(), 3)                 # line 49: maf tipm nmf
    cfg["maf"]  = mf_vals[0] if mf_vals       else 7.0
    cfg["tipm"] = mf_vals[1] if len(mf_vals)>1 else 0.1
    cfg["nmf"]  = mf_vals[2] if len(mf_vals)>2 else 0.21
    cfg["solar_flag"]      = int(_n1(nxt()))     # line 50: SOLAR
    # hourlypetflag: no dedicated CTL line — default 0
    cfg["hourly_pet_flag"] = 0
    cfg["rainfrac_flag"]   = int(_n1(nxt()))     # line 51: RAIN fraction
    sb = _getnums(nxt(), 2)                      # line 52: SUBLIMATION (0 1.00)
    cfg["subl_flag"] = int(sb[0]) if sb else 0
    cfg["subl_par"]  = sb[1] if len(sb) > 1 else 1.0
    cfg["mask_flag"]    = int(_n1(nxt()))        # line 53: MASK
    cfg["ingest_flag"]  = int(_n1(nxt()))        # line 54: INGEST
    cfg["num_ingest"]   = int(_n1(nxt()))        # line 55: number of ingest pairs
    nxt()                            # ingest file header comment
    ingest_pairs = []
    for _ in range(cfg["num_ingest"]):
        raw = nxt()
        digest = raw[:60].strip()
        ingest = raw[60:].strip()
        ingest_pairs.append((digest, ingest))
    cfg["ingest_pairs"] = ingest_pairs

    nxt()                            # output-flag section comment
    cfg["aet_flag"]  = int(_n1(nxt()))
    cfg["cld_flag"]  = int(_n1(nxt()))
    cfg["cwd_flag"]  = int(_n1(nxt()))
    cfg["evp_flag"]  = int(_n1(nxt()))
    cfg["exc_flag"]  = int(_n1(nxt()))
    cfg["mlt_flag"]  = int(_n1(nxt()))
    cfg["pet_flag"]  = int(_n1(nxt()))
    cfg["rad_flag"]  = int(_n1(nxt()))
    cfg["rch_flag"]  = int(_n1(nxt()))
    cfg["run_flag"]  = int(_n1(nxt()))
    cfg["sbl_flag"]  = int(_n1(nxt()))
    cfg["snw_flag"]  = int(_n1(nxt()))
    cfg["smd_flag"]  = int(_n1(nxt()))
    cfg["smr_flag"]  = int(_n1(nxt()))
    cfg["rvs_flag"]  = int(_n1(nxt()))

    nxt()                            # geology lookup comment
    cfg["n_geol"] = int(_scalar(nxt()))
    nxt(); nxt()                     # two header rows
    rockks = {}
    for _ in range(cfg["n_geol"]):
        parts = nxt().split()
        rockks[int(parts[0])] = float(parts[1])
    cfg["rockks_table"] = rockks

    nxt()                            # vegetation lookup comment
    cfg["n_veg"] = int(_scalar(nxt()))
    nxt(); nxt(); nxt()              # three header rows
    veg_params = {}
    for _ in range(cfg["n_veg"]):
        parts = nxt().split()
        vid = int(parts[0])
        veg_params[vid] = {
            "init_lai":   float(parts[1]),
            "up_limit":   float(parts[2]),
            "dn_limit":   float(parts[3]),
            "up_rate":    float(parts[4]),
            "dn_rate":    float(parts[5]),
            "root_depth": float(parts[6]),
            "kfactor":    [float(parts[7+m]) for m in range(12)],   # jan..dec
        }
    cfg["veg_params"] = veg_params

    # Monthly average rad and ppt files (12 each) – read until the
    # section exists; swallow any remaining content gracefully
    rad_files = []
    ppt_files = []
    try:
        nxt()                        # section header
        for _ in range(12):
            rad_files.append(_token(nxt()))
        for _ in range(12):
            ppt_files.append(_token(nxt()))
    except Exception:
        pass
    cfg["rad_ave_files"] = rad_files
    cfg["ppt_ave_files"] = ppt_files

    return cfg


# ---------------------------------------------------------------------------
# ESRI ASCII raster reader / writer
# ---------------------------------------------------------------------------
_HDR_KEYS = ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "NODATA_value")

def read_asc(path: str) -> Tuple[dict, np.ndarray]:
    """Read an ESRI ASCII raster.  Returns (header_dict, 2-D float64 array)."""
    hdr = {}
    with open(path, "r", encoding="latin-1") as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
        data_lines = f.readlines()
    arr = np.array(
        [row.split() for row in data_lines if row.strip()],
        dtype=np.float64
    )
    return hdr, arr


def write_asc(path: str, hdr: dict, arr: np.ndarray, fmt: str = "%.2f"):
    """Write a 2-D numpy array as an ESRI ASCII raster."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        f.write(f"ncols        {int(hdr['ncols']):>12d}\n")
        f.write(f"nrows        {int(hdr['nrows']):>12d}\n")
        f.write(f"xllcorner    {hdr['xllcorner']:>16.6f}\n")
        f.write(f"yllcorner    {hdr['yllcorner']:>16.6f}\n")
        f.write(f"cellsize     {hdr['cellsize']:>12.0f}\n")
        f.write(f"NODATA_value {hdr['nodata_value']:>12.1f}\n")
        np.savetxt(f, arr, fmt=fmt)


# ---------------------------------------------------------------------------
# mok_270m.inp reader
# ---------------------------------------------------------------------------
def read_inp(path: str, nrows: int, ncols: int) -> dict:
    """
    Read the terrain input file (mok_270m.inp).

    Each cell record has 46 space-delimited values (Fortran reads 45 and
    ignores a trailing extra token):
      SITE  EAST  NORTH  LAT  LON  SL  ASP  ELEV  SKY  RIDGE(1..36)  <extra>

    The Fortran RIDGE array is declared (0:36) so RIDGE(0)=0 by default;
    RIDGE(1..36) come from the file.  We store all 37 RIDGE values
    (prepending 0 for index 0) so that Fortran's SUNAZ index (0-36) maps
    directly.

    Returns a dict with keys:
      site, east, north, lat, lon, sl, asp, elev, sky, ridge
    Each scalar grid is 2-D: (nrows, ncols).
    'ridge' is 3-D: (nrows, ncols, 37) with ridge[:,:,0]=0 always.
    """
    print(f"  Reading terrain file {path} ({nrows} rows × {ncols} cols)…")
    n_total = nrows * ncols
    flat = np.zeros((n_total, 46), dtype=np.float64)   # 9 + 36 + 1 extra

    with open(path, "r", encoding="latin-1") as f:
        row_idx = 0
        for line in f:
            parts = line.split()
            if len(parts) >= 45:          # accept 45 or 46 tokens
                flat[row_idx, :min(len(parts), 46)] = [float(x) for x in parts[:46]]
                row_idx += 1
                if row_idx >= n_total:
                    break

    shape = (nrows, ncols)
    # Build RIDGE array with 37 slots (index 0..36).
    # File provides RIDGE(1..36) at flat[:,9:45]; RIDGE(0)=0.
    ridge37 = np.zeros((n_total, 37), dtype=np.float64)
    ridge37[:, 1:37] = flat[:, 9:45]    # RIDGE(1..36) from file tokens 9..44

    return {
        "site":  flat[:, 0].reshape(shape),
        "east":  flat[:, 1].reshape(shape),
        "north": flat[:, 2].reshape(shape),
        "lat":   flat[:, 3].reshape(shape),
        "lon":   flat[:, 4].reshape(shape),
        "sl":    flat[:, 5].reshape(shape),    # slope degrees
        "asp":   flat[:, 6].reshape(shape),    # aspect degrees
        "elev":  flat[:, 7].reshape(shape),    # elevation m
        "sky":   flat[:, 8].reshape(shape),    # sky-view factor
        "ridge": ridge37.reshape(nrows, ncols, 37),  # RIDGE(0..36)
    }


def write_inp(path: str, terrain: dict) -> None:
    """Write a BCM Daily terrain .inp (inverse of read_inp).

    45 tokens per cell plus a trailing 0. ridge[:,:,0] is not written.
    """
    elev = terrain["elev"]
    nrows, ncols = elev.shape
    site = terrain.get("site")
    east = terrain["east"]
    north = terrain["north"]
    lat = terrain["lat"]
    lon = terrain["lon"]
    sl = terrain["sl"]
    asp = terrain["asp"]
    sky = terrain["sky"]
    ridge = terrain["ridge"]
    n = 0
    with open(path, "w", encoding="latin-1", newline="\n") as f:
        for r in range(nrows):
            for c in range(ncols):
                n += 1
                sid = int(site[r, c]) if site is not None else n
                angs = " ".join(f"{float(ridge[r, c, k]):g}" for k in range(1, 37))
                f.write(
                    f"{sid:12d}{east[r, c]:12.1f}{north[r, c]:12.1f}"
                    f"{lat[r, c]:12.4f}{lon[r, c]:12.4f}"
                    f"{sl[r, c]:9.0f}{asp[r, c]:9.0f}"
                    f"{elev[r, c]:9.1f}{sky[r, c]:10.3f} "
                    f"{angs} 0\n"
                )


# ---------------------------------------------------------------------------
# File-naming helper (matches Fortran FORMAT 911/912/913 logic)
# ---------------------------------------------------------------------------
def day_filename(prefix: str, year: int, doy: int, suffix: str = ".asc") -> str:
    """Return a BCM-style filename such as ppt2010_001.asc."""
    return f"{prefix}{year}_{doy:03d}{suffix}"
