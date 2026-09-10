"""
BCM_Dailyv81_python.py — Python port of BCM Daily.

Copyright (c) 2026 Andrew Joros, Desert Research Institute
SPDX-License-Identifier: MIT

Usage:
    python BCM_Dailyv81_python.py

Set BCM_INDIR (input folder), BCM_CTL (control file, optional),
and BCM_OUT_DIR (output folder).
"""

from __future__ import annotations
import os
import sys
import time
import math
import numpy as np

# ── sibling modules ──────────────────────────────────────────────────────────
HERE    = os.path.dirname(os.path.abspath(__file__))
INDIR   = os.environ.get("BCM_INDIR", os.path.join(HERE, "..", "BCM_testrun_original"))
CTL     = os.environ.get("BCM_CTL", os.path.join(INDIR, "BCM_Dailyv81.ctl"))

sys.path.insert(0, HERE)
from bcm_io       import parse_ctl, read_asc, write_asc, read_inp, day_filename
from bcm_atmos_data import (PW_LAT, PW_LON, PW, BET_LAT, BET_LON, BET,
                             ALB_LAT, ALB_LON, ALB, idw_interpolate)

DR = math.pi / 180.0   # degrees → radians
RD = 180.0 / math.pi   # radians → degrees
PI = math.pi


# ============================================================================
# HELPERS
# ============================================================================
def p(path: str) -> str:
    """Resolve a filename relative to the input directory."""
    return os.path.join(
        os.environ.get("BCM_INDIR", os.path.join(HERE, "..", "BCM_testrun_original")),
        path,
    )


def is_leap(year: int) -> bool:
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def calendar_days(cfg):
    """(year, doy, leap) for each simulated day, including year wrap."""
    days = []
    for yn in range(cfg["yn1"], cfg["yn2"] + 1):
        leap = is_leap(yn)
        max_doy = 366 if leap else 365
        dn1 = cfg["dn1"] if yn == cfg["yn1"] else 1
        dn2 = cfg["dn2"] if yn == cfg["yn2"] else max_doy
        for dn in range(dn1, dn2 + 1):
            days.append((yn, dn, leap))
    return days


def parse_storms(spec: str) -> list:
    """'1-10' or '1,2,5' or '1-3,12' → storm numbers (1-based)."""
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def month_of_doy(doy: int, leap: bool) -> int:
    """Return month index (1-based) for a given day-of-year."""
    days = [31,29,31,30,31,30,31,31,30,31,30,31] if leap else \
           [31,28,31,30,31,30,31,31,30,31,30,31]
    d = 0
    for m, nd in enumerate(days, 1):
        d += nd
        if doy <= d:
            return m
    return 12


def seasonal_interp(dn: int, m1: np.ndarray, m2: np.ndarray,
                    mid1: float, span: float) -> np.ndarray:
    """Linear interpolation between two monthly values centred on mid-month."""
    frac = (mid1 - dn) / span
    return m2 - (m2 - m1) * frac


def kfactor_for_day(dn: int, kfactors: np.ndarray) -> np.ndarray:
    """
    kfactors: shape (12, nrows, ncols) — monthly Kv values (jan..dec).
    Returns kfactor interpolated to the given day-of-year, shape (nrows, ncols).
    Replicates the Fortran ELSE IF chain (non-leap year).
    """
    mid  = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
    span = [31, 28, 31,  30,  31,  30,  31,  31,  30,  31,  30,  30]
    # Month index pairs (0-based) and midpoints
    segs = [
        (dn <= 15,       0, 11, 15,  15),   # Dec→Jan early
        (dn > 349,       0, 11, 365, 15),   # Dec→Jan late
        (15 < dn <= 46,  1,  0,  46, 31),
        (46 < dn <= 74,  2,  1,  74, 28),
        (74 < dn <= 105, 3,  2, 105, 31),
        (105 < dn <= 135,4,  3, 135, 30),
        (135 < dn <= 166,5,  4, 166, 31),
        (166 < dn <= 196,6,  5, 196, 30),
        (196 < dn <= 227,7,  6, 227, 31),
        (227 < dn <= 258,8,  7, 258, 31),
        (258 < dn <= 288,9,  8, 288, 30),
        (288 < dn <= 319,10, 9, 319, 31),
        (319 < dn <= 349,11,10, 349, 30),
    ]
    for cond, ma, mb, mid_val, sp in segs:
        if cond:
            return kfactors[ma] - (kfactors[ma] - kfactors[mb]) * ((mid_val - dn) / sp)
    return kfactors[0]


# ============================================================================
# SOLAR RADIATION + PET (vectorised over the full grid per hour)
# ============================================================================
def compute_solar_pet(dn: int, hstep: int, cfg: dict, terrain: dict,
                      tmax: np.ndarray, tmin: np.ndarray,
                      bca: np.ndarray, bcb: np.ndarray, bcc: np.ndarray,
                      pt_alpha: np.ndarray, nodata: float,
                      pw_grid: np.ndarray, bet_grid: np.ndarray,
                      alb_grid: np.ndarray, oz_grid: np.ndarray) -> tuple:
    """
    Compute daily total solar radiation (MJ/m²) and PET (mm) for every cell.
    Returns (sumrad, sumtrn, petday, clouds) — all shape (nrows, ncols).
    """
    lat_rad = terrain["lat"] * DR      # shape (nrows, ncols)
    lon     = terrain["lon"]
    sl      = terrain["sl"]
    asp     = terrain["asp"]
    elev    = terrain["elev"]
    sky     = terrain["sky"]
    ridge   = terrain["ridge"]         # (nrows, ncols, 36)
    atp     = (tmax + tmin) * 0.5      # avg temp °C

    valid   = (elev != nodata) & (tmax != nodata) & (tmin != nodata)

    # ── solar geometry constants for this day ───────────────────────────────
    tau = 2 * PI * (dn - 1) / 365.0
    dec = (0.006918 - 0.399912*math.cos(tau) + 0.070257*math.sin(tau)
           - 0.006758*math.cos(2*tau) + 0.000907*math.sin(2*tau)
           - 0.002697*math.cos(3*tau) + 0.001480*math.sin(3*tau))
    ET  = (0.000075 + 0.001868*math.cos(tau) - 0.032077*math.sin(tau)
           - 0.014615*math.cos(2*tau) - 0.04089*math.sin(2*tau)) * 229.18
    EO  = 1.0 + 0.033 * math.cos(2*PI*dn/365.0)

    # Sunrise/sunset hour angles (radians) — per cell
    cos_arg  = -np.tan(lat_rad) * math.tan(dec)
    cos_arg  = np.clip(cos_arg, -1.0, 1.0)
    hasr     = np.arccos(cos_arg)              # hour angle of sunrise
    H1       = 12.0 - hasr * RD / 15.0        # sunrise local time
    H2       = 12.0 + hasr * RD / 15.0        # sunset  local time

    # Bristow-Campbell transmissivity
    tdiff    = tmax - tmin
    tdiff    = np.where(tdiff < 0.01, 0.01, tdiff)
    bctt     = bca * (1 - np.exp(-bcb * tdiff ** bcc))
    bctt     = np.clip(bctt, 0.0, 1.0)
    bcclouds = (-26.82*bctt**3 + 31.389*bctt**2 - 17.019*bctt + 12.602) / 10.0
    bcclouds = np.clip(bcclouds, 0.0, 1.0)

    # Atmospheric parameters (already interpolated per cell for this month)
    W    = pw_grid
    BETA = bet_grid
    PG   = alb_grid

    # ── hourly loop ─────────────────────────────────────────────────────────
    sum_rad = np.zeros_like(lat_rad)   # J/m²/day → accumulated
    sum_trn = np.zeros_like(lat_rad)
    sum_pet = np.zeros_like(lat_rad)

    for lst in range(1, 25, hstep):
        # Time correction (std = lon for California run)
        STD = lon
        CF  = (4 * (STD - lon) + ET) / 60.0
        T   = lst + CF
        HA  = 15.0 * (T - 12.0) * DR

        # Sun altitude
        sin_alt = (np.sin(dec) * np.sin(lat_rad)
                   + np.cos(dec) * np.cos(lat_rad) * np.cos(HA))
        sin_alt = np.clip(sin_alt, -1.0, 1.0)
        ALT     = np.arcsin(sin_alt)             # radians

        # Sun azimuth
        denom   = np.cos(ALT) * np.cos(lat_rad)
        denom   = np.where(np.abs(denom) < 1e-9, 1e-9, denom)
        az_arg  = (np.sin(ALT)*np.sin(lat_rad) - np.sin(dec)) / denom
        az_arg  = np.clip(az_arg, -1.0, 1.0)
        AZ      = np.arccos(az_arg)
        AZ      = np.where(HA < 0, PI - AZ, AZ + PI)
        ALT_deg = ALT * RD
        AZ_deg  = AZ  * RD
        ZENITH  = np.where(ALT_deg <= 0, 90.0, 90.0 - ALT_deg)

        # Atmospheric pressure at elevation
        PO   = 1013.25
        P    = np.exp(-0.0001184 * elev) * PO
        P0   = PO                             # pressure at beta measurement site

        TA   = 0.5*(tmax+tmin + (tmax-tmin)*np.cos(0.2618*(lst-15))) + 273.15
        # f90 L1643 resets P=PO before computing W, so the pressure ratio is 1.
        W_adj = W * (273.0 / TA)**0.5
        U1    = W_adj * (1.0 / (np.cos(ZENITH*DR) + 0.15*(93.885-ZENITH)**(-1.253) + 1e-9))
        U1    = np.clip(U1, 0, 500)
        WO    = 0.9
        ALPHA = 1.0

        AMSP  = (np.cos(ZENITH*DR) + 0.15*(93.885-ZENITH)**(-1.253))
        AMSP  = np.where(AMSP < 1e-9, 1e-9, AMSP)
        AMSP  = 1.0 / AMSP
        MA    = AMSP * (P / PO)
        U3    = oz_grid * AMSP   # per-cell ozone (f90 U3 = OZONE(k)*AMSP)

        TRR   = np.exp(-0.0903 * MA**0.84 * (1 + MA - MA**1.01))
        ABO   = (0.1611*U3*((1+139.48*U3)**(-0.3035))
                 - 0.002715*U3*((1+0.044*U3+0.0003*U3**2)**(-1)))
        TRO   = 1.0 - ABO
        TRG   = np.exp(-0.0127 * MA**0.26)
        ABW   = 2.4959*U1*((1+79.034*U1)**0.6828 + 6.385*U1)**(-1)
        TRW   = 1.0 - ABW
        TRA   = (0.12445*ALPHA - 0.0162) + (1.003 - 0.125*ALPHA)*np.exp(
                -BETA*ALPHA*MA*(1.089+0.5120001))
        TRAA  = 1.0 - (1-WO)*(1-MA+MA**1.06)*(1-TRA)
        TRAS  = TRA / np.where(TRAA < 1e-9, 1e-9, TRAA)

        FC    = 1.0
        ISC   = 1367.0
        cos_z = np.cos(ZENITH*DR)

        IDA   = EO*0.79*ISC*cos_z*TRO*TRG*TRW*TRAA*FC*(1-TRAS)/(1-MA+MA**1.02+1e-9)
        IDR   = EO*0.79*ISC*cos_z*TRO*TRG*TRW*TRAA*0.5*(1-TRR)/(1-MA+MA**1.02+1e-9)
        INN   = EO*0.9751*ISC*TRR*TRO*TRG*TRW*TRA
        IB    = INN * cos_z
        PAP   = 0.0685 + (1-FC)*(1-TRAS)
        IDM   = (IB+IDA+IDR) * (PG*PAP / (1-PG*PAP+1e-9))

        # Night time → zero
        night = ZENITH >= 90.0
        IB    = np.where(night, 0.0, IB)
        INN   = np.where(night, 0.0, INN)
        IDR   = np.where(night, 0.0, IDR)
        IDA   = np.where(night, 0.0, IDA)
        IDM   = np.where(night, 0.0, IDM)

        # Ridge shading of direct beam
        # SUNAZ is the azimuth bin (Fortran: 0-36); ridge is now (nrows,ncols,37)
        SUNAZ = np.floor(AZ_deg / 10.0 + 0.5).astype(int)
        SUNAZ = np.clip(SUNAZ, 0, 36)   # 37 bins: 0..36
        ridge_h = np.take_along_axis(ridge, SUNAZ[:, :, np.newaxis], axis=2)[:, :, 0]
        shade_beam  = (90 - ridge_h) <= ZENITH
        shade_diffp = (90 - ridge_h - 4.0) <= ZENITH
        IB    = np.where(shade_beam,  0.0,               IB)
        IDA   = np.where(shade_diffp, (1-1.0)*IDA,       IDA)  # csr≈1

        ID    = IDR + IDA + IDM
        ISS   = ID * sky
        IR    = (ID + IB) * PG * (1 - sky)
        IDHORIZ = ISS + IR
        IBHORIZ = IB

        # Slope incidence angle
        cos_theta = (np.cos(sl*DR)*np.cos(ZENITH*DR)
                     + np.sin(sl*DR)*np.sin(ZENITH*DR)*np.cos((AZ_deg-asp)*DR))
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        IBSLOPE   = IBHORIZ * cos_theta / (cos_z + 1e-9)
        RNSLOPE   = np.maximum(0.0, IBSLOPE + IDHORIZ)

        # Radiation J/m² per timestep
        RAD = RNSLOPE * hstep * 3600.0 * bctt

        # Net radiation for PET
        Eac    = (1 - 0.84*bcclouds)*0.0000092*TA**2 + 0.84*bcclouds
        NetShort = RAD * (1 - PG)
        NetLong  = 5.6697e-8 * (Eac - 0.98) * TA**4 * (hstep*3600.0)
        RN       = NetShort + NetLong
        # Soil heat flux: DTA = TA - TAP, TAP = tmin+273.15 fixed for the day
        # (f90 L1586/1634/1733 — TAP is NOT updated hourly).
        DTA_h    = TA - (tmin + 273.15)
        GH       = 2.1e6 * 0.01 / 24.0 * DTA_h / (hstep/24.0)

        # SSG = Δ/(Δ+γ), the 6th-order polynomial in TA (Kelvin) exactly as the
        # .f90 source (L1738). Verified against the recompiled source-matched exe:
        # it is well-behaved and positive at California temps (≈0.56 at 283 K) and
        # reproduces the exe's PET. (The earlier quadratic + ALPHA 0.95 was a
        # reverse-engineering guess made before we had a source-matched exe; it
        # cancelled against a broken soil-heat-flux term and only matched on the
        # basin mean, not per-cell — see CHANGELOG 2026-07-23.)
        SSG = (7.1300275e-12*TA**6 - 1.117951e-8*TA**5 + 7.2538824e-6*TA**4
               - 2.4939885e-3*TA**3 + 0.47946782*TA**2 - 48.898875*TA + 2067.9946)
        LHV = (2500.8 - 2.36*atp + 0.0016*atp**2 - 0.00006*atp**3) * 1e3
        # Priestley-Taylor alpha = 1.26 (f90 hardcodes varalpha(k)=1.26 at L1740,
        # overriding the PT_Alpha map). Overridable via BCM_PT_ALPHA for sweeps.
        ALPHA_PT = cfg.get("pt_alpha", 1.26)
        PET_hr = ALPHA_PT * SSG * (RN - GH) / np.where(LHV < 1e3, 1e3, LHV)
        PET_hr = np.maximum(0.0, PET_hr)

        sum_rad += RAD
        sum_trn += RN
        sum_pet += PET_hr

    # Set nodata
    nd_mask = ~valid
    sum_rad[nd_mask] = nodata
    sum_pet[nd_mask] = nodata
    bcclouds[nd_mask] = nodata

    # Convert J/m² → MJ/m²
    sum_rad_mj = np.where(sum_rad != nodata, sum_rad / 1e6, nodata)

    return sum_rad_mj, sum_trn / 1e6, sum_pet, bcclouds


# ============================================================================
# SNOW MODEL (vectorised)
# ============================================================================
def snow_step(dn: int, ppt: np.ndarray, tmax: np.ndarray, tmin: np.ndarray,
              petday: np.ndarray, elev: np.ndarray,
              prior_pack: np.ndarray, prior_ati: np.ndarray,
              prior_hdi: np.ndarray, prior_mwt: np.ndarray,
              snowaccum_t: np.ndarray, mfmax: np.ndarray, mfmin: np.ndarray,
              cfg: dict, sumrad_mj: np.ndarray, nodata: float) -> dict:
    """Vectorised Snow-17 style snowmelt for one day."""
    tipm     = cfg["tipm"]
    nmf      = cfg["nmf"]
    maf      = cfg["maf"]
    solar_flag = cfg["solar_flag"]
    plwhc    = 0.4     # liquid water holding capacity fraction

    tave  = 0.5*(tmax + tmin) - snowaccum_t
    tdiff = tmax - tmin
    tdiff = np.where(tdiff < 0.01, 0.01, tdiff)

    # Rain/snow partitioning — exact Fortran logic (f90 lines 1756-1765):
    #   IF tmin<=snowaccum AND tmax<=snowaccum → all snow
    #   ELIF tmin<=snowaccum AND tmax>snowaccum → mixed
    #   ELSE (both above snowaccum) → all rain
    all_snow  = ((tmin - snowaccum_t) <= 0) & ((tmax - snowaccum_t) <= 0)
    mixed     = ((tmin - snowaccum_t) <= 0) & ((tmax - snowaccum_t) >  0)
    # all_rain = neither all_snow nor mixed (tmin > snowaccum)
    snow_frac = np.where(all_snow, 1.0,
                 np.where(mixed,   -(tmin - snowaccum_t) / tdiff, 0.0))
    rain_frac = 1.0 - snow_frac

    snowfall = ppt * snow_frac
    pack_day = snowfall + prior_pack

    # Sublimation (equation-based)
    subl_par = (0.1309 * petday * 30 + 3.9669) / 30.0
    subl_day = np.zeros_like(pack_day)
    no_snow_today = snowfall == 0.0
    enough_pack   = pack_day >= subl_par
    some_pack     = (pack_day > 0) & (pack_day < subl_par)
    subl_day = np.where(no_snow_today & enough_pack, subl_par,
               np.where(no_snow_today & some_pack,  pack_day, 0.0))
    pack_day = np.where(no_snow_today & enough_pack, pack_day - subl_par,
               np.where(no_snow_today & some_pack,  0.0,        pack_day))

    # Melt factor (sine function)
    Sv       = 0.5 * np.sin((dn - 80) * 2*PI/366) + 0.5
    mf       = Sv * (mfmax - mfmin) + mfmin
    nmfac    = 24/6 * (mf / np.where(mfmax < 1e-9, 1e-9, mfmax)) * nmf

    # Solar temperature enhancement
    if solar_flag:
        bins = [0, 7.5, 10.5, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0]
        vals = [0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50]
        sol_tmp = np.zeros_like(sumrad_mj)
        for i in range(len(bins)-1):
            sol_tmp = np.where((sumrad_mj >= bins[i]) & (sumrad_mj < bins[i+1]),
                               vals[i], sol_tmp)
        sol_tmp = np.where(sumrad_mj >= 30.0, 2.50, sol_tmp)
    else:
        sol_tmp = np.zeros_like(pack_day)

    # ATI and HDI update
    ati    = prior_ati.copy()
    hdi    = prior_hdi.copy()
    mwat   = prior_mwt.copy()
    melt_day = np.zeros_like(pack_day)

    has_pack = pack_day > 0

    # ATI update
    many_snow  = has_pack & (snowfall >= 36)
    few_snow   = has_pack & (snowfall <  36)
    ati = np.where(many_snow, ati + 4*tave, ati)
    ati = np.where(few_snow,  ati + tipm*((tave+sol_tmp) - ati)*4, ati)
    ati = np.minimum(ati, 0.0)

    # HDI update
    delhdi = nmfac * (ati - (tave + sol_tmp))
    hdi    = np.where(has_pack, hdi + delhdi, hdi)
    hdi    = np.maximum(hdi, 0.0)

    # Melt
    tppt = np.maximum(tave, 0.0)
    esat = 2.7489e8 * np.exp(-4278.63 / (tave + 242.792 + 1e-9))
    pamb = 1013.25 * (1 - 2.25577e-5 * np.maximum(elev, 0.0))**5.25588
    uadj = 0.0864

    ros = rain_frac > cfg.get("rvs_percent", 0.3)
    # Non-rain-on-snow melt
    melt_nros = mf * (tave + sol_tmp) * 4 + 0.0125 * ppt * rain_frac * tppt
    # Rain-on-snow melt (Anderson eq 5)
    melt_ros  = (6.12e-10*24*((tave+273)**4 - 273**4)
                 + 0.0125*ppt*rain_frac*tppt
                 + 8.5*uadj*4*(0.9*esat - 6.11)
                 + 0.00057*pamb*tave)

    melt_day = np.where(has_pack & (hdi <= 0),
                        np.where(ros, melt_ros, melt_nros), 0.0)
    melt_day = np.maximum(melt_day, 0.0)

    # Liquid water retention
    liqwmx = plwhc * pack_day
    mwat   = np.where(has_pack, mwat + melt_day + ppt*rain_frac, mwat)
    # If mwat <= liqwmx, no release yet
    release = np.maximum(mwat - liqwmx, 0.0)
    melt_day = np.where(has_pack, release, melt_day)
    mwat   = np.where(has_pack, np.minimum(mwat, liqwmx), mwat)

    # Update pack
    pack_day = pack_day - melt_day
    too_low  = pack_day < 0
    melt_day = np.where(too_low, melt_day + pack_day, melt_day)
    pack_day = np.maximum(pack_day, 0.0)

    # If pack gone, reset state
    gone = pack_day == 0
    mwat = np.where(gone, 0.0, mwat)
    hdi  = np.where(gone, 0.0, hdi)
    ati  = np.where(gone, 0.0, ati)

    # Refreezing melt water
    refreeze  = mwat >= hdi
    pack_day  = np.where(refreeze,  pack_day + hdi, pack_day + mwat)
    mwat      = np.where(refreeze,  mwat - hdi,     0.0)
    hdi       = np.where(refreeze,  0.0,             np.maximum(hdi - mwat, 0.0))

    # rvsday OUTPUT = SNOW fraction, to match the compiled .exe reference, which
    # writes the snow fraction. (The .f90 source we hold instead assigns
    # rvsday = rainfrac — a source/exe conflict confirmed by the validation
    # harness: py rain_frac == 1 - exe_rvs exactly on clear days. rain_frac is
    # still used above for the melt/rain-on-snow equations; only the written
    # diagnostic grid changes here.)
    rvs_snow = snow_frac.copy()

    # Nodata propagation
    nd = (ppt == nodata) | (tmax == nodata) | (tmin == nodata) | (elev == nodata)
    for arr in [snowfall, pack_day, melt_day, subl_day, rvs_snow, ati, hdi, mwat]:
        arr[:] = np.where(nd, nodata, arr)

    return {
        "snowfall": snowfall, "packday": pack_day, "meltday": melt_day,
        "sublday":  subl_day, "rvsday":  rvs_snow,
        "ati": ati, "hdi": hdi, "mwatday": mwat,
    }


# ============================================================================
# SOIL WATER BALANCE (vectorised)
# ============================================================================
def water_balance(ppt: np.ndarray, melt: np.ndarray, snow: np.ndarray,
                  petday: np.ndarray, prior_str: np.ndarray,
                  soild_eff: np.ndarray, wpmm: np.ndarray, fcmm: np.ndarray,
                  pormm: np.ndarray, geolks: np.ndarray, kfac: np.ndarray,
                  nodata: float, rchrun_flag: int,
                  rchrun_limit: float, drydown: np.ndarray = None,
                  soil_nd: np.ndarray = None,
                  imperv: np.ndarray = None,
                  zero_ks: np.ndarray = None,
                  water_veg: np.ndarray = None,
                  imperv_release: np.ndarray = None) -> dict:
    """Vectorised soil water balance for one day.

    Replicates the Fortran "big IF" storage/AET branch logic exactly
    (BCM_Dailyv81.f90 lines 1981-2042).  The key correction vs the earlier
    Python port: when storage falls below wilting point the Fortran does NOT
    raise it back to wpmm — it lets the soil dry below wp (down to 0 via the
    drydown term).  The old `stor_day = max(stor_day, wpmm)` clamp inflated
    dry, high-capacity (deep-root) cells and was the main source of the
    ~106-130 mm str error.
    """
    pptmm    = ppt.copy()
    stor_in  = prior_str + pptmm + melt - snow

    # Runoff above porosity
    runoff   = np.where(stor_in > pormm, stor_in - pormm, 0.0)
    stor_in  = np.minimum(stor_in, pormm)

    # Potential veg AET (LAI=1 for now; full LAI in a later sprint)
    vegaet0  = petday * kfac
    vegaet0  = np.minimum(vegaet0, petday)
    vegaet0  = np.maximum(vegaet0, 0.0)

    if drydown is None:
        drydown = np.zeros_like(stor_in)   # drydownflag=0 → no extra dry-down

    # ── Fortran branch masks (exactly one holds per cell) ──────────────────
    A   = stor_in >= fcmm
    A1  = A & ((stor_in - wpmm) >= vegaet0)
    A2  = A & ~((stor_in - wpmm) >= vegaet0)
    B   = ~A
    B1  = B & (stor_in >= wpmm)
    B1a = B1 & ((stor_in - wpmm) >= vegaet0)
    B1b = B1 & ~((stor_in - wpmm) >= vegaet0)
    B2  = B & (stor_in < wpmm)

    # Branch A1: enough water for full AET; recharge excess above fc
    stor1_A1   = fcmm - vegaet0
    cl_A1      = stor1_A1 < wpmm
    vegaet_A1  = np.where(cl_A1, vegaet0 - (wpmm - stor1_A1), vegaet0)
    stor_A1    = np.where(cl_A1, wpmm, stor1_A1)
    rech_A1    = stor_in - fcmm

    # Branch A2: at/above fc but not enough headroom for full AET → stor=wpmm
    vegaet_A2  = stor_in - wpmm
    stor_A2    = np.full_like(stor_in, 0.0) + wpmm

    # Branch B1a: below fc, above wp, enough water for full AET
    stor1_B1a  = stor_in - vegaet0
    cl_B1a     = stor1_B1a < wpmm
    vegaet_B1a = np.where(cl_B1a, vegaet0 - (wpmm - stor1_B1a), vegaet0)
    stor_B1a   = np.where(cl_B1a, wpmm, stor1_B1a)

    # Branch B1b: below fc, above wp, not enough water → AET draws to wp
    vegaet_B1b = stor_in - wpmm
    stor_B1b   = np.full_like(stor_in, 0.0) + wpmm

    # Branch B2: below wp — soil dries beyond wp (drydown), NOT clamped to wp
    vegaet_B2  = vegaet0
    stor_B2    = stor_in - drydown * petday
    stor_B2    = np.where(stor_B2 < 0.0, 0.0, stor_B2)

    conds      = [A1, A2, B1a, B1b, B2]
    stor_day   = np.select(conds, [stor_A1, stor_A2, stor_B1a, stor_B1b, stor_B2],
                           default=stor_in)
    vegaet     = np.select(conds, [vegaet_A1, vegaet_A2, vegaet_B1a, vegaet_B1b, vegaet_B2],
                           default=vegaet0)
    # Full-precision carryover. The shipped exe floors str to an integer each
    # day (an artifact; Michelle: do not keep that). We do not copy it.
    recharge   = np.where(A1, rech_A1, 0.0)

    cwd       = petday - vegaet
    cwd       = np.maximum(cwd, 0.0)

    # Recharge capped by bedrock Ks
    rch3      = np.minimum(recharge, geolks)
    run2      = runoff + np.maximum(recharge - geolks, 0.0)

    # Apply rchrun scaler flag
    if rchrun_flag:
        lim = rchrun_limit
        shift  = np.where(run2 > 0, np.minimum(run2, lim), 0.0)
        rch3  += shift
        run2  -= shift
        run2   = np.maximum(run2, 0.0)

    # No-soil cells (soild_eff==0)
    no_soil   = soild_eff == 0.0
    runoff_ns = pptmm + melt - snow
    run2  = np.where(no_soil, runoff_ns, run2)
    vegaet = np.where(no_soil, 0.0, vegaet)
    rch3  = np.where(no_soil, 0.0, rch3)

    exc    = np.maximum(pptmm - vegaet, 0.0)
    smd    = fcmm - stor_day
    smr    = pormm - stor_day

    # Impervious cells: wp==0 (valid soil depth) → the exe zeroes the soil
    # (f90 L1932: soild=wp=fc=por=0) so porosity capacity is 0 and storday-pormm
    # flushes ALL water (incl. any prior/initial storage) to runoff, str→0,
    # no AET/recharge. BUT the exe then rewrites soild = root_depth + 0
    # (L1938-1942) and the L2035 override `IF(soild==0) runoff=ppt+melt-snow`
    # fires ONLY when root_depth==0. Hence, verified across all 265 such cells:
    #   root_depth  > 0 → soild≠0 → override skipped → run = prior_str+ppt+melt-snow
    #   root_depth == 0 → soild==0 → override fires   → run = ppt+melt-snow
    # (releasing prior_str conserves mass; the root==0 case does not — an exe
    # quirk we reproduce for fidelity.) `imperv_release` = imperv & root_depth>0.
    if imperv is not None:
        base_run = pptmm + melt - snow
        if imperv_release is not None:
            run_imp = np.where(imperv_release, prior_str + base_run, base_run)
        else:
            run_imp = base_run
        stor_day = np.where(imperv, 0.0, stor_day)
        vegaet   = np.where(imperv, 0.0, vegaet)
        rch3     = np.where(imperv, 0.0, rch3)
        run2     = np.where(imperv, run_imp, run2)
        cwd      = np.where(imperv, petday, cwd)
        exc      = np.where(imperv, np.maximum(pptmm, 0.0), exc)
        smd      = np.where(imperv, 0.0, smd)
        smr      = np.where(imperv, 0.0, smr)

    # ── Water / no-geology masking (mirror BCM_Dailyv81.f90) ─────────────────
    # (1) geolks==0 → geology "Water" (58) / "No Geology" (54). The Fortran sets
    #     soild=0 (L1929) and, since maxrch=geolks=0, forces run/rch/aet=0
    #     (L2067-2072). soild_eff is already 0 here (set in main), so storday,
    #     smd, smr collapse to 0 via the no-soil path; we additionally clamp
    #     run2/rch3/vegaet to 0 so rejected recharge never becomes runoff.
    if zero_ks is not None:
        stor_day = np.where(zero_ks, 0.0, stor_day)
        run2     = np.where(zero_ks, 0.0, run2)
        rch3     = np.where(zero_ks, 0.0, rch3)
        vegaet   = np.where(zero_ks, 0.0, vegaet)
        smd      = np.where(zero_ks, 0.0, smd)
        smr      = np.where(zero_ks, 0.0, smr)
    # (2) WHR vegetation "Water" (id 57): Fortran zeroes every soil term
    #     (L2073-2085); snow terms are zeroed in the driver.
    if water_veg is not None:
        stor_day = np.where(water_veg, 0.0, stor_day)
        run2     = np.where(water_veg, 0.0, run2)
        rch3     = np.where(water_veg, 0.0, rch3)
        vegaet   = np.where(water_veg, 0.0, vegaet)
        cwd      = np.where(water_veg, 0.0, cwd)
        smd      = np.where(water_veg, 0.0, smd)
        smr      = np.where(water_veg, 0.0, smr)

    # Nodata mask. The Fortran skips any cell missing a required input
    # (source line 1907), including soil depth — so water bodies / no-soil
    # cells (e.g. geology "Water", rock ID 58, which also have soild=nodata)
    # are written as nodata, NOT as zeros. Mirror that with soil_nd here.
    nd = ((ppt == nodata) | (petday == nodata) |
          (prior_str == nodata) | (geolks == nodata))
    if soil_nd is not None:
        nd = nd | soil_nd
    for arr in [stor_day, run2, rch3, exc, vegaet, cwd, smd, smr]:
        arr[:] = np.where(nd, nodata, arr)

    return {
        "storday": stor_day, "run2day": run2, "rch3day": rch3,
        "exc1day": exc,      "vegaet":  vegaet, "cwd1day": cwd,
        "smdday":  smd,      "smrday":  smr,
    }


# ============================================================================
# BASIN AGGREGATION
# ============================================================================
def basin_stats(basin_grid: np.ndarray, elev: np.ndarray,
                nodata: float, cellsize: float, **fields) -> dict:
    """
    Compute the per-basin summary averages, faithfully mirroring the Fortran
    (BCM_Dailyv81.f90 L2160-2377):

      * ``totcells`` = number of basin cells with valid elevation (and area),
        i.e. the DO-3001 accumulation guard ``elev != nodata`` (L2161).
        Used for ``area``/``rchacft``/``runacft`` (L2375-2377).
      * Every field's SUM ignores nodata cells, and each field is divided by
        **its own valid-cell count** (``pptcells``, ``tmaxcells``,
        ``tmincells``, ... ; L2352-2373) — NOT by a single shared denominator.
      * The lone exception is soil storage (``str``): the Fortran sums valid
        storage but divides by ``totcells`` (``avgstor/totcells``, L2369), so
        cells with nodata storage still dilute the basin mean.

    Prior versions divided every field (and area) by the count of ALL cells
    tagged with the basin id — including nodata-elevation cells. For basins
    whose footprint spans a large nodata region (e.g. basin 20 here, ~2.77x
    more tagged cells than valid) that diluted tmnc/tmxc/pptmm/str and inflated
    aream; basins with no interior nodata were unaffected.

    fields: name → 2-D array
    """
    basin_ids = np.unique(basin_grid)
    basin_ids = basin_ids[basin_ids != nodata]
    stats = {}
    for bid in basin_ids:
        bid       = int(bid)
        all_mask  = basin_grid == bid            # all cells tagged this basin
        good_mask = all_mask & (elev != nodata)  # Fortran totcells set (L2161)
        totcells  = int(good_mask.sum())
        if totcells == 0:
            continue
        row = {"n_cells": totcells,
               "area_m2": totcells * cellsize * cellsize}
        for name, arr in fields.items():
            valid = good_mask & (arr != nodata)
            # str divides by totcells (avgstor/totcells); all others by their
            # own valid-cell count (avg<field>/<field>cells).
            denom = totcells if name == "str" else int(valid.sum())
            row[name] = (arr[valid].sum() / denom) if denom > 0 else 0.0
        stats[bid] = row
    return stats


# ============================================================================
# MAIN DRIVER
# ============================================================================
def main():
    t0 = time.time()
    # Diagnostic hooks (no effect unless env vars are set):
    #   BCM_RESEED_STR_DIR : each day, overwrite prior_str with the Fortran
    #                        previous-day str grid from this directory. Isolates
    #                        whether the soil physics tracks once the storage
    #                        initial/spatial field matches the Fortran exe.
    #   BCM_OUT_DIR        : write outputs here instead of next to this script.
    RESEED_DIR = os.environ.get("BCM_RESEED_STR_DIR")
    OUT_DIR    = os.environ.get("BCM_OUT_DIR", HERE)
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 60)
    print("BCM Daily v8.1 — Python port")
    if RESEED_DIR:
        print(f"  [DIAG] reseeding prior_str daily from: {RESEED_DIR}")
    print("=" * 60)

    # ── parse control file ──────────────────────────────────────────────────
    ctl = os.environ.get(
        "BCM_CTL",
        os.path.join(
            os.environ.get("BCM_INDIR", os.path.join(HERE, "..", "BCM_testrun_original")),
            "BCM_Dailyv81.ctl",
        ),
    )
    indir = os.environ.get("BCM_INDIR", os.path.join(HERE, "..", "BCM_testrun_original"))
    print(f"\nReading control file: {ctl}")
    cfg = parse_ctl(ctl)
    # ── Calibration knobs (defaults preserve current behaviour) ──────────────
    # Exposed for the validation harness to sweep the open-items (see
    # Calibration_OpenItems_Memo_2026-07-09.md). None of these change defaults.
    cfg["pt_alpha"] = float(os.environ.get("BCM_PT_ALPHA", "1.26"))   # PET alpha (f90 varalpha=1.26)
    if os.environ.get("BCM_RCHRUN_DISABLE", "0") == "1":             # RCH/RUN scaler
        cfg["rchrun_flag"] = 0
    # SNOW=0: use snowaccum/mfmax/mfmin maps (Fortran .f90 always reads the maps).
    cfg["snow_flag"] = 0
    if cfg["pt_alpha"] != 0.95 or not cfg["rchrun_flag"]:
        print(f"  [DIAG] pt_alpha={cfg['pt_alpha']}, rchrun_flag={cfg['rchrun_flag']}")
    print(f"  snow_flag={cfg['snow_flag']} (0=maps)")
    print(f"  Simulation: DOY {cfg['dn1']}-{cfg['dn2']} of {cfg['yn1']}-{cfg['yn2']}")
    print(f"  Inputs -> {indir}")
    print(f"  Output -> {OUT_DIR}")

    # ── domain from DEM ─────────────────────────────────────────────────────
    print("\nLoading static grids…")
    hdr, dem  = read_asc(p(cfg["demfile"]))
    nodata    = hdr["nodata_value"]
    nrows     = int(hdr["nrows"])
    ncols     = int(hdr["ncols"])
    cellsize  = hdr["cellsize"]
    print(f"  Grid: {nrows} rows × {ncols} cols  ({nrows*ncols:,} cells)")

    _, soild    = read_asc(p(cfg["soildfile"]))
    _, wp       = read_asc(p(cfg["wpfile"]))
    _, fc       = read_asc(p(cfg["fcfile"]))
    _, por      = read_asc(p(cfg["porfile"]))
    _, geolid   = read_asc(p(cfg["geolfile"]))
    _, veg_grid = read_asc(p(cfg["vegfile"]))
    _, bca      = read_asc(p(cfg["bc_a_file"]))
    _, bcb      = read_asc(p(cfg["bc_b_file"]))
    _, bcc_arr  = read_asc(p(cfg["bc_c_file"]))
    _, pt_alpha = read_asc(p(cfg["pt_alpha_file"]))
    _, snowaccum_map = read_asc(p(cfg["snowaccumfile"]))
    _, mfmax_map     = read_asc(p(cfg["mfmaxfile"]))
    _, mfmin_map     = read_asc(p(cfg["mfminfile"]))
    _, basin_grid    = read_asc(p(cfg["areafile"]))
    print("  Static grids loaded.")

    # snow_flag=1 would replace maps with CTL scalars (1.5 / 1.8 / 0.1).
    # Forced off above so we stay on maps, same as Fortran.
    if cfg["snow_flag"]:
        snowaccum_map = np.full_like(dem, cfg["snowaccum_t"])
        mfmax_map     = np.full_like(dem, cfg["maxmf"])
        mfmin_map     = np.full_like(dem, cfg["minmf"])

    # ── bedrock Ks from lookup table ────────────────────────────────────────
    rktbl  = cfg["rockks_table"]
    geolks = np.full_like(dem, nodata)
    valid_geol = (geolid != nodata)
    for gid, ks in rktbl.items():
        geolks = np.where((geolid == gid) & valid_geol, ks, geolks)
    # Zero-Ks geology = "Water" (id 58) / "No Geology" (id 54). The Fortran
    # (BCM_Dailyv81.f90 L1929-1931) forces soild=0 on these cells and, because
    # maxrch=geolks=0, then zeroes run/rch/aet (L2067-2072). Without this guard
    # the rejected-recharge term (recharge - geolks) dumps all would-be recharge
    # into runoff, producing SPURIOUS RUNOFF over water/no-geology cells.
    zero_ks = (geolks == 0.0) & valid_geol

    # ── vegetation kfactor lookup ───────────────────────────────────────────
    vp = cfg["veg_params"]
    # kfactor: shape (12, nrows, ncols)
    kfactor_monthly = np.full((12, nrows, ncols), nodata)
    for vid, vdata in vp.items():
        mask = veg_grid == vid
        for m in range(12):
            kfactor_monthly[m] = np.where(mask, vdata["kfactor"][m],
                                           kfactor_monthly[m])

    # ── terrain file ────────────────────────────────────────────────────────
    # BCM_INP_FILE overrides the .ctl's inpfile (default-preserving) so alternate
    # terrain decks (e.g. a Skyview-regenerated .inp) can be tested without
    # touching the reference inputs.
    inp_path = os.environ.get("BCM_INP_FILE", p(cfg["inpfile"]))
    terrain = read_inp(inp_path, nrows, ncols)

    # ── pre-compute per-cell monthly atmospheric parameters ─────────────────
    # IDW interpolation for all 12 months at once (shape 12 × nrows × ncols)
    print("  Pre-computing atmospheric grids (IDW)…")
    pw_monthly  = np.stack([idw_interpolate(terrain["lat"], terrain["lon"],
                                             PW_LAT,  PW_LON,  PW[m])
                             for m in range(12)])
    bet_monthly = np.stack([idw_interpolate(terrain["lat"], terrain["lon"],
                                             BET_LAT, BET_LON, BET[m])
                             for m in range(12)])
    alb_monthly = np.stack([idw_interpolate(terrain["lat"], terrain["lon"],
                                             ALB_LAT, ALB_LON, ALB[m])
                             for m in range(12)])
    # CSR correction
    csr_monthly = 0.0264*pw_monthly**3 - 0.1238*pw_monthly**2 - 0.0347*pw_monthly + 1
    # Per-cell monthly ozone from the lat (degrees) regressions (f90 L1379-1390).
    _lat_deg = terrain["lat"]
    _oz_coef = [(-5e-7, 5e-5,  .0007, .2182), (-8e-7, 8e-5,  .0003, .2231),
                (-9e-7, 1e-4, -.0007, .2323), (-7e-7, 7e-5,  .0008, .2311),
                (-6e-7, 6e-5,  .0009, .2414), (-5e-7, 4e-5,  .0012, .2368),
                (-5e-7, 4e-5,  .0008, .2380), (-2e-7, 6e-6,  .0018, .2259),
                (-2e-7, 8e-6,  .0014, .2283), (-1e-7, -2e-7, .0015, .2185),
                (-2e-7, 1e-5,  .0014, .2182), (-2e-7, 1e-5,  .0015, .2174)]
    oz_monthly = np.stack([a*_lat_deg**3 + b*_lat_deg**2 + c*_lat_deg + d
                           for (a, b, c, d) in _oz_coef])
    print("  Atmospheric grids ready.")

    # ── soil properties (mm) ─────────────────────────────────────────────────
    veg_root = np.zeros_like(dem)
    for vid, vdata in vp.items():
        veg_root = np.where(veg_grid == vid, vdata["root_depth"], veg_root)

    # The April .f90 source adds vegsoil to soild inside the daily column loop
    # (v81 line 1940 / v86 line 1748) before computing fcmm/wpmm/pormm.
    # However, the July 2025 .exe appears NOT to augment storage CAPACITY this
    # way: a reseed diagnostic (2026-06-11) showed Python holds ~2x the water of
    # the exe specifically in high veg-root cells (35% of cells, 98.4% of the
    # total str error), while aet/rch/run match. Use BCM_SOIL_AUG to test:
    #   "full" (default) : soild_eff = soild + veg_root   (matches the .f90 source)
    #   "none"           : soild_eff = soild              (raw depth; test vs exe)
    SOIL_AUG = os.environ.get("BCM_SOIL_AUG", "full").lower()
    if SOIL_AUG == "none":
        soild_eff = soild.copy()
    else:
        soild_eff = soild + veg_root               # augmented (soild + vegsoil)
    if SOIL_AUG != "full":
        print(f"  [DIAG] soil-depth augmentation = {SOIL_AUG}")
    soild_eff = np.where(soild_eff < 0, 0.0, soild_eff)
    # Fortran L1929-1931: bedrock Ks==0 (water / no-geology) forces soil depth
    # to 0 so the cell carries no soil storage, AET, or recharge.
    soild_eff = np.where(zero_ks, 0.0, soild_eff)
    # Cells with no soil-depth data (e.g. water bodies) — Fortran writes these
    # as nodata in every water-balance grid; we blank them the same way.
    soil_nd   = (soild == nodata)
    # Cells with wilting point == 0 (but valid soil depth) are treated by the
    # exe as impervious / zero-storage (str=aet=rch=smd=smr=0, all precip→runoff).
    imperv    = (wp == 0.0) & (soild != nodata)
    # Among impervious cells, only those with veg root depth > 0 release their
    # prior/initial storage as runoff (see water_balance for the exe mechanism).
    imperv_release = imperv & (veg_root > 0.0)
    # WHR vegetation "Water" (id 57): Fortran zeroes every soil/snow term for
    # these cells (BCM_Dailyv81.f90 L2073-2085), a second water mask independent
    # of the geology map.
    water_veg = (veg_grid == 57)
    wpmm      = wp  * soild_eff * 1000.0
    fcmm      = fc  * soild_eff * 1000.0
    pormm     = por * soild_eff * 1000.0

    # ── initialize state ────────────────────────────────────────────────────
    if cfg["prior"] == 0:
        prior_pack = np.zeros_like(dem)
        prior_ati  = np.zeros_like(dem)
        prior_hdi  = np.zeros_like(dem)
        prior_mwt  = np.zeros_like(dem)
        # Initial storage. BCM_STR_INIT selects the formula:
        #   "v81src" (default): istor*fc*(soild+vegsoil)*1000 (BCM_Dailyv81 line
        #                       1312, the folder .f90 / recompiled exe). The
        #                       corrected model WITH veg-root depth in capacity.
        #   "v86"             : istor*(fc-wp)*soild_RAW*1000  (BCM_Dailyv86 line
        #                       1136 — no veg-root in the init). Domain-avg ≈ 207 mm,
        #                       reproduces the OLDER shipped exe (no veg root depth).
        #   "legacy"          : istor*(fcmm_aug - wpmm_aug)  (old port; ≈ 355 mm).
        STR_INIT  = os.environ.get("BCM_STR_INIT", "v81src").lower()
        soild_raw = np.where(soild < 0, 0.0, soild)
        if STR_INIT == "legacy":
            init_str = cfg["istor"] * (fcmm - wpmm)
        elif STR_INIT == "v81src":
            init_str = cfg["istor"] * fc * soild_eff * 1000.0
        else:
            init_str = cfg["istor"] * (fc - wp) * soild_raw * 1000.0
        if STR_INIT != "legacy":
            print(f"  [DIAG] str init formula = {STR_INIT}")
        prior_str = np.where(dem != nodata, init_str, nodata)
        prior_lai  = np.ones_like(dem)
    else:
        yn_p = cfg["yn1"]; dn_p = cfg["dn1"] - 1
        if dn_p == 0:
            yn_p -= 1; dn_p = 366 if is_leap(yn_p) else 365
        _, prior_pack = read_asc(p(day_filename("pck", yn_p, dn_p)))
        _, prior_ati  = read_asc(p(day_filename("ati", yn_p, dn_p)))
        _, prior_hdi  = read_asc(p(day_filename("hdi", yn_p, dn_p)))
        _, prior_str  = read_asc(p(day_filename("str", yn_p, dn_p)))
        _, prior_mwt  = read_asc(p(day_filename("mwt", yn_p, dn_p)))
        _, prior_lai  = read_asc(p(day_filename("lai", yn_p, dn_p)))

    print(f"\nActive cells: {int((dem != nodata).sum()):,}")
    print(f"Prior mode: {'restart' if cfg['prior'] else 'fresh start'}")

    ic_pack = prior_pack.copy()
    ic_ati  = prior_ati.copy()
    ic_hdi  = prior_hdi.copy()
    ic_mwt  = prior_mwt.copy()
    ic_str  = prior_str.copy()
    ic_lai  = prior_lai.copy()
    cal = calendar_days(cfg)
    quiet = os.environ.get("BCM_QUIET", "0") == "1"
    skip_pet = os.environ.get("BCM_SKIP_PET_READ", "0") == "1"
    write_state = os.environ.get("BCM_WRITE_STATE_MAPS", "1") != "0"
    out_hdr = (
        f"{'Year':>6}{'DOY':>5}{'Basin':>8}"
        f"{'pptmm':>10}{'petmm':>10}{'tmxC':>10}{'tmnC':>10}{'tavgC':>10}"
        f"{'snwmm':>10}{'mltmm':>10}{'sblmm':>10}{'pckmm':>10}{'excmm':>10}"
        f"{'cwdmm':>10}{'aetmm':>10}{'strmm':>10}{'rchmm':>10}{'runmm':>10}"
        f"{'rchacft':>17}{'runacft':>17}{'aream':>20}"
        f"{'mwtmm':>10}{'smdmm':>11}{'smrmm':>11}{'rvspct':>11}\n"
    )

    def run_period(out_dir, outfile_name, ppt_list=None, tmin_list=None, tmax_list=None):
        prior_pack = ic_pack.copy()
        prior_ati  = ic_ati.copy()
        prior_hdi  = ic_hdi.copy()
        prior_mwt  = ic_mwt.copy()
        prior_str  = ic_str.copy()
        prior_lai  = ic_lai.copy()
        os.makedirs(out_dir, exist_ok=True)
        outpath = os.path.join(out_dir, outfile_name)
        fout = open(outpath, "w")
        fout.write(out_hdr)
        for i, (yn, dn, leap) in enumerate(cal):
            if not quiet:
                print(f"  Day {dn:3d} / {yn}")
            if ppt_list is not None:
                ppt = ppt_list[i]
            else:
                _, ppt = read_asc(p(day_filename("ppt", yn, dn)))
            if tmax_list is not None:
                tmax = tmax_list[i].copy()
                tmin = tmin_list[i].copy()
            else:
                _, tmax = read_asc(p(day_filename("tmx", yn, dn)))
                _, tmin = read_asc(p(day_filename("tmn", yn, dn)))
            if not skip_pet:
                read_asc(p(day_filename("pet", yn, dn)))
            swap = (tmin != nodata) & (tmax != nodata) & (tmin >= tmax)
            tmax = np.where(swap, tmin + 5.0, tmax)
            if RESEED_DIR:
                pyn, pdn = yn, dn - 1
                if pdn == 0:
                    pyn -= 1
                    pdn = 366 if is_leap(pyn) else 365
                rf = os.path.join(RESEED_DIR, day_filename("str", pyn, pdn))
                if os.path.exists(rf):
                    _, prior_str = read_asc(rf)
            pw_grid  = kfactor_for_day(dn, pw_monthly)
            bet_grid = kfactor_for_day(dn, bet_monthly)
            alb_grid = kfactor_for_day(dn, alb_monthly)
            oz_grid  = kfactor_for_day(dn, oz_monthly)
            kfac_day = kfactor_for_day(dn, kfactor_monthly)
            kfac_day = np.where(kfac_day == nodata, 0.0, kfac_day)
            sumrad, sumtrn, petday, clouds = compute_solar_pet(
                dn, cfg["hstep"], cfg, terrain,
                tmax, tmin, bca, bcb, bcc_arr, pt_alpha, nodata,
                pw_grid, bet_grid, alb_grid, oz_grid
            )
            snow_out = snow_step(
                dn, ppt, tmax, tmin, petday, terrain["elev"],
                prior_pack, prior_ati, prior_hdi, prior_mwt,
                snowaccum_map, mfmax_map, mfmin_map, cfg, sumrad, nodata
            )
            wb = water_balance(
                ppt, snow_out["meltday"], snow_out["snowfall"],
                petday, prior_str, soild_eff, wpmm, fcmm, pormm, geolks,
                kfac_day, nodata, cfg["rchrun_flag"], cfg["rchrun_limit"],
                soil_nd=soil_nd, imperv=imperv,
                zero_ks=zero_ks, water_veg=water_veg,
                imperv_release=imperv_release,
            )
            if water_veg.any():
                for _k in ("snowfall", "meltday", "sublday"):
                    snow_out[_k] = np.where(water_veg, 0.0, snow_out[_k])

            # July f90 L1871: soil waterbal (not snow). ks=0 cells forced to 0 there.
            _ok = ((ppt != nodata) & (wb["storday"] != nodata)
                   & (prior_str != nodata) & (geolks != nodata) & (geolks != 0)
                   & (~water_veg))
            _wbal = (ppt - snow_out["snowfall"] + snow_out["meltday"]
                     - wb["rch3day"] - wb["run2day"]
                     - (wb["storday"] - prior_str) - wb["vegaet"])
            _w = _wbal[_ok]
            print(
                f"  waterbal {yn}_{dn:03d}  n={_w.size}  mean={_w.mean():+.4f}  "
                f"min={_w.min():+.4f}  max={_w.max():+.4f}  "
                f"maxabs={np.abs(_w).max():.4f}  "
                f"%|w|<0.01={100.0 * (np.abs(_w) < 0.01).mean():.1f}",
                flush=True,
            )

            def wr(prefix, arr, fmt="%.2f", _yn=yn, _dn=dn):
                fp = os.path.join(out_dir, day_filename(prefix, _yn, _dn))
                write_asc(fp, hdr, arr, fmt=fmt)

            if cfg["aet_flag"]:  wr("aet", wb["vegaet"])
            if cfg["cwd_flag"]:  wr("cwd", wb["cwd1day"])
            if cfg["pet_flag"]:  wr("pet", petday)
            if cfg["rad_flag"]:  wr("rad", sumrad)
            if cfg["rch_flag"]:  wr("rch", wb["rch3day"])
            if cfg["run_flag"]:  wr("run", wb["run2day"])
            if cfg["snw_flag"]:  wr("snw", snow_out["snowfall"])
            if cfg["smd_flag"]:  wr("smd", wb["smdday"])
            if cfg["smr_flag"]:  wr("smr", wb["smrday"])
            if cfg["rvs_flag"]:  wr("rvs", snow_out["rvsday"])
            if cfg["mlt_flag"]:  wr("mlt", snow_out["meltday"])
            if cfg["sbl_flag"]:  wr("sbl", snow_out["sublday"])
            if cfg["exc_flag"]:  wr("exc", wb["exc1day"])
            if write_state:
                wr("pck", snow_out["packday"])
                wr("ati", snow_out["ati"])
                wr("hdi", snow_out["hdi"])
                wr("str", wb["storday"])
                wr("mwt", snow_out["mwatday"])
                wr("lai", prior_lai)
            stats = basin_stats(
                basin_grid, terrain["elev"], nodata, cellsize,
                ppt=ppt, pet=petday, tmax=tmax, tmin=tmin,
                snw=snow_out["snowfall"], mlt=snow_out["meltday"],
                sbl=snow_out["sublday"],  pck=snow_out["packday"],
                exc=wb["exc1day"],        cwd=wb["cwd1day"],
                aet=wb["vegaet"],         str=wb["storday"],
                rch=wb["rch3day"],        run=wb["run2day"],
                mwt=snow_out["mwatday"],  smd=wb["smdday"],
                smr=wb["smrday"],         rvs=snow_out["rvsday"],
            )
            for bid, s in sorted(stats.items()):
                area = s["area_m2"]
                tavg = (s.get("tmax", 0) + s.get("tmin", 0)) * 0.5
                rchacft = (s.get("rch", 0) / 1000.0) * area / 1233.482
                runacft = (s.get("run", 0) / 1000.0) * area / 1233.482
                fout.write(
                    f"  {yn:4d}{dn:5d}{bid:8d}"
                    f"  {s.get('ppt',0):8.1f}  {s.get('pet',0):8.2f}"
                    f"  {s.get('tmax',0):8.1f}  {s.get('tmin',0):8.1f}  {tavg:8.1f}"
                    f"  {s.get('snw',0):8.1f}  {s.get('mlt',0):8.1f}"
                    f"  {s.get('sbl',0):8.1f}  {s.get('pck',0):8.1f}"
                    f"  {s.get('exc',0):8.1f}  {s.get('cwd',0):8.2f}"
                    f"  {s.get('aet',0):8.2f}  {s.get('str',0):8.1f}"
                    f"  {s.get('rch',0):8.1f}  {s.get('run',0):8.1f}"
                    f"  {rchacft:15.1f}  {runacft:15.1f}  {area:18.0f}"
                    f"  {s.get('mwt',0):8.1f}"
                    f"  {s.get('smd',0):9.2f}  {s.get('smr',0):9.2f}"
                    f"  {s.get('rvs',0):9.2f}\n"
                )
            prior_pack = snow_out["packday"]
            prior_ati  = snow_out["ati"]
            prior_hdi  = snow_out["hdi"]
            prior_mwt  = snow_out["mwatday"]
            prior_str  = wb["storday"]
        fout.close()
        return outpath

    npy_path = os.environ.get("BCM_PPT_NPY")
    storm_spec = os.environ.get("BCM_STORMS")
    if npy_path and storm_spec:
        storms = parse_storms(storm_spec)
        if not quiet:
            print(f"  Batch: {len(storms)} storms from {npy_path}")
        tmax_days, tmin_days = [], []
        for yn, dn, _leap in cal:
            _, tmax = read_asc(p(day_filename("tmx", yn, dn)))
            _, tmin = read_asc(p(day_filename("tmn", yn, dn)))
            swap = (tmin != nodata) & (tmax != nodata) & (tmin >= tmax)
            tmax_days.append(np.where(swap, tmin + 5.0, tmax))
            tmin_days.append(tmin)
        cube = np.load(npy_path, mmap_mode="r")
        run_root = os.environ.get("BCM_RUN_ROOT", OUT_DIR)
        last = None
        for s in storms:
            if not quiet:
                print(f"  Storm {s}")
            nd = len(cal)
            ppt_list = [np.asarray(cube[s - 1, i], dtype=np.float64) for i in range(nd)]
            dest = os.path.join(run_root, f"storm{s:03d}")
            last = run_period(dest, f"storm{s:03d}.out", ppt_list, tmin_days, tmax_days)
        elapsed = time.time() - t0
        print(f"\nDone. Elapsed: {elapsed:.1f} s")
        print(f"Summary output: {last}")
        return

    outpath = run_period(OUT_DIR, cfg["outfile"])
    elapsed = time.time() - t0
    print(f"\nDone. Elapsed: {elapsed:.1f} s")
    print(f"Summary output: {outpath}")


if __name__ == "__main__":
    main()
