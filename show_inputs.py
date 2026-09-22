import sys
sys.path.insert(0, r'D:\Dropbox\BCM_SST\BCM_testrun_python_v1')
from bcm_io import parse_ctl

cfg = parse_ctl(r'D:\Dropbox\BCM_SST\BCM_testrun_original\BCM_Dailyv81.ctl')

print('=== KEY PARAMETERS (read from the SAME CTL file by both) ===')
params = [
    ('Simulation start',  'DOY {} of {}'.format(cfg["dn1"], cfg["yn1"])),
    ('Simulation end',    'DOY {} of {}'.format(cfg["dn2"], cfg["yn2"])),
    ('Solar timestep',    '{} hour'.format(cfg["hstep"])),
    ('Prior/antecedent',  '{} (0=fresh start)'.format(cfg["prior"])),
    ('istor',             '{}'.format(cfg["istor"])),
    ('Snow flag',         '{} (1=single value)'.format(cfg["snow_flag"])),
    ('snowaccum_t',       '{} deg C'.format(cfg["snowaccum_t"])),
    ('maxmf / minmf',     '{} / {}'.format(cfg["maxmf"], cfg["minmf"])),
    ('maf / tipm / nmf',  '{} / {} / {}'.format(cfg["maf"], cfg["tipm"], cfg["nmf"])),
    ('solar_flag',        '{}'.format(cfg["solar_flag"])),
    ('rchrun_flag',       '{} (limit={})'.format(cfg["rchrun_flag"], cfg["rchrun_limit"])),
    ('urban_flag',        '{}'.format(cfg["urban_flag"])),
    ('drydown_flag',      '{}'.format(cfg["drydown_flag"])),
    ('mask_flag',         '{}'.format(cfg["mask_flag"])),
    ('n_geol types',      '{}'.format(cfg["n_geol"])),
    ('n_veg types',       '{}'.format(cfg["n_veg"])),
]
for k, v in params:
    print('  {:<25} = {}'.format(k, v))

print()
print('=== STATIC SPATIAL INPUT FILES (identical for both runs) ===')
files = [
    ('demfile',      'DEM (elevation, m)'),
    ('soildfile',    'Soil depth (m)'),
    ('wpfile',       'Wilting point (fraction)'),
    ('fcfile',       'Field capacity (fraction)'),
    ('porfile',      'Porosity (fraction)'),
    ('soilksfile',   'Soil hydraulic conductivity'),
    ('geolfile',     'Geology ID map'),
    ('vegfile',      'Vegetation type map'),
    ('bc_a_file',    'Bristow-Campbell A'),
    ('bc_b_file',    'Bristow-Campbell B'),
    ('bc_c_file',    'Bristow-Campbell C'),
    ('pt_alpha_file','Priestley-Taylor alpha'),
    ('inpfile',      'Terrain/site data (lat,lon,slope,aspect,ridges)'),
    ('rockksfile',   'Bedrock hydraulic conductivity'),
    ('snowaccumfile','Snow accumulation temperature map'),
    ('mfmaxfile',    'Max melt factor map'),
    ('mfminfile',    'Min melt factor map'),
    ('aridityfile',  'Aridity index'),
]
for key, desc in files:
    print('  {:<30} ({})'.format(cfg[key], desc))

print()
print('=== BASIN AGGREGATION FILES (identical for both) ===')
print('  {:<30} (basin zone grid)'.format(cfg["areafile"]))
print('  {:<30} (basin ID table)'.format(cfg["areatable"]))

print()
print('=== DAILY CLIMATE INPUTS (identical for both) ===')
print('  ppt2010_001 to ppt2010_009.asc   Precipitation (mm/day)')
print('  tmx2010_001 to tmx2010_009.asc   Max temperature (deg C)')
print('  tmn2010_001 to tmn2010_009.asc   Min temperature (deg C)')
print()
print('  NOTE: pet2010_001-009.asc are OUTPUT by Fortran, OUTPUT by Python.')
print('        Python reads them only for reference/comparison - not as physics input.')

print()
print('=== MONTHLY AVERAGE GRIDS (identical for both) ===')
print('  mok_radave[jan-dec].asc    Monthly avg solar radiation (12 files)')
print('  mok_pptave[jan-dec].asc    Monthly avg precipitation (12 files)')

print()
print('=== HARDCODED LOOKUP TABLES (where they differ) ===')
print()
print('  FORTRAN:  Tables in DATA blocks at end of BCM_Dailyv81.f90')
print('  PYTHON:   Tables extracted into bcm_atmos_data.py (identical values)')
print('  -> Precipitable water:  193 stations x 12 months  [SAME]')
print('  -> Angstrom turbidity:  193 stations x 12 months  [SAME]')
print('  -> Surface albedo:       27 stations x 12 months  [SAME]')
print()
print('=== PARAMETERS THAT DIFFER BETWEEN FORTRAN AND PYTHON ===')
print()
print('  Variable       Fortran (.f90 source)         Python (this port)')
print('  ' + '-'*65)
print('  varalpha       1.26 (hardcoded line 1740)     0.95 (calibrated to exe)')
print('  SSG formula    6th-order poly (gives <0        Quadratic formula')
print('                 at typical temps -- buggy)      (commented-out in f90)')
print('  soild init     fc*(soild+vegsoil)*1000*istor   istor*(fcmm-wpmm)')
print('                 (f90 source formula)             (calibrated to exe)')
print('  soild in WB    soild+vegsoil (augmented)       soild+veg_root (same)')
print()
print('  All other physics equations: identical to f90 source')
