# OSSE 4D-Var Setup Summary

## Overview
Created OSSE (Observing System Simulation Experiment) versions of the 4D-Var optimization system using synthetic OCO-2 observations from ORCHIDEE-ECCO2.

## Files Created

### 1. `4dvar_optimizer.osse.run`
OSSE version of the main 4D-Var optimizer script.

**Key Differences from `4dvar_optimizer.run`:**
- **Job name**: `gchp_4dvar_osse` (vs `gchp_4dvar`)
- **WORK_DIR**: `../4dvar_design_osse` (vs `../4dvar_design`)
- **Simulation period**: 2016-01-01 to 2016-01-31 (vs 2019-01-01 to 2019-01-31)
- **Output files** (all with `_osse` suffix to avoid conflicts):
  - `sigma_co2_osse.nc4` (vs `sigma_co2.nc4`)
  - `4dvar_state_osse.npz` (vs `4dvar_state.npz`)
  - `forcing_files_osse/` (vs `forcing_files/`)
  - `4dvar_output_osse/` (vs `4dvar_output/`)
  - `4dvar_optimizer_osse.log` (vs `4dvar_optimizer.log`)
  - `plots_osse/` (vs `plots/`)
- **Python script**: Uses `co2_adjoint_forcing_osse.py` instead of `co2_adjoint_forcing.py`

### 2. `co2_adjoint_forcing_osse.py`
OSSE version of the adjoint forcing computation script.

**Key Differences from `co2_adjoint_forcing.py`:**
- **Import**: Uses `convert_satelite_tracks_osse.read_oco_daily_osse` instead of `co2_sat_compare_monthly.read_oco_monthly`
- **Data path**: 
  ```python
  BASE_FOLD_OBS = '/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/'
  ```
- **New function**: `read_oco_monthly_osse(year, month, base_fold)` - aggregates daily OSSE observations into monthly dataset
- **Observation modes**: Only LG and LN (no OG for OSSE data)
- **Matching**: Uses improved time + lat/lon matching (same as updated `co2_adjoint_forcing.py`)

## Data Sources

### OSSE Observations
- **Location**: `/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/`
- **Structure**:
  - `OCO2b10_29Aug2020_LG/YYYY/MM/DD.nc` (Land Glint)
  - `OCO2b10_29Aug2020_LN/YYYY/MM/DD.nc` (Land Nadir)
  - No Ocean Glint (OG) data available
- **Time period**: 2016 data available

## File Separation Strategy

All OSSE files use distinct names/paths to prevent conflicts:

| Real Data | OSSE Data |
|-----------|-----------|
| `4dvar_optimizer.run` | `4dvar_optimizer.osse.run` |
| `co2_adjoint_forcing.py` | `co2_adjoint_forcing_osse.py` |
| `../4dvar_design/` | `../4dvar_design_osse/` |
| `sigma_co2.nc4` | `sigma_co2_osse.nc4` |
| `4dvar_state.npz` | `4dvar_state_osse.npz` |
| `forcing_files/` | `forcing_files_osse/` |
| `4dvar_output/` | `4dvar_output_osse/` |
| `4dvar_optimizer.log` | `4dvar_optimizer_osse.log` |
| `plots/` | `plots_osse/` |

## Usage

### Submit OSSE 4D-Var Job
```bash
cd gchp_4Dvar
qsub 4dvar_optimizer.osse.run
```

### Submit with custom parameters
```bash
qsub -v MAX_ITER=10,SIGMA_B=0.3 4dvar_optimizer.osse.run
```

### Resume from previous run
```bash
qsub -v RESTART=true 4dvar_optimizer.osse.run
```

## Prerequisites

Before running, ensure:
1. `../4dvar_design_osse/` directory exists with:
   - `forward_run/` containing GCHP forward model
   - `adjoint/` containing GCHP adjoint model
   - Both with `gchp` executables and `gchp.env`

2. Control files in `control_files/`:
   - `forward/` with forward run configuration
   - `adjoint/` with adjoint run configuration

3. Python environment with required packages:
   - numpy, pandas, xarray, netCDF4, metpy

4. HEMCO and ExtData.rc configured to read `../sigma_co2_osse.nc4`

## Observation Matching

Both `co2_adjoint_forcing.py` and `co2_adjoint_forcing_osse.py` now use improved matching:

1. **Time matching**: Within 60 seconds
2. **Spatial matching**: 
   - Latitude: Within 0.1°
   - Longitude: Within 0.1° (with wrapping)
3. **Algorithm**:
   - First filter by time tolerance
   - Among time-matched obs, find closest in lat/lon
   - Reject if spatial distance exceeds tolerances

This ensures each GCHP profile is matched to the correct observation location, critical when combining LG and LN modes.

## Notes

- The OSSE run is completely independent of the real data run
- Both can run simultaneously without conflicts
- Output files are clearly separated by `_osse` suffix
- Same optimization parameters (NLAT, NLON, SIGMA_B, etc.) can be used for both
