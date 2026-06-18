# ExtData.rc Updates for OSSE CO2 Emissions

## Summary
Updated `ExtData.rc` in `osse_files/` to match the new CO2 emission configuration from `sim.template` and align with `HEMCO_Config_osse.rc`.

## Changes Made

### 1. Fossil Fuel Emissions (Line ~399)
**Changed from:**
```
FOSSILCO2_ODIAC  kg/m2/s N Y F%y4-%m2-%d2T%h2:00:00 none none CO2_Flux /nobackupp27/hnesser/CO2_inversion/priors/ODIAC/HEMCO/ODIAC_%y4%m2.nc
```

**Changed to:**
```
# --- CO2 Fossil Fuel (FF_regrid2 from sim.template) ---
FOSSILCO2_ODIAC  kg/m2/s N Y F%y4-%m2-%d2T%h2:00:00 none none CO2 /nobackupp17/jliu7/INVENTORY/V4/Fossilfuel/FF_regrid2/FF_%y4%m2%d2%h2.nc
```

**Changes:**
- Path: `/nobackupp27/hnesser/...` → `/nobackupp17/jliu7/INVENTORY/V4/...`
- Model: ODIAC → FF_regrid2
- Variable name: `CO2_Flux` → `CO2`
- Temporal: Hourly (unchanged)

---

### 2. Ocean Exchange - Disabled old (Line ~402)
**Changed from:**
```
# --- CO2 Ocean Exchange ---
OCEANCO2_SCALED_MONTHLY kg/m2/s N  Y F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc
```

**Changed to:**
```
# --- CO2 Ocean Exchange (OLD - not used in OSSE) ---
#OCEANCO2_SCALED_MONTHLY kg/m2/s N  Y F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc
```

**Action:** Commented out (not used in OSSE)

---

### 3. Ocean Exchange - ECCO-Darwin (Line ~405)
**Changed from:**
```
# --- CO2 Ocean Exchange ECCO DARWIN ---
OCEANS  kg/m2/s N Y F%y4-%m2-%d2T%00:00:00 none none CO2_Flux /nobackupp27/hnesser/CO2_inversion/priors/ECCO_Darwin/HEMCO/ECCO_Darwin_%y4.nc
```

**Changed to:**
```
# --- CO2 Ocean Exchange ECCO DARWIN (ECCO-Darwin-MON-v05 from sim.template) ---
OCEAN_ECCO_DARWIN  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2 /nobackupp17/jliu7/INVENTORY/V4/Ocean/ECCO-Darwin-MON-v05/ECCO_Darwin_%y4%m2.nc
```

**Changes:**
- Variable name: `OCEANS` → `OCEAN_ECCO_DARWIN`
- Path: `/nobackupp27/hnesser/...` → `/nobackupp17/jliu7/INVENTORY/V4/...`
- Model version: Added `-MON-v05` specification
- Temporal: Daily → Monthly (`-%d2T` removed, `%h2` removed)
- Field name: `CO2_Flux` → `CO2`

---

### 4. Net Biosphere Exchange (Line ~414-415)
**Changed from:**
```
# --- CARDAMOM prior fluxes ---
#NBE kg/m2/s N Y F%y4-%m2-%d2T%03:00:00 none none NBE /nobackupp27/hnesser/CO2_inversion/priors/CARDAMOM/HEMCO/CARDAMOM_%y4.nc
NBE kg/m2/s N Y F0;013000 none none NBE /nobackupp27/hnesser/CO2_inversion/priors/CARDAMOM/HEMCO/CARDAMOM_%y4.nc
```

**Changed to:**
```
# --- CARDAMOM prior fluxes (CARDAMOM-ECCO from sim.template) ---
NBE_CARDAMOM  kg/m2/s N Y F%y4-%m2-%d2T%h2:00:00 none none CO2 /nobackupp17/jliu7/INVENTORY/V4/Balbio/CARDAMOM-ECCO/CARDAMOM_%y4%m2%d2%h2.nc
```

**Changes:**
- Variable name: `NBE` → `NBE_CARDAMOM`
- Path: `/nobackupp27/hnesser/...` → `/nobackupp17/jliu7/INVENTORY/V4/...`
- Model: CARDAMOM → CARDAMOM-ECCO
- Temporal: Fixed offset `F0;013000` → Standard hourly `F%y4-%m2-%d2T%h2:00:00`
- Field name: `NBE` → `CO2`
- File naming: `CARDAMOM_%y4.nc` → `CARDAMOM_%y4%m2%d2%h2.nc` (hourly files)

---

### 5. GPP - NEW (After line ~415)
```
# --- GPP (Gross Primary Production from ORCHIDEE) ---
GPP_ORCHIDEE  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2 /nobackupp17/jliu7/INVENTORY/V4/GPP/ORCHIDEE/ORCHIDEE_GPP_%y4%m2.nc
```

**Details:**
- Variable name: `GPP_ORCHIDEE`
- Model: ORCHIDEE
- Temporal: Monthly
- Field name: `CO2`
- Units: `kg/m2/s`

---

### 6. TER - NEW (After GPP)
```
# --- TER (Terrestrial Ecosystem Respiration from CLASS-CTEM) ---
TER_CLASS_CTEM  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2 /nobackupp17/jliu7/INVENTORY/V4/TER2/CLASS-CTEM/CLASS_CTEM_TER_%y4%m2.nc
```

**Details:**
- Variable name: `TER_CLASS_CTEM`
- Model: CLASS-CTEM
- Temporal: Monthly
- Field name: `CO2`
- Units: `kg/m2/s`
- Note: This is the adjoint-optimized component (adjId=1 in sim.template)

---

## Final CO2 ExtData.rc Configuration

| Variable Name | Model | Temporal Resolution | Field Name | Path |
|--------------|-------|---------------------|------------|------|
| `FOSSILCO2_ODIAC` | FF_regrid2 | Hourly | `CO2` | `/nobackupp17/jliu7/INVENTORY/V4/Fossilfuel/FF_regrid2/` |
| `OCEAN_ECCO_DARWIN` | ECCO-Darwin-MON-v05 | Monthly | `CO2` | `/nobackupp17/jliu7/INVENTORY/V4/Ocean/ECCO-Darwin-MON-v05/` |
| `NBE_CARDAMOM` | CARDAMOM-ECCO | Hourly | `CO2` | `/nobackupp17/jliu7/INVENTORY/V4/Balbio/CARDAMOM-ECCO/` |
| `GPP_ORCHIDEE` | ORCHIDEE | Monthly | `CO2` | `/nobackupp17/jliu7/INVENTORY/V4/GPP/ORCHIDEE/` |
| `TER_CLASS_CTEM` | CLASS-CTEM | Monthly | `CO2` | `/nobackupp17/jliu7/INVENTORY/V4/TER2/CLASS-CTEM/` |

---

## Important Notes

1. **All paths updated** from `/nobackupp27/hnesser/` to `/nobackupp17/jliu7/INVENTORY/V4/`

2. **Variable names must match** between `ExtData.rc` and `HEMCO_Config_osse.rc`

3. **Field names** in NetCDF files:
   - All now use `CO2` (not `CO2_Flux` or `NBE`)
   - Make sure the actual NetCDF files use this variable name

4. **Temporal resolution**:
   - Hourly: `F%y4-%m2-%d2T%h2:00:00`
   - Monthly: `F%y4-%m2-01T00:00:00`

5. **File naming conventions**:
   - Hourly files: `FILENAME_%y4%m2%d2%h2.nc`
   - Monthly files: `FILENAME_%y4%m2.nc`

---

## Verification Checklist

- [ ] All paths point to `/nobackupp17/jliu7/INVENTORY/V4/`
- [ ] Variable names in ExtData.rc match those in HEMCO_Config_osse.rc
- [ ] Field names (`CO2`) match actual variable names in NetCDF files
- [ ] Temporal resolution matches file availability
- [ ] File naming patterns match actual files on disk
- [ ] Old unused entries are commented out
