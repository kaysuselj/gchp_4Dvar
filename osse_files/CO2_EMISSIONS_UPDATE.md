# CO2 Emissions Update for HEMCO_Config_osse.rc

## Summary
Updated `HEMCO_Config_osse.rc` to match the CO2 emission configuration from `sim.template`.

## Changes Made

### 1. Extension Switches (Lines ~97-107)
Added two new switches and updated comments:

```
    --> FOSSIL_ODIAC           :       true     # 2000-2018 (Fossilfuel in sim.template)
    --> OCEAN_ECCO_DARWIN      :       true     # 2014-2023 (Ocean biosphere in sim.template)
    --> NBE_CARDAMOM           :       true     # 2014-2023 (Biosphere CARDAMOM-ECCO in sim.template)
    --> GPP_ORCHIDEE           :       true     # 2014-2023 (GPP ORCHIDEE in sim.template) [NEW]
    --> TER_CLASS_CTEM         :       true     # 2014-2023 (TER CLASS-CTEM in sim.template) [NEW]
```

### 2. Fossil Fuel Emissions (Line ~874)
**Changed from:**
```
0 FOSSILCO2_ODIAC  /nobackupp27/hnesser/CO2_inversion/priors/ODIAC/HEMCO/ODIAC_$YYYY$MM.nc  CO2_Flux 2014-2023/1-12/1-31/0-23 C xy kg/m2/s CO2 750 1 500
```

**Changed to:**
```
# Fossil fuel emissions from FF_regrid2 (sim.template line 79: hourly)
0 FOSSILCO2_ODIAC  /nobackupp11/jliu7/INVENTORY/V4/Fossilfuel/FF_regrid2/FF_$YYYY$MM$DD$HH.nc  CO2 2014-2023/1-12/1-31/0-23 C xy kg/m2/s CO2 750 1 500
```

### 3. Ocean Exchange (Line ~896)
**Changed from:**
```
0 OCEANS  /nobackupp27/hnesser/CO2_inversion/priors/ECCO_Darwin/HEMCO/ECCO_Darwin_$YYYY.nc  CO2_Flux 2014-2023/1-12/1-31/0 C xy kg/m2/s CO2 750 2 500
```

**Changed to:**
```
# Ocean biosphere from ECCO-Darwin-MON-v05 (sim.template line 84: monthly)
0 OCEAN_ECCO_DARWIN  /nobackupp11/jliu7/INVENTORY/V4/Ocean/ECCO-Darwin-MON-v05/ECCO_Darwin_$YYYY$MM.nc  CO2 2014-2023/1-12/1/0 C xy kg/m2/s CO2 750 2 500
```

### 4. Net Biosphere Exchange (Line ~932)
**Changed from:**
```
0 NBE  /nobackupp27/hnesser/CO2_inversion/priors/CARDAMOM/HEMCO/CARDAMOM_$YYYY.nc  NBE 2014-2023/1-12/1-31/0-23/+90minute C xy kg/m2/s CO2 750 3 500
```

**Changed to:**
```
# Net Biosphere Exchange from CARDAMOM-ECCO (Biosphere in sim.template line 83)
0 NBE_CARDAMOM  /nobackupp11/jliu7/INVENTORY/V4/Balbio/CARDAMOM-ECCO/CARDAMOM_$YYYY$MM$DD$HH.nc  CO2 2014-2023/1-12/1-31/0-23 C xy kg/m2/s CO2 750 3 500
```

### 5. GPP (Gross Primary Production) - NEW (After line ~933)
```
#==============================================================================
# --- CO2: GPP (Gross Primary Production) from ORCHIDEE ---
#
# From sim.template line 84: GPP ORCHIDEE monthly
#==============================================================================
(((GPP_ORCHIDEE
0 GPP_ORCHIDEE  /nobackupp11/jliu7/INVENTORY/V4/GPP/ORCHIDEE/ORCHIDEE_GPP_$YYYY$MM.nc  CO2 2014-2023/1-12/1/0 C xy kg/m2/s CO2 750 4 500
)))GPP_ORCHIDEE
```

### 6. TER (Terrestrial Ecosystem Respiration) - NEW (After GPP section)
```
#==============================================================================
# --- CO2: TER (Terrestrial Ecosystem Respiration) from CLASS-CTEM ---
#
# From sim.template line 85: TER CLASS-CTEM monthly
# This is the adjoint-optimized component (adjId=1)
#==============================================================================
(((TER_CLASS_CTEM
0 TER_CLASS_CTEM  /nobackupp11/jliu7/INVENTORY/V4/TER2/CLASS-CTEM/CLASS_CTEM_TER_$YYYY$MM.nc  CO2 2014-2023/1-12/1/0 C xy kg/m2/s CO2 750 5 500
)))TER_CLASS_CTEM
```

## Final CO2 Emission Configuration

### Active CO2 Emissions (5 sources):

| # | Source | Model | Temporal Resolution | Hier | Cat | File Path |
|---|--------|-------|---------------------|------|-----|-----------|
| 1 | Fossil Fuel | FF_regrid2 | Hourly | 500 | 1 | `/nobackupp11/jliu7/INVENTORY/V4/Fossilfuel/FF_regrid2/FF_$YYYY$MM$DD$HH.nc` |
| 2 | Ocean | ECCO-Darwin-MON-v05 | Monthly | 500 | 2 | `/nobackupp11/jliu7/INVENTORY/V4/Ocean/ECCO-Darwin-MON-v05/ECCO_Darwin_$YYYY$MM.nc` |
| 3 | Net Biosphere | CARDAMOM-ECCO | Hourly | 500 | 3 | `/nobackupp11/jliu7/INVENTORY/V4/Balbio/CARDAMOM-ECCO/CARDAMOM_$YYYY$MM$DD$HH.nc` |
| 4 | GPP | ORCHIDEE | Monthly | 500 | 4 | `/nobackupp11/jliu7/INVENTORY/V4/GPP/ORCHIDEE/ORCHIDEE_GPP_$YYYY$MM.nc` |
| 5 | TER | CLASS-CTEM | Monthly | 500 | 5 | `/nobackupp11/jliu7/INVENTORY/V4/TER2/CLASS-CTEM/CLASS_CTEM_TER_$YYYY$MM.nc` |

All sources:
- Use **Hier=500** (high hierarchy for optimization)
- Available for **2014-2023**
- Units: **kg/m2/s**

## Mapping to sim.template

| sim.template Entry | HEMCO Entry | Status |
|-------------------|-------------|---------|
| 1. Fossilfuel (FF_regrid2) | FOSSILCO2_ODIAC | ✅ Updated |
| 2. Biofuel (CASA-GFED3-FUEL) | - | ❌ Not enabled (fwd=F) |
| 3. Biomass Burning (GFED4) | - | ❌ Not enabled (fwd=F) |
| 4. Biosphere (CARDAMOM-ECCO) | NBE_CARDAMOM | ✅ Updated |
| 5. Ocean biosphere (ECCO-Darwin-MON-v05) | OCEAN_ECCO_DARWIN | ✅ Updated |
| 6. GPP (ORCHIDEE) | GPP_ORCHIDEE | ✅ Added |
| 7. TER (CLASS-CTEM) | TER_CLASS_CTEM | ✅ Added |

**Note:** Biofuel and Biomass Burning are disabled in sim.template (fwd=F), so they were not added to HEMCO_Config_osse.rc.
