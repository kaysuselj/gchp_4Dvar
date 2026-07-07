import xarray as xr
import numpy as np
import calendar
from datetime import datetime, timedelta
import os

# Read processed OCO-2 OSSE daily data from Junjie
def read_oco_daily_osse(date, base_fold):
    #date='20150101'
    #base_fold='/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/'
    base_date = datetime.strptime(date, '%Y%m%d')
    date_str = base_date.strftime('%Y/%m/%d')
    daily_datasets = []
    # Only LG and LN exist for OSSE (no OG)
    sub_folds = {
        'LG': 'OCO2b10_29Aug2020_LG',
        'LN': 'OCO2b10_29Aug2020_LN'
    }
    for fold_name in sub_folds.values():
        file_path = os.path.join(base_fold, fold_name, date_str + '.nc')
        if os.path.exists(file_path):
            ds = xr.open_dataset(file_path)
            daily_datasets.append(ds)
    if daily_datasets:
        # Combine datasets along a dimension 'nSamples'
        ds = xr.concat(daily_datasets, dim='nSamples')
        time_hours = ds['time'].values #hours
        # Convert hours to datetime
        time_dt = np.array([base_date + timedelta(hours=float(h)) for h in time_hours])
        #assign it back to ds['time']
        ds['time'] = ('nSamples', time_dt)
        ds['longitude'] = (ds['longitude'] % 360)

        # Extract variables needed for both track file and adjoint forcing
        ds_clean = xr.Dataset({
            'latitude':             (['nSamples'],             ds['latitude'].values),
            'longitude':            (['nSamples'],             ds['longitude'].values),
            'time':                 (['nSamples'],             time_dt),
            'pressure':             (['nSamples', 'maxLevels'], ds['pressure'].values),
            'xCO2':                 (['nSamples'],             ds['xCO2'].values),
            'xCO2-apriori':         (['nSamples'],             ds['xCO2-apriori'].values),
            'xCO2-uncertainty':     (['nSamples'],             ds['xCO2-uncertainty'].values),
            'xCO2-averagingKernel': (['nSamples', 'maxLevels'], ds['xCO2-averagingKernel'].values),
            'CO2-apriori':          (['nSamples', 'maxLevels'], ds['CO2-apriori'].values),
        })

        ds_clean = ds_clean.sortby('time')
        ds_clean = ds_clean.set_coords("time")
        ds_clean = ds_clean.swap_dims({"nSamples": "time"})
        return ds_clean
    else:
        return None


def convert_track_yearly(year, base_fold):
    """Convert daily OSSE satellite tracks to a single yearly file."""
    datasets = []
    for month in range(1, 13):
        num_days = calendar.monthrange(year, month)[1]
        n_ok = 0
        for day in range(1, num_days + 1):
            try:
                date = f'{year}{month:02d}{day:02d}'
                ds = read_oco_daily_osse(date, base_fold)
                if ds is not None:
                    datasets.append(ds)
                    n_ok += 1
            except Exception as e:
                print(f"  Failed {year}-{month:02d}-{day:02d}: {e}")
        print(f"  {year}-{month:02d}: {n_ok}/{num_days} days with observations")

    if datasets:
        print(f"Combining {len(datasets)} daily datasets for year {year}...")
        combined = xr.concat(datasets, dim='time')

        # Create dataset matching EXACT format of working file
        # Variables: latitude (float32), longitude (float32)
        # Coordinate: time (float64)
        ds_final = xr.Dataset({
            'latitude': (['time'], combined['latitude'].values.astype('float32')),
            'longitude': (['time'], combined['longitude'].values.astype('float32')),
        }, coords={
            'time': combined['time'].values
        })

        # Set attributes to match the working file format exactly
        ds_final.attrs = {}  # No global attributes

        # Clear all variable attributes (will set via encoding instead)
        ds_final['latitude'].attrs = {}
        ds_final['longitude'].attrs = {}
        ds_final['time'].attrs = {}

        # Define encoding to match working file
        # _FillValue and other attributes are set via encoding, not attrs
        encoding = {
            'latitude': {'dtype': 'float32', '_FillValue': np.float32(np.nan)},
            'longitude': {'dtype': 'float32', '_FillValue': np.float32(np.nan)},
            'time': {
                'dtype': 'float64',
                '_FillValue': np.float64(np.nan),
                'units': 'hours since 1900-01-01',
                'calendar': 'proleptic_gregorian'
            }
        }

        os.makedirs('sat_track', exist_ok=True)
        output_file = f'sat_track/track_file_osse_{year}.nc'
        ds_final.to_netcdf(output_file, encoding=encoding)
        print(f"Saved: {output_file} ({len(ds_final.time)} observations)")
    else:
        print(f"No data found for year {year}")


start_year = 2016
end_year = 2016
base_fold = '/nobackupp17/jliu7/OCO2/L2#/OCO2-B10/PSEUDO/TRENDY/ORCHIDEE-ECCO2/'

for year in range(start_year, end_year + 1):
    print(f"\n{'='*60}")
    print(f"Processing year {year}")
    print(f"{'='*60}")
    convert_track_yearly(year, base_fold)
