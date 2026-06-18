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
        else:
            print(f'File not found: {file_path}')
    if daily_datasets:
        # Combine datasets along a dimension 'nSamples'
        ds = xr.concat(daily_datasets, dim='nSamples')
        time_hours = ds['time'].values #hours
        # Convert hours to datetime
        time_dt = np.array([base_date + timedelta(hours=float(h)) for h in time_hours])
        #assign it back to ds['time']
        ds['time'] = ('nSamples', time_dt)
        ds['longitude'] = (ds['longitude'] % 360)
        ds = ds[['latitude', 'longitude', 'time']]
        ds = ds.sortby('time')
        ds = ds.set_coords("time")
        ds = ds.swap_dims({"nSamples": "time"})
        return ds
    else:
        print(f'No data found for date: {date_str}')
        return None


def convert_track_yearly(year, base_fold):
    """Convert daily OSSE satellite tracks to a single yearly file."""
    datasets = []
    for month in range(1, 13):
        num_days = calendar.monthrange(year, month)[1]
        for day in range(1, num_days + 1):
            try:
                date = f'{year}{month:02d}{day:02d}'
                print(date)
                ds = read_oco_daily_osse(date, base_fold)
                if ds is not None:
                    datasets.append(ds)
            except Exception as e:
                print(f"Failed to process {year}-{month:02d}-{day:02d}: {e}")

    if datasets:
        print(f"Combining {len(datasets)} daily datasets for year {year}...")
        combined = xr.concat(datasets, dim='time')

        # Save to new netCDF file
        # Define encoding for time
        time_encoding = {
            'time': {
                'units': 'hours since 1900-01-01 00:00:00',
                'calendar': 'proleptic_gregorian',
                'dtype': 'float64',
                '_FillValue': np.nan,
            }
        }

        output_file = f'sat_track/track_file_osse_{year}.nc'
        combined.to_netcdf(output_file, encoding=time_encoding)
        print(f"Saved: {output_file} ({len(combined.time)} observations)")
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
