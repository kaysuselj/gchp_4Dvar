#Program for comparing GCHP with satellite data (OCO-2/OCO-3), satellite tracks outputted from GCHP
import xarray as xr
import numpy as np
import glob
import os
from datetime import datetime, timedelta
import pandas as pd
from scipy.spatial import cKDTree
#from interpolation import VerticalGrid  #interpolation.py in GOOPy (unused)
from metpy.interpolate import interpolate_1d
from metpy.units import units
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import linregress
import pickle
import calendar
from dateutil.relativedelta import relativedelta

#Read processed OCO-2/OCO-3 daily data from Junjie
def read_oco_daily(date,main_fold):
    #date='20150101'
    #main_fold='/nobackupp17/jliu7/OCO2/BAKER/B11.2_L2#/'
    base_date = datetime.strptime(date, '%Y%m%d')
    date_str = base_date.strftime('%Y/%m/%d')
    daily_datasets = []
    sub_folds=['LG','LN','OG']
    for fold in sub_folds:
        file_path = os.path.join(main_fold,fold, date_str + '.nc')
        if os.path.exists(file_path):
            ds = xr.open_dataset(file_path) #xr.open_mfdataset crashes.
            # add source label variable
            ds = ds.assign(
                source=("nSamples", [fold] * ds.sizes["nSamples"])
            )
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
        ds = ds.sortby('time')
        ds = ds.set_coords("time")
        ds = ds.swap_dims({"nSamples": "time"})
        return ds
    else:
        print(f'File not found: {date_str}')

def read_oco_monthly(year,month,fold):
    #year=2014
    #month=11
    num_days = calendar.monthrange(year, month)[1]
    datasets = []
    for day in range(1, num_days + 1):
        try:
            date=f'{year}{month:02d}{day:02d}'
            print(date)
            ds = read_oco_daily(date,fold)
            if ds is not None:
                datasets.append(ds)
        except Exception as e:
            print(f"Failed to process {year}-{month:02d}-{day:02d}: {e}")
    if datasets:
        #print(datasets)
        combined = xr.concat(datasets, dim='time')
    # Save to new netCDF file
    # Define encoding for time
    time_encoding = {
        'time': {
            'units': 'hours since 1900-01-01 00:00:00',
            'calendar': 'proleptic_gregorian',
            'dtype': 'float64',  # or 'double'
            '_FillValue': np.nan,
        }
    }
    combined['time'] = pd.to_datetime(combined['time'].values).round('s')
    return combined

def read_gchp_monthly(year,month,expname):
    #GCHP
    date = datetime(year, month, 1)
    prev_date = date - relativedelta(months=1)
    # Extract new year and month
    prev_year = prev_date.year
    prev_month = prev_date.month
    print(prev_year)
    print(prev_month)
    if expname in ['CATRINE_C90_CO2','CATRINE_C180_CO2']:
        fname_gchp=glob.glob('/nobackupp17/lwu5/GCHP/GCHP_run/CATRINE_CO2/'+expname+'/OutputDir/GEOSChem.sat_track.'+f"{prev_year}{prev_month:02d}"+'*z.nc4')
    elif expname=='CATRINE_C24_CO2':
        fname_gchp=glob.glob('/nobackupp17/lwu5/GCHP/GCHP_run/CATRINE_CO2/C24_CO2/OutputDir/GEOSChem.sat_track.'+f"{prev_year}{prev_month:02d}"+'*z.nc4')
    elif expname=='CATRINE_C360_CO2':
        fname_gchp=glob.glob('/nobackupp17/lwu5/GCHP/GCHP_run/CATRINE_CO2/C360_CO2/OutputDir/GEOSChem.sat_track.'+f"{prev_year}{prev_month:02d}"+'*z.nc4')  
    else:
        fname_gchp=glob.glob('/nobackupp17/lwu5/GCHP/GCHP_run/gchp_merra2_carbon_CO2/'+expname+'/OutputDir/GEOSChem.sat_track.'+f"{prev_year}{prev_month:02d}"+'*z.nc4')
    print(fname_gchp)
    ds_gchp=xr.open_dataset(fname_gchp[0])
    ds_gchp['time'] = pd.to_datetime(ds_gchp['time'].values).round('s')
    return ds_gchp

#Read satellite and collocated GCHP xCO2 data
def read_xco2_monthly(year,month,expname):
    #Observations
    fold_obs='/nobackupp17/jliu7/OCO2/BAKER/B11.2_L2#/'
    ds_obs=read_oco_monthly(year,month,fold_obs)
    #GCHP
    ds_mod=read_gchp_monthly(year,month,expname)
    co2_mod=ds_mod['SpeciesConcVV_CO2'].transpose("time", "lev")
    #prs_mod=ds_mod['Met_PMID'].transpose("time", "lev")
    prs_mod=ds_mod['Met_PMIDDRY'].transpose("time", "lev")
    if ds_obs is not None:
        prs_mod=prs_mod*units.hPa
        co2_mod=co2_mod*units.K #fake unit
        #ps_mod=prs_mod[:,0] #GC surface pressure 
        prs_obs=ds_obs['pressure']*units.hPa
        #ps_obs=prs_obs[:,0]*units.hPa
        #Linear interpolation
        nobs,nlev=prs_obs.shape
        co2_interp = np.empty((nobs, nlev), dtype=np.float32)
        # Loop over each profile
        for i in range(nobs):
            co2_interp[i] = interpolate_1d(prs_obs[i], prs_mod[i], co2_mod[i])
            # Replace NaN at the first interpolated level
            if np.isnan(co2_interp[i, 0]):
                co2_interp[i, 0] = co2_mod[i, 0]
        #--------------------------------------------------------------
        # Apply GOS observation operator
        #
        #   x_hat = x_a + A_k ( x_m - x_a ) 
        #  
        #  where  
        #    x_hat = GC modeled column as seen by TES [vmr]
        #    x_a   = GOS apriori column               [vmr]
        #    x_m   = GC modeled column                [vmr]
        #    A_k   = GOS averaging kernel 
        #--------------------------------------------------------------
        co2_pert=co2_interp-ds_obs['CO2-apriori']
        xco2_hat=ds_obs['xCO2-apriori']+np.sum(ds_obs['xCO2-averagingKernel']*co2_pert,axis=1)
        xco2_obs=ds_obs['xCO2']
        lat=ds_obs['latitude']
        lon=ds_obs['longitude']
        time=ds_obs['time']
        return xco2_obs,xco2_hat,lat,lon,time
    else:
        return None, None, None, None, None

def read_xco2_yearly(year, expname):
    xco2_obs_list = []
    xco2_hat_list = []
    lat_list = []
    lon_list = []
    time_list = []

    for month in range(1, 13):
        xco2_obs, xco2_hat, lat, lon, time = read_xco2_monthly(year, month, expname)
        if xco2_obs is not None:
            xco2_obs_list.append(xco2_obs)
            xco2_hat_list.append(xco2_hat)
            lat_list.append(lat)
            lon_list.append(lon)
            time_list.append(time)

    if xco2_obs_list:
        xco2_obs_all = np.concatenate(xco2_obs_list)*1e6 #ppm
        xco2_hat_all = np.concatenate(xco2_hat_list)*1e6 #ppm
        lat_all = np.concatenate(lat_list)
        lon_all = np.concatenate(lon_list)
        time_all = np.concatenate(time_list)
        
        return xco2_obs_all, xco2_hat_all, lat_all, lon_all, time_all
    else:
        return None, None, None, None, None


def calculate_error(obs,gchp):
    # Compute RMS and bias
    bias = np.mean(gchp - obs)
    rms = np.sqrt(np.mean((gchp - obs) ** 2))
    return bias,rms

def plot_hist2d(expname,region,year,obs,gchp,bias,rms):
    #############
    # Hist 2D
    #############
    # Define histogram binning
    #bins = 100
    #vmax=480
    #vmin=380
    #bins=45
    #vmax=430
    #vmin=385
    #bins=460
    #vmax=800
    #vmin=340
    bins=100
    #vmax=obs.max()
    #vmin=obs.min()
    vmax=430
    vmin=400
    range_ = [[vmin, vmax], [vmin, vmax]]
    # Compute 2D histogram
    hist, xedges, yedges = np.histogram2d(obs, gchp, bins=bins, range=range_)
    # Convert to percentage
    hist_percent = (hist / hist.sum()) * 100
    # Plot
    plt.figure(figsize=(8, 6))
    mesh = plt.pcolormesh(xedges, yedges, hist_percent.T, norm=LogNorm(vmin=1e-3, vmax=hist_percent.max()), cmap='viridis')
    plt.xlabel("Observed xCO₂ (ppm)",fontsize=16)
    plt.ylabel(expname+" xCO₂ (ppm)",fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    #plt.title(f"{year}{month:02d}", fontsize=16)
    plt.title(f"{expname}, {region}", fontsize=16)
    plt.plot([vmin, vmax], [vmin, vmax], 'r--', label='1:1 line')
    #plt.colorbar(mesh, label="Percentage (%)")
    cbar = plt.colorbar(mesh)
    cbar.set_label("Percentage (%)", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    # Fit a line: obs = a * gchp + b
    fit_mask = (gchp >= vmin) & (gchp <= vmax) & (obs >= vmin) & (obs <= vmax)
    #a, b = np.polyfit(gchp[fit_mask], obs[fit_mask], 1)
    a, b, r_value, p_value, std_err = linregress(gchp[fit_mask], obs[fit_mask])
    print(a,b)
    # Add fitted line to plot
    xfit = np.linspace(vmin, vmax, 100)
    yfit = a * xfit + b
    plt.plot(xfit, yfit, 'r-', label='Fit', linewidth=2)
    plt.xlim(vmin,vmax)
    plt.ylim(vmin,vmax)
    # Annotate RMS, bias, and fit equation
    plt.text(
        0.05, 0.95,
        f"Bias: {float(bias):.2f} ppm\nRMSE: {float(rms):.2f} ppm\nOBS = {a:.2f} × GCHP + {b:.2f}",
        fontsize=14,
        verticalalignment='top',
        transform=plt.gca().transAxes,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray')
    )
    ## Annotate RMS and bias
    #plt.text(
    #    0.05, 0.95,  # x, y in axes-relative coordinates (top-left)
    #    f"Bias: {float(bias):.2f} ppmv\nRMSE: {float(rms):.2f} ppmv",
    #    fontsize=14,
    #    verticalalignment='top',
    #    transform=plt.gca().transAxes,
    #    bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray')
    #)
    plt.tight_layout()
    #plt.savefig("Figure/CO2_2dhist_below2km_"+expname+"_"+f"{year}{month:02d}"+".png", dpi=300,bbox_inches='tight')
    #plt.savefig("Figure/CO2_2dhist_"+expname+"_"+f"{year}"+".png", dpi=300,bbox_inches='tight')
    plt.savefig("Figure/CO2_2dhist_OCO2_"+expname+"_"+region+"_"+f"{year}"+".png", dpi=300,bbox_inches='tight')
    plt.close()


def plot_hist2d_1exp(expname,region,year):
    #year=2021

    #date_start='20150106'
    #date_end='20150107'

    #expname='CO2C24_GridFED_smoothed'
    xco2_obs,xco2_hat,lat,lon,time=read_xco2_yearly(year,expname)
    #with open('data/XCO2_obs_'+str(year)+'_'+expname+'.pkl', 'wb') as file:
    with open('data/XCO2_obs_dry_sattype_'+str(year)+'_'+expname+'.pkl', 'wb') as file:
        pickle.dump((xco2_obs,xco2_hat,lat,lon,time), file)
    #with open('data/XCO2_obs_2021_'+expname+'.pkl', 'rb') as file:
    #    xco2_obs,xco2_hat,lat,lon,time=pickle.load(file)
    
    if region=='Tropics':
        indx=np.where(abs(lat)<30)[0]
        xco2_obs=xco2_obs[indx]
        xco2_hat=xco2_hat[indx]
    elif region=='SH':
        indx=np.where(lat<-30)[0]
        xco2_obs=xco2_obs[indx]
        xco2_hat=xco2_hat[indx]
    elif region=='NH':
        indx=np.where(lat>30)[0]
        xco2_obs=xco2_obs[indx]
        xco2_hat=xco2_hat[indx]
    elif region=='EastAsia':
        indx = np.where((lat > 10) & (lat < 50) & (lon > 100) & (lon < 140))[0]
        xco2_obs=xco2_obs[indx]
        xco2_hat=xco2_hat[indx]

    bias,rmse=calculate_error(xco2_obs,xco2_hat)
    plot_hist2d(expname,region,year,xco2_obs,xco2_hat,bias,rmse)



if __name__ == '__main__':
    #expnames=['CO2C24_GridFED_smoothed','CO2C90_GridFED_smoothed','CO2C360_GridFED_smoothed','CO2C24_GridFED','CO2C90_GridFED','CO2C360_GridFED']
    #regions=['Global','Tropics','NH','SH']
    #expnames=['C24_GridFed','C90_GridFed','C360_GridFed','C24_GridFed_smooth','C90_GridFed_smooth','C360_GridFed_smooth']
    #expnames=['C90_4x5']
    #expnames=['CATRINE_C90_CO2','CATRINE_C180_CO2','CATRINE_C24_CO2']
    expnames=['CATRINE_C360_CO2']
    regions=['Global']
    year=2023
    for expname in expnames:
        for region in regions:
            plot_hist2d_1exp(expname,region,year)



