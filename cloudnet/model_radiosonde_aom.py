#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday 29 Jan 2024

@author: heather guy

Convert the AOM radiosonde data to "model" format 
required by the cloudnepy alogirthm. 
This includes calculation of atmospheric attenuation 
using the mwps-0.8.1 package.
"""

from glob import glob
import pandas as pd
import numpy  as np
import datetime as dt
import os
import sys
import xarray as xr
import subprocess
sys.path.append(os.getcwd())

from model_utils import *

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    Input: date, infile
    date: in format str YYMMDD
    infile: Path to radiosonde input data
    """    
    date = args_in[1]
    infile = args_in[2]
    try:
        d=dt.datetime.strptime(date,'%y%m%d')
        print('Processing %s'%d)
    except:
        print('Incorrect date format, should be YYMMDD')
        return None
    try:
        radiosondes = xr.open_dataset(infile)
    except:
        print('Cannot open radiosonde file, check path.')
        return None
    return date, radiosondes

def main():
    if get_args(sys.argv) == None:
        print('Input Error')
        sys.exit()
    else:
        date, radiosondes = get_args(sys.argv)

    daily_sonds = radiosondes.sel(time=slice(dt.datetime.strptime(date,'%y%m%d')-dt.timedelta(days=1),dt.datetime.strptime(date,'%y%m%d')+dt.timedelta(days=2)))
    
    # interpolate radiosonde data to 1-hourly resolution then crop to one day 
    daily_sonds_hourly=daily_sonds.resample(time='1H').interpolate('linear').sel(time=slice(dt.datetime.strptime(date,'%y%m%d'),dt.datetime.strptime(date,'%y%m%d')+dt.timedelta(days=1)))

    # rename existing variables
    renamed=daily_sonds_hourly.rename({'u':'uwind','v':'vwind','sphum':'q','height':'level'})

    # add new coords
    starttime=renamed.time.data[0]
    renamed = renamed.assign_coords(frequency=[35,94])
    renamed['time'],start_date = get_time_hours(renamed.time.data)

    # height needs to be 2D
    renamed['height'] = (('time','level'), np.tile(renamed.level, (len(renamed.time), 1)))

    # temp in kelvin
    renamed['temperature']= renamed['temperature']+273.15

    # pressure in pa
    renamed['pressure']=renamed.pressure*100

    # Convert rh
    # Needs to be with respect to ice below 273.16
    renamed['rh']=renamed['RH'].where(renamed['temperature']>273.16,other=renamed['RHice'])/100

    frequencies = np.asarray([35,94])
    temperatures = renamed['temperature'].to_numpy()
    pressures = renamed['pressure'].to_numpy()
    rhs = renamed['rh'].to_numpy()
    save_loc='./data/temp/'
    gen_prop_input(frequencies,temperatures,pressures,rhs,save_loc,overwrite=True)

    # Run the model
    # Run propagation ccal
    # propagation: reads frequency, temperature, pressure and 
    # humidity and returns propagation parameters

    command_loc='./mwps-0.8.1/'
    [gas_atten_1way,liq_atten_1way,k2] = run_prop_cal(save_loc+'atmos.temp',command_loc)
    [gas_dry_1way,liq_dry_1way,k2_dry] = run_prop_cal(save_loc+'atmos_dry.temp',command_loc)
    [gas_sat_1way,liq_sat_1way,k2_sat] = run_prop_cal(save_loc+'atmos_sat.temp',command_loc)

    # First reshape the output back onto the time-height-frequency demain
    specific_gas_1way=np.reshape(gas_atten_1way,[len(frequencies),len(renamed.time),len(renamed.level)])
    specific_liq_1way=np.reshape(liq_atten_1way,[len(frequencies),len(renamed.time),len(renamed.level)])

    gas_attenuation=calc_2way_gas_atten(specific_gas_1way,frequencies,renamed.time.to_numpy(),renamed.level.to_numpy())
    liquid_attenuation=calc_2way_gas_atten(specific_liq_1way,frequencies,renamed.time.to_numpy(),renamed.level.to_numpy())
    
    # Add new variables 

    renamed['gas_atten']=(('frequency','time', 'level'), gas_attenuation)
    renamed.gas_atten.attrs['missing_value'] = -999.9
    renamed.gas_atten.attrs['long_name'] = "Two-way attenuation from the ground due to atmospheric gases"
    renamed.gas_atten.attrs['units'] = "dB"
    renamed.gas_atten.attrs['_FillValue'] = -999.9

    renamed['K2']=(('frequency','time', 'level'), np.reshape(k2,[len(frequencies),len(renamed.time),len(renamed.level)]))
    renamed.K2.attrs['units_html'] = "db km<sup>-1</sup>"
    renamed.K2.attrs['missing_value'] = -999.9
    renamed.K2.attrs['_FillValue'] = -999.9
    renamed.K2.attrs['long_name'] = "Dielectric parameter (|K|^2) of liquid water"
    renamed.K2.attrs['units'] = "dB km-1"

    renamed['specific_dry_gas_atten']=(('frequency','time', 'level'), np.reshape(gas_dry_1way,[len(frequencies),len(renamed.time),len(renamed.level)]))
    renamed.specific_dry_gas_atten.attrs['units_html'] = "db km<sup>-1</sup>"
    renamed.specific_dry_gas_atten.attrs['missing_value'] = -999.9
    renamed.specific_dry_gas_atten.attrs['long_name'] = "Specific one-way attenuation due to atmospheric gases for dry air (no water vapour)"
    renamed.specific_dry_gas_atten.attrs['units'] = "dB km-1"
    renamed.specific_dry_gas_atten.attrs['_FillValue'] = -999.9

    renamed['specific_gas_atten']=(('frequency','time', 'level'), specific_gas_1way)
    renamed.specific_gas_atten.attrs['long_name'] = "Specific one-way attenuation due to atmospheric gases"
    renamed.specific_gas_atten.attrs['units'] = "dB km-1"
    renamed.specific_gas_atten.attrs['units_html'] = "db km<sup>-1</sup>"
    renamed.specific_gas_atten.attrs['missing_value'] = -999.9
    renamed.specific_gas_atten.attrs['_FillValue'] = -999.9

    renamed['specific_liquid_atten']=(('frequency','time', 'level'), specific_liq_1way)
    renamed.specific_liquid_atten.attrs['units_html'] = "(dB km<sup>-1</sup>)/(g m<sup>-3</sup>)"
    renamed.specific_liquid_atten.attrs['missing_value'] = -999.9
    renamed.specific_liquid_atten.attrs['units'] = "(dB km-1)/(g m-3)"
    renamed.specific_liquid_atten.attrs['_FillValue'] = -999.9
    renamed.specific_liquid_atten.attrs['long_name'] = "Specific one-way attenuation due to liquid water, per unit liquid water content"

    renamed['specific_saturated_gas_atten']=(('frequency','time', 'level'), np.reshape(gas_sat_1way,[len(frequencies),len(renamed.time),len(renamed.level)]))
    renamed.specific_saturated_gas_atten.attrs['long_name'] = "Specific one-way attenuation due to atmospheric gases for saturated air (saturated with respect to ice below 0 degrees C)"
    renamed.specific_saturated_gas_atten.attrs['units'] = "dB km-1"
    renamed.specific_saturated_gas_atten.attrs['units_html'] = "db km<sup>-1</sup>"
    renamed.specific_saturated_gas_atten.attrs['missing_value'] = -999.9
    renamed.specific_saturated_gas_atten.attrs['_FillValue'] = -999.9

    # Add attributes
    # label as cloudnet model file
    renamed.attrs['cloudnet_file_type'] = 'model'
    renamed.attrs['year'] = dt.datetime.strptime(date,'%y%m%d').year
    renamed.attrs['month'] = dt.datetime.strptime(date,'%y%m%d').month
    renamed.attrs['day'] = dt.datetime.strptime(date,'%y%m%d').day
    renamed.attrs['location'] = 'AOM Ice Station'

    renamed.height.attrs['units']='m'
    renamed.height.attrs['_FillValue'] = -999.0
    renamed.height.attrs['long_name'] = "Height above ground"
    renamed.height.attrs['units'] = "m"
    renamed.height.attrs['missing_value'] = -999.0
    renamed.height.attrs['standard_name'] = "height"
    renamed.height.attrs['positive'] = "up"

    renamed.time.attrs['long_name'] = "Hours UTC";
    renamed.time.attrs['units'] = "hours since %s"%(dt.datetime.strftime(start_date,'%Y-%m-%d %H:%M:%S +00:00'))
    renamed.time.attrs['standard_name'] = "time"
    renamed.time.attrs['axis'] = "T"
    renamed.time.attrs['calendar'] = "standard"

    renamed.temperature.attrs['units'] = "K"
    renamed.temperature.attrs['long_name'] = "Temperature"

    renamed.pressure.attrs['long_name'] = "pressure"
    renamed.pressure.attrs['units'] = "Pa"

    renamed.rh.attrs['long_name'] = "Relative humidity"
    renamed.rh.attrs['units'] = "1"
    renamed.rh.attrs['missing_value'] = -999.0
    renamed.rh.attrs['comment'] = "With respect to liquid above 0 degrees C and with respect to ice below 0 degrees C."
    renamed.rh.attrs['standard_name'] = "relative_humidity"
    renamed['relative_humidity']= (('time', 'level'),renamed.rh.to_numpy())

    # Save test file
    renamed.to_netcdf('./data/model/%s_obs_sonde.nc'%date)
    print('Saved file %s_obs_sonde.nc'%date)

    return
            
if __name__ == '__main__':   
    main()        