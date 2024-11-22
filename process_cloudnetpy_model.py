#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import sys
from glob import glob
import metpy.calc as mpcalc
from metpy.units import units
from model_utils import *

####### INPUTS #######
# Data location:
#in_loc = path
#mwps = path
#save = path

#start='yymmdd'
#stop='yymmdd'

#Example usage: 
# python process_cloudnetpy_model.py $in_loc $mwps $save $start $stop
#############################################################

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    try:
        dpath = args_in[1]
        mwps_path = args_in[2]
        save = args_in[3]
        start = dt.datetime.strptime(args_in[4],'%y%m%d')
        stop = dt.datetime.strptime(args_in[5],'%y%m%d')
        calc_atten=args_in[6]
    except:
        print('Input error')
        sys.exit()
    # return values:
    return dpath,start,stop,save,mwps_path,calc_atten

def decimaldayofyear(date):
    '''
    this function returns the decimal day of year.
    Input is a date in datetime64 format, changed into datetime object here.
    The minimum time is microseconds.
    eg. date = numpy.datetime64('2023-05-09T11:33:02.639000000'), output: 129.481385 (with 6 digit precision),
    output given in fractional day of year, with respect to 2023 (1 Jan 2023 = day 1)
    '''
    seconds_in_a_day=24*60*60
    total_seconds = (pd.Timestamp(date).hour*60*60) + (pd.Timestamp(date).minute*60) + pd.Timestamp(date).second + pd.Timestamp(date).microsecond*10**(-6)
    fractional_day_of_year =np.round(pd.Timestamp(date).dayofyear + total_seconds/seconds_in_a_day,6)
    return fractional_day_of_year


def main():
    """
    Gets corrected hatpro temperature and humidity data
    Winds from interpolated radiosondes
    Calculates radar attenuation characteristics using mwps (optional)
    Generates the 'model' file for use in cloudnetpy. 
    Adds additional metadata and formatting
    
    Parameters:
        dpath: Path to input data
        mwps_path: Path to mwps model
        save:  Output directory
        start: Start datetime for processing
        stop:  Stop datetime for processing
        calc_atten: True or false, calculate radar attenuation characteristics using mwps?
        
    Returns:
        None
    
    """
    
    # check / get args:
    dpath,start,stop,outdir,mwps_path,calc_atten = get_args(sys.argv)
    all_dates = pd.date_range(start,stop)

    # Get mwr data
    model_in_fil = dpath+'mwr/Hatpro_profiles_ARTofMELT_20230507_20230614_v01.nc'
    in_fil_name = model_in_fil.split('/')[-1]
    print(in_fil_name)
    hatpro=xr.open_dataset(model_in_fil)

    # Get radiosonde data (for pressure and winds)
    radiosondes = xr.open_dataset(dpath+'soundings/radiosonde_ARTofMELT_20230509_20230613_OnCommonGrid_v01.nc')
    radiosonde_fil = 'radiosonde_ARTofMELT_20230509_20230613_OnCommonGrid_v01.nc'
    print(radiosonde_fil)

    # Process one day at a time (easier for mwps)
    for date in all_dates:
        day_str = dt.datetime.strftime(date,'%y%m%d')
        print(day_str)
        t_corrected = hatpro['air_temperature_corrected'].sel(time=slice(date,date+dt.timedelta(hours=23,minutes=59,seconds=59)))
        rh_corrected = hatpro['relative_humidity_corrected'].sel(time=slice(date,date+dt.timedelta(hours=23,minutes=59,seconds=59)))
        
        # Convert rh
        # Needs to be with respect to ice below 273.16
        # See eq 2.16 in Rodgers & Yau, T is in K
        rh = rh_corrected.where(t_corrected>273.16,other=((rh_corrected/100) * ((273 / t_corrected)**2.66)*100))
         
        # interpolate radiosonde data to 1-minutely resolution
        pressure = radiosondes['air_pressure'].sel(time=slice(date-dt.timedelta(hours=3),date+dt.timedelta(hours=29))).resample(time='1min').interpolate('linear')
        u = radiosondes['eastward_wind'].sel(time=slice(date-dt.timedelta(hours=3),date+dt.timedelta(hours=29))).resample(time='1min').interpolate('linear')
        v = radiosondes['northward_wind'].sel(time=slice(date-dt.timedelta(hours=3),date+dt.timedelta(hours=29))).resample(time='1min').interpolate('linear')
        
        # reindex pressure to match hatpro time/height. 
        pressure = pressure.reindex(time=t_corrected.time,height=t_corrected.height,method='nearest')
        u = u.reindex(time=t_corrected.time,height=t_corrected.height,method='nearest')
        v = v.reindex(time=t_corrected.time,height=t_corrected.height,method='nearest')
        
        # Calculate the specific humidity from the hatpro files, use units kg/kg
        #metpy.calc.dewpoint_from_relative_humidity(temperature, relative_humidity)
        dewpoint = mpcalc.dewpoint_from_relative_humidity(t_corrected.to_numpy()* units.K, rh_corrected.to_numpy()* units.percent)
        #metpy.calc.specific_humidity_from_dewpoint(pressure, dewpoint)
        sphum= mpcalc.specific_humidity_from_dewpoint(pressure * units.hPa, dewpoint).to('kg/kg').magnitude
        
        # start the new dataset
        model = xr.Dataset(coords={'time':t_corrected.time,'level':t_corrected.height.to_numpy(),'frequency':[35,94]})
        
        # height needs to be 2D
        # temperautre in kelvin
        # pressure in Pa
        model = model.assign(variables={'height':(('time','level'), np.tile(t_corrected.height, (len(t_corrected.time), 1))),
                       'temperature':(('time','level'),t_corrected.to_numpy()),
                            'pressure':(('time','level'),pressure.to_numpy()*100),
                               'rh':(('time','level'),rh.to_numpy()/100),
                               'uwind':(('time','level'),u.to_numpy()),
                               'vwind':(('time','level'),v.to_numpy()),
                               'q':(('time','level'),sphum),})

        if calc_atten==True: 
            print('Calculating attenuations..')
            # Run the model to calculate radar attenuation
            # mwps-0.8.1
            # Run propagation ccal
            # propagation: reads frequency, temperature, pressure and 
            # humidity and returns propagation parameters
            save_loc='./temp/'
            frequencies = np.asarray([35,94])
            temperatures = model['temperature'].to_numpy()
            pressures = model['pressure'].to_numpy()
            rhs = model['rh'].to_numpy()
            print('Generating input for mwps..')
            gen_prop_input(frequencies,temperatures,pressures,rhs,save_loc,overwrite=True)
            print('Running mwps')
            [gas_atten_1way,liq_atten_1way,k2] = run_prop_cal(save_loc+'atmos.temp',mwps_path)
            [gas_dry_1way,liq_dry_1way,k2_dry] = run_prop_cal(save_loc+'atmos_dry.temp',mwps_path)
            [gas_sat_1way,liq_sat_1way,k2_sat] = run_prop_cal(save_loc+'atmos_sat.temp',mwps_path)
            
            # First reshape the output back onto the time-height-frequency demain
            specific_gas_1way=np.reshape(gas_atten_1way,[len(frequencies),len(model.time),len(model.level)])
            specific_liq_1way=np.reshape(liq_atten_1way,[len(frequencies),len(model.time),len(model.level)])
            
            # Fill the highest level with a duplication 
            specific_gas_1way[:,:,-1] = specific_gas_1way[:,:,-2]
            specific_liq_1way[:,:,-1] = specific_liq_1way[:,:,-2]
            
            gas_attenuation=calc_2way_gas_atten(specific_gas_1way,frequencies,model.time.to_numpy(),model.level.to_numpy())
            liquid_attenuation=calc_2way_gas_atten(specific_liq_1way,frequencies,model.time.to_numpy(),model.level.to_numpy())
            
            # Add new variables 
            model['gas_atten']=(('frequency','time', 'level'), gas_attenuation)
            model.gas_atten.attrs['missing_value'] = np.nan
            model.gas_atten.attrs['long_name'] = "Two-way attenuation from the ground due to atmospheric gases"
            model.gas_atten.attrs['units'] = "dB"
            model.gas_atten.attrs['_FillValue'] = np.nan
            
            model['K2']=(('frequency','time', 'level'), np.reshape(k2,[len(frequencies),len(model.time),len(model.level)]))
            model.K2.attrs['units_html'] = "db km<sup>-1</sup>"
            model.K2.attrs['missing_value'] = np.nan
            model.K2.attrs['_FillValue'] = np.nan
            model.K2.attrs['long_name'] = "Dielectric parameter (|K|^2) of liquid water"
            model.K2.attrs['units'] = "dB km-1"
            
            model['specific_dry_gas_atten']=(('frequency','time', 'level'), np.reshape(gas_dry_1way,[len(frequencies),len(model.time),len(model.level)]))
            model.specific_dry_gas_atten.attrs['units_html'] = "db km<sup>-1</sup>"
            model.specific_dry_gas_atten.attrs['missing_value'] = np.nan
            model.specific_dry_gas_atten.attrs['long_name'] = "Specific one-way attenuation due to atmospheric gases for dry air (no water vapour)"
            model.specific_dry_gas_atten.attrs['units'] = "dB km-1"
            model.specific_dry_gas_atten.attrs['_FillValue'] = np.nan
            
            model['specific_gas_atten']=(('frequency','time', 'level'), specific_gas_1way)
            model.specific_gas_atten.attrs['long_name'] = "Specific one-way attenuation due to atmospheric gases"
            model.specific_gas_atten.attrs['units'] = "dB km-1"
            model.specific_gas_atten.attrs['units_html'] = "db km<sup>-1</sup>"
            model.specific_gas_atten.attrs['missing_value'] = np.nan
            model.specific_gas_atten.attrs['_FillValue'] = np.nan
            
            model['specific_liquid_atten']=(('frequency','time', 'level'), specific_liq_1way)
            model.specific_liquid_atten.attrs['units_html'] = "(dB km<sup>-1</sup>)/(g m<sup>-3</sup>)"
            model.specific_liquid_atten.attrs['missing_value'] = np.nan
            model.specific_liquid_atten.attrs['units'] = "(dB km-1)/(g m-3)"
            model.specific_liquid_atten.attrs['_FillValue'] = np.nan
            model.specific_liquid_atten.attrs['long_name'] = "Specific one-way attenuation due to liquid water, per unit liquid water content"
            
            model['specific_saturated_gas_atten']=(('frequency','time', 'level'), np.reshape(gas_sat_1way,[len(frequencies),len(model.time),len(model.level)]))
            model.specific_saturated_gas_atten.attrs['long_name'] = "Specific one-way attenuation due to atmospheric gases for saturated air (saturated with respect to ice below 0 degrees C)"
            model.specific_saturated_gas_atten.attrs['units'] = "dB km-1"
            model.specific_saturated_gas_atten.attrs['units_html'] = "db km<sup>-1</sup>"
            model.specific_saturated_gas_atten.attrs['missing_value'] = np.nan
            model.specific_saturated_gas_atten.attrs['_FillValue'] = np.nan
        
        # Add attributes
        # label as cloudnet model file
        model.attrs['cloudnet_file_type'] = 'obs'
        model.attrs['year'] = date.year
        model.attrs['month'] = date.month
        model.attrs['day'] = date.day

        model.height.attrs['units']='m'
        model.height.attrs['_FillValue'] = np.nan
        model.height.attrs['long_name'] = "Height above ground"
        model.height.attrs['units'] = "m"
        model.height.attrs['missing_value'] = np.nan
        model.height.attrs['standard_name'] = "height"
        model.height.attrs['positive'] = "up"
        
        model=model.assign(time = np.array([(pd.to_datetime(d.data)-date).total_seconds()/60/60 for d in model.time]))
        model.time.attrs['long_name'] = "Hours UTC";
        model.time.attrs['units'] = "hours since %s"%(dt.datetime.strftime(date,'%Y-%m-%d %H:%M:%S +00:00'))
        model.time.attrs['standard_name'] = "time"
        model.time.attrs['axis'] = "T"
        model.time.attrs['calendar'] = "standard"
        
        model.temperature.attrs['units'] = "K"
        model.temperature.attrs['long_name'] = "Temperature"
        model.temperature.attrs['_FillValue'] = np.nan
        model.temperature.attrs['missing_value'] = np.nan
        
        model.pressure.attrs['long_name'] = "pressure"
        model.pressure.attrs['units'] = "Pa"
        model.pressure.attrs['_FillValue'] = np.nan
        model.pressure.attrs['missing_value'] = np.nan
        
        model.rh.attrs['long_name'] = "Relative humidity"
        model.rh.attrs['units'] = "1"
        model.rh.attrs['missing_value'] = np.nan
        model.rh.attrs['_FillValue'] = np.nan
        model.rh.attrs['comment'] = "With respect to liquid above 0 degrees C and with respect to ice below 0 degrees C."
        model.rh.attrs['standard_name'] = "relative_humidity"
        model['relative_humidity']= (('time', 'level'),model.rh.to_numpy())
        model.relative_humidity.attrs['_FillValue'] = np.nan
        model.relative_humidity.attrs['missing_value'] = np.nan
        
        model.uwind.attrs['long_name']='east wind component'
        model.uwind.attrs['units']='m/s'
        model.uwind.attrs['_FillValue'] = np.nan
        model.uwind.attrs['missing_value'] = np.nan
        
        model.vwind.attrs['long_name']='north wind component'
        model.vwind.attrs['units']='m/s'
        model.vwind.attrs['_FillValue'] = np.nan
        model.vwind.attrs['missing_value'] = np.nan
        
        model.q.attrs['long_name']='specific humidity'
        model.q.attrs['units']='kg/kg'
        model.q.attrs['_FillValue'] = np.nan
        model.q.attrs['missing_value'] = np.nan

        # Fill all missing values with -999
        model = model.fillna(np.nan)
        
        # Add some global attributes, source, location, ect. 
        model.attrs['creator_name'] =	'Heather Guy'
        model.attrs['creator_email'] =	'heather.guy@ncas.ac.uk'
        model.attrs['creator_url'] =	'https://orcid.org/0000-0003-3525-0766'
        model.attrs['platform'] =	'Swedish Icebreaker Oden'
        model.attrs['platform_type'] =	'moving_platform'
        model.attrs['deployment_mode'] =	'ship'
        model.attrs['title'] =	'Atmospheric temperature and humidity profiles from ARTofMELT 2023, for input into cloudnetpy'
        model.attrs['featureType'] =	'timeSeries'
        model.attrs['time_coverage_start'] = dt.datetime.strftime(date + dt.timedelta(hours=float(model.time.min().data)),'%d %b %Y, %H:%M UTC')
        model.attrs['time_coverage_end'] = dt.datetime.strftime(date + dt.timedelta(hours=float(model.time.max().data)),'%d %b %Y, %H:%M UTC')
        model.attrs['geospatial_bounds'] = hatpro.geospatial_bounds 
        model.attrs['platform_altitude'] = "Oden 4th deck, ~20 m a.s.l"
        model.attrs['location_keywords'] = "Oden, Arctic Ocean, Fram Strait, atmosphere, profile, on the ship"
        model.attrs['date_created'] = dt.datetime.strftime(dt.datetime.now(),'%d %b %Y %H:%M')
        model.attrs['institution']="Stockholm University and the University of Leeds"
        model.attrs['processing_software'] = "mwps-0.8.1; https://github.com/heatherguy/AOM-WP1-BL/blob/main/process_cloudnetpy_model.py"
        model.attrs['sampling_interval'] = "1 minute"
        model.attrs['product_version'] = "v01"
        model.attrs['project'] = "ARTofMELT"
        model.attrs['project_principal_investigator'] = "Michael Tjernström"
        model.attrs['project_principal_investigator_email'] = "michaelt@misu.su.se"
        model.attrs['project_principal_investigator_url'] = "https://orcid.org/0000-0002-6908-7410"
        model.attrs['additional_creator_name'] = "Sonja Murto"
        model.attrs['additional_creator_email'] = "sonja.murto@misu.su.se "
        model.attrs['additional_creator_url'] = "https://orcid.org/0000-0002-4966-9077"
        model.attrs['input_file_reference'] = 'Michael Tjernström, Sonja Murto, Michail Karalis, John Prytherch (2024) Temperature and humidity profiles from microwave radiometer during expedition ARTofMELT, Arctic Ocean, 2023. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-microwave-profiles-1\nSonja Murto, Michael Tjernström, Ian Brooks, Heather Guy, Michail Karalis, Timo Vihma, Gabin Urbancic, John Prytherch (2024) Radiosonde profiles from expedition ARTofMELT, Arctic Ocean, 2023. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-radiosonde-1'
        if calc_atten==True:
            model.attrs['history'] = "Input file %s created on %s. Input file %s created on %s. Radar attenuation parameters calculated using mwps-0.8.1. Names modified to match the correct format for cloudnetpy input."%(in_fil_name,hatpro.last_revised_date, radiosonde_fil,radiosondes.last_revised_date)
            model.attrs['comment'] = 'This file consists of vertical atmospheric profiles of air temperature and humidity as well as radar attenuation characteristics from the location of the Icebreaker Oden during the ARTofMELT campaign. The file has been specifically formatted as an input *model* file for the cloudnetpy retrieval algorithm. The data are derived from the temperature and humidity retrievals from an RPG-HATPRO radiometer that has had a bias correction applied by interpolating 6-hourly radiosonde profiles. See the *input_file_reference* for more information regarding the bias correction. Wind and pressure data included in the file are just from the interpolated radiosonde data, and the attenuation parameters were calculated from the temperature and humidity profiles using mwps version 0.8.1 by Robin Hogan, 2002.'
        else:
            model.attrs['history'] = "Input file %s created on %s. Input file %s created on %s. Names modified to match the correct format for cloudnetpy input."%(in_fil_name,hatpro.last_revised_date, radiosonde_fil,radiosondes.last_revised_date)
            model.attrs['comment'] = 'This file consists of vertical atmospheric profiles of air temperature and humidity from the location of the Icebreaker Oden during the ARTofMELT campaign. The file has been specifically formatted as an input *model* file for the cloudnetpy retrieval algorithm. The data are derived from the temperature and humidity retrievals from an RPG-HATPRO radiometer that has had a bias correction applied by interpolating 6-hourly radiosonde profiles. See the *input_file_reference* for more information regarding the bias correction. Wind and pressure data included in the file are just from the interpolated radiosonde data.'            
        
        model.attrs['product_name'] = "RPG-HATPRO radiometer temperature and humidity profiles with (6-hourly) radiosonde bias correction"
        model.attrs['source'] = "RPG-HATPRO radiometer; 6-hourly Vaisala Digicora 5.4.0 radiosounding system with RS41 sondes"
        model.attrs['location'] = "Arctic Ocean; Fram Strait"
        
        # add decimal doy
        model = model.assign(day_of_year = (['time'],[decimaldayofyear(t) for t in model.time.data], {'units':"1",'long_name':"Day of Year",'description':"time as decimal day of year"}))
        model['time'].attrs['_FillValue']=False
        
        print('Saving file cloudnet_profile_obs_ARTofMELT_%s_V01.nc'%dt.datetime.strftime(date,'%y%m%d'))
        model.to_netcdf(outdir+'cloudnet_profile_obs_ARTofMELT_%s_V01.nc'%dt.datetime.strftime(date,'%y%m%d'))
   
    return

if __name__ == '__main__':
    main()          
