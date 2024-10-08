#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 2023

@author: heather

Process MRR raw spectra with IMProToo
Generate netcdf file

Inputs: 
    in_dir: Raw spectra directory
    in_dat: Date of raw input file, format 'MMDD'
    out_dir: Output location

Example usage: 
   cd /Users/heather/IMProToo/
   source venv/bin/activate
   python MRR_process.py $in_dir $in_dat $out_dir >> /Users/heather/IMProToo/output/logs/${in_dat}_log.txt 2>&1
"""

import sys,os
import IMProToo
import matplotlib.pyplot as plt
import xarray as xr
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
import pandas as pd
#from matplotlib.colors import DivergingNorm

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
    in_dir=sys.argv[1]
    in_dat=sys.argv[2]
    out_dir=sys.argv[3]
    
    print('Loading data...')
    rawData = IMProToo.mrrRawData(str(in_dir)+str(in_dat) + '.raw')
    processedSpec = IMProToo.MrrZe(rawData)
    
    processedSpec.co["ncCreator"] = "Heather Guy, University of Leeds, heather.guy@ncas.ac.uk"
    processedSpec.co["ncDescription"] = 'Micro rain radar data processed with IMProToo from expedition ARTofMELT, Arctic Ocean, 2023'
    processedSpec.co["ncInstitution"] = "University of Leeds and the National Centre for Atmospheric Science"
    processedSpec.co["ncLocation"] = "Ice breaker Oden foredeck, during the ARTofMELT field campaign. 78.13 to 80.52 degrees N, -3.87 to 14.14 degrees E."

    processedSpec.co["dealiaseSpectrum"] = True
    processedSpec.rawToSnow()
    
    print('Writing data...')
    # Write processed netcdf file
    YYYYMMDD='2023%s'%(in_dat)
    processedSpec.writeNetCDF(out_dir+'ncas-mrr-1_ARTofMELT_%s_IMProToo_v0.nc'%(YYYYMMDD),ncForm="NETCDF3_CLASSIC")

    # Modify netcdf attributes and resave
    print('Updating attributes..')
    radar = xr.open_dataset(out_dir+'ncas-mrr-1_ARTofMELT_%s_IMProToo_v0.nc'%(YYYYMMDD))

    # Add latitude and longitude (lets pull them from the radar files?)
    RPG_radar = xr.open_dataset('/Users/heather/Desktop/ARTofMELT/data/radar/'+'RPGFMCW94_ARTofMELT_%s_radar_v1.nc'%(YYYYMMDD[2:]))
    radar['latitude'] = RPG_radar.latitude
    radar['latitude'] = RPG_radar.longitude   

    # Add altitude 
    radar = radar.assign(altitude =np.array(20.))
    radar['altitude'].attrs = {'units':"m",'long_name':"Altitude of MRR above mean sea level",'standard_name':"altitude"}

    # add meta data
    #radar.attrs['comment'] = 'Radar artifacts have been removed by visual inspection. See quality_flag_artifact. Contact Heather Guy for further information'
    radar.attrs['creator_name'] =	'Heather Guy'
    radar.attrs['creator_email'] =	'heather.guy@ncas.ac.uk'
    radar.attrs['creator_url'] =	'https://orcid.org/0000-0003-3525-0766'
    radar.attrs['platform'] =	'Swedish Icebreaker Oden'
    radar.attrs['platform_type'] =	'moving_platform'
    radar.attrs['deployment_mode'] =	'ship'
    radar.attrs['featureType'] =	'timeSeries'
    radar.attrs['time_coverage_start'] = str(radar.time.min().data)
    radar.attrs['time_coverage_end'] = str(radar.time.max().data)
    radar.attrs['geospatial_bounds'] = "%sN, %sE, %sN, %sE"%(str(radar.latitude.max().data),str(radar.longitude.min().data),str(radar.latitude.min().data),str(radar.longitude.max().data))
    radar.attrs['platform_altitude'] = "Oden foredeck, ~10 m a.s.l"
    radar.attrs['location_keywords'] = "Oden, Arctic Ocean, Fram Strait, atmosphere, profile, on the ship, radar, snowfall, MRR"
    radar.attrs['date_created'] = dt.datetime.strftime(dt.datetime.now(),'%d %b %Y %H:%M')
    radar.attrs['processing_software'] = "IMProToo v0.107"
    radar.attrs['sampling_interval'] = "10s"
    radar.attrs['product_version'] = "v01"
    radar.attrs['project'] = "ARTofMELT"
    radar.attrs['project_principal_investigator'] = "Michael Tjernström"
    radar.attrs['project_principal_investigator_email'] = "michaelt@misu.su.se"
    radar.attrs['project_principal_investigator_url'] = "https://orcid.org/0000-0002-6908-7410"
    radar.attrs['Conventions'] =	'CF-1.8'
    radar.attrs['year']=str(2023)
    radar.attrs['month']=str(pd.to_datetime(mrr.time[0].data).month).zfill(2)
    radar.attrs['day']=str(pd.to_datetime(mrr.time[0].data).day).zfill(2)
    
    # add decimal doy
    radar = radar.assign(day_of_year = (['time'],[decimaldayofyear(t) for t in radar.time.data], {'units':"1",'long_name':"Day of Year",'description':"time as decimal day of year"}))

    # Edit existing attributes
    del radar.attrs['contact_person']
    del radar.attrs['author']
    # rename 'properties' to 'algorithm_settings'
    radar.attrs['IMProToo_settings'] = mrr.attrs.pop('properties')
    radar.attrs['source'] = "NCAS AMOF METEK MRR2 Micro Rain Radar (ncas-mrr-1)"

    # Save
    print('Saving masked netcdf..')
    radar.to_netcdf(out_dir+'ncas-mrr-1_ARTofMELT_%s_IMProToo_v1.nc'%(YYYYMMDD))

if __name__ == '__main__':
    main()  