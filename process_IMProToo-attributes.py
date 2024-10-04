#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 2023

@author: heather

Update attributes in IMProToo output file

Inputs: 
    in_dir: IMProToo directory
    in_dat: Date 'MMDD'
    out_dir: Output location

"""

import sys,os
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
    
    YYYYMMDD='2023%s'%(in_dat)
    radar = xr.open_dataset(out_dir+'ncas-mrr-1_ARTofMELT_%s_IMProToo_v0.nc'%(YYYYMMDD))

    # Add latitude and longitude (lets pull them from the radar files?)
    try: 
        RPG_radar = xr.open_dataset('/Users/heather/Desktop/ARTofMELT/data/radar/'+'RPGFMCW94_ARTofMELT_%s_radar_v1.nc'%(YYYYMMDD[2:]))
        radar['latitude'] = RPG_radar.latitude
        radar['longitude'] = RPG_radar.longitude  
    except:
        print('Cant find radar file to update latitude and longitude')
        return
        
    # Add altitude 
    radar = radar.assign(altitude =np.array(10.))
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
    radar.attrs['month']=str(pd.to_datetime(radar.time[0].data).month).zfill(2)
    radar.attrs['day']=str(pd.to_datetime(radar.time[0].data).day).zfill(2)
    
    # add decimal doy
    radar = radar.assign(day_of_year = (['time'],[decimaldayofyear(t) for t in radar.time.data], {'units':"1",'long_name':"Day of Year",'description':"time as decimal day of year"}))

    # Edit existing attributes
    del radar.attrs['contact_person']
    del radar.attrs['author']
    # rename 'properties' to 'algorithm_settings'
    radar.attrs['IMProToo_settings'] = radar.attrs.pop('properties')
    radar.attrs['source'] = "NCAS AMOF METEK MRR2 Micro Rain Radar (ncas-mrr-1)"

    # Edit quality description
    radar['quality'].attrs = {'units':"bin",'long_name':"Quality control flag (binary), indicating the status of the IMProToo algoirthm.",'description':"Quality flag consists of 18-bits, flag 1 is the least significant. Flags 1-5 can usually be ignored. Flags 8-12 indicate reasons why a spectrum does not contain a peak. Flags 16-18 are serious errors, do not use these data. See flag_meanings for the definition of each bit.",'flag_meanings':'1: spectrum interpolated around 0 and 12 m s^-1 \n2: peak streches over interpolated part  \n3: peak is dealiased  \n4: first Algorithm to determine peak failed, used backup  \n5: dealiasing went wrong, but is corrected  \n6: Not used \n7: Not used \n8: spectrum was incompletely recorded  \n9: the variance test indicated no peak  \n10: spectrum is not processed due to according setting  \n11: peak removed since not wide enough  \n12: peak removed, because too few neighbours show signal, too  \n13: Not used \n14: Not used \n15: Not used \n16: peak is at the very border to bad data  \n17: in this area there are still strong velocity jumps, indicates failed dealiasing  \n18: during dealiasing, a warning was triggered, applied to whole columm \n'}

    # Save
    print('Saving netcdf v1..')
    radar.to_netcdf(out_dir+'ncas-mrr-1_ARTofMELT_%s_IMProToo_v1.nc'%(YYYYMMDD))

if __name__ == '__main__':
    main()  