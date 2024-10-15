#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import sys

####### INPUTS #######
# Data location:
#in_loc = path
#save = path
#qcfil = filepath to qc file

#start='202305050000'
#stop='202306140000'

#Example usage: 
# python process_radar.py $in_loc $start $stop $save $qcfil
#############################################################

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

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    try:
        dpath = args_in[1]
        save = args_in[4]
        start = dt.datetime.strptime(args_in[2],'%Y%m%d%H%M')
        stop = dt.datetime.strptime(args_in[3],'%Y%m%d%H%M')
        artifacts = pd.read_csv(args_in[5],parse_dates=[0,1],date_format='%Y-%m-%d %H:%M:%S')
    except:
        print('Input error')
        sys.exit()
    # return values:
    return dpath,start,stop,save,artifacts


def main():
    """
    Applies artifact removal and flagging to cloudnetpy netcdfs of  
    RPG-FMCW-94-SP/DP 94 GHz W-band Cloud Doppler Radar data. 
    Adds additional metadata and formatting
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Cloudnetpy radar filepath
        save:  Output directory (optional)
        
    Returns:
        Modified radar xarray dataframe
    
    """
    
    # check / get args:
    dpath,start,stop,save,artifacts = get_args(sys.argv)
    mask_vars=['Zh', 'v','width','zdr']
    all_dates = pd.date_range(start,stop)

    # Generate the artifact mask and mask the processed_nc files
    for i in range(0,len(all_dates)):
        day = all_dates[i].to_pydatetime()
        print(day)
        #radar = all_radar.sel(time=slice(day,day+dt.timedelta(hours=23,minutes=59,seconds=59,microseconds=999999)))
        radar = xr.open_dataset(dpath + '%s_radar.nc'%dt.datetime.strftime(day,'%y%m%d'))
    
        # find artifacts for this day
        artifacts_subset=artifacts.loc[(artifacts['date1']>=day) & (artifacts['date2']<=day+dt.timedelta(days=1))]
    
        # generate the mask (ones are good, twos are artifacts)
        radar = radar.assign(quality_flag_artifact = (['time','range'],np.ones(np.shape(radar['Zh']))))
    
        # add meta data
        radar.attrs['comment'] = 'Radar artifacts have been removed by visual inspection. See quality_flag_artifact. Contact Heather Guy for further information'
        radar.attrs['creator_name'] =	'Heather Guy'
        radar.attrs['creator_email'] =	'heather.guy@ncas.ac.uk'
        radar.attrs['creator_url'] =	'https://orcid.org/0000-0003-3525-0766'
        radar.attrs['platform'] =	'Swedish Icebreaker Oden'
        radar.attrs['platform_type'] =	'moving_platform'
        radar.attrs['deployment_mode'] =	'ship'
        radar.attrs['title'] =	'RPG-FMCW-94 cloud radar from ARTofMELT 2023'
        radar.attrs['featureType'] =	'timeSeries'

        radar.attrs['time_coverage_start'] = str(radar.time.min().data)
        radar.attrs['time_coverage_end'] = str(radar.time.max().data)
        radar.attrs['geospatial_bounds'] = "%sN, %sE, %sN, %sE"%(str(radar.latitude.max().data),str(radar.longitude.min().data),str(radar.latitude.min().data),str(radar.longitude.max().data))
        radar.attrs['platform_altitude'] = "Oden 4th deck, ~20 m a.s.l"

        # Make sure the height is correct assuming a height amsl of 20m on the 4th deck
        radar.altitude.data=np.array(20.)
        radar.height.data = radar.range.data + 20.0
        
        radar.attrs['location_keywords'] = "Oden, Arctic Ocean, Fram Strait, atmosphere, profile, on the ship"
        radar.attrs['date_created'] = dt.datetime.strftime(dt.datetime.now(),'%d %b %Y %H:%M')
        radar.attrs['institution']="Stockholm University and the University of Leeds"
        radar.attrs['processing_software'] = "cloudnetpy version 1.46.4"
        radar.attrs['sampling_interval'] = "5 s"
        radar.attrs['product_version'] = "v01"
        radar.attrs['project'] = "ARTofMELT"
        radar.attrs['project_principal_investigator'] = "Michael Tjernström"
        radar.attrs['project_principal_investigator_email'] = "michaelt@misu.su.se"
        radar.attrs['project_principal_investigator_url'] = "https://orcid.org/0000-0002-6908-7410"
        radar.attrs['additional_creator_name'] = "John Prytherch"
        radar.attrs['additional_creator_email'] = "john.prytherch@misu.su.se"
        radar.attrs['additional_creator_url'] = "https://orcid.org/0000-0003-1209-289X"

        # add decimal doy
        radar = radar.assign(day_of_year = (['time'],[decimaldayofyear(t) for t in radar.time.data], {'units':"1",'long_name':"Day of Year",'description':"time as decimal day of year"}))

        # Update latitude and longitude using ship data
        # Get ship data
        ship_path='/Users/heather/Desktop/ARTofMELT/data/shipdata/'
        fils = glob.glob(ship_path+'*.csv')
        fils.sort()
        ship_data = pd.concat([pd.read_csv(fil,parse_dates=[0],index_col=0) for fil in fils])
        lats = ship_data[' Oden.Ship.LatitudeDegreesFixed'].reindex(radar.time,method='nearest',tolerance='5s')
        lons = ship_data[' Oden.Ship.LongitudeDegreesFixed'].reindex(radar.time,method='nearest',tolerance='5s')
        radar = radar.assign(latitude = (['time'],lats.to_numpy(), {'units':"degree_north",'long_name':"Latitude of site",'standard_name':"latitude"}))
        radar = radar.assign(longitude = (['time'],lons.to_numpy(), {'units':"degree_east",'long_name':"Longitude of site",'standard_name':"longitude"}))

        
        if len(artifacts_subset)>0:
            print('Masking radar artifacts...')
            for j in range(0,len(artifacts_subset)):
                d1=artifacts_subset.iloc[j]['date1']
                d2=artifacts_subset.iloc[j]['date2']
                h1=artifacts_subset.iloc[j]['height1']
                h2=artifacts_subset.iloc[j]['height2']
                for var in mask_vars:
                    var_atts = radar[var].attrs
                    bad_inds=np.where(radar['quality_flag_artifact'].where((((radar['time']>=d1) & (radar['time']<=d2)) & ((radar['range']>=h1)&(radar['range']<=h2))),other=0).to_numpy()>0)
                    new_var = radar[var].to_numpy()
                    new_var[bad_inds]=np.nan
                    radar[var]=(['time','range'],new_var)
                    # write attributes back in
                    radar[var].attrs = var_atts
                    radar[var].attrs['_FillValue'] = -999.0
                    radar[var].attrs['missing_value'] = -999.0
                    new_qc=radar['quality_flag_artifact'].to_numpy()
                    new_qc[bad_inds]=2
                    radar['quality_flag_artifact']=(['time','range'],new_qc)
    
        radar['quality_flag_artifact'].attrs = {'units' :1,
                                            'long_name' :'Quality flag indicating where radar artifacts have been removed',
                                            'definition':'flag value 1 = good data, flag value 2 = artifacts have been removed'}
        radar['time'].attrs['_FillValue']=False
     
        # Save the masked files
        print('Saving masked netcdf..')
        radar.to_netcdf(save+'RPGFMCW94_ARTofMELT_%s_radar_v1.nc'%(dt.datetime.strftime(day,'%y%m%d')))
    return radar



    
if __name__ == '__main__':
    main()          

