#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs November 28 17:22:21 2023

@author: Heather Guy
"""

import sys,os
sys.path.append(os.getcwd())
import numpy as np      
import datetime as dt
import pandas as pd
import os
import sys

from NC_functions_v1 import *
from surface_met_functions import *
from netCDF4 import Dataset, date2num
import glob

import warnings
warnings.filterwarnings("ignore")

####### INPUTS #######
# Data location:
#in_loc = '/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/'
#out_loc = '/gws/nopw/j04/ncas_radar_vol1/heather/AoM2023/ice-station-data/final_nc/'
# start_1 = '202305160000'
# stop_1 = '202305220000'
# start_2 = '202305290000'
# stop_2 = '202306120000'

#Example usage: 
# python parse_surface-met.py $in_loc $out_loc $start_1 $stop_1
#############################################################

def valminmax(ncf,varn):  
    if isinstance(ncf.variables[varn].valid_min,str):
        ncf.variables[varn].valid_min=np.nanmin(ncf.variables[varn][:])
        ncf.variables[varn].valid_max=np.nanmax(ncf.variables[varn][:])
    if np.isnan(ncf.variables[varn].valid_min):
        ncf.variables[varn].valid_min=np.nanmin(ncf.variables[varn][:])
        ncf.variables[varn].valid_max=np.nanmax(ncf.variables[varn][:])
    if np.nanmin(arr)< ncf.variables[varn].valid_min:
        ncf.variables[varn].valid_min=np.nanmin(ncf.variables[varn][:])
    try:
        if max(arr)> ncf.variables[varn].valid_max:
            ncf.variables[varn].valid_max=np.nanmax(ncf.variables[varn][:])
    except:
        if max([arr])> ncf.variables[varn].valid_max:
            ncf.variables[varn].valid_max=np.nanmax([ncf.variables[varn][:]])   

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    # get values:
    in_loc = args_in[1]
    out_loc = args_in[2]

    # return values:
    return in_loc,out_loc

def main():
    """
    main function generates netCDF and stores in out_loc
    """

    # check / get args:
    try:
        in_loc,out_loc,start,stop = get_args(sys.argv)
    except:
        print('Input error')
        sys.exit()

    # start and end time
    time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='1min')[:]

    # Global attributes
    meta_f = os.getcwd()+'/meta-data/surface-met-metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = os.getcwd()+'/meta-data/surface-met.xlsx'
    var = pd.read_excel(var_f)
    

    # Get GPS data

    print('Extracting GPS data...')
    gps_fils1 = glob.glob(in_loc+'*.GPS')
    gps_fils2 = glob.glob(in_loc+'radiometer_stand/*.gps')
    gps1,met_latlon = get_gps(glob.glob(in_loc+'*.GPS'))
    gps2,rad_latlon = get_gps(glob.glob(in_loc+'radiometer_stand/*.gps'))
    # Apply correction for step change 
    met_corrected=correct_gps(met_latlon,['longitude','latitude']).reindex(time_list,method='nearest',tolerance='1min')
    rad_corrected = correct_gps(rad_latlon,['longitude','latitude']).reindex(time_list,method='nearest',tolerance='1min')


    # Get air pressure (from licor)

    print('Getting pressure form licor...')
    lifils = glob.glob(in_loc+'licor_processed/licor_*')
    licor=get_licor(lifils)
    pressure = licor['P'].resample('1min').mean().reindex(time_list,method='nearest',tolerance='1min')
    # QC pressure
    qc_pressure = np.ones(len(pressure))
    # Sanity check 
    qc_pressure[np.where(np.asarray(pressure)<90000)[0]]=2
    qc_pressure[np.where(np.asarray(pressure)>110000)[0]]=2  


    # Get HMP110 data

    print('Extracting HMP110 data...')
    hmp_fils = glob.glob(in_loc+'*.HMP110')
    trh=(get_hmp(hmp_fils)).reindex(time_list,method='nearest',tolerance='1min')
    #qc_flag_relative_humidity
    qc_flag_relative_humidity=np.ones(len(trh))
    qc_flag_relative_humidity[np.where(np.asarray(trh['qc'])==3)[0]]=2
    #qc_flag_temperature
    qc_flag_temperature=np.ones(len(trh))
    qc_flag_temperature[np.where(np.asarray(trh['qc'])==2)[0]]=2


    # Get 3D sonic data
    # Processed 3D sonic data should already be corrected to a 
    # right-handed coordinate sytesm
    # +ve z is upwards, +ve x is 'northwards', +ve y is 'eastward'

    print('Extracting metek and calculating winds...')
    metekfils = glob.glob(in_loc+'metek_processed/metek_*')
    metek=get_metek(metekfils)
    # Do the 1 minute mean of the vector wind components
    # get the 1 minute mean metek data first
    m_vector=metek[['x','y','z']].resample('1min').mean().reindex(time_list,method='nearest',tolerance='1min')

    # Calculate windspeed and direction from 3D sonic data
    m_vector['fwd_azimuth']=np.nan
    m_vector['gps_distance']=np.nan
    m_vector['wind_speed']=np.nan
    m_vector['wind_from_direction']=np.nan
    for i in range(0,len(m_vector)):
        datetime = m_vector.index[i]

        # Derive heading correction
        # Metek 'North' axis is oriented along the direction from the met mast
        # to the radiometer stand. 
        # GPS ellipsoid is WGS84
        # units are degrees and m
        fwd_azimuth,distance = get_sonic_direction(datetime,met_corrected,rad_corrected)

        # flag as bad if horizontal distance > 100 m 
        if distance < 100 : 
            m_vector.loc[datetime,'fwd_azimuth']=fwd_azimuth
            m_vector.loc[datetime,'gps_distance']=distance
        else: 
            m_vector.loc[datetime,'fwd_azimuth']=np.nan
            m_vector.loc[datetime,'gps_distance']=np.nan
        
        # fwd_azimuth is the direction the x-axis of the sonic is pointing (clockwise from north, in degrees)
        # So to correct, we need to rotate anticlockwise by the fwd_azimuth value
        vector = [m_vector.loc[datetime]['x'],m_vector.loc[datetime]['y'],m_vector.loc[datetime]['z']]
        #   rotation correction by rotating around z
        rot = rotate_around(vector,'z',fwd_azimuth)

        # calculate windspeed
        m_vector.loc[datetime,'wind_speed']=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)

        # Calculate wind direction from u,v
        # remember we want wind 'from' direction (function accounts for this)
        wdir_relative = get_windd(rot[0],rot[1])
        m_vector.loc[datetime,'wind_from_direction']=wdir_relative

    # Wind qc flags
    # 0 not_used
    # 1 good_data 
    # 2 suspect_data_measured_wind_speed_==_0_m_s-1
    # 3 bad_data_wind_speed_outside_sensor_operational_range
    # 4suspect_data_time_stamp_error
    qc_flag_wind_speed = np.ones(len(m_vector))
    qc_flag_wind_from_direction = np.ones(len(m_vector))
    qc_flag_wind_speed[np.where(np.asarray(m_vector['wind_speed'])==0)[0]]=2
    qc_flag_wind_from_direction[np.where(np.asarray(m_vector['wind_speed'])==0)[0]]=2
    qc_flag_wind_speed[np.where(np.asarray(m_vector['wind_speed'])>50)[0]]=3
    qc_flag_wind_from_direction[np.where(np.asarray(m_vector['wind_from_direction'])<0)[0]]=3
    qc_flag_wind_from_direction[np.where(np.asarray(m_vector['wind_from_direction'])>360)[0]]=3


    # Get radiation data

    print('Extracting radiometer data...')
    rad_fils= glob.glob(in_loc+'radiometer_stand/*_CR1000_Kipp-and-Zonen_radiation_KnZRad.dat')
    qcfil = 'radiometer_qc.txt'
    radiation = get_rad(rad_fils,qcfil)
    rad_one_min = radiation.resample('1min').mean().reindex(time_list,method='nearest',tolerance='1min')
    downwelling_longwave_flux_in_air = rad_one_min['CGR_up']
    downwelling_shortwave_flux_in_air = rad_one_min['CMP_up']
    downwelling_total_irradiance = rad_one_min['CGR_up']+rad_one_min['CMP_up']
    upwelling_longwave_flux_in_air = rad_one_min['CGR_dn']
    upwelling_shortwave_flux_in_air = rad_one_min['CMP_dn']
    upwelling_total_irradiance = rad_one_min['CGR_dn'] + rad_one_min['CMP_dn']
    # +ve towards the surface 
    net_total_irradiance  = downwelling_total_irradiance - upwelling_total_irradiance# +ve towards the surface 
    # QC
    #1:good_data 
    #2:bad_data_longwave_radiation_outside_sensor_operational_range
    #3:bad_data_shortwave_radiation_outside_sensor_operational_range
    #4:suspect_data_instrument_log 
    #5:suspect_data_time_stamp_error
    qc_flag_downwelling_radiation = np.ones(len(rad_one_min))
    qc_flag_upwelling_radiation = np.ones(len(rad_one_min))
    # if not 1 from get_rad, then the qc is from the instrument log. 
    qc_flag_downwelling_radiation[np.where(np.asarray(rad_one_min['qc'])!=1)[0]] = 3                       
    qc_flag_upwelling_radiation[np.where(np.asarray(rad_one_min['qc'])!=1)[0]] = 3                       


    # Get HFP and snow temperature profile

    print('Extracting HFP data...')
    hfpfils=[dloc+'ice-station-1/IceStation1_CR1000_Ice_T-Profile_TC_HFP.dat',dloc+'ice-station-2/IceStation2_CR1000_Ice_T-Profile_TC_HFP.dat']
    ice_to_snow_heat_flux,qc_ice_to_snow_heat_flux,thermistor_string = get_hf(hfpfils)


    # Get KT15 data
    # K
    # skin_temperature_1 = KT15 (met)
    # skin_temperature_2 = KT15-1 (radiometer)
    # qc_flag_skin_temperature

    print('Extracting KT15 data...')
    kt1_fils = glob.glob(in_loc+'*.KT15')
    kt2_fils = glob.glob(in_loc+'radiometer_stand/*.KT15-1')
    kt1_fils.sort()
    kt2_fils.sort()
    kt1,kt1_amb,kt1_qc = get_kt15(kt1_fils,1)
    kt2,kt2_amb,kt2_qc = get_kt15(kt2_fils,2)      
 

    # Set up netcdf files.

    print('Writing netcdf file.')
            
    f1 = 'icestation'    #instrument
    f2 = 'ARTofMELT' #platform name  
    f3 = '%s_%s'%(dt.datetime.strftime(start,'%Y%m%d'),dt.datetime.strftime(stop,'%Y%m%d'))
    f5 = "v1" #version number
    f6 = ".nc"
    fn = out_loc + 'surface-met/' +f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'surface-met' + chr(95) + f5 + f6
    nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

    NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
    NC_Dimensions(nc, len(time_list), index=6)  
    NC_CommonVariables(nc, time_list, np)
    NC_SpecificVariables(nc, var, np)

    # Write in data

    nc.variables['longitude'][:]=met_corrected.to_numpy()
    nc.variables['latitude'][:]=met_corrected.to_numpy()
    nc.variables['air_pressure'][:]=pressure.to_numpy()/100 # hPa
    nc.variables['air_temperature'][:]=trh['temperature'].to_numpy()+273.15 # Kelvin
    nc.variables['relative_humidity'][:]==trh['rh'].to_numpy()
    nc.variables['wind_speed'][:]=m_vector['wind_speed'].to_numpy()
    nc.variables['wind_from_direction'][:]=m_vector['wind_from_direction'].to_numpy()
    nc.variables['downwelling_longwave_flux_in_air'][:]=downwelling_longwave_flux_in_air.to_numpy()
    nc.variables['downwelling_shortwave_flux_in_air'][:]=downwelling_shortwave_flux_in_air.to_numpy()
    nc.variables['downwelling_total_irradiance'][:]=downwelling_total_irradiance.to_numpy()
    nc.variables['upwelling_longwave_flux_in_air'][:]=upwelling_longwave_flux_in_air.to_numpy()
    nc.variables['upwelling_shortwave_flux_in_air'][:]=upwelling_shortwave_flux_in_air.to_numpy()
    nc.variables['upwelling_total_irradiance'][:]=upwelling_total_irradiance.to_numpy()
    nc.variables['net_total_irradiance'][:]=net_total_irradiance.to_numpy()
    nc.variables['ice_to_snow_heat_flux'][:]=ice_to_snow_heat_flux.reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['skin_temperature_1'][:]=kt1.reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['skin_temperature_2'][:]=kt2.reindex(time_list,method='nearest',tolerance='1min').to_numpy()

    nc.variables['snow_temperature'][:,0]=thermistor_string[-0.02].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['snow_temperature'][:,1]=thermistor_string[-0.04].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['snow_temperature'][:,2]=thermistor_string[-0.08].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['snow_temperature'][:,3]=thermistor_string[-0.16].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['snow_temperature'][:,4]=thermistor_string[-0.32].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['snow_temperature'][:,5]=thermistor_string[-0.64].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['height_relative_to_snow_surface'][:]=np.asarray([-0.02,-0.04,-0.08,-0.16,-0.32,-0.64])

    nc.variables['qc_flag_temperature'][:]=qc_flag_temperature
    nc.variables['qc_flag_relative_humidity'][:]=qc_flag_relative_humidity
    nc.variables['qc_flag_pressure'][:]=qc_pressure
    nc.variables['qc_flag_wind_speed'][:]=qc_flag_wind_speed
    nc.variables['qc_flag_wind_from_direction'][:]=qc_flag_wind_from_direction
    nc.variables['qc_flag_downwelling_radiation'][:]=qc_flag_downwelling_radiation
    nc.variables['qc_flag_upwelling_radiation'][:]=qc_flag_upwelling_radiation
    nc.variables['qc_ice_to_snow_heat_flux'][:]=pd.DataFrame(index=ice_to_snow_heat_flux.index,data=qc_ice_to_snow_heat_flux).reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['qc_snow_temperature'][:]=thermistor_string['qc'].reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['qc_flag_skin_temperature_1'][:]=pd.DataFrame(index=kt1.index,data=kt1_qc).reindex(time_list,method='nearest',tolerance='1min').to_numpy()
    nc.variables['qc_flag_skin_temperature_2'][:]=pd.DataFrame(index=kt2.index,data=kt2_qc).reindex(time_list,method='nearest',tolerance='1min').to_numpy()
                 
    # Calculation valid max and min
    valminmax(nc,'longitude')
    valminmax(nc,'latitude')
    valminmax(nc,'air_pressure')
    valminmax(nc,'air_temperature')
    valminmax(nc,'relative_humidity')
    valminmax(nc,'wind_speed')
    valminmax(nc,'wind_from_direction')
    valminmax(nc,'downwelling_longwave_flux_in_air')
    valminmax(nc,'downwelling_shortwave_flux_in_air')
    valminmax(nc,'downwelling_total_irradiance')
    valminmax(nc,'upwelling_longwave_flux_in_air')
    valminmax(nc,'upwelling_shortwave_flux_in_air')
    valminmax(nc,'upwelling_total_irradiance')
    valminmax(nc,'net_total_irradiance')
    valminmax(nc,'ice_to_snow_heat_flux')
    valminmax(nc,'snow_temperature')
    valminmax(nc,'height_relative_to_snow_surface')
    valminmax(nc,'skin_temperature_1')
    valminmax(nc,'skin_temperature_2')

    # Optional add a comment
    #base_str = 'Example string'
    #nc.setncattr('comment', base_str)

    # Close netcdf file
    nc.close()
    
if __name__ == '__main__':
    main()
