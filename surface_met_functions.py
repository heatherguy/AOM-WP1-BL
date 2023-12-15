#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs November 28 17:22:21 2023

@author: heather

Various functions for processing/ calculating icestation surface met data from ARTofMELT
"""

import sys,os
import numpy as np      
import datetime as dt
import pandas as pd
import xarray as xr
from netCDF4 import Dataset, date2num
import glob
import pyproj

def dms2dd(degrees, minutes, seconds, direction):
    # Converts degrees minutes to decimal degrees
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
    if direction == 'S' or direction == 'W':
        dd *= -1
    return dd

def get_gps(infils):
    # Extract GPS data from files and apply basic qc. 
    infils.sort()
    
    header=[0,
            'week',
            'seconds',
            'fix_date_ddmmyy',
            'fix_time_hhmmss',
            'leap_second_count',
            'latitude',
            'lat_hem',
            'longitude',
            'lon_hem',
            'mode',
            'fix_type',
            'speed_kmh',
            'course_deg',
            'position_dilution_of_precision',
            'time_dilution_of_prescision']
    
    all_pdfs=[]
    for fil in infils:
        all_pdfs.append(pd.read_csv(fil,delim_whitespace=True,parse_dates=[[0,1,2,3,4,5]],index_col=0, date_format='%Y %m %d %H %M %S.%f',header=None,dtype='str'))
    
    gps = pd.concat(all_pdfs)
    gps.columns=header
    # basic qc
    gps['qc']=np.ones(len(gps))
    # must have a fix
    gps['fix_type'] = pd.to_numeric(gps['fix_type'],errors='coerce')
    gps.loc[gps['fix_type']==0,'qc']=2
    gps.loc[np.isnan(gps['fix_type'].to_numpy()),'qc']=2
    
    # convert lat/lon
    lats=np.empty(len(gps))*np.nan
    lons = np.empty(len(gps))*np.nan
    for i in range(0,len(gps)):
        # ddmm.mmmm
        try:
            dlat=int(gps['latitude'].iloc[i][0:2])
            mlat=float(gps['latitude'].iloc[i][2:]) 
            lats[i]=dms2dd(dlat,mlat,0,gps['lat_hem'].iloc[i])
        except:
            lats[i]=np.nan
            continue
        
        # dddmm.mmmm
        try:
            dlon=int(gps['longitude'].iloc[i][0:3])
            mlon=float(gps['longitude'].iloc[i][3:]) 
            lons[i]=dms2dd(dlon,mlon,0,gps['lon_hem'].iloc[i])
        except:
            lons[i]=np.nan
            continue
    
    latlon=pd.DataFrame(index=gps.index,data={'latitude':lats,'longitude':lons})
    latlon['qc']=gps['qc']
    
    # Sanity check Latitude > 40, Longitude < 90
    latlon.loc[latlon['latitude']<40,'qc']=3
    latlon.loc[latlon['longitude']>90,'qc']=3

    # 1 min mean, ignoring bad qc
    latlon.loc[latlon['qc']!=1,['latitude','longitude']]=np.nan
    latlon=latlon.resample('1min').mean()
    
    return gps,latlon

def get_hmp(infils):
    # Extract HMP temperature and humidity data from files
    # and apply a basic qc. 
    
    infils.sort()
    all_pdfs=[]
    for fil in infils:
        try:
            all_pdfs.append(pd.read_csv(fil,delim_whitespace=True,parse_dates=[[0,1,2,3,4,5]],index_col=0, date_format='%Y %m %d %H %M %S.%f',header=None))
        except:
            print('Could not parse %s'%fil)
        
    hmp = pd.concat(all_pdfs)
    trh=pd.DataFrame({'temperature':hmp[7],'rh':hmp[10]})
    trh['qc']=np.ones(len(trh))
    
    # qc for realistic values
    trh.loc[trh['temperature']<-50,'qc']=2
    trh.loc[trh['temperature']>30,'qc']=2  
    trh.loc[trh['rh']<0,'qc']=3
    trh.loc[trh['rh']>110,'qc']=3 
    
    return trh

def get_rad(infils,qcfil):
    # Extract radiometer data from Campbell logger
    # files and apply a base qc. 
    
    headers = ["TIMESTAMP","RECORD","BattV","PanelT","CMP_up","CMP_dn","CGR_up","CGR_dn","CMP_up_TsK","CMP_dn_TsK","CGR_up_TsK","CGR_dn_TsK"]

    all_pdfs=[]
    for fil in infils:
        all_pdfs.append(pd.read_csv(fil,skiprows=4,names=headers,index_col=0,parse_dates=[0]))

    rad_dat = pd.concat(all_pdfs)
    rad_dat.sort_index(inplace=True)
    rad_dat['CMP_up'] = rad_dat['CMP_up'].astype(float)
    rad_dat['CMP_dn'] = rad_dat['CMP_dn'].astype(float)
    rad_dat['CGR_up'] = rad_dat['CGR_up'].astype(float)
    rad_dat['CGR_dn'] = rad_dat['CGR_dn'].astype(float)
    
    # QC for cleaning times
    rad_dat['qc']=np.ones(len(rad_dat))
    qc_fil=pd.read_csv(qcfil,parse_dates=[0,1])

    for i in range(0,len(qc_fil)):
        rad_dat.loc[qc_fil.iloc[i]['qc_start']:qc_fil.iloc[i]['qc_stop'],'qc']=2

    return rad_dat

def get_licor(infils):
    # Load licor data from pre-processed licor data files. 
    
    all_pdfs=[]
    for fil in infils:
        all_pdfs.append(pd.read_csv(fil, index_col=0, parse_dates=[0]))
    
    licor = pd.concat(all_pdfs)
    licor.sort_index(inplace=True)
    
    return licor

def get_metek(infils):
    # Load metek data from pre-processed metek data files. 
    
    all_pdfs=[]
    for fil in infils:
        all_pdfs.append(pd.read_csv(fil, index_col=0, parse_dates=[0]))
    
    metek = pd.concat(all_pdfs)
    metek.sort_index(inplace=True)
    
    return metek

def get_sonic_direction(datetime,met_gps,rad_gps):
    # sonic +ve x-axis (North) is oriented along the heading from met_gps 
    # to rad_gps. 
    # This function calculates that heading (fwd_azimuth from met_gps to
    # rad gps) and returns the fwd_azimuth and the horizontal distance 
    # (in m) between met_gps and rad_gps
    # GPS ellipsoid is WGS84
   
    geodesic = pyproj.Geod(ellps='WGS84')
    # Check we have both GPS points or exit with error
    try: 
        met_gps.loc[datetime]['longitude']
    except:
        #print('Missing met gps')
        return np.nan, np.nan
    try: 
        rad_gps.loc[datetime]['longitude']
    except:
        #print('Missing rad gps')
        return np.nan, np.nan
    
    # units are degrees and m
    fwd_azimuth,back_azimuth,distance = geodesic.inv(met_gps.loc[datetime]['longitude'], met_gps.loc[datetime]['latitude'], rad_gps.loc[datetime]['longitude'], rad_gps.loc[datetime]['latitude'])
    
    return fwd_azimuth,distance

def get_windd(uwind,vwind):
    """
    Converts u and v wind components into meteorological degrees.
    meteorological degrees = clockwise, north =0, direction wind is coming from.
    
    Parameters:
        uwind: wind in u direction (north)
        vwind: wind in v direction (east)
        
    Returns:
        dir: Direction wind is comping from in degrees
    
    """   
    dir = 180 + ((np.arctan2(uwind/(np.sqrt(uwind**2 + vwind**2)), vwind/(np.sqrt(uwind**2 + vwind**2)))) * 180/np.pi) # degrees
    
    return dir

def deg_rot(orig,diff):
    """
    Adds a fixed number of degrees to a wind direction. 
    To account for sonics not being oriented to north. 
    
    Parameters:
        orig: list or series of wind directions in degrees. 
        diff: number of degrees clockwise to rotate (can be negative).
        
    Returns:
        new: Corrected direction wind is comping from in degrees
    
    """   
    new = orig
    new = new + diff
    new[new<0] = 360 + (orig[new<0] + diff)
    new[new>360] = diff - (360 - orig[new>360])
    new[new==360] = 0

    return new

def rotate_around(vector,axis,theta):
    # Using RHR convention, and rotates anticlockwise
    # convert theta to radians
    theta = np.deg2rad(theta)
    
    # Rx - rotates anticlockwise around the x axis by angle theta
    rx =np.asarray([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])

    # Ry - rotates anticlockwise around the y axis by angle theta
    ry=np.asarray([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])

    # Rz - rotates anticlockwise around the z axis by angle theta
    rz=np.asarray([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])

    if axis=='x':
        rotated = np.dot(rx, vector)
    elif axis=='y':
        rotated = np.dot(ry, vector)
    elif axis=='z':
        rotated = np.dot(rz, vector)
    else: 
        print('incorrect axis')
        rotated=np.nan
    
    return rotated

def get_hf(infils):
    # Extract heat flux plate and thermistor string data from campbell logger file
    # Apply qc specific to ARTofMELT instrument log (hard coded)
    # Note using the factory calibrated HFP data for consistency (self calibrations in field were inconcsistent)
    # return:  ice_to_snow_heat_flux, qc_ice_to_snow_heat_flux, thermistor_string
    
    headers = ["TIMESTAMP","RECORD","HFP_factory_cal","HFP_self_cal","HeaterFlag","TC(1)","TC(2)","TC(3)","TC(4)","TC(5)","TC(6)"]
    snow_depth=[-0.02,-0.04,-0.08,-0.16,-0.32,-0.64]

    data=[]
    for fil in infils:
        therm_dat = pd.read_csv(fil,skiprows=4,names=headers,index_col=0,parse_dates=[0])
        therm_dat.sort_index(inplace=True)
        therm_dat = therm_dat[~therm_dat.index.duplicated()]
        data.append(therm_dat)

    therm_dat = pd.concat(data)

    # Get HFP data
    # Clean self calibration times
    cal_times = [dt.datetime(2023,6,3,0),dt.datetime(2023,6,5,0),dt.datetime(2023,6,7,0),dt.datetime(2023,6,9,0),dt.datetime(2023,6,11,0)]
    for c in cal_times:
        therm_dat.loc[c:c+dt.timedelta(minutes=60),['HFP_factory_cal', 'HFP_self_cal']]=np.nan

    # Use factory calibration for consistency (self calibrations in field were inconcsistent)
    ice_to_snow_heat_flux = -pd.to_numeric(therm_dat['HFP_factory_cal'],errors='coerce')
    # Resample to 1 min mean
    ice_to_snow_heat_flux = ice_to_snow_heat_flux.resample('1min').mean()

    # QC for known bad times
    #1:good_data 
    #2:bad_data_sensor_operational_range
    #3:suspect_data_instrument_log 
    #4:suspect_data_time_stamp_error
    qc_ice_to_snow_heat_flux=np.ones(len(ice_to_snow_heat_flux))
    qc_ice_to_snow_heat_flux[np.where(np.abs(np.asarray(ice_to_snow_heat_flux))>100)[0]]=2

    dt_index = np.asarray(ice_to_snow_heat_flux.index.to_pydatetime())
    qc_ice_to_snow_heat_flux[np.where(dt_index<dt.datetime(2023,5,17,13))[0]]=3 # Installation
    qc_ice_to_snow_heat_flux[np.where((dt_index>dt.datetime(2023,5,18,2,10))&(dt_index<dt.datetime(2023,5,18,2,30)))[0]]=3 # Installation

    # Get thermistor string data
    thermistor_string=pd.DataFrame(index=therm_dat.index,columns=snow_depth)
    thermistor_string[-0.02][:] = pd.to_numeric(therm_dat['TC(1)'],errors='coerce')+273.15
    thermistor_string[-0.04][:] = pd.to_numeric(therm_dat['TC(2)'],errors='coerce')+273.15
    thermistor_string[-0.08][:] = pd.to_numeric(therm_dat['TC(3)'],errors='coerce')+273.15
    thermistor_string[-0.16][:] = pd.to_numeric(therm_dat['TC(4)'],errors='coerce')+273.15
    thermistor_string[-0.32][:] = pd.to_numeric(therm_dat['TC(5)'],errors='coerce')+273.15
    thermistor_string[-0.64][:] = pd.to_numeric(therm_dat['TC(6)'],errors='coerce')+273.15
    thermistor_string['qc']=np.ones(len(thermistor_string))

    # sanity check
    thermistor_string.loc[(thermistor_string>283).any(axis=1),'qc']=2   
    thermistor_string.loc[((thermistor_string<223)&(thermistor_string>3)).any(axis=1),'qc']=2

    # QC for instrument log
    # polar bear destruction
    thermistor_string.loc[dt.datetime(2023,6,10,21,15):dt.datetime(2023,6,11,17),'qc']=3
    # suspect data 
    thermistor_string.loc[dt.datetime(2023,5,18,4):dt.datetime(2023,5,18,4),'qc']=3
    
    thermistor_string = thermistor_string.resample('1min').mean()
    
    return ice_to_snow_heat_flux,qc_ice_to_snow_heat_flux,thermistor_string

def get_kt15(fils,qcf):
    # Get KT15 data from files
    # qc for sanity values and instrument log
    # qcf = 1 : use the log from the met mast to apply qc
    # qcf = 2 : use the log from the radiometer mast to apply qc
    # return skin temperature, ambient temperature, and skin temperature qc
    amb=[]
    skin=[]
    for fil in fils: 
        try:
            #print(fil)
            dat = pd.read_csv(fil,parse_dates=[[0,1,2,3,4,5]],header=None,delim_whitespace='True',index_col=0, date_format='%Y %m %d %H %M %S.%f',names=[0,1,2,3,4,5,'skin_t','skin_unit','amb_t','amb_unit'])
        except:
            print('Could not parse %s'%fil)
            continue

        amb.append(pd.to_numeric(dat[dat['skin_t']=='AMB']['amb_t'],errors='coerce') + 273.15)
        skin.append(pd.to_numeric(dat[dat['skin_t']!='AMB']['skin_t'],errors='coerce') + 273.15)

    amb = pd.concat(amb)
    skin = pd.concat(skin)

    # resample to 1 minute mean: 
    amb = amb.resample('1min').mean()
    skin = skin.resample('1min').mean()

    # basic qc for log and unrealistic data
    # qc for realistic values
    #1: good_data 
    #2: suspect_data_detailed_in_log 
    #3: bad_data_temperature_outside_sensor_operational_range
    skin_qc=skin.copy()*np.zeros(len(skin))+1
    skin_qc.loc[skin<223]=3
    skin_qc.loc[skin>300]=3  

    # qc based on instrument log
    if qcf==1: 
        # Suspect times, kt1: 
        skin_qc.loc[dt.datetime(2023,5,16,19):dt.datetime(2023,5,16,20,10)]=2
        skin_qc.loc[dt.datetime(2023,5,20,10,25):dt.datetime(2023,5,20,10,30)]=2
        skin_qc.loc[dt.datetime(2023,5,31,3):dt.datetime(2023,5,31,16)]=2
        skin_qc.loc[dt.datetime(2023,6,1,13,25):dt.datetime(2023,6,1,16,15)]=2
        skin_qc.loc[dt.datetime(2023,6,5,16,36):dt.datetime(2023,6,5,16,38)]=2
        skin_qc.loc[dt.datetime(2023,6,8,13,35):dt.datetime(2023,6,8,13,40)]=2
    elif qcf==2:
        # Suspect times, kt2: 
        skin_qc.loc[dt.datetime(2023,5,21,18,46):dt.datetime(2023,5,21,18,50)]=2
        skin_qc.loc[dt.datetime(2023,6,2,10,30):dt.datetime(2023,6,2,10,35)]=2

    return skin,amb,skin_qc