#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday 30 October 2024

@author: heather guy

Generates categorize and classification files
using the cloudnetpy algorithm

"""

import sys
import datetime as dt
from cloudnetpy.categorize import generate_categorize
from cloudnetpy.products import generate_classification
import xarray as xr
import glob
import numpy as np
import pandas as pd
import dask
from scipy import stats
from voodoonet import version

# ST thresholds
#Temperature
h_t_f=-0.3	+273.15	#Temperature distinguishing generic liquid and ice, K
th_t_f_dr=1.0+273.15		#Temperature above which everything must be liquid phase (buffer due to T interpolation), K
th_t_liq=-25.0	+273.15	#Lowest temperature at which all-liquid clouds can exist, K
th_t_hom=-40.0	+273.15	#Homogeneous freezing temperature, no liquid colder, K
#Height
th_z_i=4.0		#Climatological height at T of th_ti, km
th_z_f=1.0		#Climatological height at T of th_tf, km
#Reflectivity
th_db_low=-90.0		#Lowest acceptible reflectivity for hydrometeor return, dBZ
th_db_snow=5.0		#Reflectivity threshold between ice cloud and snow, dBZ
th_db_liq=-17.0 	#Reflectivity threshold between liquid cloud and drizzle/mixed, dBZ
th_db_rain=5		#Lowest acceptible reflectivity for rain, dBZ
#Velocity
th_ve_liq=1.0		#Velocity threshold above which there is no cloud liquid cloud, m/s
th_ve_rain=2.5		#Velocity threshold above which rain occurs, m/s
th_ve_snow=1.0		#Velocity threshold above which snow occurs, m/s
#Spectrum Width
th_sw_liq=0.1		#Spectrum width threshold above which liquid typically occurs, m/s
#if radarstream eq 'mos' or radarstream eq 'nsakazr' or radarstream eq 'oli' or radarstream eq 'olige' then th_sw_liq=0.25
# Using 0.1 for RPG radar based on distribution from ARTofMELT

#LWP (MWR or other)
th_lwp=20.0 /1000.  #LWP above which there must be liquid classified in the column, kg/m2
th_lwp_low = 5.0 /1000. # LWP below which there should be no liquid in the column. 
#if keyword_set(lwpthresh) then th_lwp=lwpthresh   ;Modify the lwp threshold if specified by keyword. 


def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    Input arguments:
    1) date: in format str YYMMDD
    2) radar input data directory
    3) lidar input data directory
    4) lwp input data directory
    5) atmospheric profile input data directory
    6) present weather sensor input data directory
    7) cloudnetpy output directory
    8) artifacts: directory for radar artifacts
    9) Use voodoo? optional LV0 directory
    """   
    date = args_in[1]
    try:
        day=dt.datetime.strptime(date,'%y%m%d')
        print('Processing %s'%day)
    except:
        print('Incorrect date format, should be YYMMDD')
        return None

    radar_dir = args_in[2]
    lidar_dir = args_in[3]
    lwp_dir = args_in[4]
    model_dir = args_in[5]
    pw_dir = args_in[6]
    output_dir = args_in[7]
    artifacts_dir=args_in[8]
    if len(args_in)>9:
        voodoo_dir = args_in[9]
    else:
        voodoo_dir=False
    return date, radar_dir,lidar_dir,lwp_dir,model_dir,pw_dir,output_dir,artifacts_dir,voodoo_dir

#def correct_voodoo_artifacts(cat,outfil):
#    # Input is a cloudnet cateogorize file opened using xarray    
#    # and a string to save the post-processed file. 
#    # This function looks for times where radar artifacts have been
#    # manually filtered but not from the voodoo algorithm, resulting
#    # in spurious liquid water detection from radar artifacts. 
#    # We remove the spurious liquid and resave the categorise file. 
#    
#    # liquid identified from voodoo when liquid_prob > 0.55
#    # Falling hydrometeors are radar signals that are
#    # a) not insects b) not clutter. Furthermore, falling hydrometeors
#    # are strong lidar pixels excluding liquid layers (thus these pixels
#    # are ice or rain). They are also weak radar signals in very cold
#    # temperatures.
#
#    # So, if liquid_prob > 0.55, and falling=0, We know that there is a radar signal
#    # from voodoo, but no good radar signal from the reflectivity. 
#    # Therefore these are the artificats. 
#    # in this case if liq=1, we need to set liq=0
#
#    # decode
#    liq,falling,cold,melting,aerosol,insects=decode_cat(cat['category_bits'].to_numpy())
#
#    # find bad liquid signals from radar artifacts
#    condition = (liq==1) & (falling==0) & (cat['liquid_prob'].to_numpy()>0.55)
#
#    # modify liq
#    mod_liq = liq.copy()
#    mod_liq[condition]=False
#
#    # recode
#    new_cat_bits = recode_cat(mod_liq,falling,cold,melting,aerosol,insects)
#
#    # get original attrs
#    atts = cat['category_bits'].attrs
#    # add a note
#    atts['comment1']='Cleaned for bad voodoo results due to radar artifacts. If liquid_prob>0.55 and falling=0, any liq flags have been removed.'
#
#    # write new cat bits
#    cat['category_bits'] = (('time','height'),new_cat_bits)
#
#    # add attributes
#    cat['category_bits'].attrs.update(atts)
#
#    # save...
#    #cat.to_netcdf(outfil)
#    return cat

def correct_voodoo_artifacts(cat,qcfil):
    # Input is a cloudnet cateogorize file opened using xarray    
    # Also the radar data qc file for masking. 
    # This function looks for times where radar artifacts have been
    # manually filtered but not from the voodoo algorithm, resulting
    # in spurious liquid water detection from radar artifacts. 
    # We remove the spurious liquid and return the updated categorise file. 

    # get radar artifacts
    artifacts = pd.read_csv(qcfil,parse_dates=[0,1],date_format='%Y-%m-%d %H:%M:%S')
    
    # find artifacts for this day
    day=cat.time[0].data
    artifacts_subset=artifacts.loc[(artifacts['date1']>=(day-np.timedelta64(1, 'h'))) & (artifacts['date2']<=day+ np.timedelta64(1, 'D'))]
    
    # decode
    liq,falling,cold,melting,aerosol,insects=decode_cat(cat['category_bits'].to_numpy())

    # find and remove bad liquid signals from radar artifacts
    mod_liq = liq.copy()
    mod_falling=falling.copy()
    if len(artifacts_subset)>0:
        print('Masking radar artifacts...')
        for j in range(0,len(artifacts_subset)):
            d1=artifacts_subset.iloc[j]['date1']
            d2=artifacts_subset.iloc[j]['date2']
            h1=artifacts_subset.iloc[j]['height1']
            h2=artifacts_subset.iloc[j]['height2']
            bad_inds=np.where(cat['category_bits'].where((((cat['time']>=d1) & (cat['time']<=d2)) & ((cat['height']>=h1)&(cat['height']<=h2))),other=0).to_numpy()>0)
            mod_liq[bad_inds]=False
            mod_falling[bad_inds]=False

    # recode
    new_cat_bits = recode_cat(mod_liq,mod_falling,cold,melting,aerosol,insects)

    # get original attrs
    atts = cat['category_bits'].attrs
    # add a note
    atts['comment1']='Cleaned for bad voodoo results due to radar artifacts.'

    # write new cat bits
    cat['category_bits'] = (('time','height'),new_cat_bits)

    # add attributes
    cat['category_bits'].attrs.update(atts)
    return cat

def replace_nan(data, mask):
    result = np.full(mask.shape, np.nan, dtype=float)
    result[~mask] = data
    return result

def decode_cat(arr):
    # Bit 0: Small liquid droplets are present.
    # Bit 1: Falling hydrometeors are present; if Bit 2 is set then these are most
    # likely ice particles, otherwise they are drizzle or rain drops.
    # Bit 2: Wet-bulb temperature is less than 0 degrees C, implying the phase of
    # Bit-1 particles.
    # Bit 3: Melting ice particles are present.
    # Bit 4: Aerosol particles are present and visible to the lidar.
    # Bit 5: Insects are present and visible to the radar.
    
    nan_mask = np.isnan(arr)
    arr_nonan = arr[~nan_mask].astype(int)  # Convert to int after NaN removal

    # Preallocate output arrays
    liq = np.zeros_like(arr, dtype=float)
    falling = np.zeros_like(arr, dtype=float)
    cold = np.zeros_like(arr, dtype=float)
    melting = np.zeros_like(arr, dtype=float)
    aerosol = np.zeros_like(arr, dtype=float)
    insects = np.zeros_like(arr, dtype=float)

    # Extract bits using bitwise operations
    liq[~nan_mask] = (arr_nonan & 0b000001) > 0
    falling[~nan_mask] = (arr_nonan & 0b000010) > 0
    cold[~nan_mask] = (arr_nonan & 0b000100) > 0
    melting[~nan_mask] = (arr_nonan & 0b001000) > 0
    aerosol[~nan_mask] = (arr_nonan & 0b010000) > 0
    insects[~nan_mask] = (arr_nonan & 0b100000) > 0

    # Restore NaNs
    liq[nan_mask] = np.nan
    falling[nan_mask] = np.nan
    cold[nan_mask] = np.nan
    melting[nan_mask] = np.nan
    aerosol[nan_mask] = np.nan
    insects[nan_mask] = np.nan

    return liq.astype('bool'), falling.astype('bool'), cold.astype('bool'), melting.astype('bool'), aerosol.astype('bool'), insects.astype('bool')


def recode_cat(liq, falling, cold, melting, aerosol, insects):
    # Ensure inputs are NumPy arrays
    liq, falling, cold, melting, aerosol, insects = map(np.asarray, (liq, falling, cold, melting, aerosol, insects))

    # Construct integer values using bitwise operations
    cat_bits = ((insects.astype(int) << 5) | 
                (aerosol.astype(int) << 4) | 
                (melting.astype(int) << 3) | 
                (cold.astype(int) << 2) | 
                (falling.astype(int) << 1) | 
                (liq.astype(int)))

    return cat_bits


def decode_quality(arr):
    # Bit 0: An echo is detected by the radar.
    # Bit 1: An echo is detected by the lidar.
    # Bit 2: The apparent echo detected by the radar is ground clutter or some
    #        other non-atmospheric artifact.
    # Bit 3: The lidar echo is due to clear-air molecular scattering.
    # Bit 4: Liquid water cloud, rainfall or melting ice below this pixel
    #        will have caused radar and lidar attenuation; if bit 5 is set then
    #        a correction for the radar attenuation has been performed;
    #        otherwise do not trust the absolute values of reflectivity factor.
    #        No correction is performed for lidar attenuation.
    # Bit 5: Radar reflectivity has been corrected for liquid-water attenuation
    #        using the microwave radiometer measurements of liquid water path
    #        and the lidar estimation of the location of liquid water cloud;
    #        be aware that errors in reflectivity may result.

    nan_mask = np.isnan(arr)
    arr_nonan = arr[~nan_mask].astype(int)  # Convert to int after removing NaNs

    # Preallocate output arrays
    radar = np.zeros_like(arr, dtype=float)
    lidar = np.zeros_like(arr, dtype=float)
    radar_clutter = np.zeros_like(arr, dtype=float)
    lidar_clearair = np.zeros_like(arr, dtype=float)
    attenuated = np.zeros_like(arr, dtype=float)
    atten_corrected = np.zeros_like(arr, dtype=float)

    # Extract bits using bitwise operations
    radar[~nan_mask] = (arr_nonan & 0b000001) > 0
    lidar[~nan_mask] = (arr_nonan & 0b000010) > 0
    radar_clutter[~nan_mask] = (arr_nonan & 0b000100) > 0
    lidar_clearair[~nan_mask] = (arr_nonan & 0b001000) > 0
    attenuated[~nan_mask] = (arr_nonan & 0b010000) > 0
    atten_corrected[~nan_mask] = (arr_nonan & 0b100000) > 0

    # Restore NaNs
    radar[nan_mask] = np.nan
    lidar[nan_mask] = np.nan
    radar_clutter[nan_mask] = np.nan
    lidar_clearair[nan_mask] = np.nan
    attenuated[nan_mask] = np.nan
    atten_corrected[nan_mask] = np.nan

    return radar.astype('bool'), lidar.astype('bool'), radar_clutter.astype('bool'), lidar_clearair.astype('bool'), attenuated.astype('bool'), atten_corrected.astype('bool')

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
    date, radar_dir,lidar_dir,lwp_dir,model_dir,pw_dir,output_dir,artifacts_dir,voodoo_dir = get_args(sys.argv)
    r_fil='RPGFMCW94_ARTofMELT_%s_radar_v1.nc'%date
    l_fil='CL31_cloudnet_ARTofMELT_%s_v01.nc'%date
    mwr_fil='%s_hatpro.nc'%date
    model_fil='cloudnet_profile_obs_ARTofMELT_%s_V01.nc'%date#'%s_obs_sonde.nc'%date
    
    if voodoo_dir: 
        print('Applying voodoo...')
        input_files = {'radar': radar_dir+r_fil,
                    'lidar': lidar_dir+l_fil,
                    'mwr': lwp_dir+mwr_fil,
                    'model': model_dir+model_fil,
                    'lv0_files':glob.glob(voodoo_dir+'%s/%s_*_P06_ZEN.LV0'%(date,date))}
    else: 
        input_files = {'radar': radar_dir+r_fil,
                    'lidar': lidar_dir+l_fil,
                    'mwr': lwp_dir+mwr_fil,
                    'model': model_dir+model_fil} 

    print(input_files)

    cat_fil = output_dir+'categorize/%s_categorize.nc'%date
    cat_uuid = generate_categorize(input_files, cat_fil+'_tmp')
    
    ## Correct voodoo radar artifacts
    if voodoo_dir: 
        print('Correcting voodoo artifacts...')
        cat = xr.open_dataset(cat_fil+'_tmp')
        cat = correct_voodoo_artifacts(cat,artifacts_dir)
    else:
        cat = xr.open_dataset(cat_fil+'_tmp')
    
    # Sort variables and metdata ect
    # Get the cloudnetpy version 
    # Add heathers interation version and comments. 
    
    # Add in present weather sensor visibility (10 min average) 
    # Get present weather sensor data
    vis=xr.open_dataset(pw_dir+'ACAS_AoM2023_PWD_30s_v1_0.nc')
    
    # Adjust insect flag based on visibility 
    print('Adjusting category bits for fog...')
    liq,falling,cold,melting,aerosol,insects=decode_cat(cat['category_bits'].to_numpy())
    
    fog = np.transpose(np.tile((vis['vis_10min'].reindex(time=cat.time,method='nearest',tolerance='45s')<1000).to_numpy(),(np.shape(insects)[1],1)))
    # Turn off insects during fog
    insects[np.where(fog)]=False

    # Add in a horizontal visibility field
    cat['horizontal_visibility'] = (('time'),vis['vis_10min'].reindex(time=cat.time,method='nearest',tolerance='45s').data)
    hv_atts={'long_name': 'Visibility, 10 minute average', 
             'units': 'm',
             'source_file':'ACAS_AoM2023_PWD_30s_v1_0.nc',
             'doi': 'https://doi.org/10.17043/oden-artofmelt-2023-present-weather-1'}
    cat['horizontal_visibility'].attrs.update(hv_atts)
    
    ## Add in a ST algorithm liquid flag 
    #print('Generating a ST liquid flag...')
    #radar,lidar,radar_clutter,lidar_clearair,attenuated,atten_corrected=decode_quality(cat['quality_bits'].to_numpy())
    #ST_liq=np.zeros((len(cat.time),len(cat.height)))
    
    # When liquid is detected by the lidar, set liquid flag to one. 
    #lidar_liq_condition= (((lidar==1) & liq==1))
    #radar_only_condition = (((lidar==0) & (radar==1)))
    #radar_liq_condition= (((lidar==0) & (radar==1)) & (liq==1))
    #ST_liq[lidar_liq_condition]=1
    
    # If wbt > t_thres, force liquid
    #ST_liq[(radar_only_condition & ((cat['Tw'].data>th_t_f_dr)))]=1
    
    # For below freezing: 
    #below_freezing=((cat['Tw'].data<th_t_f_dr))
    
    # liquid
    # spec with is greater than than threshold
    # reflectivity is less than that of ice cloud
    # velocity is small
    #wh= ((cat['width'].data > th_sw_liq) & ((cat['Z'].data < th_db_liq) & (-cat['v'] < th_ve_liq)))
    #ST_liq[below_freezing & wh]=1
    
    # mixed phase
    # spec width greater than threshold, 
    # reflectivity greater than liquid cloud only threshold
    # or velocity greater than liquid cloud only threshold
    # Note shupes sign convention for velocity is different (positive downards)
    #whm = ((cat['width'].data > th_sw_liq) & (((cat['Z'].data>th_db_liq) | (-cat['v'].data>th_ve_liq))))
    #ST_liq[below_freezing & whm]=1

    # Turn off ST liquid flag if there's no LWP detected. 
    #ST_liq[cat['lwp'].to_numpy() < th_lwp_low]=0
    #ST_liq[np.transpose(np.tile((cat['lwp'].to_numpy() < 0.005), (np.shape(cat['Z'])[1],1)))]=0
    
    # Adjust droplet flag:
    # If there's no droplets in the column and lwp>20 gm-2 and ST liquid flag is set in the column
    # Change droplet where ST liquid flag is set to true
    # Using xx pixels as a minimum number for 'droplets in the column' for now - might want to check sensitivity for this and see what's reasonable. 
    #use_st_liq=(((np.sum(liq,axis=1)<3) & (np.sum(ST_liq,axis=1)>3) ) & (cat['lwp'].to_numpy()>th_lwp))
    #liq[use_st_liq,:]=ST_liq[use_st_liq,:]
    #using_st_liq_flag = np.zeros(len(use_st_liq))
    #using_st_liq_flag[use_st_liq]=1
    
    # # Add in a ST algorithm frozen precipitation flag ## TO DO
    # Could also think about including present weather data in this, or MRR data
    #ST_frozen_precip=np.zeros((len(cat.time),len(cat.height)))
    # For temperautre less than zero
    # Use reflectivity thresold for snow
    #ST_frozen_precip[(((cat['Tw'].data<th_t_f_dr) & (cat['Z'].data>th_db_snow)))] = 1
    
    # Recode categorisation 
    new_cat_bits = recode_cat(liq,falling,cold,melting,aerosol,insects)

    # To do: Add in a SLDR threshold?
    
    # write new cat bits
    # get original attrs
    atts = cat['category_bits'].attrs
    cat['category_bits'] = (('time','height'),new_cat_bits)
    
    # add attributes
    # add a note
    atts['comment2']='Removed false insect detection during fog events.'
    #atts['comment3']='Using Shupe-Turner algorithm liquid identification when LWP>%.1f g/m2 but otherwise no liquid is identified in the column. See variables ST_drop_flag for where this has been implemented.'%(th_lwp*1000.)
    cat['category_bits'].attrs.update(atts)
    
    # Add in some kind of warning / uncertainty flag for when classifying liquid is difficult. ## To Do.....
    
    # Add new variables into the categorisation file
    #cat = cat.assign(ST_drop_flag=(('time'),
    #                         using_st_liq_flag,
    #                         {'long_name':'Flag indicating where the droplet bit is derived from the Shupe-Turner algorithm',
    #                          'description':'1: droplet bit derived from Shupe-Turner algorithm, 0: droplet bit dervied using cloudnet defaults (including voodoo if voodoo is used)',
    #                          'comment':'Using Shupe-Turner algorithm for liquid identification when LWP>%.1f g/m2 but otherwise no liquid is identified in the column. The thresholds used for the ST algorithm are: spectral width %.1f m/s, reflectivity %.1f dBz, and vertical velocity %.1f m/s'%(th_lwp*1000.,th_sw_liq,th_db_liq,th_ve_liq),
    #                          'units':1,}))
    #cat = cat.assign(ST_drop=(('time','height'),
    #                    ST_liq,
    #                    {'long_name':'Pixels idenfied as containing liquid or mixed phase particles using the Shupe-Turner algorithm.',
    #                     'description':'1: pixel identified as containing liquid or mixed-phase particles, 0: pixel does not contain liquid or mixed-phase particles',
    #                     'comment':'The thresholds used for the ST algorithm to identify liquid or mixed-phase pixels from radar data are: spectral width %.1f m/s, reflectivity %.1f dBz, and vertical velocity %.1f m/s, see doi:10.1029/2007GL031008 for further details. Note that this flag is automatically set to 1 when lidar data is available and liquid is detected from the lidar using the cloudnet algorithm'%(th_sw_liq,th_db_liq,th_ve_liq),
    #                          'units':1,
    #                    }))
    #cat = cat.assign(ST_frozen_precip=(('time','height'),
    #                             ST_frozen_precip,
    #                             {'long_name':'Pixels idenfied as containing snow or falling frozen precipitation using the Shupe-Turner algorithm.',
    #                              'description':'1: pixel identified as containing frozen precipitation, 0: pixel does not contain frozen precipitation',
    #                              'comment':'The reflectivity threshold used by the ST algorithm to seperate cloud ice from frozen precipitation is >%.1f dBZ. See doi:10.1029/2007GL031008 for further details.'%th_db_snow,
    #                              'units':1,}))
    
    
    # Save and close the updated cat file temporarily before generating the classification
    cat.to_netcdf(cat_fil)
    cat.close()
    
    
    # Generate classification
    uuid = generate_classification(cat_fil, output_dir+'classification/%s_classification.nc'%date)
    
    # Add in metadata ect. 
    classi=xr.open_dataset(output_dir+'classification/%s_classification.nc'%date)
    cat=xr.open_dataset(cat_fil)
    
    for cnp_fil in [cat,classi]:
        # add decimal doy
        cnp_fil = cnp_fil.assign(day_of_year = (['time'],[decimaldayofyear(t) for t in cnp_fil.time.data], {'units':"1",'long_name':"Day of Year",'description':"time as decimal day of year"}))
        #cnp_fil['time'].attrs['_FillValue']=False
        
        # delete obsolete data
        del cnp_fil.attrs['source_file_uuids']
        cnp_fil.attrs['location'] = 'ARTofMELT'
        cnp_fil.attrs['creator_name'] =	'Heather Guy'
        cnp_fil.attrs['creator_email'] =	'heather.guy@ncas.ac.uk'
        cnp_fil.attrs['creator_url'] =	'https://orcid.org/0000-0003-3525-0766'
        cnp_fil.attrs['platform'] =	'Swedish Icebreaker Oden'
        cnp_fil.attrs['platform_type'] =	'moving_platform'
        cnp_fil.attrs['deployment_mode'] =	'ship'
        cnp_fil.attrs['featureType'] =	'timeSeries'
        #cnp_fil.attrs['time_coverage_start'] = dt.datetime.strftime(dt.datetime.strptime(date,'%y%m%d') + dt.timedelta(hours=float(cnp_fil.time.min().data)),'%d %b %Y, %H:%M UTC')
        #cnp_fil.attrs['time_coverage_end'] = dt.datetime.strftime(dt.datetime.strptime(date,'%y%m%d') + dt.timedelta(hours=float(cnp_fil.time.max().data)),'%d %b %Y, %H:%M UTC')
        cnp_fil.attrs['geospatial_bounds'] = "%sN, %sE, %sN, %sE"%(str(cnp_fil.latitude.max().data),str(cnp_fil.longitude.min().data),str(cnp_fil.latitude.min().data),str(cnp_fil.longitude.max().data))
        cnp_fil.attrs['platform_altitude'] = "Oden 4th deck, ~20 m a.s.l"
        cnp_fil.attrs['location_keywords'] = "Oden, Arctic Ocean, Fram Strait, atmosphere, profile, on the ship"
        cnp_fil.attrs['date_created'] = dt.datetime.strftime(dt.datetime.now(),'%d %b %Y %H:%M')
        cnp_fil.attrs['institution']="Stockholm University and the University of Leeds"
        cnp_fil.attrs['sampling_interval'] = "30 s"
        cnp_fil.attrs['product_version'] = "v01"
        cnp_fil.attrs['project'] = "ARTofMELT"
        cnp_fil.attrs['project_principal_investigator'] = "Michael Tjernström"
        cnp_fil.attrs['project_principal_investigator_email'] = "michaelt@misu.su.se"
        cnp_fil.attrs['project_principal_investigator_url'] = "https://orcid.org/0000-0002-6908-7410"
        cnp_fil.attrs['input_file_references'] = 'Radar: Heather Guy, Ian Brooks, Sonja Murto, Michael Tjernström, Michail Karalis, John Prytherch (2024) Cloud radar data from expedition ARTofMELT, Arctic Ocean, 2023. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-cloud-radar-1,\n Lidar: Heather Guy, Sonja Murto, Michael Tjernström, Michail Karalis, John Prytherch, Ian Brooks (2024) Cloud base heights and atmospheric backscatter observations from expedition ARTofMELT, Arctic Ocean, 2023 -- cloudnetpy input. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/lzv64l-1\n, LWP: Michael Tjernström, Sonja Murto, Michail Karalis, John Prytherch (2024) Integrated column water vapor and liquid water paths from expedition ARTofMELT, Arctic Ocean, 2023 — raw data. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-microwave-column-water-raw-1\n, thermodynamic profile: Michael Tjernström, Sonja Murto, Michail Karalis, John Prytherch (2024) Temperature and humidity profiles from microwave radiometer during expedition ARTofMELT, Arctic Ocean, 2023. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-microwave-profiles-1'
        cnp_fil.attrs['input_file_names'] = 'Radar: %s,\n Lidar: %s\n, LWP: %s\n, thermodynamic profile: %s'%(r_fil,l_fil,mwr_fil,model_fil)
        if voodoo_dir:
            cnp_fil.attrs['voodoonet_version']='%s'%(version.__version__)
            #cnp_fil.attrs['voodoonet_model']='RPG_merged_cloudy_trained-model.pt'
            cnp_fil.attrs['voodoonet_model']='RPG_merged_cloudy_resampledv_trained-model.pt'
        if cnp_fil.cloudnet_file_type == "categorize":
            print('Saving categorize')
            # Mask out radar artifacts from voodoo net results
            if 'liquid_prob' in cnp_fil.variables:
                print('correcting artifacts in voodoo data...')
                artifacts = pd.read_csv(artifacts_dir,parse_dates=[0,1],date_format='%Y-%m-%d %H:%M:%S')
                var_atts = cnp_fil.variables['liquid_prob'].attrs
                new_prob= cnp_fil.variables['liquid_prob'].to_numpy()
                for j in range(0,len(artifacts)):
                    d1=artifacts.iloc[j]['date1']
                    d2=artifacts.iloc[j]['date2']
                    h1=artifacts.iloc[j]['height1']
                    h2=artifacts.iloc[j]['height2']
                    bad_inds=np.where(cnp_fil.variables['liquid_prob'].where((((cnp_fil['time']>=d1) & (cnp_fil['time']<=d2)) & ((cnp_fil['height']>=h1)&(cnp_fil['height']<=h2))),other=-999).to_numpy()!=-999)
                    new_prob[bad_inds]=0
                    
                cnp_fil = cnp_fil.drop_vars('liquid_prob')
                cnp_fil = cnp_fil.assign(liquid_prob = (['time','height'],new_prob, var_atts))
            
            cnp_fil.attrs['title'] ='Cloud categorization products from ARTofMELT 2023'
            cnp_fil['model_height'].attrs['comment']='Note that while model_time and model_height are used in the cloudnet files, the data corresponding to these dimensions is observational data (see thermodynamic profile in input_file_references for the data source)'
            cnp_fil['model_time'].attrs['comment']='Note that while model_time and model_height are used in the cloudnet files, the data corresponding to these dimensions is observational data (see thermodynamic profile in input_file_references for the data source)'
            
            # Add in details about changes made to default cloudnet algorithm ## TO DO 
            cnp_fil.attrs['comment1']='Temperature threshold for insects has been increased relative to the default algorithm to 274 K'
            cnp_fil.attrs['comment2']='A horizontal visibility field has been added based on 10 minute average data from the present weather sensor (doi: 10.17043/oden-artofmelt-2023-present-weather-1), insect classifications have been removed where visibility is < 1000 m per Achtert et al., 2020 (DOI: 10.5194/acp-20-14983-2020)'
            #cnp.fil.attrs['comment3']='Where LWP>0 but no liquid is detected by the voodoo algorithm, the droplet bit is set for liquid pixels detected by the ST algorithm (ST_liq_flag).'
            
            cnp_fil.to_netcdf(output_dir+'categorize/cloudnet_categorize_ARTofMELT_%s_V01.nc'%date)
            cnp_fil.close()
        
        elif cnp_fil.cloudnet_file_type == "classification":
            cnp_fil.attrs['title'] ='Cloud classification products from ARTofMELT 2023'
            
            print('saving classification')
            cnp_fil.to_netcdf(output_dir+'classification/cloudnet_classification_ARTofMELT_%s_V01.nc'%date)
            cnp_fil.close()
        else:
            print('here')



if __name__ == '__main__':   
    main() 
