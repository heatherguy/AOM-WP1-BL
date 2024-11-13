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

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    Input arguments:
    1) date: in format str YYMMDD
    2) radar input data directory
    3) lidar input data directory
    4) lwp input data directory
    5) atmospheric profile input data directory
    6) cloudnetpy output directory
    7) Use voodoo? optional LV0 directory
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
    output_dir = args_in[6]
    if len(args_in)>7:
        voodoo_dir = args_in[7]
    else:
        voodoo_dir=False
    return date, radar_dir,lidar_dir,lwp_dir,model_dir,output_dir,voodoo_dir

def correct_voodoo_artifacts(cat,outfil):
    # Input is a cloudnet cateogorize file opened using xarray
    # and a string to save the post-processed file. 
    # This function looks for times where radar artifacts have been
    # manually filtered but not from the voodoo algorithm, resulting
    # in spurious liquid water detection from radar artifacts. 
    # We remove the spurious liquid and resave the categorise file. 
    
    # liquid identified from voodoo when liquid_prob > 0.55
    # Falling hydrometeors are radar signals that are
    # a) not insects b) not clutter. Furthermore, falling hydrometeors
    # are strong lidar pixels excluding liquid layers (thus these pixels
    # are ice or rain). They are also weak radar signals in very cold
    # temperatures.

    # So, if liquid_prob > 0.55, and falling=0, We know that there is a radar signal
    # from voodoo, but no good radar signal from the reflectivity. 
    # Therefore these are the artificats. 
    # in this case if liq=1, we need to set liq=0

    # decode
    liq,falling,cold,melting,aerosol,insects=decode_cat(cat['category_bits'].to_numpy())

    # find bad liquid signals from radar artifacts
    condition = (liq==1) & (falling==0) & (cat['liquid_prob'].to_numpy()>0.55)

    # modify liq
    mod_liq = liq.copy()
    mod_liq[condition]=False

    # recode
    new_cat_bits = recode_cat(mod_liq,falling,cold,melting,aerosol,insects)

    # get original attrs
    atts = cat['category_bits'].attrs
    # add a note
    atts['comment1']='Cleaned for bad voodoo results due to radar artifacts. If liquid_prob>0.55 and falling=0, any liq flags have been removed.'

    # write new cat bits
    cat['category_bits'] = (('time','height'),new_cat_bits)

    # add attributes
    cat['category_bits'].attrs.update(atts)

    # save...
    cat.to_netcdf(outfil)
    return

def replace_nan(nonan,mask):
    #if non nan is a flat array with nans removed
    # and mask is the orginal masked array where nans are
    # add the nans back into nonan and return. 
    c =  np.empty_like(mask)
    c.fill(np.nan)
    c[~mask] = nonan
    return c

def decode_cat(arr):
    # Bit 0: Small liquid droplets are present.
    # Bit 1: Falling hydrometeors are present; if Bit 2 is set then these are most
    # likely ice particles, otherwise they are drizzle or rain drops.
    # Bit 2: Wet-bulb temperature is less than 0 degrees C, implying the phase of
    # Bit-1 particles.
    # Bit 3: Melting ice particles are present.
    # Bit 4: Aerosol particles are present and visible to the lidar.
    # Bit 5: Insects are present and visible to the radar.
    arr_flat = arr.flatten()
    nan_mask = np.isnan(arr_flat)
    arr_nonan = arr_flat[~nan_mask]
    arr_decode=[bin(int(c))[2:].zfill(6) for c in arr_nonan]
    liq = np.asarray([int(decode[-1]) for decode in arr_decode ])
    falling = np.asarray([int(decode[-2]) for decode in arr_decode ])
    cold = np.asarray([int(decode[-3]) for decode in arr_decode ])
    melting = np.asarray([int(decode[-4]) for decode in arr_decode ])
    aerosol = np.asarray([int(decode[-5]) for decode in arr_decode ])
    insects = np.asarray([int(decode[-6]) for decode in arr_decode ])

    # reinsert nans and reshape
    liq = replace_nan(liq,nan_mask).reshape(np.shape(arr))
    falling = replace_nan(falling,nan_mask).reshape(np.shape(arr))
    cold = replace_nan(cold,nan_mask).reshape(np.shape(arr))
    melting = replace_nan(melting,nan_mask).reshape(np.shape(arr))
    aerosol = replace_nan(aerosol,nan_mask).reshape(np.shape(arr))
    insects = replace_nan(insects,nan_mask).reshape(np.shape(arr))

    return liq,falling,cold,melting,aerosol,insects

def recode_cat(liq,falling,cold,melting,aerosol,insects):
    # flatten the input arrays
    l=liq.flatten()
    f=falling.flatten()
    c=cold.flatten()
    m=melting.flatten()
    a=aerosol.flatten()
    ins=insects.flatten()

    # Make the binary strings
    bi_str=['%i%i%i%i%i%i'%(ins[i],a[i],m[i],c[i],f[i],l[i]) for i in range(0,len(l)) ]
    # cal the integers
    ints = [int(s,2) for s in bi_str]
    # reshape
    cat_bits=np.asarray(ints).reshape(np.shape(liq))
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
    arr_flat = arr.flatten()
    nan_mask = np.isnan(arr_flat)
    arr_nonan = arr_flat[~nan_mask]
    arr_decode=[bin(int(c))[2:].zfill(6) for c in arr_nonan]
    radar = np.asarray([int(decode[-1]) for decode in arr_decode ])
    lidar = np.asarray([int(decode[-2]) for decode in arr_decode ])
    radar_clutter = np.asarray([int(decode[-3]) for decode in arr_decode ])
    lidar_clearair = np.asarray([int(decode[-4]) for decode in arr_decode ])
    attenuated = np.asarray([int(decode[-5]) for decode in arr_decode ])
    atten_corrected = np.asarray([int(decode[-6]) for decode in arr_decode ])

    # reinsert nans and reshape
    radar = replace_nan(radar,nan_mask).reshape(np.shape(arr))
    lidar = replace_nan(lidar,nan_mask).reshape(np.shape(arr))
    radar_clutter = replace_nan(radar_clutter,nan_mask).reshape(np.shape(arr))
    lidar_clearair = replace_nan(lidar_clearair,nan_mask).reshape(np.shape(arr))
    attenuated = replace_nan(attenuated,nan_mask).reshape(np.shape(arr))
    atten_corrected = replace_nan(atten_corrected,nan_mask).reshape(np.shape(arr))

    return radar,lidar,radar_clutter,lidar_clearair,attenuated,atten_corrected

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
    date, radar_dir,lidar_dir,lwp_dir,model_dir,output_dir,voodoo_dir = get_args(sys.argv)
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
    
    # Correct voodoo radar artifacts
    if voodoo_dir: 
        print('Correcting voodoo artifacts...')
        cat = xr.open_dataset(cat_fil+'_tmp')
        cat = correct_voodoo_artifacts(cat,cat_fil)
    else:
        cat = xr.open_dataset(cat_fil+'_tmp')

        # Sort variables and metdata ect
        # Get the cloudnetpy version 
        # Add heathers interation version and comments. 
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
        
        # delte obsolete data
        del cnp_fil.attrs['source_file_uuids']
        del cnp_fil.attrs['location']
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
        cnp_fil.attrs['additional_creator_name'] = "John Prytherch"
        cnp_fil.attrs['additional_creator_email'] = "john.prytherch@misu.su.se"
        cnp_fil.attrs['additional_creator_url'] = "https://orcid.org/0000-0003-1209-289X"
        cnp_fil.attrs['additional_creator_name'] = "Sonja Murto"
        cnp_fil.attrs['additional_creator_email'] = "sonja.murto@misu.su.se "
        cnp_fil.attrs['additional_creator_url'] = "https://orcid.org/0000-0002-4966-9077"
        cnp_fil.attrs['input_file_references'] = 'Radar: Heather Guy, Ian Brooks, Sonja Murto, Michael Tjernström, Michail Karalis, John Prytherch (2024) Cloud radar data from expedition ARTofMELT, Arctic Ocean, 2023. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-cloud-radar-1,\n Lidar: Heather Guy, Sonja Murto, Michael Tjernström, Michail Karalis, John Prytherch, Ian Brooks (2024) Cloud base heights and atmospheric backscatter observations from expedition ARTofMELT, Arctic Ocean, 2023 -- cloudnetpy input. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/lzv64l-1\n, LWP: ......\n, thermodynamic profile: .....'
        cnp_fil.attrs['input_file_names'] = 'Radar: %s,\n Lidar: %s\n, LWP: %s\n, thermodynamic profile: %s'%(r_fil,l_fil,mwr_fil,model_fil)
        if voodoo_dir:
            cnp_fil.attrs['voodoonet_version']='%s'%(version.__version__)
        if cnp_fil.cloudnet_file_type == "categorize":
            print('Saving categorize')
            cnp_fil.attrs['title'] ='Cloud categorization products from ARTofMELT 2023'
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
