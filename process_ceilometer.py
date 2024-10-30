#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import datetime as dt
import xarray as xr
import sys
import glob
from scipy.ndimage import gaussian_filter
from numpy import ma
from numpy.lib.stride_tricks import sliding_window_view

####### INPUTS #######
# Data location:
#in_loc = path
#save = path

#start='202305050000'
#stop='202306140000'

#Example usage: 
# python process_ceilometer.py $in_loc $start $stop $save
#############################################################

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    try:
        dpath = args_in[1]
        save = args_in[4]
        start = dt.datetime.strptime(args_in[2],'%Y%m%d%H%M')
        stop = dt.datetime.strptime(args_in[3],'%Y%m%d%H%M')
    except:
        print('Input error')
        sys.exit()
    # return values:
    return dpath,start,stop,save

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

def mdiff(array: np.ndarray) -> float:
    """Returns median difference of 1-D array."""
    return float(ma.median(ma.diff(array)))
    
def calc_sigma_units(
    time_vector: np.ndarray, range_los: np.ndarray
) -> tuple[float, float]:
    """Calculates Gaussian peak std parameters.

    The amount of smoothing is hard coded. This function calculates
    how many steps in time and height corresponds to this smoothing.

    Args:
        time_vector: 1D vector (fraction hour).
        range_los: 1D vector (m).

    Returns:
        tuple: Two element tuple containing number of steps in time and height to
            achieve wanted smoothing.

    """
    #if len(time_vector) == 0 or np.max(time_vector) > 24:
    #    raise ValueError("Invalid time vector")
    minutes_in_hour = 60
    sigma_minutes = 2
    sigma_metres = 5
    time_step = mdiff(time_vector) * minutes_in_hour
    alt_step = mdiff(range_los)
    x_std = sigma_minutes / time_step
    y_std = sigma_metres / alt_step
    return x_std, y_std

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def convert_times_hours(times):
    frac_hours=np.zeros(len(times))*np.nan
    # Calculate fractional hours since 1st may 2023
    for i in range(0,len(times)):
        tss = pd.to_datetime(times[i])
        fractional_hours = (tss - pd.Timestamp(dt.datetime(2023,5,1))) / pd.Timedelta(hours=1)
        frac_hours[i]=fractional_hours

    return frac_hours

def ceilometer_noise_filter(beta_raw,heights,time,snr_limit,gates,noise_min,range_corrected=1):
    '''
    # Noise filters raw ceilometer data based on a snr limit
    # Following methodology described in Kotthaus et al., 2016, doi:10.5194/amt-9-3769-2016
    #
    # Inputs : 
    # beta_raw: raw backscatter in sr-1 m-1 (numpy 2D array)
    # heights: range to center of bin in m (numpy 1D array)
    # time: time should be in fractional hours  (numpy 1D array)
    # range_corrected: 1 or 0, if 1, input data is already range corrected, zero if not. 
    # snr_limit: snr limit to filter noise
    # gates: index of the top gates to use for the noise filter, int or tuple, if int, use the top n gates. 
    #        if it'ls a tuple, use gate a to gate b. 
    # noise_min: minimum expected noise level for a particular ceilometer (use 1e-9 for CT25K)
    #
    # Outputs: 
    # beta_raw_filtered: Raw backscatter noise filtered
    # beta_smooth_filtered: Smooth backscatter noise filtered
    # noise_floor: 1D array noise floor for each time step
    # SNR: 2D array signal to noise ration 
    '''

    # Get the non-ranged corrected data to calculate the noise floor
    if range_corrected==1:
        non_range_corrected = beta_raw.copy() / (heights**2)
    else: 
        non_range_corrected = beta_raw.copy()
        beta_raw = beta_raw * (heights**2)

    # Calculated the smoothed backscatter using the same gaussian filter used in cloudnet
    sigma = calc_sigma_units(time, heights)
    beta_smooth_nrc = gaussian_filter(non_range_corrected, sigma)
    beta_smooth = beta_smooth_nrc * (heights**2)

    # Calculate the noise floor:  the mean β plus standard deviation σβ of the attenuated backscatter β 
    #    (i.e. before range correction) across a certain number of gates from the top of the profile. 
    #    Suggest highest 300 m. In cloudnetpy, they use the top 0.1 gates. 
    #    note that for the summit, the top 5 gates are zero anyway, so use to 7215 to 7515 (gates -16 to -5)
    if type(gates) ==int: 
        beta_noise_mean = np.mean(non_range_corrected[:,-gates:],axis=1)
        beta_noise_sd = np.std(non_range_corrected[:,-gates:],axis=1)
    elif type(gates)==tuple:
        beta_noise_mean = np.mean(non_range_corrected[:,gates[0]:gates[1]],axis=1)
        beta_noise_sd = np.std(non_range_corrected[:,gates[0]:gates[1]],axis=1)
    else: 
        print('Error: gates either needs to be an integer (top x gates), or a tuple (gate a to gate b)')
        return

    noise_floor = beta_noise_mean + beta_noise_sd
    noise_floor1=noise_floor.copy()

    # Check the profiles used to calculate the noise floor don't contain any cirruse clouds
    #      Calculate the relative variance (RV) or coefficient of variation)
    #      The relative variance is the ration of the standard deviation to the mean, with statistics applied
    #      over moving winds along range and time (here wr = wt = 3) - that's +/-, so the window size is 6. 
    #      if RV is sufficiently small, the backscatter is interpreted as a cloud and is not used to estimate the noise floor
    #      RV should exceen 1, indicating that the variability exceeds the mean signal 

    relative_mean = np.mean(sliding_window_view(non_range_corrected,[6,6]),axis=(2,3))
    relative_std = np.std(sliding_window_view(non_range_corrected,[6,6]),axis=(2,3))
    relative_variance = (relative_std / relative_mean)**2
    # pad the sides of the relative variance by duplicating the nearest value
    relative_variance = np.insert(relative_variance,0,relative_variance[0,:],axis=0)
    relative_variance = np.insert(relative_variance,0,relative_variance[0,:],axis=0)
    relative_variance = np.insert(relative_variance,0,relative_variance[0,:],axis=0)
    relative_variance = np.insert(relative_variance,-1,relative_variance[-1,:],axis=0)
    relative_variance = np.insert(relative_variance,-1,relative_variance[-1,:],axis=0)

    # nan noise where RV exceeds one
    if type(gates) ==int: 
        noise_floor[np.min(relative_variance[:,-gates:],axis=1)<1]=np.nan
    elif type(gates)==tuple:
        noise_floor[np.min(relative_variance[:,gates[0]:gates[1]],axis=1)<1]=np.nan
        
    # replace negatives & set minimum 
    noise_floor[noise_floor<noise_min]=noise_min

    # Interpolate across any gaps in the noise floor. 
    nans, x= nan_helper(noise_floor)
    noise_floor[nans]= np.interp(x(nans), x(~nans), noise_floor[~nans])

    # calculate the SNR = smoothed non-range corrected beta / noise_floor

    SNR = beta_smooth_nrc / np.reshape(np.repeat(noise_floor,np.shape(non_range_corrected)[1]),np.shape(non_range_corrected))

    #  Mask any points where SNR < optimal SNR
    beta_raw_filtered =  beta_raw.copy() 
    beta_raw_filtered[((SNR<snr_limit)|(np.isnan(SNR)))]=np.nan

    beta_smooth_filtered = beta_smooth.copy()
    beta_smooth_filtered[((SNR<snr_limit)|(np.isnan(SNR)))]=np.nan

    return beta_raw_filtered, beta_smooth_filtered,noise_floor,noise_floor1, SNR, beta_smooth_nrc



def main():
    """
    Applies a noise filtering procedure to the raw ceilometer backscatter
    and generates files in the correct format for input into cloudnetpy.
    Adds additional metadata and formatting
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Path to input ceilometer data
        save:  Output directory (optional)
        
    Returns:
        None
    
    """
    
    # check / get args:
    dpath,start,stop,save = get_args(sys.argv)
    all_dates = pd.date_range(start,stop)

    cl_all = xr.open_dataset(dpath,decode_times=True)

    # Update latitude and longitude using ship data
    # Get ship data
    ship_path='/Users/heather/Desktop/ARTofMELT/data/shipdata/'
    fils = glob.glob(ship_path+'*.csv')
    fils.sort()
    ship_data = pd.concat([pd.read_csv(fil,parse_dates=[0],index_col=0) for fil in fils])

    for i in range(0,len(all_dates)):
        day = all_dates[i].to_pydatetime()
        print(day)
        cl_in = cl_all.sel(time=slice(day,day+dt.timedelta(hours=24,minutes=60,seconds=59))).copy()
        
        if cl_in.time.size == 0:
            print('No data found, skipping ')
            continue
            
        # Use bs_prof, and convert units from sr-1 km-1 to sr-1 m-1 (/1000)
        #beta_raw (time, range)
        #units :
        #sr-1 m-1
        #long_name :
        #Attenuated backscatter coefficient
        #comment :
        #Non-screened attenuated backscatter coefficient.
        cl_in['beta_raw']=cl_in['backscatter_profile']/1000.
        cl_in['beta_raw'].attrs['long_name'] = 'Attenuated backscatter coefficient'
        cl_in['beta_raw'].attrs['units'] = 'sr-1 m-1'
        cl_in['beta_raw'].attrs['comment'] = 'Non-screened attenuated backscatter coefficient.'
        
        # Add wavelength = 910.
        #wavelength
        #()
        #float32
        #units :
        #nm
        #long_name :
        #laser wavelength
        cl_in['wavelength']=910.0
        cl_in['wavelength'].attrs['units'] = 'nm'
        cl_in['wavelength'].attrs['long_name'] = 'laser wavelength'

        # rename ceil_range as height (accounts for height of ceilometer above sea level)
        # height
        #(range)
        #float32
        #units :
        #m
        #long_name :
        #Height above mean sea level
        #standard_name :
        #height_above_mean_sea_level
        cl_in['height']=cl_in['ceilometer_range']
        cl_in['height'].attrs['units'] = 'm'
        cl_in['height'].attrs['long_name'] = 'Height above mean sea level'
        cl_in['height'].attrs['standard_name'] = 'height_above_mean_sea_level'

        # Set equal to 1
        # calibration_factor (one value, no coords)
        #float32
        #units :
        #1
        #long_name :
        #Attenuated backscatter calibration factor
        #comment :
        #Calibration factor applied.
        cl_in['calibration_factor']=1.
        cl_in['calibration_factor'].attrs['units'] = '1'
        cl_in['calibration_factor'].attrs['long_name'] = 'Attenuated backscatter calibration factor'
        cl_in['calibration_factor'].attrs['comment'] = 'Calibration factor applied (Note that if calibration factor=1, this is equivalent to no calibration applied).'
        
        # Set equial to 1
        #zenith_angle (one value, no coords)
        #float32
        #units :
        #degree
        #long_name :
        #Zenith angle
        #standard_name :
        #zenith_angle
        #comment :
        #Angle to the local vertical. A value of zero is directly overhead.
        # Looking at the raw data seems to vary between 0-4 
        cl_in['zenith_angle']=1.
        cl_in['zenith_angle'].attrs['units'] = 'degree'
        cl_in['zenith_angle'].attrs['long_name'] = 'Zenith angle'
        cl_in['zenith_angle'].attrs['standard_name'] = 'zenith_angle'
        cl_in['zenith_angle'].attrs['comment'] = 'This is set at a single value for the purpose of the cloudnetpy algorithm and corresponds to how the ceilometer was configured assuming a level platform. In reality the zenith angle varied between 0-4 degrees when the ship was in motion. This will have essentially no impact on the data for the majority of uses.' 

        cl_in.attrs['software_version'] = 202

        # modify range
        cl_in['range']=cl_in['height'] - 25.
        cl_in['range'].attrs['units'] = 'm'
        cl_in['range'].attrs['long_name'] = 'Range from instrument'
        cl_in['range'].attrs['standard_name'] = 'Distance from instrument to centre of each range bin.'
        cl_in['range'].attrs['axis'] ='Z'
        cl_in = cl_in.swap_dims({'range_levels':'range'})

        # do the noise filtering 
        # noise filtering options: 
        
        # The snr limit is calculated by looking at the welch t-test acceptance plot (like fig 9 in kotthause et al., 2016)
        # when more than 90% of values are significantly (p<0.01) from the noise. This is about 10. 
        snr_limit=10
        #  for the summit data I'm looking at, the top 5 gates are zero anyway, so use to 7215 to 7515 (gates -16 to -5)
        # for aom top 300 m is the top 31 gates
        gates = 31 
        # To get this value, I'm using the median value of the noise floor for the AOM campaign
        noise_min = 4.509295348092409e-14 
        range_corrected=1
        frac_hours = convert_times_hours(cl_in['time'].to_numpy())
        
        # Do the noise filtering
        beta_raw = cl_in['beta_raw'].to_numpy()
        beta_raw_filtered, beta_smooth_filtered,noise_floor,noise_floor1, SNR, beta_smooth_nrc = ceilometer_noise_filter(beta_raw,cl_in['range'].to_numpy(),frac_hours,snr_limit,gates,noise_min,range_corrected)
        
        # add to dataframe
        cl_in['beta_smooth']=xr.DataArray(beta_smooth_filtered, dims=['time', 'range'],
                                        coords={'time':cl_in.time,'range':(cl_in['height'] - 25.).data},
                                   attrs={'long_name': 'Attenuated backscatter coefficient',
                                          'units' :"sr-1 m-1",
                                          'comment': "SNR-screened attenuated backscatter coefficient.\nWeak background smoothed using Gaussian 2D-kernel. SNR threshold applied: %s"%snr_limit})# ,
                                          #  '_FillValue': 9.96921E36,
                                          #  '_ChunkSizes': [1440, 385]})
        
        
        cl_in['beta_raw']=xr.DataArray(beta_raw, dims=['time', 'range'], 
            attrs={'long_name': 'Attenuated backscatter coefficient',
            'units' :"sr-1 m-1",
            'comment': "Non-screened attenuated backscatter coefficient." ,
              '_FillValue': 9.96921E36,
              '_ChunkSizes': [1440, 385]})
        
        cl_in['beta']=xr.DataArray(beta_raw_filtered, dims=['time', 'range'], 
            attrs={'long_name': 'Attenuated backscatter coefficient',
            'units' :"sr-1 m-1",
            'comment': "SNR-screened attenuated backscatter coefficient. SNR threshold applied: %s"%snr_limit ,
              '_FillValue': 9.96921E36,
              '_ChunkSizes': [1440, 385]})
        
        # Update latitude and longitude using ship data
        # Get ship data
        lats = ship_data[' Oden.Ship.LatitudeDegreesFixed'].reindex(cl_in.time,method='nearest',tolerance='5s')
        lons = ship_data[' Oden.Ship.LongitudeDegreesFixed'].reindex(cl_in.time,method='nearest',tolerance='5s')
        cl_in = cl_in.assign(latitude = (['time'],lats.to_numpy(), {'units':"degree_north",'long_name':"Latitude of site",'standard_name':"latitude"}))
        cl_in = cl_in.assign(longitude = (['time'],lons.to_numpy(), {'units':"degree_east",'long_name':"Longitude of site",'standard_name':"longitude"}))

        # Sort out the rest of attributes and variables. 
        # sort times
        frac_hours=np.zeros(len(cl_in.time))*np.nan
        for i in range(0,len(cl_in.time)):
            tss = pd.to_datetime(cl_in['time'][i].to_numpy())
            fractional_hours = (tss - day) / pd.Timedelta(hours=1)                     
            frac_hours[i]=fractional_hours

        cl_in['time']=frac_hours
        cl_in['time'].attrs['units'] = "hours since %s 00:00:00 +00:00"%dt.datetime.strftime(day,'%Y-%m-%d')
        cl_in['time'].attrs['long_name'] = "Time UTC"
        cl_in['time'].attrs['standard_name'] = 'time'
        cl_in['time'].attrs['calendar'] ='standard'
        tunits="hours since %s 00:00:00 +00:00"%dt.datetime.strftime(day,'%Y-%m-%d')
        cl_in.time.encoding['units'] = tunits
        cl_in.time.encoding['calendar'] = 'standard'
        cl_in.attrs['cloudnet_file_type'] = "lidar"
        cl_in.attrs['year'] = str(list(set(cl_in.year.to_numpy()))[0])
        cl_in.attrs['month'] = str(list(set(cl_in.month.to_numpy()))[0])    
        cl_in.attrs['day'] = str(list(set(cl_in.day.to_numpy()))[0])    

        # Drop obsolete variables: 
        cl_in = cl_in.drop_vars(['backscatter_profile', 'ceilometer_range','range_levels'])

        # add SNR
        cl_in = cl_in.assign(SNR = (['time','range'],SNR, {'units':"1",'long_name':"Signal to noise ratio",'description':"Signal to noise ratio calculated following the method of Kotthaus et al., 2016 (doi:10.5194/amt-9-3769-2016) (Eq. 15). Minimum noise level is the median value of the nosie floor throughout the ARTofMELT campaign: %.2e. "%noise_min}))
        cl_in['time'].attrs['_FillValue']=False
        
        # add meta data
        del cl_in.attrs['acknowledgement']
        del cl_in.attrs['comments']
        cl_in.attrs['comment'] = 'Noise filtering has been applied following the method of Kotthaus et al., 2016 (doi:10.5194/amt-9-3769-2016). Minimum noise level (the median value of the nosie floor throughout the ARTofMELT campaign): %.2e. SNR limit (calculated using welch t-tests as in fig. 9, kotthause et al., 2016): %s'%(noise_min,snr_limit)
        cl_in.attrs['creator_name'] =	'Heather Guy'
        cl_in.attrs['creator_email'] =	'heather.guy@ncas.ac.uk'
        cl_in.attrs['creator_url'] =	'https://orcid.org/0000-0003-3525-0766'
        cl_in.attrs['platform'] =	'Swedish Icebreaker Oden'
        cl_in.attrs['platform_type'] =	'moving_platform'
        cl_in.attrs['deployment_mode'] =	'ship'
        cl_in.attrs['title'] =	'CL31 ceilometer data from ARTofMELT 2023'
        cl_in.attrs['featureType'] =	'timeSeries'

        cl_in.attrs['time_coverage_start'] = dt.datetime.strftime(day + dt.timedelta(hours=float(cl_in.time.min().data)),'%d %b %Y, %H:%M UTC')
        cl_in.attrs['time_coverage_end'] = dt.datetime.strftime(day + dt.timedelta(hours=float(cl_in.time.max().data)),'%d %b %Y, %H:%M UTC')
        cl_in.attrs['geospatial_bounds'] = "%sN, %sE, %sN, %sE"%(str(cl_in.latitude.max().data),str(cl_in.longitude.min().data),str(cl_in.latitude.min().data),str(cl_in.longitude.max().data))
        cl_in.attrs['platform_altitude'] = "Oden 7th deck, ~25 m a.s.l"

        # Make sure the height is correct assuming a height amsl of 25m on the 4th deck
        cl_in = cl_in.assign(altitude = np.array(25.))
        
        cl_in.attrs['location_keywords'] = "Oden, Arctic Ocean, Fram Strait, atmosphere, profile, on the ship"
        cl_in.attrs['date_created'] = dt.datetime.strftime(dt.datetime.now(),'%d %b %Y %H:%M')
        cl_in.attrs['institution']="Stockholm University and the University of Leeds"
        cl_in.attrs['processing_software'] = "https://github.com/heatherguy/AOM-WP1-BL/blob/main/process_ceilometer.py"
        cl_in.attrs['sampling_interval'] = "30 s"
        cl_in.attrs['product_version'] = "v01"
        cl_in.attrs['project'] = "ARTofMELT"
        cl_in.attrs['project_principal_investigator'] = "Michael Tjernström"
        cl_in.attrs['project_principal_investigator_email'] = "michaelt@misu.su.se"
        cl_in.attrs['project_principal_investigator_url'] = "https://orcid.org/0000-0002-6908-7410"
        cl_in.attrs['additional_creator_name'] = "John Prytherch"
        cl_in.attrs['additional_creator_email'] = "john.prytherch@misu.su.se"
        cl_in.attrs['additional_creator_url'] = "https://orcid.org/0000-0003-1209-289X"
        cl_in.attrs['additional_creator_name'] = "Sonja Murto"
        cl_in.attrs['additional_creator_email'] = "sonja.murto@misu.su.se "
        cl_in.attrs['additional_creator_url'] = "https://orcid.org/0000-0002-4966-9077"
        cl_in.attrs['input_file_reference'] = 'Sonja Murto, Michael Tjernström, Michail Karalis, John Prytherch (2024) Cloud base heights and atmospheric backscatter observations from expedition ARTofMELT, Arctic Ocean, 2023. Dataset version 1. Bolin Centre Database. https://doi.org/10.17043/oden-artofmelt-2023-ceilometer-1'
        cl_in.attrs['history'] = "Noise filtering procedure applied and variable names modified to match the correct format for cloudnetpy input, applied to input file: %s"%dpath[-54:]

        # add decimal doy
        cl_in = cl_in.assign(day_of_year = (['time'],[decimaldayofyear(t) for t in cl_in.time.data], {'units':"1",'long_name':"Day of Year",'description':"time as decimal day of year"}))
        cl_in['time'].attrs['_FillValue']=False
     
        # Save the masked files
        print('Saving masked netcdf..')
        cl_in.to_netcdf(save+'CL31_cloudnet_ARTofMELT_%s_v01.nc'%(dt.datetime.strftime(day,'%y%m%d')))

    return



    
if __name__ == '__main__':
    main()          

