#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:06:27 2020

@author: heather

Various functions for processing/ calculating turbulent fluxes
"""

import numpy as np      
import datetime as dt
from scipy.signal import medfilt, detrend
import pandas as pd

def replace_outliers(var,sd):
    """
    Replace outliers with median filter. This function uses a window 
    size of 11 data points for the median filter and defines outliers as
    data that are greater than 3 standard deviations away from the median.
    Interpolates over Nan's.
   
    Parameters:
        var: input variable
        sd:  standard deviation
        
    Returns:
        var with median filter applied and nans interpolated over.
    
    """
    # Replace outliers with median filter. This function uses a window 
    # size of 11 data points for the median filter and defines outliers as
    # data that are greater than 3 standard deviations away from the median. 
    
    var=var.astype(float)
    jj = ~np.isnan(var)                # Ignore nans
    temp = var[jj]
    mf = medfilt(temp,11)             # Get median filter.
    ii = np.abs(temp - mf) > 3*sd     # Get outliers.
    temp[ii] = mf[ii]                 # Replace these outliers with the median
    var_clean = var
    var_clean[jj] = temp              # Put back into orginal array
    
    # Interpolate over any Nans (default method = linear)
    var_clean = var_clean.interpolate()
    
    return var_clean


def clean_metek(m_in):
    """
    Clean Metek (3D sonic) data
    Based on clean_metek.m by IAM, 25/7/2014
    Calls replace_outliers
    wind compents and sonic temperature are cleaned up:
    - the time series are filtered for outliers, these are replaced with  
       median filtered values
    - missing values from error messages (NaNs) are then interpolated over
    
    Parameters:
        m_in: metek  - data strucure
        
    Returns:
        m_out: cleaned data structure. 
    
    """
    m_out = m_in.copy()
    jj_z = ~np.isnan(m_in.z)         # Not nan indices - use to calculate standard deviation.
    jj_T = ~np.isnan(m_in['T'])
    m_sd_z = np.std(m_in.z[jj_z])    # standard deviation of vertical wind component. 
    m_sd_T = np.std(m_in['T'][jj_T]) # standard deviation of Temperature. 
    m_out['x']=replace_outliers(m_in['x'],m_sd_z)
    m_out['y']=replace_outliers(m_in['y'],m_sd_z)
    m_out['z']=replace_outliers(m_in['z'],m_sd_z)   
    m_out['T']=replace_outliers(m_in['T'],m_sd_T)

    return m_out


def Ts_sidewind_correction(Ts,u,v,w):
    """
    Do cross wind temperature correction
    Adapted from Ts_sidewind_correction.m by IMB July 2007
    Correct contamination of sonic temperature measurement for lengthening of
    sonic path by sidewind. 
    
    Correction follows van Dijk et al. 2004: The
    principles of surface flux physics: theory, practice and description of
    the ECPACK library (www.met.wau.nl/projects/jep). See also: Liu et al.
    2001: BLM, 100, 459-468, and Schotanus et al. 1983: BLM, 26, 81-93.
    
    Parameters:
        Ts    : sonic temperature (K)
        u,v,w : wind components in sonic measurement frame (before any rotations,
                motion correction, etc) (m/s)
        
    Returns:
        T     : corrected temperature (K)
    
    """
    vn2 = (3/4)*(u**2 + v**2) + 0.5*w**2
    T = Ts + vn2/403
    return T


def rotate_to_run(m,avp):
    """
    Correct tilt and align with horizontal streamline over a single run
    Adapted from rotate_to_run.m by IMB July 2006
    references:
    Wilczak et al. 2001: sonic anemometer tilt corrections. BLM, 99, 127-150
    (but beware typos in equations)
    Kaimal & Finnigan, 1994: Atmospheric Boundary Layer Flows: Their
    Structure and Measurement. Oxford University Press
    van Dijk et al. 2004: The principles of surface flux physics: theory,
    practice and description of the ECPACK library
    www.met.wau.nl/projects/jep
    
    Parameters:
        m     : Unrotated metek data structure.
        avp   : Length of a single run (minutes). i.e. averaging period. 
        
    Returns:
        m_out     : Wind components in streamline oriented reference frame
    
    """
    
    # First rotate to align x-axis with mean wind direction in sonic's
    # reference frame
    
    m_g = m.groupby(pd.Grouper(freq='%sMin'%avp))     # Split into single runs (of length avp minutes) 
    m_out = pd.DataFrame(columns=['x','y','z','T','u','v','w','theta','phi'])
    
    # Loop through each run to perform correction: 
    for group in m_g:
        m = group[1]
        
        # First rotate to align x-axis with mean wind direction in sonic's
        # reference frame
        
        theta=np.arctan2(np.mean(m['y']),np.mean(m['x']))
        u1 = m['x']*np.cos(theta) + m['y']*np.sin(theta)
        v1 = -m['x']*np.sin(theta) + m['y']*np.cos(theta)
        w1 = m['z']

        # Next rotate u and w so that x-axis lies along mean streamline and 
        # mean(w) is zero
        
        phi = np.arctan2(np.mean(w1),np.mean(u1))
        m['u'] = u1*np.cos(phi) + w1*np.sin(phi)
        m['v'] = v1
        m['w']=  -u1*np.sin(phi) + w1*np.cos(phi)
        
        # Theta is angle of rotation um-to-vm (anticlockwise or righthanded)
        # to align u with mean wind (degrees)
        
        m['theta'] = theta*180/np.pi

        # phi is tilt angle (+ve tilts x-axis upwards) to align x-axis with
        # mean streamline and force <w>=0
        
        m['phi'] = phi*180/np.pi

        m_out = pd.concat([m_out, m])
    return m_out


def eddyflux(x,y):
    """
    Function to calculate instantaneous flux 
    x is time series of vertical velocity
    y is parameter of interest.
           
    Returns:
        flux: Instantaneous flux
        std: Standard deviation of instantaneous flux
    
    """   
    xprime = detrend(x)
    yprime = detrend(y)
    flux = np.mean(xprime * yprime) # Instantaneous flux
    std = np.std(xprime * yprime)   # Standard deviation of instantaneous flux

    return flux,std


def licor_mmr(T,P,Nconc):
    """
    Calculate mass mixing ratio from LiCOR. 
    Requires P from Licor, Nconc from Licor, T from HMP155
    Converts Nconc (molar number concentration) to Mass mixing ration q (kg/kg)
    
    Parameters:
        T: Temperature form HMP155 sensor
        P: Air pressure from licor
        Nconc: Molar number concentration H2O from licor
        
    Returns:
        q:      Specific humidity (kg/kg)
        PPwet:  Partial pressure of water vapor
        PPdry:  Partial pressure of dry air
        mmr:    Mass mixing ratio
    
    """
   
    Ma = 28.96        # Molar mass of dry air (g/mol)
    Mh = 18.01528     # Molar mass of H2O (g/mol)
    R = 8.314         # Universal gas constant (J/K/mol)
    MW = 0.01802;     # Molecular weight of water, kg mol-1
    Rw = R/MW         # Gas constand for water vapor (J/K/kg)
    
    # Join and interpolate temperature (linear interpolation)
    T = T.reindex(P.index,method='nearest',tolerance='1min')
    temp_df = pd.concat([T, P, Nconc],axis=1,keys=['T','P','Nconc'])
    temp_df['T'] = temp_df['T'].interpolate(limit=200)
    
    # Calculate mass mixing ratio (q)
    mass_conc_dryair = (Ma*temp_df['P'])/(R*temp_df['T']) / 1000    # kg dry air/ m3
    mass_conc_h2o = (Mh * temp_df['Nconc']) / 1000                  # kg water / m3
    mmr = mass_conc_h2o / mass_conc_dryair                          # Mass mixing ratio   
    q = mass_conc_h2o / (mass_conc_h2o + mass_conc_dryair)
    q = q.rename('q')
    
    # Wet air partial pressure (Pa)
    PPwet = mass_conc_h2o * temp_df['T'] * Rw
    
    # Dry air partial pressure (pa)
    PPdry = temp_df['P'] - PPwet
    
    return q, PPwet, PPdry, mmr


def rho_m(T,P,q):
    """
    Estimate air density using time varying T, P and q.
    Assume ideal gas behaviour
    
    Parameters:
        P: static pressure (mb) from licor
        T: air temperature (K) from HMP1
        q: specific humidity (kg/kg) from licor.
        
    Returns:
        rho: air density (kg m-2)
    
    """

    Rd = 287 # J K-1 kg-1, gas constant for dry air.     

    temp_df = pd.concat([T, P, q],axis=1,keys=['T','P','q'])
    temp_df['T'] = temp_df['T'].interpolate(limit=200) # Interpolate temperature
    
    rho_m =  1/ ((Rd * temp_df['T']) *(1+0.61*temp_df['q'])) * temp_df['P']
    
    return rho_m


def get_windd(uwind,vwind):
    """
    Converts u and v wind components into meteorological degrees.
    meteorologica degrees = clockwise, north =0, direction wind is coming from.
    
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

def stationarity(x, y):
    """
    Test for stationarity: 
    Divide averaging period into smaller segments say 5 minutes each. 
    Calculate the covariance of each interval, then calculate the mean. 
    Also calculate the covariance of the whole interval. 
    If there is a difference of less than 30% between the covariances, 
    then the measurement is considered to be stationary. [foken and wichura 1995]
        
    Input:
        x, y: variables to test
    Output:
        fs: list of stationarity parameter from each period, defined above. 
        p_cov: list of covariance between variables for each period. 
        rol_cov: list of 5 min rolling covariances between variables for each period. 
    """    
    xw=detrend(x)
    yw=detrend(y)
 
    df = pd.DataFrame(list(zip(xw, yw)), 
               columns =['x', 'y'])

    # Calculate covariance over the entire period. 
        
    p_cov = df['x'].cov(df['y'])
        
    # Calculate a rolling covariance, window size = 5 mins (3000 data points).
        
    rol_cov = df.rolling(3000).cov().unstack()['x']['y']

    # Calculate stationarity parameter. 
    
    mean_segs = np.mean(rol_cov)
    fs = np.abs((mean_segs - p_cov) / p_cov)
        
    return fs, p_cov, rol_cov


def skew(x):
    """
    Skew: calculate skew of turbulent variable
    The direction and magnitude of a datasets deviation from normal.
    sk = (1/N) * sum((x-xbar)**3) / sigma **3
    N = length, xbar = mean, sigma = std
    Abs sk > 2 = bad
    Abs sk > 1 = OK
    Otherwise = good. 
    
    Input:
        x: variable ot test
    Output:
        sk: skew of x
    """      
    N = len(x)
    xbar = np.mean(x)
    sigma = np.std(x)
    sk = (1.0/N) * sum((x-xbar)**3) / (sigma**3)
    
    return sk

def kurtosis(x):
    """
    Kurtosis: calculate kurtosis of turbulent variable
    Kurtosis is a measure of ouliers/ dispersion of data
    Pearson kurtossis
    Larger than 3=sharper than gaussian, smaller than 3 = flatter than gaussian
    kurt = (1/N) * sum((x-xbar)**4) / sigma**4
    kurt < 1 or kurt > 8 = bad
    kurt <2 or >5 = OK
    otherwise = good.  
    
    Input:
        x: variable ot test
    Output:
        kurt: Kurtosis of x
    """        
    N = len(x)
    xbar = np.mean(x)
    sigma = np.std(x)
    kurt = (1.0/N) * sum((x-xbar)**4) / (sigma**4)
    
    return kurt

def flux_devel_test(itc_w,sst):
    """
    Assign QC flag for flux development based on the 
    combination of the integral scale test result and 
    the test for stationarity. 
   
    Input:
        itc_w: Result of itc test
        sst:   Result of stationarity test
    Output:
        QC: QC flag associated with flux development. 
            0 = bad data
            1 = good data fundamental research
            2 = good data general use
            3 = Highly suspect data
    """    
    # Check for nan's
    if np.isnan(sst) or np.isnan(itc_w):
        QC = np.nan
    elif itc_w < 0 or sst < 0:
        QC = 0
    elif itc_w <=75 and sst <=30:
        QC = 1
    elif itc_w <=100 and sst <=100:
        QC = 2
    elif itc_w <=1000 and sst <=1000:
        QC = 3
    else:
        QC = 0     
    return QC

def skew_flag(skew):
    """
    Assign QC flag for skew
   
    Input:
        skew: Variable skew
    Output:
        QC: QC flag associated with skew. 
            0 = bad data
            1 = good data fundamental research
            2 = good data general use
    """    
    if np.isnan(skew):
        flag = 0
    elif np.abs(skew)>2:
        flag = 0
    elif np.abs(skew)>1:
        flag = 2
    else:
        flag = 1
    return flag
    
def kurt_flag(kurt):
    """
    Assign QC flag for kurtosis
   
    Input:
        kurt: Variable kurtosis
    Output:
        QC: QC flag associated with kurtosis. 
            0 = bad data
            1 = good data fundamental research
            2 = good data general use
    """    
    if np.isnan(kurt):
        flag = 0
    elif kurt<1 or kurt>8:
        flag = 0
    elif kurt <2 or kurt>5:
        flag = 2
    else:
        flag = 1
    return flag
