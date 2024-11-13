#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday 29 Jan 2024

@author: heather guy

Functions for converting data to the input "model"
data required by the CLOUDNETpy algorithm. 
"""

from glob import glob
import pandas as pd
import numpy  as np
import datetime as dt
import os
import sys
import xarray as xr
import subprocess

def get_unixtime(dt64):
    # Converts datetime64 to unix timestamp in units: 
    # number of seconds that have elapsed since January 1, 1970
    #  (the Unix epoch) at 00:00:00 UTC
    return dt64.astype('datetime64[s]').astype('int')

def get_time_hours(time_array):
    # This converts a numpy array of datetime64's to
    # hours since the start of the file.
    sd = get_unixtime(time_array[0])
    # standard utctime from timestamp
    sd = dt.datetime.utcfromtimestamp(sd)
    # get start of day
    sd = sd.replace(hour=0,minute=0,second=0)
    # convert back to numpy.datetime64:
    start_date = np.datetime64(sd)
    return (get_unixtime(time_array) - get_unixtime(start_date))/60/60, pd.to_datetime(start_date)

def gen_prop_input(frequencies,temperatures,pressures,rhs,save_loc,overwrite=True):
    # Generates input files for E O'Connors mwps-0.8.1 program to 
    # calculate microwave propagation parameters
    #
    # Input: 
    # Microwave frequencies (GHz), 1D array
    # temperautres (K), array-like, shape is [time, height]
    # pressures (Pa), array-like, shape is [time, height]
    # rh's (unitless) (0 to 1), with respect to ice below 273.16 K, array-like, shape is [time, height]
    # save_loc: Directory for output
    # overwrite: If true deletes any existing input files in save_loc, else
    # appends to existing files. Caution using False here - may lead to unexpected
    # results if the propagtion output aren't interpreted correctly. 
    #
    # Output: 
    # Saves three input files in save_loc:
    # 'atmos.temp'     Using all inputs
    # 'atmos_sat.temp' RH hardcoded to 1
    # 'atmos_dry.temp' RH hardcoded to 0
    # Format of the files is as follows: 
    # fprintf(fid,'%6.2f %12.8f %12.8f %12.8f\n',[freq(:) temperature(:) pressure(:) rh(:)]');
    
    fid=save_loc+'atmos.temp'
    fid_sat=save_loc+'atmos_sat.temp'
    fid_dry=save_loc+'atmos_dry.temp'

    if overwrite==True:
        # Check if these files already exist, and delete them if they do
        # Check if the file exists
        for f in [fid,fid_sat,fid_dry]:
            if os.path.exists(f):
                # Delete the file
                os.remove(f)
                print(f"The file '{f}' has been deleted.")
    
    # Check flattend list length 
    ll = len(frequencies)*np.shape(rhs)[0]*np.shape(rhs)[1]

    # Create flattened lists
    # columns = freq,temp,pres,rh
    c=0
    for i in range(0,len(frequencies)):
        for j in range(0,np.shape(rhs)[0]):
            for k in range(0,np.shape(rhs)[1]):
                line = '%6.2f %12.8f %12.8f %12.8f\n'%(frequencies[i],temperatures[j,k],pressures[j,k],rhs[j,k])
                line_sat = '%6.2f %12.8f %12.8f %12.8f\n'%(frequencies[i],temperatures[j,k],pressures[j,k],np.ones(np.shape(rhs))[j,k])
                line_dry = '%6.2f %12.8f %12.8f %12.8f\n'%(frequencies[i],temperatures[j,k],pressures[j,k],np.zeros(np.shape(rhs))[j,k])       
                #print(line)
                # Write line to file
                # Open the file in append mode ('a')
                with open(fid, 'a') as file:
                    file.write(line)
                with open(fid_sat, 'a') as file:
                    file.write(line_sat)
                with open(fid_dry, 'a') as file:
                    file.write(line_dry)
                c=c+1

    print('Number of lines: %s'%c)
    if c==ll:
        print('Number of lines is as expected')
    else:
        print('Warning! Expected number of lines should be %s'%ll)
    return

def run_prop_cal(infil,command_loc):
    # Runs E O'Connors mwps-0.8.1 program to 
    # calculate microwave propagation parameters
    # propagation: reads frequency, temperature, pressure and 
    # humidity and returns propagation parameters from input files
    # generated by gen_prop_input. 
    # Then goes through and unpacks the output. 
    #
    # Input:
    # infil: path of standard input file generated by gen_prop_input
    # command_loc: Location of the command for mwps-0.8.1 propogation calc
    #
    # Output:
    # Returns 1D arrays of [gas_atten_1way,liq_atten_1way,k2]
    # gas_atten_1way: 1-way specific gaseous attenuation (dB/km)
    # liq_atten_1way: 1-way liquid water specific attenuation (dB/km)/(g/m3)
    # k2: dielectric parameter |K|^2
    
    command=command_loc + 'propagation ' + infil
    output=subprocess.run([command,'-verbose'],capture_output=True,text=True,shell=True)
    print('Error:',output.stderr)
    print('Return code: ',output.returncode)

    # parse output
    outstr = [line.split() for line in output.stdout.split('\n')]
    outfloat = [[float(i) for i in lst] for lst in outstr]
    # remove empty lines
    outfloat = list(filter(None, outfloat))
    # conver to array
    outarr = np.asarray(outfloat)
    gas_atten_1way=outarr[:,0]   # 1-way specific gaseous attenuation (dB/km)
    liq_dc_real=outarr[:,1]      # real part of the liquid water dielectric constant
    liq_dc_img =outarr[:,2]      # imaginary part of the liquid water dielectric constant\n"
    k2=outarr[:,3]               # dielectric parameter |K|^2
    liq_atten_1way = outarr[:,4] # 1-way liquid water specific attenuation (dB/km)/(g/m3)

    return [gas_atten_1way,liq_atten_1way,k2]

def calc_2way_gas_atten(specific_atten,frequencies, times,heights):
    # "Two-way attenuation from the ground due to atmospheric gases, dB"
    # To calculatue the 2-way attenuation, need to do a cumulative
    # sum with height and multiply by 2 for there and back. 
    #
    # Inputs: 
    # Specific gas attenuation in db/km with shape (frequencies,times,heights)
    # frequencies
    # times
    # heights in (m)
    # 
    # Outputs: 
    # array of shape (frequencies,times,heights) of two-way attenuation from 
    # ground in db. 

    gas_atten=np.empty(np.shape(specific_atten))*np.nan
    
    # Difference in km between height levels in km
    dh=np.diff(np.insert(heights, 0, 0))/1000. 
    
    for i in range(0,len(frequencies)):
        for j in range(0,len(times)):
            arr = specific_atten[i,j,:]
            twoway_cum = np.nancumsum(arr * dh) * 2
            gas_atten[i,j,:]=twoway_cum

    return gas_atten

def wind_to_uv(wspd,wdir):
    """
    calculated the u and v wind components from wind speed and direction
    Input:
        wspd: wind speed
        wdir: wind direction
    Output:
        u: u wind component
        v: v wind component
    """    
   
    rad = 4.0*np.arctan(1)/180.
    u = -wspd*np.sin(rad*wdir)
    v = -wspd*np.cos(rad*wdir)

    return u,v