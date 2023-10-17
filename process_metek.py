#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs May 25 23:55 2023

@author: Heather Guy
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import glob
import matplotlib.dates as md
import matplotlib.colors as colors
from flux_functions import *
flatten = lambda l: [item for sublist in l for item in sublist]
import os
import io
import sys


####### INPUTS #######
# Data location:
#in_loc = '/Users/heather/Desktop/ARTofMELT/ice-station-data/'
#save = 'Desktop/ARTofMELT/ice-station-data/metek_processed/'

#start='20230516'
#stop='20230522'

#Example usage: 
# python process_metek.py $in_loc $start $stop $save
#############################################################


def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    try:
        dpath = args_in[1]
        save = args_in[4]
        start = dt.datetime.strptime(args_in[2],'%Y%m%d')
        stop = dt.datetime.strptime(args_in[3],'%Y%m%d')
    except:
        print('Input error')
        sys.exit()
    # return values:
    return dpath,start,stop,save


def main():
    """
    Extracts 3D sonic data from raw output. 
    Cleans bad data.
    Converts to SI units.
    NOTE: THIS IS FOR USA METEK-100 CORRECTS Left-handed to right-handed coordinate (-ve y)
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
        
    Returns:
        m1, m2: Clean dataframe for lower 3D sonic (m1) and upper (m2)
    
    """
    
    # check / get args:
    
    dpath,start,stop,save = get_args(sys.argv)
    #except:
    #    print('Input error')
    #    sys.exit()
    
    # raw data format:
    #2019 04 02 16 23 41.734 M:x =    14 y =    -1 z =    12 t =  2357
    # M = measured data heater off
    # H = measured data heater on
    # D = measured data heater defect
    # x,y,z componenets of wind in cm/s
    # t = acoustic temperature 2357 = 23.57C
    
    os.chdir(dpath)                    # Change directory to where the data is
    all_files = glob.glob('*.metek')  # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]
    dfs.sort()

    # Initialise empty data frames
    m1 = pd.DataFrame()

    # Get data.
    for f in dfs: 
        
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue
        
        # Filter out bad data files
        try:
            fd = io.open(f,"r",errors='replace')
            f_dat = fd.readlines()
        except:
            print('Data error with %s'%f)
            continue   

        # Iterate through all data lines. 
        # Delete duplicated error lines (E)
        # Replace bad data lines with nans. 
        # Iterate backwards so you can delete in place.
        # Replace any lines of non standard length with nan line.

        nan_line = 'x =   nan y =   nan z =   nan t = nan  \n'
        for i in range(len(f_dat) -1 , -1, -1):
            if len(f_dat[i])!= 66:
                f_dat[i] = f_dat[i][0:26] + nan_line
            if f_dat[i][24]=='D':
                f_dat[i] = f_dat[i][0:26] + nan_line
            elif f_dat[i][24]=='E':
                f_dat[i-1] = f_dat[i-1][0:26] + nan_line
                del f_dat[i]

        # Filter out incomplete data files
        if len(f_dat) < 72000:
            print('Incomplete file %s'%f)
            continue
        
        # Ignore extra data points
        f_dat = f_dat[0:72000]
        
        # Now split up strings and put in pandas df.
        try:
            pdf = pd.DataFrame(f_dat)
            pdf[1] = pdf[0].str.split()
            df = pd.DataFrame(pdf[1].values.tolist(), columns=['Year','Month','Day','Hour','Minute','Second','Status','junk1','x','junk2','junk3','y','junk4','junk5','z','junk6','junk7','T'])
            df['Logger_Date'] = pd.to_datetime(df[['Year','Month','Day','Hour','Minute','Second']])
            del df['junk1'],df['junk2'],df['junk3'],df['junk4'],df['junk5'],df['junk6'],df['junk7'],df['Year'],df['Month'],df['Day'],df['Hour'],df['Minute'],df['Second']
        except:
            print('Column error with %s'%f)
            continue

        # Make a standard 20Hz time series. 
        # Start time
        st = pd.to_datetime(f[0:16],format='%y%m%d_%H%M%S.%f')
        # End time
        et = st + pd.Timedelta(hours=1)
        et = et - pd.Timedelta(seconds=0.05)
        # 10 Hz time series (72000 data points in one hour)
        df['Date'] = pd.date_range(st,et,periods=72000)
        
        # Tidy and sort units.
        df = df.set_index('Date')
        df['T'] = pd.to_numeric(df['T'], errors='coerce')
        df['T']=df['T'].astype(float)
        df['x'] = pd.to_numeric(df['x'], errors='coerce')
        df['x']=df['x'].astype(float)
        df['y'] = pd.to_numeric(df['y'], errors='coerce')
        df['y']=df['y'].astype(float)
        df['z'] = pd.to_numeric(df['z'], errors='coerce')
        df['z']=df['z'].astype(float)

        df['T']=df['T']/100.0

        df = df.sort_values('Date')
        df.index = pd.DatetimeIndex(df.index)
        #df = df[~df.index.duplicated()]

        # metek usa100 coordinate system: 
        # it's left-handed. If you just multiply the y-axis velocity by -1 
        # that sets it all back to something sensible.

        # metek.x = data(:,7)/100; % m/s +ve to North
        # metek.y = -data(:,8)/100; % m/s +ve to East
        #metek.z = data(:,9)/100; % m/s +ve up
        # metek.sontemp = data(:,10)/100; % degC
        
        # Change units from cm/s to m/s, (* 0.01), and T to kelvin
        
        df['x']=df['x']*0.01
        df['y']= - df['y']*0.01
        df['z']=df['z']*0.01
        df['T']= df['T']+ 273.15

        # Append to relavent dataframe
        #try             
        m1 = pd.concat([m1, df])
        #except:
        #    print('Data error with %s'%f)
        #    continue
            
    # crop data for date/time
    m1=m1[start:stop]
        
    if save: 
        m1.to_csv(save+'metek_%s'%str(start.date()))
   
    return m1 

if __name__ == '__main__':
    main()          



