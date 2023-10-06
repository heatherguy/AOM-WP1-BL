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
#save = '/Users/heather/Desktop/ARTofMELT/ice-station-data/licor_processed/'

#start='202305160000'
#stop='202305220000'

#Example usage: 
# python process_licor.py $in_loc $start $stop $save
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


def main():
    """
    Extracts licor data from raw output. 
    Corrects for data transmission delay. 
    QC's licor data. 
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
        
    Returns:
        Clean licor dataframe
    
    """
    
    # check / get args:
    
    dpath,start,stop,save = get_args(sys.argv)
    #except:
    #    print('Input error')
    #    sys.exit()

    
    # Raw data format: 
    #2019 04 03 11 11 56.453 89	189	0.16469	35.4518	0.04404	297.105	20.74	99.0	1.5224
    # Date, Ndx, DiagVal, CO2R, CO2D, H2OR, H2OD, T, P, cooler

    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*.licor')  # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]
    dfs.sort()

    # Initialise empty data frames
    licor = pd.DataFrame()

    # Function to convert to float otherwise write nan
    def convert_float(x):
        try:
            return np.float(x)
        except:
            return np.nan

    # Extract the data
    for f in dfs: 

        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue 
        
        try:
            headers = ['Year','Month','Day','Hour','Minute','Second','DiagV','CO2R','CO2D','H2OR','H2OD','T','P']
            li = pd.read_csv(f, names=headers, delim_whitespace=True,error_bad_lines=False)
        except:
            print('Data Error with: '+f+'\n')
            continue        

        # Filter out incomplete data files
        if len(li) < 71995:
            print('Incomplete file %s'%f)
            continue
        
        # Ignore extra data points
        li = li[0:72000]
 
        # Sort out the date referencing
        li['Logger_Date']=pd.to_datetime(li[['Year','Month','Day','Hour','Minute','Second']])
        li.index=li['Logger_Date']- pd.Timedelta(seconds=0.2)
        li = li[~li.index.duplicated()]

        # Make a standard 20Hz time series. 

        # DONT FORGET, LICOR DATA TRANSMISSION IS DELAYED BY 200MS (TWO DATA POINTS)
        # EACH DATA POINT SHOULD BE SHIFTED BACK IN TIME BY THIS AMOUNT. ••

        # Start time
        st = pd.to_datetime(f[0:16],format='%y%m%d_%H%M%S.%f') - pd.Timedelta(seconds=0.2) # •• Implemented here. 
        # End time
        et = st + pd.Timedelta(hours=1)
        et = et - pd.Timedelta(seconds=0.05)
        # 20 Hz time series (72000 data points in one hour)
        #li['Date'] = pd.date_range(st,et,periods=72000)
        
        # new datafrom to append to 
        li_reindex = pd.DataFrame(index=pd.date_range(st,et,periods=72000),columns=['DiagV','CO2R','CO2D','H2OR','H2OD','T','P'])
        
        # Sort the rest of the columns
        li_reindex['DiagV']=li['DiagV'].astype('int').reindex(li_reindex.index,method='nearest',tolerance='0.05s')
        li_reindex['CO2R'] = li['CO2R'].reindex(li_reindex.index,method='nearest',tolerance='0.05s')
        li_reindex['CO2R']=li_reindex['CO2R'].apply(convert_float)
        li_reindex['CO2D'] = li['CO2D'].reindex(li_reindex.index,method='nearest',tolerance='0.05s')
        li_reindex['CO2D']=li_reindex['CO2D'].apply(convert_float)/1000
        li_reindex['H2OR'] = li['H2OR'].reindex(li_reindex.index,method='nearest',tolerance='0.05s')
        li_reindex['H2OR']=li_reindex['H2OD'].apply(convert_float)
        li_reindex['H2OD'] = li['H2OD'].reindex(li_reindex.index,method='nearest',tolerance='0.05s')
        li_reindex['H2OD']=li_reindex['H2OD'].apply(convert_float)/1000 # mol/m3
        li_reindex['T'] = (li['T'].astype('float')+273.15).reindex(li_reindex.index,method='nearest',tolerance='0.05s')      # K
        li_reindex['P'] = (li['P'].astype('float')*1000 ).reindex(li_reindex.index,method='nearest',tolerance='0.05s')      # Pa
         
        licor = licor.append(li_reindex)
    
    # OC licor data
    # Interp diagnostic bitmap
    #The cell diagnostic value (diag) is a 1 byte unsigned integer (value between 0 and 255) with the
    #following bit map:
    #bit 7 bit 6 bit 5 bit 4 bit 3 bit 2 bit 1 bit 0
    #Chopper Detector PLL Sync <------------------------ ----AGC / 6.25 ------------------------>
    #1=ok 1=ok 1=ok 1=ok
    #Example: a value is 125 (01111101) indicates Chopper not ok, and AGC = 81% (1101 is 13,
    #times 6.25)
    
    #   For now 1 is good, zero is bad, 2=no met data
    licor['QC']=np.ones(len(licor))
    diag_list = licor['DiagV'].to_list()
    chopper = []
    detector = []
    pll = []
    sync = []
    agc = []
    
    for i in range(0,len(diag_list)):
        try:
            binar = format(int(diag_list[i]), "008b")
        except:
            chopper.append(np.nan)
            detector.append(np.nan)
            pll.append(np.nan)
            sync.append(np.nan)
            agc_temp = np.nan
            agc.append(np.nan)
            continue
        
        chopper.append(int(binar[0]))
        detector.append(int(binar[1]))
        pll.append(int(binar[2]))
        sync.append(int(binar[3]))
        agc_temp = binar[4:]
        agc.append(int(agc_temp, 2)* 6.25)
    
    licor['chopper']=chopper
    licor['detector']=detector
    licor['pll']=pll
    licor['sync']=sync
    licor['agc']=agc

    # Filter out bad values
    #licor['QC'][licor['chopper']==0] = 0 
    #licor['QC'][licor['detector']==0] = 0
    licor.loc[licor['pll']==0,'QC'] = 0
    licor.loc[licor['sync']==0,'QC'] = 0
    licor.loc[licor['agc']>90,'QC'] = 0
    licor.loc[np.isnan(licor['pll']),'QC'] = 0
    licor.loc[np.isnan(licor['sync']),'QC'] = 0
    licor.loc[np.isnan(licor['agc']),'QC'] = 0

    # Clean pressure for values over 80000 or under 50000
    # Clean temperature for values under 200 or over 300
    licor.loc[licor['P']<50000,'QC']=0
    licor.loc[licor['P']>110000,'QC']=0
    licor.loc[licor['T']<200,'QC']=0
    licor.loc[licor['T']>300,'QC']=0

    # Filter clear outliers.
    jj = ~np.isnan(licor['H2OD']) # Not nan indices
    sd = np.std(licor['H2OD'][jj]) # standard deviation 
    licor['H2OD']=replace_outliers(licor['H2OD'],sd)
        
    # Remove negative concentrations
    licor.loc[licor['H2OD']<0,'QC']=0
    licor.loc[licor['H2OR']<0,'QC']=0

    if save: 
        licor.to_csv(save+'licor_%s'%str(start.date()))

    return licor    

    
if __name__ == '__main__':
    main()          

