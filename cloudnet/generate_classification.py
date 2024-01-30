
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday 29 Jan 2024

@author: heather guy

Generates categorize and classification files 
using the cloudnetpy algorithm
"""
import sys
import datetime as dt
from cloudnetpy.categorize import generate_categorize
from cloudnetpy.products import generate_classification

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    Input: date, data_dir
    date: in format str YYMMDD
    data_dir: Path to cloudnetpy data directory
    """    
    date = args_in[1]
    data_dir = args_in[2]
    try:
        d=dt.datetime.strptime(date,'%y%m%d')
        print('Processing %s'%d)
    except:
        print('Incorrect date format, should be YYMMDD')
        return None
    return date, data_dir

def main():
    if get_args(sys.argv) == None:
        print('Input Error')
        sys.exit()
    else:
        date, data_dir = get_args(sys.argv)
        #print(date)
        input_files = {'radar': data_dir+'radar/processed_nc/%s_radar.nc'%date,
                        'lidar': data_dir+'lidar/processed_nc/%s_CL.nc'%date,
                        'mwr': data_dir+'hatpro/processed_nc/%s_hatpro.nc'%date,
                        'model': data_dir+'model/%s_obs_sonde.nc'%date}
        
        cat_uuid = generate_categorize(input_files, data_dir+'categorize/%s_categorize.nc'%date)
        uuid = generate_classification(data_dir+'categorize/%s_categorize.nc'%date, data_dir+'classification/%s_classification.nc'%date)
    return
            
if __name__ == '__main__':   
    main()        