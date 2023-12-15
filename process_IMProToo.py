#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 2023

@author: heather

Process MRR raw spectra with IMProToo
Generate netcdf file

Inputs: 
    in_dir: Raw spectra directory
    in_dat: Date of raw input file, format 'MMDD'
    out_dir: Output location

Example usage: 
   cd /Users/heather/IMProToo/
   source venv/bin/activate
   python MRR_process.py $in_dir $in_dat $out_dir >> /Users/heather/IMProToo/output/logs/${in_dat}_log.txt 2>&1
"""

import sys,os
import IMProToo
import matplotlib.pyplot as plt
import xarray as xr
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
#from matplotlib.colors import DivergingNorm

def main():
    in_dir=sys.argv[1]
    in_dat=sys.argv[2]
    out_dir=sys.argv[3]
    
    print('Loading data...')
    rawData = IMProToo.mrrRawData(str(in_dir)+str(in_dat) + '.raw')
    processedSpec = IMProToo.MrrZe(rawData)
    
    processedSpec.co["ncCreator"] = "Heather Guy, University of Leeds, heather.guy@ncas.ac.uk"
    processedSpec.co["ncDescription"] = "MRR data processed with IMProToo"
    processedSpec.co["ncInstitution"] = "National Centre for Atmospheric Science, University of Leeds"
    processedSpec.co["ncLocation"] = "Ice breaker Oden foredeck, during the ARTofMELT field campaign. 78.13 to 80.52 degrees N, -3.87 to 14.14 degrees E."

    processedSpec.co["dealiaseSpectrum"] = True
    processedSpec.rawToSnow()
    
    print('Writing data...')
    # Write processed netcdf file
    YYYYMMDD='2023%s'%(in_dat)
    processedSpec.writeNetCDF(out_dir+'ncas-mrr-1_ARTofMELT_%s_IMProToo_v1.nc'%(YYYYMMDD),ncForm="NETCDF3_CLASSIC")
      
if __name__ == '__main__':
    main()  