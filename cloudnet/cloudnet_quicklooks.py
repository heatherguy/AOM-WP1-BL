#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 13 01:06 2023

@author: heather

Generate and save cloudnet quicklook plots

Inputs: 
    in_dir: data input directory
    in_dat: Date to process, string 'YYMMDD'
    out_dir: Quicklooks save directory

Example usage: 
   cd /Users/heather/cloudnetpy/
   source venv/bin/activate
   in_dir='/Users/heather/cloudnetpy/data/'
   in_dat='0510'
   out_dir='/Users/heather/cloudnetpy/data/quicklooks/'
   
   python cloudnet_quicklooks.py $in_dir $in_dat $out_dir
"""

from cl_nport_dat import *
import sys

# ARTofMELT instruments: 
from cloudnetpy.instruments import rpg2nc
from cloudnetpy.instruments import ceilo2nc
from cloudnetpy.instruments import hatpro2nc

from cloudnetpy.categorize import generate_categorize
from cloudnetpy.products import generate_classification

# plotting
from cloudnetpy.plotting import generate_figure

def main():
    in_dir=sys.argv[1]
    in_dat=sys.argv[2]
    out_dir=sys.argv[3]
    
    site_meta={'name': 'AOM23', 'altitude': 50}
    
    # ArtOfMelt HATPro
    print('processing hatpro...')
    uuid = hatpro2nc(in_dir+'hatpro/%s/'%in_dat,in_dir+'hatpro/processed_nc/%s_hatpro.nc'%in_dat, site_meta)
    
    # ArtOfMelt cloud radar
    # RPG-FMCW-94
    # instruments.rpg2nc(path_to_l1_files: str, output_file: str, site_meta: dict, uuid: Optional[str] = None, date: Optional[str] = None)→ tuple[str, 
    print('processing radar...')
    uuid = rpg2nc(in_dir+'radar/%s/'%in_dat,in_dir+'radar/processed_nc/%s_radar.nc'%in_dat, site_meta)
    
    # ArtofMELT ceilometer: Vaisala CL31
    print('processing lidar...')
    in_cl=in_dir+'lidar/%s/'%in_dat
    out_cl=in_dir+'lidar/'
    # preprocess
    cl_nport_dat(in_dat,in_cl,out_cl)
    
    uuid=ceilo2nc(in_dir+'lidar/%s_CL.DAT'%in_dat,in_dir+ 'lidar/processed_nc/%s_CL.nc'%in_dat, site_meta)

    input_files = {
    'radar': in_dir + 'radar/processed_nc/%s_radar.nc'%in_dat,
    'lidar': in_dir +'lidar/processed_nc/%s_CL.nc'%in_dat,
    'mwr': in_dir + 'hatpro/processed_nc/%s_hatpro.nc'%in_dat
    }
    
    
    # Make figures
    print('Generating quicklooks...')
    fig=generate_figure(input_files['radar'], ['Zh', 'v', 'width'],image_name=out_dir+in_dat+'_radar.png')
    
    fig=generate_figure(input_files['lidar'], ['beta'],max_y=6,image_name=out_dir+in_dat+'_CL.png')
    
    fig=generate_figure(input_files['mwr'], ['lwp','iwv'],image_name=out_dir+in_dat+'_mwr.png')

    
if __name__ == '__main__':
    main()  