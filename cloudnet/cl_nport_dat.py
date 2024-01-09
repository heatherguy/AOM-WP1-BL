#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 22:28 2023

@author: heather

Preprocessing function to transform nport logged ceilometer data
into .dat format used by cloudnetpy
"""

import numpy as np
import re
import glob
import datetime as dt

flatten = lambda l: [item for sublist in l for item in sublist]

def cl_nport_dat(date,inloc,outloc):
    """
    Parameters:
        date: String 'YYMMDD'
        inloc: input data directory
        outloc: output data directory
        
    Returns:
        Saves .dat file in output directory
    
    """
    
    # Get daily files
    fils=glob.glob(inloc+date+'*.CL')
    fils.sort()
    
    # Join one day of data
    dat=[]
    for fil in fils: 
        # Modify the file to make input for cl2nc
        with open(fil,mode='r') as file: 
            try:
                X = file.readlines()
                dat.append([X[n] for n in range(0,len(X))])
            except:
                print('Data issue with %s'%fil)
                continue  
            
    new_dat = flatten(dat)
    ll = len(new_dat)
    
    # find the start of a message cycle
    cycle=[]
    for n in range(0,ll):
        temp=new_dat[n]
        if np.char.find(temp,'CL')!=-1:      
            cycle.append(n)
    
    # Write the output file
    with open(outloc+'%s_CL.DAT'%date,mode='w') as file: 
        for i in range(0,len(cycle)-1):
            d = dt.datetime.strptime(new_dat[cycle[i]][0:23],'%Y %m %d %H %M %S.%f')
            d_str = dt.datetime.strftime(d,'-%Y-%m-%d %H:%M:%S\n')
            file.write(d_str)
            file.writelines(nd[24:] for nd in new_dat[cycle[i]:cycle[i+1]])
            file.write('\n')
    
    return
