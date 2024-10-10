#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs March 12 17:22:21 2020

@author: Heather Guy

Functions to generate netCDF files. 

"""
import pandas as pd
import numpy as np

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

def NC_SpecificVariables(fn_nc, var, np):
    """
    Collects specific variables from standard excel file and
    writes to file. 
    
    Parameters:
        fn_nc : NetCDF file to write on. 
        var :   Specific variables DataFrame from standard
                excel file. Read by pd.read_excel()
        np:     numpy
    
    """
    type_dict = {'float32': np.float32,
             'float64': np.float64,
             'int32'  : np.int32,
             'byte'   : np.byte}
    
    for i in range(0,len(var['Variable'].dropna())):
        if i+1 == len(var['Variable'].dropna()):
            att_end = var.index[-1]
        else: 
            att_end = var['Variable'].dropna().index[i+1] 
           
        att_list = var['Attribute'][var['Variable'].dropna().index[i]:att_end+1].dropna()
        val_list = var['Value'][att_list.index]       
    
        varname = var['Variable'].dropna().iloc[i]
        datatype = type_dict.get(val_list[att_list[att_list=='type'].index].to_numpy()[0])
        dimensions = dimensions=val_list.iloc[1].split(', ')
    
        if len(att_list[att_list=='_FillValue'])==1:
            fv=float(val_list[att_list[att_list=='_FillValue'].index])
        else:
            fv=None
       
        varn = fn_nc.createVariable(varname,datatype,dimensions,fill_value=fv)
        
        for j in range(0,len(att_list)):
            if att_list.iloc[j][0]=='_':
                continue
            else:
                if att_list.iloc[j]=='units':
                    varn.setncattr(att_list.iloc[j],str(val_list.iloc[j]))
                else:
                    varn.setncattr(att_list.iloc[j],val_list.iloc[j])
                    #print(att_list.iloc[j],val_list.iloc[j])
    return


def NC_CommonVariables(fn_nc, time_list,lat_list,lon_list, np):
    """
    Writes common variables to netCDF file. 
    
    Parameters:
        fn_nc :     NetCDF file to write on. 
        time_list:  Pandas date_range for the time dimension
        np:         numpy
    
    """
    import netCDF4
    #time
    times = fn_nc.createVariable('time', np.float64, ('time',))
    #variable attributes
    times.type = 'float64'
    times.units = 'seconds since 1970-01-01 00:00:00 UTC'
    times.standard_name = 'time'
    times.long_name = 'Time (seconds since 1970-01-01 00:00:00)'
    times.axis = 'T'
    times.valid_min = np.float64(netCDF4.date2num(min(time_list),units=times.units))
    times.valid_max = np.float64(netCDF4.date2num(max(time_list),units=times.units))
    times.calendar = 'standard'
    #write data
    times[:] = np.float64(netCDF4.date2num(time_list.to_list(),units=times.units))
 
    #year
    years = fn_nc.createVariable('year', np.int32, ('time',))
    #variable attributes
    years.type = 'int32'
    years.units = '1'
    years.long_name = 'Year'
    years.valid_min = np.int32(min(time_list).year)
    years.valid_max = np.int32(max(time_list).year) 
    #write data
    years[:] = np.int32(time_list.year.to_numpy())

    #month
    months = fn_nc.createVariable('month', np.int32, ('time',))
    #variable attributes
    months.type = 'int32'
    months.units = '1'
    months.long_name = 'Month'
    months.valid_min = np.int32(min(time_list).month)
    months.valid_max = np.int32(max(time_list).month) 
    #write data
    months[:] = np.int32(time_list.month.to_numpy())
   
    #day
    days = fn_nc.createVariable('day', np.int32, ('time',))
    #variable attributes
    days.type = 'int32'
    days.units = '1'
    days.long_name = 'Day'
    days.valid_min = np.int32(min([d.day for d in time_list]))
    days.valid_max = np.int32(max([d.day for d in time_list]))
    #write data
    days[:] = np.int32(np.int32(time_list.day.to_numpy()))
   
    #hour
    hours = fn_nc.createVariable('hour', np.int32, ('time',))
    #variable attributes
    hours.type = 'int32'
    hours.units = '1'
    hours.long_name = 'Hour'
    hours.valid_min = np.int32(min([d.hour for d in time_list]))
    hours.valid_max = np.int32(max([d.hour for d in time_list] ))
    #write data
    hours[:] = np.int32(time_list.hour.to_numpy())
   
    #minute
    minutes = fn_nc.createVariable('minute', np.int32, ('time',))
    #variable attributes
    minutes.type = 'int32'
    minutes.units = '1'
    minutes.long_name = 'Minute'
    minutes.valid_min = np.int32(min([d.minute for d in time_list]))
    minutes.valid_max = np.int32(max([d.minute for d in time_list])) 
    #write data
    minutes[:] = np.int32(time_list.minute.to_numpy())
   
    #second
    seconds = fn_nc.createVariable('second', np.float32, ('time',))
    #variable attributes
    seconds.type = 'float32'
    seconds.units = '1'
    seconds.long_name = 'Second'
    seconds.valid_min = np.float32(min([d.second for d in time_list]))
    seconds.valid_max = np.float32(max([d.second for d in time_list]))
    #write data
    seconds[:] = np.int32(time_list.second.to_numpy())
   
    #doy
    doys = fn_nc.createVariable('day_of_year', np.float32, ('time',))
    all_doys = np.float32(np.asarray([decimaldayofyear(time_list[i]) for i in range(0,len(time_list))]))
    #variable attributes
    doys.type = 'float32'
    doys.units = '1'
    doys.long_name = 'Day of Year'
    doys.description = "time as decimal day of year"
    doys.valid_min = np.float32(min(all_doys))
    doys.valid_max = np.float32(max(all_doys))
    #write data
    doys[:] = np.float32(np.asarray([decimaldayofyear(time_list[i]) for i in range(0,len(time_list))]))
    
    lats = fn_nc.createVariable('latitude', np.float32, ('time',))
    #variable attributes
    lats.type = 'float32'
    lats.units = 'degree_north'
    lats.long_name = 'Latitude'
    lats.cell_method= 'time: point'
    lats[:] = lat_list
    lats.valid_min= np.float32(np.nanmin(lat_list))
    lats.valid_max= np.float32(np.nanmax(lat_list))
   
    lons = fn_nc.createVariable('longitude', np.float32, ('time',))
    #variable attributes
    lons.type = 'float32'
    lons.units = 'degree_east'
    lons.long_name = 'Longitude'
    lons.cell_method= 'time: point'
    lons[:] = lon_list
    lons.valid_min= np.float32(np.nanmin(lon_list))
    lons.valid_max= np.float32(np.nanmax(lon_list))
   
   
    return
      
   
def NC_Dimensions(fn_nc, len_time,index=False):
    """
    Writes dimensions to netCDF file
    
    Parameters:
        fn_nc :     NetCDF file to write on. 
        len_time :  Length of time dimension
        index:      length of optional index dimension
    
    """
    fn_nc.createDimension('time', len_time )
    #fn_nc.createDimension('latitude', len_time)
    #fn_nc.createDimension('longitude', len_time) 
    if index:
        index = fn_nc.createDimension('index', index)
    
    return


def NC_Global_Attributes(fn_nc, meta, start_date,end_date):
    """
    Writes global attributes to netCDF file
    
    Parameters:
        fn_nc :     NetCDF file to write on.
        meta :      DataFrame of global attributes from standard
                    excel file. Read by pd.read_excel()
        start_date: Start datetime of data in file
        end_date:   end datetime of data in file

    """
    from datetime import datetime
    import numpy as np
    name = meta.loc[:, 'Name':'Name':1].values
    exp = meta.loc[:, 'Example':'Example':1].values
    pos = exp[31]
    pos = str(pos[0])
    ix1 = pos.find('N')
    ix2 = pos.find(' ')
    ix3 = pos.find('E')
    try:
        lat = np.float32(pos[0:ix1])
        lon = np.float32(pos[ix2+1:ix3])
    except:
        lat = np.nan
        lon=np.nan
    pos = exp[32]
    pos = pos[0]
    ix1 = pos.find('m')
    base_height = np.float32(pos[0:ix1])
   
    for i in range(0,len(meta)):
       msg1 = np.array(name[i])
       msg2 = np.array(exp[i])
       fn_nc.setncattr(msg1[0], str(msg2[0]))
   
    fn_nc.last_revised_date = datetime.utcnow().isoformat()  
    fn_nc.time_coverage_start = start_date.isoformat()
    fn_nc.time_coverage_end = end_date.isoformat()
   
    return

