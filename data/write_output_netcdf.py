#!/usr/bin/env py36
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 11:18:47 2017

write MEMLS output into netcdf

@author: Clara Burgard
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import datetime
import glob

##########################################################################

#went once through all experiments
ee='75N00W-p4'
#ee='NorthPole-p4'
ee2=ee.split("-")[0]

inputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/INPUT/original_files/'
outputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/INPUT/netcdf_files/'
inputpath2 = '/work/mh0033/m300411/SatSim/MEMLS_exp/OUTPUT/original_files/'
outputpath2 = '/work/mh0033/m300411/SatSim/MEMLS_exp/OUTPUT/netcdf_files/'

#################################################################

sens_exp = ['complex','simpletemp','simplesalfunc','simpleallfunc','simplesalconst','simpleallconst']



#time goes from July 2005 to December 2009
#1120176000 : 01.07.2005
#1262217600 : 31.12.2009
# create timestamps in between, with 2 steps per day (12 hours steps)
base0 = datetime.datetime(2005,7,1,0,0,0)
numdays = 1643
date_list0 = [base0 + datetime.timedelta(days=x) for x in range(0, numdays)]
base1 = datetime.datetime(2005,7,1,12,0,0)      
date_list1 = [base1 + datetime.timedelta(days=x) for x in range(0, numdays)]    

idate = []
for jj in range(1643):
  idate.append(date_list0[jj])
  if jj < 1642:
    idate.append(date_list1[jj])


#read in the brightness temperatures
#ds_output = None
if icepack == 'no':
    for sens in sens_exp:
        print(sens)
        for snow_prop in ['snowno','snowyes']:
        	ds2=None
        	for i in range(0, 3285):
          		#print(i)
          		ff = os.path.join(inputpath2, 
                            'output_timestep_{0:04d}_{1}_{2}.dat'.format(i,sens,snow_prop))
          		temp_data = pd.read_table(
              				ff,
              				delimiter=" ",
              				names=['frequency','emis_h', 'emis_v',
                     		'tb_h', 'tb_v',
                     		'teff_h', 'teff_v',
                     		'pend_h', 'pend_v'
                     		],
              				)
          		temp_data = temp_data.set_index(['frequency'])
          		temp_ds = xr.Dataset.from_dataframe(temp_data)
          		temp_ds = temp_ds.expand_dims('time')
          		temp_ds['time'] = np.array((idate[i], ))
          		if ds2 is None:
          			ds2 = temp_ds.copy()
          		else:
          			ds2= xr.concat([ds2, temp_ds], dim='time')
          		print(i, sens)  
        	ds2.to_netcdf(outputpath2+'outputMEMLS_'+sens+'_'+ee2+'_'+snow_prop+'.nc',mode='w')
        

    dss2 = {}
    ds_output = None
    for sens in sens_exp:
            dss2[sens] = xr.open_dataset(outputpath2+'outputMEMLS_'+sens+'_'+ee2+'_snowno.nc')
            dss2[sens] = dss2[sens].expand_dims('sens_exp')
            dss2[sens]['sens_exp'] = np.array((sens, ))
            if ds_output is None:
              ds_output = dss2[sens].copy()
            else:
              ds_output = xr.concat([ds_output, dss2[sens]], dim='sens_exp')
        #      print(sens)
    ds_output.to_netcdf(outputpath2+'outputMEMLS_'+ee2+'_snowno.nc',mode='w')
    
    dss2 = {}
    ds_output = None
    for sens in sens_exp:
            dss2[sens] = xr.open_dataset(outputpath2+'outputMEMLS_'+sens+'_'+ee2+'_snowyes.nc')
            dss2[sens] = dss2[sens].expand_dims('sens_exp')
            dss2[sens]['sens_exp'] = np.array((sens, ))
            if ds_output is None:
              ds_output = dss2[sens].copy()
            else:
              ds_output = xr.concat([ds_output, dss2[sens]], dim='sens_exp')
        #      print(sens)
    ds_output.to_netcdf(outputpath2+'outputMEMLS_'+ee2+'_snowyes.nc',mode='w')


