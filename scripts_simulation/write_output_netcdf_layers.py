#!/usr/bin/env py36
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 11:18:47 2017

write MEMLS output into netcdf taking into account that the input was on different amount
of layers

@author: Clara Burgad
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import datetime

##########################################################################

ee2='75N00W'
#ee2='NorthPole'

if ee2 == '75N00W':
    folder = 'exp053_75N00W_20191106'
elif ee2 == 'NorthPole':
    folder = 'exp054_NorthPole_20191106'


inputpath = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/INPUT/original_files/'
outputpath = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/INPUT/netcdf_files/'
inputpath2 = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/OUTPUT/original_files/'
outputpath2 = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/OUTPUT/netcdf_files/'

#################################################################

sens_exp = ['complex','simpleallfunc','simpleallconst']

#which amount of layers were computed?
layers = [3,5,7,10]


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

for ll in layers:
        #read in the brightness temperatures
        #ds_output = None
        for sens in sens_exp:
            print(sens,ll)
            ds2=None
            for i in range(0, 3285):
              #print(i)
              ff = os.path.join(inputpath2, 
                                'output_timestep_{0:04d}_{2:02d}_{1}.dat'.format(i,sens,ll))
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
              print(i, sens, ll)  
            ds2.to_netcdf(outputpath2+'outputMEMLS_'+sens+'_'+str(ll).zfill(2)+'_'+ee2+'.nc',mode='w')
          
        dss2 = {}
        ds_output = None
        for sens in sens_exp:
                dss2[sens] = xr.open_dataset(outputpath2+'outputMEMLS_'+sens+'_'+str(ll).zfill(2)+'_'+ee2+'.nc')
                dss2[sens] = dss2[sens].expand_dims('sens_exp')
                dss2[sens]['sens_exp'] = np.array((sens, ))
                if ds_output is None:
                  ds_output = dss2[sens].copy()
                else:
                  ds_output = xr.concat([ds_output, dss2[sens]], dim='sens_exp')
            #      print(sens)
        ds_output.to_netcdf(outputpath2+'outputMEMLS_'+str(ll).zfill(2)+'_'+ee2+'.nc',mode='w')

