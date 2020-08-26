#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 11:11:53 2017

write MEMLS input into netcdf, once with normal layers and once with interpolations to constant
layer thicknesses

updated to python3 on 11.05.2018

@author: m300411
"""

import numpy as np
import pandas as pd
import xarray as xr
import os
import datetime

##########################################################################


ee2='75N00W'
#ee='NorthPole-p4'

if ee2 == '75N00W':
    folder = 'exp053_75N00W_20191106'
elif ee2 == 'NorthPole':
    folder = 'exp054_NorthPole_20191106'


inputpath = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/INPUT/original_files/'
outputpath = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/INPUT/netcdf_files/'

#################################################################

sens_exp = ['complex','simpleallfunc','simpleallconst']
sens_exp_simp = ['simpleallfunc','simpleallconst']


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

###############################################
#OPTIONS
###############################################

rewrite_orig = 'yes'
#rewrite_orig = 'no'
#rewrite_interp = 'yes'
rewrite_interp = 'no'

#layertype = 'same'
layertype = 'diff'

layers = [3,5,7,10]

#############################################
#original layer thicknesses
#############################################

for ll in layers:

    if rewrite_orig=='yes':
      #read in the input
      for sens in sens_exp:
        ds1=None
        for i in range(0, 3285):
          ff = os.path.join(inputpath, 
                            'timestep{0:04d}_{1}_{2:02d}layers.dat'.format(i,sens,ll))
          temp_data = pd.read_table(
              ff,
              delimiter=" ",
              names=['number','temperature', 'wetness',
                     'density', 'correlation_length',
                     'exp_correlation_length', 'thickness',
                     'salinity','type','snow_ice', 'liquid_water_fraction'
                 ],
          ).drop(['number'], axis=1)
          temp_data.index += 102-len(temp_data)
          temp_data['depth'] = np.cumsum(temp_data['thickness'][::-1])
          temp_data = temp_data.reindex(np.arange(101, 0, -1))
          temp_data.index.name = 'top_100_to_bottom_0'
          temp_ds = xr.Dataset.from_dataframe(temp_data)
          temp_ds = temp_ds.expand_dims('time')
          temp_ds['time'] = np.array((idate[i], ))
          if ds1 is None:
            ds1 = temp_ds.copy()
          else:
            ds1 = xr.concat([ds1, temp_ds], dim='time')
          print(i, sens, ll,'orig')
        ds1.to_netcdf(outputpath+'inputMEMLS_'+sens+'_'+str(ll).zfill(2)+'layers_'+ee2+'.nc',mode='w')
          
    #collect them all in one netcdf file
      if layertype == 'diff':
          dss = {}
          ds_input = None
          for sens in sens_exp_simp:
                print(sens)
                dss[sens] = xr.open_dataset(outputpath+'inputMEMLS_'+sens+'_'+str(ll).zfill(2)+'layers_'+ee2+'.nc')
                dss[sens] = dss[sens].expand_dims('sens_exp')
                dss[sens]['sens_exp'] = np.array((sens, ))
                if ds_input is None:
                  ds_input = dss[sens].copy()
                else:
                  ds_input = xr.concat([ds_input, dss[sens]], dim='sens_exp')
            #      print sens 
          ds_input.to_netcdf(outputpath+'inputMEMLS_'+str(ll).zfill(2)+'layers_'+ee2+'_simpall.nc',mode='w')
        
        
    
    ###############################################################################
    #depth interpolated
    ###############################################################################
      
    if rewrite_interp=='yes':
      #read in the input
      #write different netcdfs for what is needed
      for sens in sens_exp:
      #for sens in ['complex']:
        ds1=None
        for tt in range(0, 3285):
        #for tt in range(120, 150):
          ff = os.path.join(inputpath, 
                            'timestep{0:04d}_{1}_{2:02d}layers.dat'.format(tt,sens,ll))
          temp_data = pd.read_table(
              ff,
              delimiter=" ",
              names=['number','temperature', 'wetness',
                     'density', 'correlation_length',
                     'exp_correlation_length', 'thickness',
                     'salinity','type','snow_ice', 'liquid_water_fraction'
                 ],
          ).drop(['number'], axis=1)
          temp_data.index += 102-len(temp_data)
          temp_data['depth'] = np.cumsum(temp_data['thickness'][::-1])
          temp_data = temp_data.reindex(np.arange(101, 0, -1))
          temp_data.index.name = 'top_100_to_bottom_0'
          new_var = pd.DataFrame()
          dd = temp_data['depth'].max()
          if np.isnan(dd):
              dd = 0.0
          new_var['depth']=np.arange(0,dd+0.011,0.01)
          new_var.index = new_var['depth']
          variables = ['temperature','wetness','density','correlation_length','exp_correlation_length','thickness','salinity','type','snow_ice']
          for vv in variables:
            new_var[vv] = np.interp(new_var['depth'],temp_data['depth'],temp_data[vv]) 
          #new_var['liquid_water_fraction'] = memls.Vb(new_var['temperature'].values,new_var['salinity'].values)  
          temp_ds = xr.Dataset.from_dataframe(new_var)
          temp_ds = temp_ds.expand_dims('time')
          temp_ds['time'] = np.array((idate[tt], ))
          if ds1 is None:
            ds1 = temp_ds.copy()
          else:
            ds1 = xr.concat([ds1, temp_ds], dim='time')
          print(tt, sens, ll, 'interp')
        ds1.to_netcdf(outputpath+'inputMEMLS_depthidx_'+sens+'_'+str(ll).zfill(2)+'layers_'+ee2+'.nc',mode='w')
          
      #collect them all in one netcdf file
      dss = {}
      ds_input = None
      for sens in sens_exp:
            print(sens)
            dss[sens] = xr.open_dataset(outputpath+'inputMEMLS_depthidx_'+sens+'_'+str(ll).zfill(2)+'layers_'+ee2+'.nc')
            dss[sens] = dss[sens].expand_dims('sens_exp')
            dss[sens]['sens_exp'] = np.array((sens, ))
            if ds_input is None:
              ds_input = dss[sens].copy()
            else:
              ds_input = xr.concat([ds_input, dss[sens]], dim='sens_exp')
        #      print sens 
      ds_input.to_netcdf(outputpath+'inputMEMLS_depthidx_'+str(ll).zfill(2)+'layers_'+ee2+'.nc',mode='w')

  
