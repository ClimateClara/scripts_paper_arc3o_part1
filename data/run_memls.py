#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:19:41 2017

run MEMLS from netcdf with the functions being in memls_functions.py

!BEWARE: we tried to compute the penetration depth but the formula was not verified and the result was
not used in our study! We do not recommend using this function!

@author: Clara Burgard
"""

import numpy as np
import pandas as pd
import glob
import xarray as xr
import memls_functions as memls


##########################################################################

#went once through all experiments
ee='75N00W-p4'
#ee='NorthPole-p4'
ee2=ee.split("-")[0]

inputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/INPUT/netcdf_files/'
outputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/OUTPUT/original_files'


####################################################

sens_exp = ['complex','simpleallfunc','simpletemp','simplesalfunc','simplesalconst','simpleallconst'] 

for sens in sens_exp:
  print(sens+' started')
  netcdf_input = xr.open_dataset(inputpath+'inputMEMLS_'+sens+'_'+ee2+'.nc') #open the input netcdf
  for tt,timet in enumerate(netcdf_input['time']):
    print(str(tt)+' '+sens+' '+ee)
    
    for snow_prop in ['snowyes','snowno']: #run two times: once with snow and once without
        if snow_prop == 'snowyes':
            workdata=netcdf_input.isel(time=tt)
        elif snow_prop == 'snowno':
            workdata0=netcdf_input.isel(time=tt)
            workdata = workdata0.where(workdata0['snow_ice'] == 1).dropna('top_100_to_bottom_0')
        
        num = range(1,len(workdata['temperature'].dropna('top_100_to_bottom_0').values)+1)
            
        if not num:
          frequency = [1.4, 6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0]
          eh = np.empty(len(frequency))
          eh[:] = np.nan
          ev = np.empty(len(frequency))
          ev[:] = np.nan
          TBH = np.empty(len(frequency))
          TBH[:] = np.nan
          TBV = np.empty(len(frequency))
          TBV[:] = np.nan
          TeffH = np.empty(len(frequency))
          TeffH[:] = np.nan
          TeffV = np.empty(len(frequency))
          TeffV[:] = np.nan
          pend_h = np.empty(len(frequency))
          pend_h[:] = np.nan
          pend_v = np.empty(len(frequency))
          pend_v[:] = np.nan
        else:
          #Ti = workdata['temperature'].dropna('top_100_to_bottom_0').values.copy()
          Ti = np.flipud(workdata['temperature'][0:num[-1]].values.copy())
          Wi = np.flipud(workdata['wetness'][0:num[-1]].values.copy())
          roi = np.flipud(workdata['density'][0:num[-1]].values.copy())
          pci = np.flipud(workdata['correlation_length'][0:num[-1]].values.copy())
          epci = np.flipud(workdata['exp_correlation_length'][0:num[-1]].values.copy())
          di = np.flipud(workdata['thickness'][0:num[-1]].values.copy())
          sal = np.flipud(workdata['salinity'][0:num[-1]].values.copy())
          sitype = np.flipud(workdata['type'][0:num[-1]].values.copy())
          si = np.flipud(workdata['snow_ice'][0:num[-1]].values.copy())
          print('read the input data')
          
          frequency, eh, ev, TBH, TBV, TeffH, TeffV, pend_h, pend_v = memls.memls_mod(np.array(num),di,Ti,Wi,roi,pci,sal,sitype,si)

          
          print('memls through')
        
        fileinput = pd.DataFrame()
        fileinput['Freq'] = frequency
        fileinput['eh'] = eh
        fileinput['ev'] = ev
        fileinput['TBH'] = TBH
        fileinput['TBV'] = TBV
        fileinput['TeffH'] = TeffH
        fileinput['TeffV'] = TeffV
        fileinput['Pen. depth H'] = pend_h
        fileinput['Pen. depth V'] = pend_v
        #fileinput['Pen depth'] = pend
        #ff = outputpath+exp+'/timestep'+str(tt).zfill(4)+'_orig2lay.dat'
        ff = outputpath+'/output_timestep_'+str(tt).zfill(4)+'_'+sens+'_'+snow_prop+'.dat'
        np.savetxt(ff, fileinput.values, fmt='%1.2f') 
    #print('file written')

    
  print(sens+' finished')
print(ee)