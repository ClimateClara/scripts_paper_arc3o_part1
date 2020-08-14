#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:19:41 2017

run MEMLS from netcdf with the functions being in memls_functions.py
THis script especially makes a loop over different amount of layers

!BEWARE: we tried to compute the penetration depth but the formula was not verified and the result was
not used in our study! We do not recommend using this function!

@author: Clara Burgard
"""

import numpy as np
import pandas as pd
import glob
import xarray as xr
import clalib.memls_functions as memls


##########################################################################

#went once through all experiments
ee='75N00W-p4'
#ee='NorthPole-p4'
ee2=ee.split("-")[0]

inputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/INPUT/netcdf_files/'
outputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/OUTPUT/original_files'


####################################################

#what layer amounts were prepared?
layers = [3,5,7,10]


sens_exp = ['complex','simpleallfunc','simpleallconst']
for ll in layers:
    for sens in sens_exp:
      print(sens+' '+str(ll)+' started')
      netcdf_input = xr.open_dataset(inputpath+'inputMEMLS_'+sens+'_'+str(ll).zfill(2)+'layers_'+ee2+'.nc')
      for tt,timet in enumerate(netcdf_input['time']):
        print(str(tt)+' '+sens+' '+ee+' '+str(ll))
        
        workdata=netcdf_input.isel(time=tt)
        
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
          #pci = zeros(len(workdata['pci']))
          #pci[:] = 0.5
          epci = np.flipud(workdata['exp_correlation_length'][0:num[-1]].values.copy())
          di = np.flipud(workdata['thickness'][0:num[-1]].values.copy())
          sal = np.flipud(workdata['salinity'][0:num[-1]].values.copy())
          sitype = np.flipud(workdata['type'][0:num[-1]].values.copy())
          #sitype = zeros(len(workdata['sitype']))
          #sitype[:] = 3.
          si = np.flipud(workdata['snow_ice'][0:num[-1]].values.copy())
          #print('read the input data')
          
          frequency, eh, ev, TBH, TBV, TeffH, TeffV, pend_h, pend_v = memls.memls_mod(np.array(num),di,Ti,Wi,roi,pci,sal,sitype,si)
          #print 'pend_h =',pend_h
          #print 'pend_v =',pend_v
          
          #print('memls through')
        
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
        ff = outputpath+'/output_timestep_'+str(tt).zfill(4)+'_'+str(ll).zfill(2)+'_'+sens+'.dat'
        np.savetxt(ff, fileinput.values, fmt='%1.2f') 
        #print('file written')
    
        
      print(sens+' finished')
    print(ee)
print(ll)