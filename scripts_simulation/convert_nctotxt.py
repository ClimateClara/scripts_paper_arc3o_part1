#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 11:47:17 2018

prepare input for SAMSIM from ERA data

@author: claraburgard
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
import scipy.io as sio
import datetime

#inputpath = "/home/mpim/m300411/SatSim/SAMSIM/input/ERA-interim/"
inputpath = "/Users/claraburgard/mistral_home/SatSim/SAMSIM/input/ERA-interim/"
outputpath = "/Users/claraburgard/mistral_home/SatSim/SAMSIM/input/ERA-interim/"


d = {}

both = 'True'
# both = False

# for ee in ['75N00W','75N180E','80N90E','NorthPole','82N50E','74N170E','80N160W','85N50W','82N120W']:
for ee in ['77N39E']:
  if both == 'True':
    file0 = glob.glob(inputpath + 'ERA_interim_' + ee + '_forecast_2005-2009.nc')
  else:
    file0 = glob.glob(inputpath + 'ERA-' + ee + '.nc')
  fid0 = sio.netcdf_file(file0[0])

  precip = fid0.variables['tp']
  d['precip_' + ee] = precip[:, 0, 0] * precip.scale_factor + precip.add_offset

  snowfall = fid0.variables['sf']
  d['snowfall_' + ee] = snowfall[:, 0, 0] * snowfall.scale_factor + snowfall.add_offset

  swdown = fid0.variables['ssrd']
  d['swdown0_' + ee] = (swdown[:, 0, 0] * swdown.scale_factor + swdown.add_offset)  # /(3*3600.)
  d['swdown_' + ee] = np.zeros(len(d['swdown0_' + ee]))
  for n, j in enumerate(range(0, len(d['swdown0_' + ee]), 4)):
    d['swdown_' + ee][j + 0] = d['swdown0_' + ee][j + 0] / (3 * 3600.)
    d['swdown_' + ee][j + 1] = d['swdown0_' + ee][j + 1] / (6 * 3600.)
    d['swdown_' + ee][j + 2] = d['swdown0_' + ee][j + 2] / (9 * 3600.)
    d['swdown_' + ee][j + 3] = d['swdown0_' + ee][j + 3] / (12 * 3600.)

  if both == 'True':
    swnet = fid0.variables['ssr']
    d['swnet0_' + ee] = (swnet[:, 0, 0] * swnet.scale_factor + swnet.add_offset)  # /(3*3600.)
    d['swnet_' + ee] = np.zeros(len(d['swnet0_' + ee]))
    for n, j in enumerate(range(0, len(d['swnet0_' + ee]), 4)):
      d['swnet_' + ee][j + 0] = d['swnet0_' + ee][j + 0] / (3 * 3600.)
      d['swnet_' + ee][j + 1] = d['swnet0_' + ee][j + 1] / (6 * 3600.)
      d['swnet_' + ee][j + 2] = d['swnet0_' + ee][j + 2] / (9 * 3600.)
      d['swnet_' + ee][j + 3] = d['swnet0_' + ee][j + 3] / (12 * 3600.)
    d['swup_' + ee] = d['swdown_' + ee] - d['swnet_' + ee]

  lwdown = fid0.variables['strd']
  d['lwdown0_' + ee] = (lwdown[:, 0, 0] * lwdown.scale_factor + lwdown.add_offset)  # /(3*3600.)
  d['lwdown_' + ee] = np.zeros(len(d['lwdown0_' + ee]))
  for n, j in enumerate(range(0, len(d['lwdown0_' + ee]), 4)):
    d['lwdown_' + ee][j + 0] = d['lwdown0_' + ee][j + 0] / (3 * 3600.)
    d['lwdown_' + ee][j + 1] = d['lwdown0_' + ee][j + 1] / (6 * 3600.)
    d['lwdown_' + ee][j + 2] = d['lwdown0_' + ee][j + 2] / (9 * 3600.)
    d['lwdown_' + ee][j + 3] = d['lwdown0_' + ee][j + 3] / (12 * 3600.)

  if both == 'True':
    lwnet = fid0.variables['str']
    d['lwnet0_' + ee] = (lwnet[:, 0, 0] * lwnet.scale_factor + lwnet.add_offset)  # /(3*3600.)
    d['lwnet_' + ee] = np.zeros(len(d['lwnet0_' + ee]))
    for n, j in enumerate(range(0, len(d['lwnet0_' + ee]), 4)):
      d['lwnet_' + ee][j + 0] = d['lwnet0_' + ee][j + 0] / (3 * 3600.)
      d['lwnet_' + ee][j + 1] = d['lwnet0_' + ee][j + 1] / (6 * 3600.)
      d['lwnet_' + ee][j + 2] = d['lwnet0_' + ee][j + 2] / (9 * 3600.)
      d['lwnet_' + ee][j + 3] = d['lwnet0_' + ee][j + 3] / (12 * 3600.)
    d['lwup_' + ee] = d['lwdown_' + ee] - d['lwnet_' + ee]

  precip = fid0.variables['tp']
  d['precip0_' + ee] = (precip[:, 0, 0] * precip.scale_factor + precip.add_offset)  # /(3*3600.)
  d['precip_' + ee] = np.zeros(len(d['precip0_' + ee]))
  for n, j in enumerate(range(0, len(d['precip0_' + ee]), 4)):
    d['precip_' + ee][j + 0] = d['precip0_' + ee][j + 0] / (3 * 3600.)
    d['precip_' + ee][j + 1] = d['precip0_' + ee][j + 1] / (6 * 3600.)
    d['precip_' + ee][j + 2] = d['precip0_' + ee][j + 2] / (9 * 3600.)
    d['precip_' + ee][j + 3] = d['precip0_' + ee][j + 3] / (12 * 3600.)

  #  lwdown=fid0.variables['strd']
  #  d["lwdown0_"+ee] = (lwdown[:,0,0]*lwdown.scale_factor + lwdown.add_offset)/(6*3600.)

  u10 = fid0.variables['u10']
  d['u10_' + ee] = u10[:, 0, 0] * u10.scale_factor + u10.add_offset

  v10 = fid0.variables['v10']
  d['v10_' + ee] = v10[:, 0, 0] * v10.scale_factor + v10.add_offset

  temp = fid0.variables['t2m']
  d['temp_' + ee] = temp[:, 0, 0] * temp.scale_factor + temp.add_offset

  dewtemp = fid0.variables['d2m']
  d['dewtemp_' + ee] = dewtemp[:, 0, 0] * dewtemp.scale_factor + dewtemp.add_offset

  d['time_' + ee] = fid0.variables['time'][:]
  base0 = datetime.datetime(1979, 1, 1, 0, 0, 0)
  numdays = int(len(d['precip0_' + ee]) / 8)
  date_list0 = [base0 + datetime.timedelta(days=x) for x in np.arange(0, numdays, 0.125)]

  # d["time_"+ee] = time[:]*time.scale_factor + time.add_offset
  # d["time_{0}".format(ee)]=fid0.variables['time'][:]

  tt = np.argsort(d['time_' + ee])
  for var in ['precip', 'snowfall', 'swdown', 'swup', 'lwdown', 'lwup', 'u10', 'v10', 'temp', 'dewtemp', 'time']:
    print(var + '_' + ee + '_sorted')
    d[var + '_' + ee + '_sorted'] = d[var + '_' + ee][tt]

  # precip0 = array(d['precip_'+ee+'_sorted'])
  # precip0[precip0<=0.0001] = 0.00
  # d['precip_'+ee+'_sorted'] = precip0

  # time
  # 3-hourly values
  timesteps = len(d['time_' + ee + '_sorted'])
  daystep = timesteps / 8.

  # create different subplots
  f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)
  plt.title(ee)

  # SHORTWAVE
  # figure()
  ax1.plot(np.arange(timesteps) / 8., d['swdown_' + ee + '_sorted'], 'r-')
  ax1.set_title('Shortwave radiation ' + ee)
  ax1.set_ylabel('Shortwave radiation [W/m$^2$]')
  ax1.grid()

  # LONGWAVE
  # figure()
  ax2.plot(np.arange(timesteps) / 8., d['lwdown_' + ee + '_sorted'], 'g-')
  ax2.set_title('Longwave radiation ' + ee)
  ax2.set_ylabel('Longwave radiation [W/m$^2$]')
  ax2.grid()

  # PRECIPITATION
  # figure()
  ax3.plot(np.arange(timesteps) / 8., d['precip_' + ee + '_sorted'], 'b-')
  ax3.set_title('Precipitation ' + ee)
  ax3.set_ylabel('Precipitation [m]')
  ax3.grid()

  # TEMPERATURE
  # figure()
  ax4.plot(np.arange(timesteps) / 8., d['temp_' + ee + '_sorted'], 'm-')
  ax4.set_title('2m temperature ' + ee)
  ax4.set_ylabel('Temperature [K]')
  ax4.grid()

  # create different subplots for the extra snow variables
  f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True)
  plt.title(ee)

  # SNOWFALL
  # figure()
  ax1.plot(np.arange(timesteps) / 8., d['snowfall_' + ee + '_sorted'], 'r-')
  ax1.set_title('Snowfall ' + ee)
  ax1.set_ylabel('Snowfall [m]')
  ax1.grid()

  # DEW POINT (to compute relative humidity)
  # figure()
  ax2.plot(np.arange(timesteps) / 8., d['dewtemp_' + ee + '_sorted'], 'g-')
  ax2.set_title('Dew-point temperature ' + ee)
  ax2.set_ylabel('Dew-point temperature [K]')
  ax2.grid()

  # U-WIND - remember to look at absolute value!
  # figure()
  ax3.plot(np.arange(timesteps) / 8., d['u10_' + ee + '_sorted'], 'b-')
  ax3.set_title('u-Wind ' + ee)
  ax3.set_ylabel('u-Wind [m/s]')
  ax3.grid()

  # V-WIND - remember to look at absolute value!
  # figure()
  ax4.plot(np.arange(timesteps) / 8., d['v10_' + ee + '_sorted'], 'm-')
  ax4.set_title('v-Wind ' + ee)
  ax4.set_ylabel('v-Wind [m/s]')
  ax4.grid()

  ####################### WRITE THIS IN A TEXT FILE FOR SAMSIM

  fileinput = pd.DataFrame()
  fileinput['Precipitation'] = d['precip_' + ee + '_sorted']
  fileinput['Longwave radiation'] = d['lwdown_' + ee + '_sorted']
  fileinput['Shortwave radiation'] = d['swdown_' + ee + '_sorted']
  fileinput['2m Temperature'] = d['temp_' + ee + '_sorted'] - 273.15

  #### 
  fileinput['Snowfall'] = d['snowfall_' + ee + '_sorted']
  fileinput['Dew-point temperature'] = d['dewtemp_' + ee + '_sorted']
  fileinput['u-wind'] = d['u10_' + ee + '_sorted']
  fileinput['v-wind'] = d['v10_' + ee + '_sorted']
  ####

  ff = outputpath + '/Clara/' + ee + '/flux_lw.txt.input'
  np.savetxt(ff, fileinput['Longwave radiation'].values, fmt='%f')

  ff = outputpath + '/Clara/' + ee + '/flux_sw.txt.input'
  np.savetxt(ff, fileinput['Shortwave radiation'].values, fmt='%f')

  ff = outputpath + '/Clara/' + ee + '/precip.txt.input'
  np.savetxt(ff, fileinput['Precipitation'].values, fmt='%10.6e')

  ff = outputpath+'/Clara/'+ ee+'/T2m.txt.input'
  np.savetxt(ff, fileinput['2m Temperature'].values, fmt='%f')
##################################

