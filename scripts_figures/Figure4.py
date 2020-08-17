#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:48:16 2018

plot with salinity profiles

@author: claraburgard
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import clalib.simplification_functions as sf
import seaborn as sns

sns.set_context('paper')

# params = {'legend.fontsize': 'x-large',
#           'figure.figsize': (10,7),
#          'axes.labelsize': 'xx-large',
#          'axes.titlesize':'xx-large',
#          'xtick.labelsize':'xx-large',
#          'ytick.labelsize':'xx-large'}
# mpl.rcParams.update(params)

home_path='/Users/claraburgard'

outputpath = home_path+'/mistral_work/SatSim/Paper_plots/'
#outputpath = '/work/mh0033/m300411/SatSim/Paper_plots/'

depth = np.arange(1.1,0,-0.01)

plt.figure(figsize=(8.27/1.5,8.27/2))
plt.axvline(x = 5.,c='orange',linestyle='--',label='FYI constant salinity')
plt.axvline(x = 1.,c='deepskyblue',linestyle='--',label='MYI constant salinity')
plt.plot(sf.sal_approx_fy(depth),-depth,c='orange',linestyle='-',label='FYI salinity function of depth')
plt.plot(sf.sal_approx_my(depth),-depth,c='deepskyblue',linestyle='-',label='MYI salinty function of depth')
plt.ylim(-1,0)
plt.xlim(0,25)
plt.legend()
plt.yticks(np.arange(-1,0.19,0.2),np.round(np.arange(1.001,0,-0.2),2))
plt.xlabel('Salinity [g/kg]')
plt.ylabel('Normalized depth')
sns.despine()
plt.tight_layout()
plt.savefig(outputpath+'salinity_profiles_new.pdf',rasterize=True,bbox_inches='tight')

### SUBPLOTS FOR SENSITIVITY STUDIES

