#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 14:48:16 2018

plot the salinity profiles functions

@author: claraburgard
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.append('../scripts_simulation/')
import simplification_functions as sf

sns.set_context('paper')

# params = {'legend.fontsize': 'x-large',
#           'figure.figsize': (10,7),
#          'axes.labelsize': 'xx-large',
#          'axes.titlesize':'xx-large',
#          'xtick.labelsize':'xx-large',
#          'ytick.labelsize':'xx-large'}
# mpl.rcParams.update(params)

outputpath = '/your/path/to/PLOTS'

depth = np.arange(1.1,0,-0.01)

plt.figure(figsize=(8.27/1.5,8.27/2))
plt.axvline(x = 5.,c='orange',linestyle='--',linewidth=2.0,label='First-year ice constant salinity')
plt.axvline(x = 1.,c='deepskyblue',linestyle='--',linewidth=2.0,label='Multiyear ice constant salinity')
plt.plot(sf.sal_approx_fy(depth),-depth,c='orange',linestyle='-',linewidth=2.0,label='First-year ice salinity function of depth')
plt.plot(sf.sal_approx_my(depth),-depth,c='deepskyblue',linestyle='-',linewidth=2.0,label='Multiyear ice salinty function of depth')
plt.ylim(-1,0)
plt.xlim(0,25)
plt.legend()
plt.yticks(np.arange(-1,0.19,0.2),np.round(np.arange(1.001,0,-0.2),2))
plt.xlabel('Salinity [g/kg]')
plt.ylabel('Normalized depth')
sns.despine()
plt.tight_layout()
plt.savefig(outputpath+'Figure5.pdf',rasterize=True,bbox_inches='tight')