#!/usr/bin/env python
# coding: utf-8


"""
Created on Wed Nov 06 15:27 2019

final figures paper 1 before submitting

@author: Clara Burgard
"""

import numpy as np
import xarray as xr
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.stats import kde
import seaborn as sns
import pandas as pd

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('./data/')

import memls_functions as memls
import simplification_functions as sf
import analysis_functions as ana



sns.set_context('paper')

####################################################
### READ IN THE DATA FIRST-YEAR ICE
####################################################


exp_nb='049' # this was the number in my experiments, you can change it if needed

#went once through all experiments
ee='75N00W-p4'
#ee='82N50E-p4'
#ee='85N50W-p4'
#ee='82N120W-p4'
#ee='80N160W-p4'
#ee='77N39E-p4'
ee2=ee.split("-")[0]

home_path='/home/mpim/m300411'

subfolder=glob.glob('/home/mpim/m300411/SatSim/git_scripts/Documentation/Experiments/README_exp'+exp_nb+'*.txt')
datet = subfolder[0].split("_")[4].split(".")[0]
inputpath2 = '/work/mh0033/m300411/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/INPUT/netcdf_files/'
inputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/OUTPUT/netcdf_files/'
outputpath = '/work/mh0033/m300411/SatSim//MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/RESULTS/'

#input interpolated on equidistant layers - good to compare by depth
memls_input_FYI = xr.open_dataset(inputpath2+'inputMEMLS_depthidx_'+ee2+'.nc')#.resample(time='1D').mean()
#input interpolated on equidistant layers - good to compare by depth - reference
memls_input_comp_FYI = xr.open_dataset(inputpath2+'inputMEMLS_depthidx_complex_'+ee2+'.nc')
#input on original layers (between 0 and 100) 
memls_input_FYI_orig = xr.open_dataset(inputpath2+'inputMEMLS_'+ee2+'.nc')#.resample(time='1D').mean()
#brightness temperatures assuming there is no snow layer
memls_output_FYI_snowno = xr.open_dataset(inputpath+'outputMEMLS_'+ee2+'_snowno.nc')#.resample(time='1D').mean()
#brightness temperatures assuming there is a snow layer
memls_output_FYI_snowyes = xr.open_dataset(inputpath+'outputMEMLS_'+ee2+'_snowyes.nc')#.resample(time='1D').mean()

timelength_FYI = memls_input_FYI['time']
ts_month_FYI = timelength_FYI['time.month']

inputpath_samsim = '/work/mh0033/m300411/SatSim/SAMSIM/'

#timeseries of snow thickness
thick_sn_file_FYI = np.loadtxt(inputpath_samsim+'75N00W-p4/dat_snow.dat') #thick_snow,T_snow,psi_l_snow,psi_s_snow
thick_sn_FYI = thick_sn_file_FYI[:,0]
thick_sn_FYI[thick_sn_FYI==0] = np.nan

#timeseries of ice thickness
thick_FYI = np.loadtxt(inputpath_samsim+"75N00W-p4/dat_thick.dat")
thick_ice_FYI = np.sum(thick_FYI,axis=1) #still contains timesteps where there is actually no ice, needs a proper mask

#add these timeseries to the input arrays
memls_input_FYI['thick_snow'] = xr.DataArray(thick_sn_FYI, coords=[timelength_FYI], dims=['time'])#.resample(time='1D').mean()
memls_input_FYI['thick_ice'] = xr.DataArray(thick_ice_FYI, coords=[timelength_FYI], dims=['time'])#.resample(time='1D').mean()

####################################################
### READ IN THE MULTIYEAR ICE DATA
####################################################


exp_nb='050'

#went once through all experiments
ee='NorthPole-p4'
ee2=ee.split("-")[0]

home_path='/home/mpim/m300411'

subfolder=glob.glob('/home/mpim/m300411/SatSim/git_scripts/Documentation/Experiments/README_exp'+exp_nb+'*.txt')
datet = subfolder[0].split("_")[4].split(".")[0]
inputpath2 = '/work/mh0033/m300411/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/INPUT/netcdf_files/'
inputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/OUTPUT/netcdf_files/'
outputpath = '/work/mh0033/m300411/SatSim//MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/RESULTS/'

#input interpolated on equidistant layers - good to compare by depth
memls_input_MYI = xr.open_dataset(inputpath2+'inputMEMLS_depthidx_'+ee2+'.nc')#.resample(time='1D').mean()
#input on original layers (between 0 and 100) 
memls_input_MYI_orig = xr.open_dataset(inputpath2+'inputMEMLS_'+ee2+'.nc')#.resample(time='1D').mean()
#brightness temperatures assuming there is no snow layer
memls_output_MYI_snowno = xr.open_dataset(inputpath+'outputMEMLS_'+ee2+'_snowno.nc')#.resample(time='1D').mean()
#brightness temperatures assuming there is a snow layer
memls_output_MYI_snowyes = xr.open_dataset(inputpath+'outputMEMLS_'+ee2+'_snowyes.nc')#.resample(time='1D').mean()

#timeseries of snow thickness
thick_sn_file_MYI = np.loadtxt(inputpath_samsim+'NorthPole-p4/dat_snow.dat') #thick_snow,T_snow,psi_l_snow,psi_s_snow
thick_sn_MYI = thick_sn_file_MYI[:,0]
thick_sn_MYI[thick_sn_MYI==0] = np.nan

#timeseries of ice thickness
thick_MYI = np.loadtxt(inputpath_samsim+"NorthPole-p4/dat_thick.dat")
thick_ice_MYI = np.sum(thick_MYI,axis=1) #still contains timesteps where there is actually no ice, needs a proper mask

#add these timeseries to the input arrays
memls_input_MYI['thick_snow'] = xr.DataArray(thick_sn_MYI, coords=[timelength_FYI], dims=['time'])#.resample(time='1D').mean()
memls_input_MYI['thick_ice'] = xr.DataArray(thick_ice_MYI, coords=[timelength_FYI], dims=['time'])#.resample(time='1D').mean()

timelength = memls_input_MYI['time']

#this experiment contains one year of first-year ice, then multiyear ice, so we filter out the FYI
timelength_MYI = timelength[memls_input_MYI['type'].sel(sens_exp='complex').median('depth')==4.].drop('sens_exp')
memls_output_MYI_snowno = memls_output_MYI_snowno.sel(time=timelength_MYI)
memls_output_MYI_snowyes = memls_output_MYI_snowyes.sel(time=timelength_MYI)
memls_input_MYI = memls_input_MYI.sel(time=timelength_MYI)
memls_input_MYI_orig = memls_input_MYI_orig.sel(time=timelength_MYI)
ts_month_MYI = timelength_MYI['time.month']

####################################################
### COMPUTE THE SURFACE LIQUID WATER FRACTION
####################################################


surf_lwf_ice0 = xr.DataArray(np.zeros((len(memls_input_FYI['time']),len(memls_input_FYI['sens_exp']))), coords=[
                                ('time',memls_input_FYI.time.values),
                                ('sens_exp',memls_input_FYI.sens_exp.values)])
surf_lwf_ice_FYI = surf_lwf_ice0*np.nan
surf_lwf_ice_MYI = surf_lwf_ice0*np.nan

#write out the brine volume fraction at the surface
for sens_exp in memls_input_FYI.sens_exp:
    print(sens_exp.values)
    #FYI
    for tt,timet in enumerate(memls_input_FYI['time']):
        #take part of the array without snow to get the top ice layer
        temporary = memls_input_FYI_orig['liquid_water_fraction'].sel(time=timet,sens_exp=sens_exp).where(memls_input_FYI_orig['snow_ice'].sel(time=timet,sens_exp=sens_exp)>0).dropna('top_100_to_bottom_0')
        if len(temporary)>0:
            surf_lwf_ice_FYI.loc[dict(time=timet,sens_exp=sens_exp)] = temporary.values[0]
        #MYI
        if timet.values in memls_input_MYI['time']:
            #take part of the array without snow to get the top ice layer
            temporary2 = memls_input_MYI_orig['liquid_water_fraction'].sel(time=timet,sens_exp=sens_exp).where(memls_input_MYI_orig['snow_ice'].sel(time=timet,sens_exp=sens_exp)>0).dropna('top_100_to_bottom_0')
            if len(temporary2)>0:
                surf_lwf_ice_MYI.loc[dict(time=timet,sens_exp=sens_exp)] = temporary2.values[0]

####################################################
#### PREPARE FIGURE SPECIFICS
####################################################

outputpath_fig = '/work/mh0033/m300411/SatSim/Paper_plots/'

moncol = ['mediumpurple','darkorchid','palevioletred','sandybrown','gold','darkorange','maroon','red','olive','deepskyblue','navy','blue']

#params = {'legend.fontsize': 'xx-large',
#          'figure.figsize': (7,7),
#         'axes.labelsize': 'xx-large',
#         'axes.titlesize':'xx-large',
#         'xtick.labelsize':'xx-large',
#         'ytick.labelsize':'xx-large'}
#mpl.rcParams.update(params)

alph = 0.2

siz = 10

xx = range(60,300) 


####################################################
### PREPARE FIGURE 3: SURFACE LIQUID WATER FRACTION AGAINST TBV
#################################################

#change these two to snowno for the right panel
TBV_comp_FYI = memls_output_FYI_snowyes['tb_v'].sel(sens_exp='complex',frequency=6.9)
TBV_simp_FYI = memls_output_FYI_snowyes['tb_v'].sel(sens_exp='simpleallfunc',frequency=6.9)

surf_lwf_FYI = surf_lwf_ice_FYI

y_actual_FYI = np.append(surf_lwf_ice_FYI.sel(sens_exp='complex'),surf_lwf_ice_FYI.sel(sens_exp='simpleallfunc'))
y_predicted_FYI = np.append(TBV_comp_FYI,TBV_simp_FYI)
ydiff_FYI = y_predicted_FYI - y_actual_FYI

corrcoeff_FYI = np.corrcoef(y_actual_FYI[~np.isnan(y_actual_FYI) & ~np.isnan(y_predicted_FYI)],y_predicted_FYI[~np.isnan(y_actual_FYI) & ~np.isnan(y_predicted_FYI)])[0,1]

TBV_comp_MYI = memls_output_MYI_snowyes['tb_v'].sel(sens_exp='complex',frequency=6.9).sel(time=memls_input_MYI['time'])
TBV_simp_MYI = memls_output_MYI_snowyes['tb_v'].sel(sens_exp='simpleallfunc',frequency=6.9).sel(time=memls_input_MYI['time'])

surf_lwf_MYI = surf_lwf_ice_MYI.sel(time=timelength_MYI)

y_actual_MYI = np.append(surf_lwf_ice_MYI.sel(sens_exp='complex').sel(time=timelength_MYI),surf_lwf_ice_MYI.sel(sens_exp='simpleallfunc').sel(time=timelength_MYI))
y_predicted_MYI = np.append(TBV_comp_MYI,TBV_simp_MYI)
ydiff_MYI = y_predicted_MYI - y_actual_MYI

corrcoeff_MYI = np.corrcoef(y_actual_MYI[~np.isnan(y_actual_MYI) & ~np.isnan(y_predicted_MYI)],y_predicted_MYI[~np.isnan(y_actual_MYI) & ~np.isnan(y_predicted_MYI)])[0,1]

####################################################
# PLOT FIGURE 3: use snowyes in the two first lines above to produce the left panel and snowno 
# to produce the right panel
####################################################

siz = 5
#mpl.rcParams.update({'font.size': 14})
#mpl.rcParams.update({'axes.labelsize': 16})

#try to plot this against the brightness temperature including snow in the calculations
f = plt.figure(dpi=200)
f.set_size_inches(8.27/2,8.27/2)

plt.scatter(surf_lwf_FYI.sel(sens_exp='complex',time=ana.is_summer(ts_month_FYI)),TBV_comp_FYI.sel(time=ana.is_summer(ts_month_FYI)),c='r',marker='o',alpha=alph,s=siz,edgecolors='None')
plt.scatter(surf_lwf_FYI.sel(sens_exp='complex',time=ana.is_winter(ts_month_FYI)),TBV_comp_FYI.sel(time=ana.is_winter(ts_month_FYI)),c='b',marker='o',alpha=alph,s=siz,edgecolors='None')
plt.scatter(surf_lwf_MYI.sel(sens_exp='complex',time=ana.is_summer(ts_month_MYI)),TBV_comp_MYI.sel(time=ana.is_summer(ts_month_MYI)),marker='P',c='r',alpha=alph,s=siz,edgecolors='None')
plt.scatter(surf_lwf_MYI.sel(sens_exp='complex',time=ana.is_winter(ts_month_MYI)),TBV_comp_MYI.sel(time=ana.is_winter(ts_month_MYI)),marker='P',c='b',alpha=alph,s=siz,edgecolors='None')
plt.xlabel('Surface liquid water fraction')
plt.ylabel('Brightness temperature V-pol 6.9 GHz [K]')
plt.text(0.6,270,'FYI r = '+str(np.round(corrcoeff_FYI,2)))
plt.text(0.6,265,'MYI r = '+str(np.round(corrcoeff_MYI,2)))
#sns.despine(offset=10)
sns.despine()
plt.tight_layout()
#f.savefig(outputpath_fig+'liquidwater_TBV_withsnow_'+ee2+'_Nov2019.pdf',rasterize=True,bbox_inches='tight')
f.savefig(outputpath_fig+'liquidwater_TBV_withoutsnow_'+ee2+'_Nov2019.pdf',rasterize=True,bbox_inches='tight')



####################################################
## SENSITIVITY STUDIES FIGURE SUMMER MYI and FYI V
## FIGURE 6
####################################################

siz = 5

corrcoeff = np.zeros((5,2))
mean_diff = np.zeros((5,2))
std_diff = np.zeros((5,2))

memls_output_FYI = memls_output_FYI_snowno.copy()
memls_output_MYI = memls_output_MYI_snowno.copy()

for i,sity in enumerate(['FYI','MYI']):
    if sity == 'FYI':
        y_actual = memls_output_FYI['tb_v'].sel(sens_exp='complex',frequency=6.9,time=ana.is_summer(ts_month_FYI))
    elif sity == 'MYI':
         y_actual = memls_output_MYI['tb_v'].sel(sens_exp='complex',frequency=6.9,time=ana.is_summer(ts_month_MYI))
        
    for j,exp in enumerate(['simpleallconst','simpletemp','simplesalconst','simpleallfunc','simplesalfunc']):
        if sity == 'FYI':
            y_predicted = memls_output_FYI['tb_v'].sel(sens_exp=exp,frequency=6.9,time=ana.is_summer(ts_month_FYI))
        elif sity == 'MYI':
             y_predicted = memls_output_MYI['tb_v'].sel(sens_exp=exp,frequency=6.9,time=ana.is_summer(ts_month_MYI))

        ydiff = abs(y_predicted - y_actual)
#        mean_diff[i,j] = np.mean(abs(ydiff))
#        std_diff[i,j] = np.std(abs(ydiff))
        mean_diff[j,i] = np.mean(ydiff)
        std_diff[j,i] = np.std(ydiff)
        corrcoeff[j,i] = np.corrcoef(y_actual[~np.isnan(y_actual) & ~np.isnan(y_predicted)],y_predicted[~np.isnan(y_actual) & ~np.isnan(y_predicted)])[0,1]

# mpl.rcParams.update({'font.size': 30})
# mpl.rcParams.update({'axes.labelsize': 14})

f, axs = plt.subplots(5, 2, sharex=True, sharey=True)

#define format     
f.set_size_inches(5.5,10)
#f.subplots_adjust(bottom=0.4)
#f.subplots_adjust(left=0.4)
# f.subplots_adjust(right=0.95)
# f.subplots_adjust(top=0.92)
# f.subplots_adjust(hspace=0.07, wspace=0.05)

y1 = memls_output_FYI['tb_v'].sel(sens_exp='complex',frequency=6.9)
y2 = memls_output_FYI['tb_v'].sel(sens_exp='simpleallfunc',frequency=6.9)
y3 = memls_output_FYI['tb_v'].sel(sens_exp='simpleallconst',frequency=6.9)
y4 = memls_output_FYI['tb_v'].sel(sens_exp='simpletemp',frequency=6.9)
y5 = memls_output_FYI['tb_v'].sel(sens_exp='simplesalfunc',frequency=6.9)
y6 = memls_output_FYI['tb_v'].sel(sens_exp='simplesalconst',frequency=6.9)

y1b = memls_output_MYI['tb_v'].sel(sens_exp='complex',frequency=6.9)
y2b = memls_output_MYI['tb_v'].sel(sens_exp='simpleallfunc',frequency=6.9)
y3b = memls_output_MYI['tb_v'].sel(sens_exp='simpleallconst',frequency=6.9)
y4b = memls_output_MYI['tb_v'].sel(sens_exp='simpletemp',frequency=6.9)
y5b = memls_output_MYI['tb_v'].sel(sens_exp='simplesalfunc',frequency=6.9)
y6b = memls_output_MYI['tb_v'].sel(sens_exp='simplesalconst',frequency=6.9)


axs[0,0].set_xlim(160,280)
axs[0,0].set_ylim(160,280)

axs[3, 0].scatter(y1.sel(time=ana.is_summer(ts_month_FYI)),y2.sel(time=ana.is_summer(ts_month_FYI)), marker='.', c='r', edgecolors='None', alpha=alph, s=siz)
axs[0, 0].scatter(y1.sel(time=ana.is_summer(ts_month_FYI)),y3.sel(time=ana.is_summer(ts_month_FYI)), marker='.', c='r', edgecolors='None', alpha=alph, s=siz)
axs[1, 0].scatter(y1.sel(time=ana.is_summer(ts_month_FYI)),y4.sel(time=ana.is_summer(ts_month_FYI)), marker='.', c='r', edgecolors='None',alpha=alph, s=siz)
axs[4, 0].scatter(y1.sel(time=ana.is_summer(ts_month_FYI)),y5.sel(time=ana.is_summer(ts_month_FYI)), marker='.', c='r', edgecolors='None',alpha=alph, s=siz)
axs[2, 0].scatter(y1.sel(time=ana.is_summer(ts_month_FYI)),y6.sel(time=ana.is_summer(ts_month_FYI)), marker='.', c='r', edgecolors='None',alpha=alph, s=siz)
  
axs[3, 1].scatter(y1b.sel(time=ana.is_summer(ts_month_MYI)),y2b.sel(time=ana.is_summer(ts_month_MYI)), marker='P', c='r', edgecolors='None',alpha=alph, s=siz)
axs[0, 1].scatter(y1b.sel(time=ana.is_summer(ts_month_MYI)),y3b.sel(time=ana.is_summer(ts_month_MYI)), marker='P', c='r', edgecolors='None',alpha=alph, s=siz)
axs[1, 1].scatter(y1b.sel(time=ana.is_summer(ts_month_MYI)),y4b.sel(time=ana.is_summer(ts_month_MYI)), marker='P', c='r', edgecolors='None',alpha=alph, s=siz)
axs[4, 1].scatter(y1b.sel(time=ana.is_summer(ts_month_MYI)),y5b.sel(time=ana.is_summer(ts_month_MYI)), marker='P', c='r', edgecolors='None',alpha=alph, s=siz)
axs[2, 1].scatter(y1b.sel(time=ana.is_summer(ts_month_MYI)),y6b.sel(time=ana.is_summer(ts_month_MYI)), marker='P', c='r', edgecolors='None',alpha=alph, s=siz)

     
for i in range(5):
    for j in range(2):
        axs[i,j].plot(xx,xx,'k--')
        #axs[i, j].text(162,272,'r = '+str(np.round(corrcoeff[i,j],2)))
        axs[i, j].text(162,272,'Mean abs diff')
        axs[i, j].text(162,260,str(round(mean_diff[i,j],1))+'$\pm$'+str(round(std_diff[i,j],1))+' K')
        #axs[i, j].text(162,265,'Mean abs $\Delta$ = '+str(round(mean_diff[i,j],2))+'$\pm$'+str(round(std_diff[i,j],2))+' K')
        #axs[i,j].tick_params(axis='both', which='major')#, labelsize=30)
   
     
plt.setp([a.get_xticklabels() for a in axs[0, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
plt.setp([a.get_xticklabels() for a in axs[:, 0]])
plt.setp([a.get_xticklabels() for a in axs[:, 1]])

f.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel(' "Real" brightness temperatures [K]')
plt.ylabel('Brightness temperatures simulated from simplified profiles [K]')

#f.text(0.5,0.01,' "Real" brightness temperatures [K]' ,ha='center')
#f.text(0.03,0.7, 'Brightness temperatures simulated from simplified profiles [K]' ,ha='center',rotation='vertical')
#f.savefig(outputpath+'sensitivitystudies_V_FYI.pdf',bbox_inches='tight',orientation='landscape')
sns.despine()
plt.tight_layout()

#f.savefig(outputpath_fig+'sensitivity_TBV_summer_wosnow.png')
f.savefig(outputpath_fig+'sensitivity_TBV_summer_withoutsnow'+ee2+'_Nov2019.pdf',rasterize=True,bbox_inches='tight')


####################################################
### SENSITIVITY STUDIES FIGURE WINTER MYI and FYI V
## FIGURE 5
####################################################

siz = 5

corrcoeff = np.zeros((5,2))
mean_diff = np.zeros((5,2))
std_diff = np.zeros((5,2))

memls_output_FYI = memls_output_FYI_snowno.copy()
memls_output_MYI = memls_output_MYI_snowno.copy()


for i,sity in enumerate(['FYI','MYI']):
    if sity == 'FYI':
        y_actual = memls_output_FYI['tb_v'].sel(sens_exp='complex',frequency=6.9,time=ana.is_winter(ts_month_FYI))
    elif sity == 'MYI':
        y_actual = memls_output_MYI['tb_v'].sel(sens_exp='complex',frequency=6.9,time=ana.is_winter(ts_month_MYI))
        
    for j,exp in enumerate(['simpleallconst','simpletemp','simplesalconst','simpleallfunc','simplesalfunc']):
        if sity == 'FYI':
            y_predicted = memls_output_FYI['tb_v'].sel(sens_exp=exp,frequency=6.9,time=ana.is_winter(ts_month_FYI))
        elif sity == 'MYI':
            y_predicted = memls_output_MYI['tb_v'].sel(sens_exp=exp,frequency=6.9,time=ana.is_winter(ts_month_MYI))

        ydiff = abs(y_predicted - y_actual)
#        mean_diff[i,j] = np.mean(abs(ydiff))
#        std_diff[i,j] = np.std(abs(ydiff))
        mean_diff[j,i] = np.mean(ydiff)
        std_diff[j,i] = np.std(ydiff)
        corrcoeff[j,i] = np.corrcoef(y_actual[~np.isnan(y_actual) & ~np.isnan(y_predicted)],y_predicted[~np.isnan(y_actual) & ~np.isnan(y_predicted)])[0,1]

# mpl.rcParams.update({'font.size': 30})
# mpl.rcParams.update({'axes.labelsize': 14})

f, axs = plt.subplots(5, 2, sharex=True, sharey=True)

#define format     
f.set_size_inches(5.5,10)
#f.subplots_adjust(bottom=0.4)
#f.subplots_adjust(left=0.4)
# f.subplots_adjust(right=0.95)
# f.subplots_adjust(top=0.92)
# f.subplots_adjust(hspace=0.07, wspace=0.05)

y1 = memls_output_FYI['tb_v'].sel(sens_exp='complex',frequency=6.9)
y2 = memls_output_FYI['tb_v'].sel(sens_exp='simpleallfunc',frequency=6.9)
y3 = memls_output_FYI['tb_v'].sel(sens_exp='simpleallconst',frequency=6.9)
y4 = memls_output_FYI['tb_v'].sel(sens_exp='simpletemp',frequency=6.9)
y5 = memls_output_FYI['tb_v'].sel(sens_exp='simplesalfunc',frequency=6.9)
y6 = memls_output_FYI['tb_v'].sel(sens_exp='simplesalconst',frequency=6.9)

y1b = memls_output_MYI['tb_v'].sel(sens_exp='complex',frequency=6.9)
y2b = memls_output_MYI['tb_v'].sel(sens_exp='simpleallfunc',frequency=6.9)
y3b = memls_output_MYI['tb_v'].sel(sens_exp='simpleallconst',frequency=6.9)
y4b = memls_output_MYI['tb_v'].sel(sens_exp='simpletemp',frequency=6.9)
y5b = memls_output_MYI['tb_v'].sel(sens_exp='simplesalfunc',frequency=6.9)
y6b = memls_output_MYI['tb_v'].sel(sens_exp='simplesalconst',frequency=6.9)


#axs[0,0].set_xlim(240,270)
#axs[0,0].set_ylim(240,270)

axs[3, 0].scatter(y1.sel(time=ana.is_winter(ts_month_FYI)),y2.sel(time=ana.is_winter(ts_month_FYI)), marker='.', c='b', edgecolors='None', alpha=alph, s=siz)
axs[0, 0].scatter(y1.sel(time=ana.is_winter(ts_month_FYI)),y3.sel(time=ana.is_winter(ts_month_FYI)), marker='.', c='b', edgecolors='None', alpha=alph, s=siz)
axs[1, 0].scatter(y1.sel(time=ana.is_winter(ts_month_FYI)),y4.sel(time=ana.is_winter(ts_month_FYI)), marker='.', c='b', edgecolors='None',alpha=alph, s=siz)
axs[4, 0].scatter(y1.sel(time=ana.is_winter(ts_month_FYI)),y5.sel(time=ana.is_winter(ts_month_FYI)), marker='.', c='b', edgecolors='None',alpha=alph, s=siz)
axs[2, 0].scatter(y1.sel(time=ana.is_winter(ts_month_FYI)),y6.sel(time=ana.is_winter(ts_month_FYI)), marker='.', c='b', edgecolors='None',alpha=alph, s=siz)
  
axs[3, 1].scatter(y1b.sel(time=ana.is_winter(ts_month_MYI)),y2b.sel(time=ana.is_winter(ts_month_MYI)), marker='P', c='b', edgecolors='None',alpha=alph, s=siz)
axs[0, 1].scatter(y1b.sel(time=ana.is_winter(ts_month_MYI)),y3b.sel(time=ana.is_winter(ts_month_MYI)), marker='P', c='b', edgecolors='None',alpha=alph, s=siz)
axs[1, 1].scatter(y1b.sel(time=ana.is_winter(ts_month_MYI)),y4b.sel(time=ana.is_winter(ts_month_MYI)), marker='P', c='b', edgecolors='None',alpha=alph, s=siz)
axs[4, 1].scatter(y1b.sel(time=ana.is_winter(ts_month_MYI)),y5b.sel(time=ana.is_winter(ts_month_MYI)), marker='P', c='b', edgecolors='None',alpha=alph, s=siz)
axs[2, 1].scatter(y1b.sel(time=ana.is_winter(ts_month_MYI)),y6b.sel(time=ana.is_winter(ts_month_MYI)), marker='P', c='b', edgecolors='None',alpha=alph, s=siz)

     
for i in range(5):
    for j in range(2):
        axs[i,j].plot(xx,xx,'k--')
        #axs[i, j].text(162,272,'r = '+str(np.round(corrcoeff[i,j],2)))
        axs[i, j].text(242,268,'Mean abs diff')
        axs[i, j].text(242,264,str(round(mean_diff[i,j],1))+'$\pm$'+str(round(std_diff[i,j],1))+' K')
        #axs[i, j].text(162,265,'Mean abs $\Delta$ = '+str(round(mean_diff[i,j],2))+'$\pm$'+str(round(std_diff[i,j],2))+' K')
        #axs[i,j].tick_params(axis='both', which='major')#, labelsize=30)
   
     
plt.setp([a.get_xticklabels() for a in axs[0, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
plt.setp([a.get_xticklabels() for a in axs[:, 0]])
plt.setp([a.get_xticklabels() for a in axs[:, 1]])

f.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel(' "Real" brightness temperatures [K]')
plt.ylabel('Brightness temperatures simulated from simplified profiles [K]')

#f.text(0.5,0.01,' "Real" brightness temperatures [K]' ,ha='center')
#f.text(0.03,0.7, 'Brightness temperatures simulated from simplified profiles [K]' ,ha='center',rotation='vertical')
#f.savefig(outputpath+'sensitivitystudies_V_FYI.pdf',bbox_inches='tight',orientation='landscape')
sns.despine()
plt.tight_layout()

#f.savefig(outputpath_fig+'sensitivity_TBV_summer_wosnow.png')
f.savefig(outputpath_fig+'sensitivity_TBV_winter_withoutsnow'+ee2+'_Nov2019.pdf',rasterize=True,bbox_inches='tight')




