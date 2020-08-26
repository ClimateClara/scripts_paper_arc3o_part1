#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:05:32 2017

This script runs through the SAMSIM output data and creates simplifications from the reference profiles

updated to python3 on 11.05.2018

@author: Clara Burgard
"""

import numpy as np
import pandas as pd
import simplification_functions as sf
import memls_functions as memls


##########################################################################

ee2='75N00W'
#ee2='NorthPole'
#######
## additional experiments
#ee2='85N50W'
#ee2='82N120W'
#ee2='80N160W'
#ee2='77N39E'
#ee2='74N170E'
#######

if ee2 == '75N00W':
    folder = 'exp049_75N00W_20191009'
elif ee2 == 'NorthPole':
    folder = 'exp050_NorthPole_20191011'
#######
## additional experiments
elif ee2 == '85N50W':
    folder = 'exp045_85N50W_20190715'
elif ee2 == '82N120W':
    folder = 'exp046_82N120W_20190716'
elif ee2 == '80N160W':
    folder = 'exp047_80N160W_20190716'
elif ee2 == '77N39E':
    folder = 'exp051_77N39E_20191030'
elif ee2 == '74N170E':
    folder = 'exp052_74N170E_20191106'
#######

#Output from SAMSIM
inputpath = '/work/mh0033/m300411/SatSim/data_repo_part1/SAMSIM_input_output/SAMSIM_OUTPUT/'+ee2
#ERA-Interim input data
inputpath3 = '/work/mh0033/m300411/SatSim/data_repo_part1/SAMSIM_input_output/SAMSIM_INPUT/ERA-Interim/'
#output path to put the resulting .dat-files
outputpath = '/work/mh0033/m300411/SatSim/data_repo_part1/MEMLS_input_output/'+folder+'/INPUT/original_files/'

###################################
## choose your salinity scheme
###################################
    
#saltype = 'func' #simplified salinity as a function of depth
saltype = 'const' #simplified salinity constant

######################################################################################
## choose which experiments to run (possibilities: simpletemp, simplesal, simpledens)
######################################################################################

# default always written out are 'complex' and 'simpleall'

#experiments = ['simpletemp','simplesal']
experiments = ['simplesal']
#experiments = ['simpletemp']
#experiments=[]

###############################################
## keep amount of layers same as SAMSIM output
###############################################

layertype = 'same' #same amount of layers as reference
# layertype = 'diff' # not used in this script but rather in the run_simplifications_layers.py
ll=10

##########################
# do we run with snow on top
##########################

# Contrary to what the name suggests, it is not random, it is the snow cover from SAMSIM! :)
random_snow = 'yes'
#random_snow = 'no'

################################################################
#### READ IN THE DATA FROM SAMSIM OUTPUT
free_flag    = 0  

var1name     = 'Bulk salinity'
var1unit     = '[g/kg]'
var1         = np.loadtxt(inputpath+ee2+"/dat_S_bu.dat")

var2name     = 'Temperature'
var2unit     = '[C]'
var2         = np.loadtxt(inputpath+ee2+"/dat_T.dat")

var3name     = 'Liquid fraction'
var3unit     = '[fraction]'
var3         = np.loadtxt(inputpath+ee2+"/dat_psi_l.dat")

thick      = np.loadtxt(inputpath+ee2+"/dat_thick.dat")
freeboard  = np.loadtxt(inputpath+ee2+"/dat_freeboard.dat")

thick_sn_file = np.loadtxt(inputpath+ee2+'/dat_snow.dat') #thick_snow,T_snow,psi_l_snow,psi_s_snow
thick_sn = thick_sn_file[:,0]
T_sn = thick_sn_file[:,1]
T2m_file = np.loadtxt(inputpath+ee2+'/dat_T2m_T_top.dat') #T2m, Ttop
T2m = T2m_file[:,0]
Ttop = T2m_file[:,1]+273.15
airtemp = np.loadtxt(inputpath3+ee2+'/T2m.txt.input')
thick_ice = np.sum(thick,axis=1) #still contains timesteps where there is actually no ice, needs a proper mask, but enough for our purpose here

T_a = airtemp[:-20:4]

#compute temperature at the top of the ice from the temperature at the top of snow, assuming that the profile in the snow is linear
Ttop_new = (2*T_sn-var2[:,0])+273.15
Ttop_new[Ttop_new>273.15] = 273.15

#
#plt.figure()
#plt.plot(T_a+273.15,label='air')
#plt.plot(T_sn+273.15,label='snow')
#plt.plot(new_Ttopice,label='new top ice')
#plt.plot(var2[:,0]+273.15,label='old top ice')
#plt.plot(T2m[400:450],label='T2m')
#plt.plot(Ttop,label='T_top')
#plt.plot(Ttop_new,label='T_top computed')
#plt.legend()
#
#plt.figure()
##plt.plot(T_a[400:450],label='air')
#plt.plot(T_sn+273.15,label='snow')
#plt.plot(new_Ttopice,label='new top ice')
#plt.plot(var2[:,0]+273.15,label='old top ice')
##plt.plot(T2m[400:450],label='T2m')
#plt.plot(Ttop,label='T_top')
#plt.legend()
#
#
#plt.figure()
#plt.plot(new_Ttopice-273.15-var2[:,0])

timesteps = len(thick)

#Setting freeboard to zero if free_flag = 0
if   free_flag == 0:
  freeboard[:] = 0.

ylen = len(thick[0,:])
xlen = len(thick[:,0])


#Restructuring the data so it represents the finite volume nature of SAMSIM
depth_step=np.hstack((thick,thick))
var1_step =np.hstack((thick,thick))
var2_step =np.hstack((thick,thick))
var3_step =np.hstack((thick,thick))

ylen = len(thick[0,:])
xlen = len(thick[:,0])
i=0
j=0
while (i<xlen):
  while (j<ylen):
    depth_step[i,2*j]     = -sum(thick[i,0:j])+freeboard[i]
    depth_step[i,2*j+1]   = -sum(thick[i,0:j])-thick[i,j]+freeboard[i]
    var1_step[i,2*j]      = var1[i,j]
    var1_step[i,2*j+1]    = var1[i,j]
    var2_step[i,2*j]      = var2[i,j]
    var2_step[i,2*j+1]    = var2[i,j]
    var3_step[i,2*j]      = var3[i,j]
    var3_step[i,2*j+1]    = var3[i,j]
    j=j+1
  i=i+1
  j=0
########################################################################
      
    
###############################################
###### CONSTRUCT WHAT IS NEEDED FOR MEMLS #####
###############################################
  
###########################################
######### WRITING THE FILES
############################################
      
timeseries_depth = []
timeseries_depth_ip = []
icetype_series = []    
#icetype_series2 = []     


#loop over each timestep
for tt in range(timesteps):
  
########## COMPLEX (=REFERENCE)
#  if 'complex' in experiments:
  print(tt, 'complex')
  [T_i, num, W_i, rho_i, pc_i, d_i, sal_i, sitype, s_i, icetype, psi_l] = sf.orig_complex_iceonly(tt,var2_step, depth_step, var1_step, var3_step,timeseries_depth)
  icetype_series.append(icetype) 

  #create dataframe to write into the file
  fileinput = pd.DataFrame()

  if random_snow=='yes' and thick_sn[tt]>0.:
      fileinput['Layer number'] = np.append(num,num[-1]+1)
      fileinput['Layer temperature'] = np.append(T_i,T_sn[tt]+273.15) #(T_i[-1]+airtemp[tt]+273.15)/2
      fileinput['Layer wetness'] = np.append(W_i,0)
      fileinput['Layer density'] = np.append(rho_i,330.)
      fileinput['Layer correlation length'] = np.append(pc_i,0.15)
      fileinput['Layer exp. correlation length'] = np.append(pc_i,0.15)
      fileinput['Layer thickness'] = np.append(d_i,thick_sn[tt])
      fileinput['Layer salinity'] = np.append(sal_i,0.)
      fileinput['Layer type'] = np.append(sitype,2)
      fileinput['Snow/ice'] = np.append(s_i,0.)
      fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

  else:
      fileinput['Layer number'] = num
      fileinput['Layer temperature'] = T_i
      fileinput['Layer wetness'] = W_i
      fileinput['Layer density'] = rho_i
      fileinput['Layer correlation length'] = pc_i
      fileinput['Layer exp. correlation length'] = pc_i
      fileinput['Layer thickness'] = d_i
      fileinput['Layer salinity'] = sal_i
      fileinput['Layer type'] = sitype
      fileinput['Snow/ice'] = s_i
      #fileinput['Liquid_water_frac'] = psi_l
      fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])
  

    
  ff = outputpath+'/timestep'+str(tt).zfill(4)+'_complex.dat'
  np.savetxt(ff, fileinput.values, fmt='%f')
##############

######### SIMPLE ALL
#  if 'simpleall' in experiments:
  print(tt, 'simple all')

#do not forget to put double_linear to yes or no
  if layertype == 'same':
      [nn, d_i2, T_i2, W_i2, pc_i2, sal_i2, rho_i2, s_i2, sitype_2, icetype_2] = sf.samelayers(T_i,d_i,W_i,pc_i,sal_i,saltype,rho_i,s_i,icetype,sitype,thick_sn[tt],thick_ice[tt],Ttop_new[tt],'yes')
  elif layertype == 'diff':
      [nn, d_i2, T_i2, W_i2, pc_i2, sal_i2, rho_i2, s_i2, sitype_2, icetype_2] = sf.difflayers(ll,T_i,d_i,W_i,pc_i,sal_i,saltype,rho_i,s_i,icetype,sitype,thick_sn[tt],thick_ice[tt],Ttop_new[tt],'yes')
   
  #create dataframe to write into the file
  fileinput = pd.DataFrame()

  if random_snow=='yes' and thick_sn[tt]>0.:
      fileinput['Layer number'] = np.arange(nn+1)+1
      fileinput['Layer temperature'] = np.append(T_i2,(T_i2[-1]+Ttop_new[tt])/2)
      fileinput['Layer wetness'] = np.append(W_i2,0)
      fileinput['Layer density'] = np.append(rho_i2,330.)
      fileinput['Layer correlation length'] = np.append(pc_i2,0.15)
      fileinput['Layer exp. correlation length'] = np.append(pc_i2,0.15)
      fileinput['Layer thickness'] = np.append(d_i2,thick_sn[tt])
      fileinput['Layer salinity'] = np.append(sal_i2,0.)
      fileinput['Layer type'] = np.append(sitype_2,2)
      fileinput['Snow/ice'] = np.append(s_i2,0.)
      fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

  else:
      fileinput['Layer number'] = np.arange(nn)+1
      fileinput['Layer temperature'] = T_i2
      fileinput['Layer wetness'] = W_i2
      fileinput['Layer density'] = rho_i2
      fileinput['Layer correlation length'] = pc_i2
      fileinput['Layer exp. correlation length'] = pc_i2
      fileinput['Layer thickness'] = d_i2
      fileinput['Layer salinity'] = sal_i2
      fileinput['Layer type'] = sitype_2
      fileinput['Snow/ice'] = s_i2
      #fileinput['Liquid_water_frac'] = psi_l
      fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])


  ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simpleall'+saltype+'.dat'
  np.savetxt(ff, fileinput.values, fmt='%f')  

########### SIMPLE TEMP
  if 'simpletemp' in experiments:
    print(tt, 'simple temp')
    
    if random_snow=='yes' and thick_sn[tt]>0.:
        
        if layertype == 'same':
            rho_i3 = sf.icerho(T_i2,sal_i)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.append(num,num[-1]+1)
            fileinput['Layer temperature'] = np.append(T_i2,(T_i2[-1]+Ttop_new[tt])/2)
            fileinput['Layer wetness'] = np.append(W_i,0)
            fileinput['Layer density'] = np.append(rho_i3,330.)
            fileinput['Layer correlation length'] = np.append(pc_i,0.15)
            fileinput['Layer exp. correlation length'] = np.append(pc_i,0.15)
            fileinput['Layer thickness'] =  np.append(d_i,thick_sn[tt])
            fileinput['Layer salinity'] = np.append(sal_i,0.)
            fileinput['Layer type'] = np.append(sitype,2)  
            fileinput['Snow/ice'] = np.append(s_i,0.)
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

            
        elif layertype == 'diff':
            sal_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),sal_i)
            rho_i3 = sf.icerho(T_i2,sal_diff)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.arange(nn+1)+1
            fileinput['Layer temperature'] = np.append(T_i2,(T_i2[-1]+Ttop_new[tt])/2)
            fileinput['Layer wetness'] = np.append(W_i2,0)
            fileinput['Layer density'] = np.append(rho_i3,330.)
            fileinput['Layer correlation length'] = np.append(pc_i2,0.15)
            fileinput['Layer exp. correlation length'] = np.append(pc_i2,0.15)
            fileinput['Layer thickness'] = np.append(d_i2,thick_sn[tt])
            fileinput['Layer salinity'] = np.append(sal_i2,0.)
            fileinput['Layer type'] = np.append(sitype_2,2)
            fileinput['Snow/ice'] = np.append(s_i2,0.)   
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])
            


    else:
        
        if layertype == 'same':
            rho_i3 = sf.icerho(T_i2,sal_i)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = num
            fileinput['Layer temperature'] = T_i2
            fileinput['Layer wetness'] = W_i
            fileinput['Layer density'] = rho_i3
            fileinput['Layer correlation length'] = pc_i
            fileinput['Layer exp. correlation length'] = pc_i
            fileinput['Layer thickness'] = d_i
            fileinput['Layer salinity'] = sal_i
            fileinput['Layer type'] = sitype  
            fileinput['Snow/ice'] = s_i 
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        elif layertype == 'diff':
            sal_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),sal_i)
            rho_i3 = sf.icerho(T_i2,sal_diff)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.arange(nn)+1
            fileinput['Layer temperature'] = T_i2
            fileinput['Layer wetness'] = W_i2
            fileinput['Layer density'] = rho_i3
            fileinput['Layer correlation length'] = pc_i2
            fileinput['Layer exp. correlation length'] = pc_i2
            fileinput['Layer thickness'] = d_i2
            fileinput['Layer salinity'] = sal_diff
            fileinput['Layer type'] = sitype_2 
            fileinput['Snow/ice'] = s_i2    
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        
    
    ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simpletemp.dat'
    np.savetxt(ff, fileinput.values, fmt='%f')   
  
########## SIMPLE SAL
  if 'simplesal' in experiments:
    print(tt, 'simple sal')
    
    if random_snow=='yes' and thick_sn[tt]>0.:
        
        if layertype == 'same':
            rho_i4 = sf.icerho(T_i,sal_i2)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.append(num,num[-1]+1)
            fileinput['Layer temperature'] = np.append(T_i,T_sn[tt]+273.15)
            fileinput['Layer wetness'] = np.append(W_i,0)
            fileinput['Layer density'] = np.append(rho_i4,330.)
            fileinput['Layer correlation length'] = np.append(pc_i,0.15)
            fileinput['Layer exp. correlation length'] = np.append(pc_i,0.15)
            fileinput['Layer thickness'] = np.append(d_i,thick_sn[tt])
            fileinput['Layer salinity'] = np.append(sal_i2,0.)
            fileinput['Layer type'] = np.append(sitype,2)   
            fileinput['Snow/ice'] = np.append(s_i,0.)
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        elif layertype == 'diff':
            T_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),T_i)
            rho_i4 = sf.icerho(T_diff,sal_i2)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.arange(nn+1)+1
            fileinput['Layer temperature'] = np.append(T_diff,T_sn[tt]+273.15)
            fileinput['Layer wetness'] = np.append(W_i2,0)
            fileinput['Layer density'] = np.append(rho_i4,330.)
            fileinput['Layer correlation length'] = np.append(pc_i2,0.15)
            fileinput['Layer exp. correlation length'] = np.append(pc_i2,0.15)
            fileinput['Layer thickness'] = np.append(d_i2,thick_sn[tt])
            fileinput['Layer salinity'] = np.append(sal_i2,0.)
            fileinput['Layer type'] = np.append(sitype_2,2)  
            fileinput['Snow/ice'] = np.append(s_i2,0.)
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

            
    
    else:
    
        if layertype == 'same':
            rho_i4 = sf.icerho(T_i,sal_i2)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = num
            fileinput['Layer temperature'] = T_i
            fileinput['Layer wetness'] = W_i
            fileinput['Layer density'] = rho_i4
            fileinput['Layer correlation length'] = pc_i
            fileinput['Layer exp. correlation length'] = pc_i
            fileinput['Layer thickness'] = d_i
            fileinput['Layer salinity'] = sal_i2
            fileinput['Layer type'] = sitype  
            fileinput['Snow/ice'] = s_i 
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        elif layertype == 'diff':
            T_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),T_i)
            rho_i4 = sf.icerho(T_diff,sal_i2)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.arange(nn)+1
            fileinput['Layer temperature'] = T_diff
            fileinput['Layer wetness'] = W_i2
            fileinput['Layer density'] = rho_i4
            fileinput['Layer correlation length'] = pc_i2
            fileinput['Layer exp. correlation length'] = pc_i2
            fileinput['Layer thickness'] = d_i2
            fileinput['Layer salinity'] = sal_i2
            fileinput['Layer type'] = sitype_2 
            fileinput['Snow/ice'] = s_i2  
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        

    ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simplesal'+saltype+'.dat'
    np.savetxt(ff, fileinput.values, fmt='%f') 

########## SIMPLE DENS
  if 'simpledens' in experiments:
    print(tt, 'simple dens')

    if random_snow=='yes' and thick_sn[tt]>0.:

        if layertype == 'same':
            rho_i2 = sf.icerho(T_i,sal_i)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.append(num,num[-1]+1)
            fileinput['Layer temperature'] = np.append(T_i,T_sn[tt])
            fileinput['Layer wetness'] = np.append(W_i,0)
            fileinput['Layer density'] =  np.append(rho_i2,330.)
            fileinput['Layer correlation length'] = np.append(pc_i,0.15)
            fileinput['Layer exp. correlation length'] = np.append(pc_i,0.15)
            fileinput['Layer thickness'] = np.append(d_i,thick_sn[tt])
            fileinput['Layer salinity'] =  np.append(sal_i,0.)
            fileinput['Layer type'] = np.append(sitype,2)    
            fileinput['Snow/ice'] = np.append(s_i,0.) 
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        elif layertype == 'diff':
            T_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),T_i)
            sal_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),sal_i)
            rho_i2 = sf.icerho(T_diff,sal_diff)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.arange(nn+1)+1
            fileinput['Layer temperature'] = np.append(T_diff,T_sn[tt])
            fileinput['Layer wetness'] = np.append(W_i2,0)
            fileinput['Layer density'] =  np.append(rho_i2,330.)
            fileinput['Layer correlation length'] = np.append(pc_i2,0.15)
            fileinput['Layer exp. correlation length'] = np.append(pc_i2,0.15)
            fileinput['Layer thickness'] =  np.append(d_i2,thick_sn[tt])
            fileinput['Layer salinity'] = np.append(sal_diff,0.)
            fileinput['Layer type'] = np.append(sitype_2,2)   
            fileinput['Snow/ice'] = np.append(s_i2,2) 
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

    
    else:

        if layertype == 'same':
            rho_i2 = sf.icerho(T_i,sal_i)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = num
            fileinput['Layer temperature'] = T_i
            fileinput['Layer wetness'] = W_i
            fileinput['Layer density'] = rho_i2
            fileinput['Layer correlation length'] = pc_i
            fileinput['Layer exp. correlation length'] = pc_i
            fileinput['Layer thickness'] = d_i
            fileinput['Layer salinity'] = sal_i
            fileinput['Layer type'] = sitype  
            fileinput['Snow/ice'] = s_i  
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

        elif layertype == 'diff':
            T_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),T_i)
            sal_diff = np.interp(np.cumsum(d_i2),np.cumsum(d_i),sal_i)
            rho_i2 = sf.icerho(T_diff,sal_diff)
            fileinput = pd.DataFrame()
            fileinput['Layer number'] = np.arange(nn)+1
            fileinput['Layer temperature'] = T_diff
            fileinput['Layer wetness'] = W_i2
            fileinput['Layer density'] = rho_i2
            fileinput['Layer correlation length'] = pc_i2
            fileinput['Layer exp. correlation length'] = pc_i2
            fileinput['Layer thickness'] = d_i2
            fileinput['Layer salinity'] = sal_diff
            fileinput['Layer type'] = sitype_2 
            fileinput['Snow/ice'] = s_i2  
            fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])



    ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simpledens.dat'
    np.savetxt(ff, fileinput.values, fmt='%f') 


      
#####################################################  
###### has to be run in the end when saltype == "func", this part computes the transition between FYI and MYI
####################################################

# the melting and freezing seasons were investigated by examining the time series of ice thickness for each location.

if ee2=='75N00W':
  meltseasons = np.array([i for j in (range(0,57),range(558,867),range(1331,1578),range(2048,2329),range(2771,2998)) for i in j])
  freezeseasons = np.array([i for j in (range(57,558),range(867,1331),range(1578,2048),range(2329,2771),range(2998,3285)) for i in j])
elif ee2=='NorthPole':
  meltseasons = np.array([i for j in (range(0,105),range(579,840),range(1328,1600),range(2078,2330),range(2807,3070)) for i in j])
  freezeseasons = np.array([i for j in (range(105,579),range(840,1328),range(1600,2078),range(2330,2807),range(3070,3285)) for i in j])
elif ee2=='85N50W':
  meltseasons = np.array([i for j in (range(0,121),range(714,875),range(1432,1666),range(2167,2349),range(2887,3118)) for i in j])
  freezeseasons = np.array([i for j in (range(121,714),range(875,1432),range(1666,2167),range(2349,2887),range(3118,3285)) for i in j])
elif ee2=='82N120W':
  meltseasons = np.array([i for j in (range(0,132),range(666,852),range(1418,1612),range(2132,2379),range(2868,3143)) for i in j])
  freezeseasons = np.array([i for j in (range(132,666),range(852,1418),range(1612,2132),range(2379,2868),range(3143,3285)) for i in j])
elif ee2=='80N160W':
  meltseasons = np.array([i for j in (range(0,126),range(661,847),range(1401,1613),range(2128,2389),range(2866,3147)) for i in j])
  freezeseasons = np.array([i for j in (range(126,661),range(847,1401),range(1613,2128),range(2389,2866),range(3147,3285)) for i in j])
elif ee2=='74N170E':
  meltseasons = np.array([i for j in (range(0,128),range(660,857),range(1385,1628),range(2128,2349),range(2868,3036)) for i in j])
  freezeseasons = np.array([i for j in (range(128,660),range(857,1385),range(1628,2128),range(2349,2868),range(3036,3285)) for i in j])
elif ee2=='77N39E':
    meltseasons = np.array([i for j in (range(0,141),range(685,889),range(1443,1588),range(2154,2349),range(2887,3089)) for i in j])
    freezeseasons = np.array([i for j in (range(141,685),range(889,1443),range(1588,2154),range(2349,2887),range(3089,3285)) for i in j])



if saltype == 'func':
  tfreeze = 0
  tmelt = 0
  for tt in range(1,timesteps):
    #where are we looking at?
    if tt > 200 and tt-1 in meltseasons and tt in freezeseasons:
      tfreeze = tt
      print(tfreeze)
    elif tt > 200 and tt-1 in freezeseasons and tt in meltseasons:
      tmelt = tt
      print(tmelt)
      
    #is there a transition going on?  
    if tfreeze > tmelt and tfreeze == tt:
      print("in the loop")
      if icetype_series[tmelt] == 'OW' and icetype_series[tfreeze] == 'FYI':
        change = True
        print('1')
      elif icetype_series[tmelt] == 'FYI' and icetype_series[tfreeze] == 'OW':
        change = True
        print('5')
      elif icetype_series[tmelt] == 'FYI' and icetype_series[tfreeze] == 'FYI':
        change = True
        print('2')
      elif icetype_series[tmelt] == 'FYI' and icetype_series[tfreeze] == 'MYI':
        change = True
        print('3')
      else:
        change = False
        print('4')
      
      if change == True:
        for t in range(tmelt,tfreeze):
  #        print t
          
          #### rewrite the simpleall output
          print('rewrite simple all from ',tmelt,' to ',tfreeze,' now at ',t)
          ff2 = outputpath+'/timestep'+str(t).zfill(4)+'_simpleall'+saltype+'.dat'
          data2 = pd.read_table(ff2, delimiter=" ", names=['Layer number','Layer temperature',\
                                                           'Layer wetness','Layer density',\
                                                           'Layer correlation length','Layer exp. correlation length',\
                                                           'Layer thickness','Layer salinity','Layer type','Snow/ice','Liquid_water_frac'])

          if random_snow=='yes' and thick_sn[t]>0.:
              d_i = data2['Layer thickness'][:-1:]
              sal_i = data2['Layer salinity'][:-1:]
              T_i = data2['Layer temperature'][:-1:]
              rho_i = data2['Layer density'][:-1:]          
              
          else:                                               
              d_i = data2['Layer thickness']
              sal_i = data2['Layer salinity']
              T_i = data2['Layer temperature']
              rho_i = data2['Layer density']

          tmf = (t-tmelt*1.0)/(len(range(tmelt,tfreeze)))
          
          if layertype == 'same':
              [nn, sal_i5] = sf.samelayers_fytomy(d_i,sal_i,tmf)
              rho_i5 = sf.icerho(T_i,sal_i5)
          elif layertype == 'diff':
              [nn, sal_i5] = sf.difflayers_fytomy(ll,d_i,sal_i,tmf)
              rho_i5 = sf.icerho(T_i,sal_i5)
            
          
          fileinput = pd.DataFrame()
          fileinput['Layer number'] = data2['Layer number']
          fileinput['Layer temperature'] = data2['Layer temperature']
          fileinput['Layer wetness'] = data2['Layer wetness']
          if random_snow=='yes' and thick_sn[t]>0.:
              fileinput['Layer density'] = np.append(rho_i5, 330.)         
          else:
              fileinput['Layer density'] = rho_i5              
          fileinput['Layer correlation length'] = data2['Layer correlation length']
          fileinput['Layer exp. correlation length'] = data2['Layer exp. correlation length']
          fileinput['Layer thickness'] = data2['Layer thickness']          
          if random_snow=='yes' and thick_sn[t]>0.:
              fileinput['Layer salinity'] = np.append(sal_i5, 0.)          
          else:
              fileinput['Layer salinity'] = sal_i5          
          fileinput['Layer type'] = data2['Layer type']
          fileinput['Snow/ice'] = data2['Snow/ice']
          fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])
        
          ff = outputpath+'/timestep'+str(t).zfill(4)+'_simpleall'+saltype+'.dat'
          np.savetxt(ff, fileinput.values, fmt='%f')  
          
          #### rewrite the simplesal output
          print('rewrite simple sal from ',tmelt,' to ',tfreeze,' now at ',t)
          ff = outputpath+'/timestep'+str(t).zfill(4)+'_complex.dat'
          data = pd.read_table(ff, delimiter=" ", names=['Layer number','Layer temperature',\
                                                       'Layer wetness','Layer density',\
                                                       'Layer correlation length','Layer exp. correlation length',\
                                                       'Layer thickness','Layer salinity','Layer type','Snow/ice','Liquid_water_frac'])

          if random_snow=='yes' and thick_sn[t]>0.:
              d_i2 = data['Layer thickness'][:-1:]
              sal_i2 = data['Layer salinity'][:-1:]
              T_i2 = data['Layer temperature'][:-1:]
              
          else:                                               
              d_i2 = data['Layer thickness']
              sal_i2 = data['Layer salinity']
              T_i2 = data['Layer temperature']



          if layertype == 'same':
              [nn, sal_i6] = sf.samelayers_fytomy(d_i2,sal_i2,tmf)
              rho_i6 = sf.icerho(T_i2,sal_i6)
              fileinput = pd.DataFrame()
              fileinput['Layer number'] = data['Layer number']
              fileinput['Layer temperature'] = data['Layer temperature']
              fileinput['Layer wetness'] = data['Layer wetness']
              if random_snow=='yes' and thick_sn[t]>0.:
                  fileinput['Layer density'] = np.append(rho_i6, 330.)         
              else:
                  fileinput['Layer density'] = rho_i6       
              fileinput['Layer correlation length'] = data['Layer correlation length']
              fileinput['Layer exp. correlation length'] = data['Layer exp. correlation length']
              fileinput['Layer thickness'] = data['Layer thickness']
              if random_snow=='yes' and thick_sn[t]>0.:
                  fileinput['Layer salinity'] = np.append(sal_i6, 0.)          
              else:
                  fileinput['Layer salinity'] = sal_i6   
              fileinput['Layer type'] = data['Layer type']
              fileinput['Snow/ice'] = data['Snow/ice']
              fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

          elif layertype == 'diff':
              [nn, sal_i6] = sf.difflayers_fytomy(ll,d_i,sal_i2,tmf)
              T_diff = np.interp(np.cumsum(d_i),np.cumsum(d_i2),T_i2)
              rho_i6 = sf.icerho(T_diff,sal_i6)
              fileinput = pd.DataFrame()
              fileinput['Layer number'] = data2['Layer number']
              if random_snow=='yes' and thick_sn[t]>0.:
                  fileinput['Layer temperature'] = np.append(T_diff,(T_diff[-1]+airtemp[t]+273.15)/2)
              else:
                  fileinput['Layer temperature'] = T_diff
              fileinput['Layer wetness'] = data2['Layer wetness']
              if random_snow=='yes' and thick_sn[t]>0.:
                  fileinput['Layer density'] = np.append(rho_i6, 330.)         
              else:
                  fileinput['Layer density'] = rho_i6  
              fileinput['Layer correlation length'] = data2['Layer correlation length']
              fileinput['Layer exp. correlation length'] = data2['Layer exp. correlation length']
              fileinput['Layer thickness'] = data2['Layer thickness']
              if random_snow=='yes' and thick_sn[t]>0.:
                  fileinput['Layer salinity'] = np.append(sal_i6, 0.)          
              else:
                  fileinput['Layer salinity'] = sal_i6   
              fileinput['Layer type'] = data2['Layer type']
              fileinput['Snow/ice'] = data2['Snow/ice']
              fileinput['Liquid_water_frac'] = memls.Vb(fileinput['Layer temperature'],fileinput['Layer salinity'])

              

        
          ff = outputpath+'/timestep'+str(t).zfill(4)+'_simplesal'+saltype+'.dat'
          np.savetxt(ff, fileinput.values, fmt='%f')  
