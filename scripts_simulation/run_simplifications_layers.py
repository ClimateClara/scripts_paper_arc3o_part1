#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:05:32 2017

run centralized_simplifications but only by changing stuff in simplification_functions

updated to python3 on 11.05.2018

@author: m300411
"""

import numpy as np
import pandas as pd
import simplification_functions as sf
import memls_functions as memls

##########################################################################

#went once through all experiments
#ee='75N00W-p4'
ee='NorthPole-p4'
#ee='82N50E-p4'
#ee='85N50W-p4'
#ee='82N120W-p4'
#ee='80N160W-p4'
#ee='77N39E-p4'
#ee='74N170E-p4'

ee2=ee.split("-")[0]

home_path= #define your home path here

#Output from SAMSIM
inputpath = home_path+'/mistral_work/SatSim/SAMSIM/'
#ERA-Interim input data
inputpath3 = home_path+'/mistral_home/SatSim/SAMSIM/input/ERA-interim/Clara/'

# change outputpath if you want to make several experiments, otherwise there is a danger of overwriting previous results

#output path to put the resulting .dat-files
outputpath = home_path+'/mistral_work/SatSim/MEMLS_exp/INPUT/original_files'


###################################
## choose your salinity scheme
###################################
    
#saltype = 'func' #simplified salinity constant
saltype = 'const' #simplified salinity constant

###################################
## choose which experiments to run
###################################

# default always written out are 'complex' and 'simpleall'

#experiments = ['simpletemp','simplesal']
#experiments = ['simplesal']
#experiments = ['simpletemp']
experiments=[]
#,'simpledens'

#################################
## define the amount of layers
#################################

# layertype = 'same' # not used in this script but rather in the run_simplifications.py
layertype = 'diff'
ll = 10

#############################
# do we have icepack data?
############################

#icepack = 'yes'
icepack = 'no'

##########################
# do we run with snow layers
##########################

random_snow = 'yes'
#random_snow = 'no'

################################################################
#### READ IN THE DATA FROM SAMSIM OUTPUT
free_flag    = 0  

var1name     = 'Bulk salinity'
var1unit     = '[g/kg]'
var1         = np.loadtxt(inputpath+ee+"/dat_S_bu.dat")

var2name     = 'Temperature'
var2unit     = '[C]'
var2         = np.loadtxt(inputpath+ee+"/dat_T.dat")

var3name     = 'Liquid fraction'
var3unit     = '[fraction]'
var3         = np.loadtxt(inputpath+ee+"/dat_psi_l.dat")

thick      = np.loadtxt(inputpath+ee+"/dat_thick.dat")
freeboard  = np.loadtxt(inputpath+ee+"/dat_freeboard.dat")

thick_sn_file = np.loadtxt(inputpath+ee+'/dat_snow.dat') #thick_snow,T_snow,psi_l_snow,psi_s_snow
thick_sn = thick_sn_file[:,0]
T_sn = thick_sn_file[:,1]
T2m_file = np.loadtxt(inputpath+ee+'/dat_T2m_T_top.dat') #T2m, Ttop
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

    
#####################################
###### WHAT IS NEEDED FOR MEMLS #####
#####################################
  
###########################################
######### WRITING THE FILES
############################################
      
timeseries_depth = []
timeseries_depth_ip = []
icetype_series = []    



#need loop over each timestep
for tt in range(timesteps):
  
########## COMPLEX  
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
  

    
  ff = outputpath+'/timestep'+str(tt).zfill(4)+'_complex_'+str(ll).zfill(2)+'layers.dat'
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


  ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simpleall'+saltype+'_'+str(ll).zfill(2)+'layers.dat'
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

        
    
    ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simpletemp_'+str(ll).zfill(2)+'layers.dat'
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

        

    ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simplesal'+saltype+'_'+str(ll).zfill(2)+'layers.dat'
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



    ff = outputpath+'/timestep'+str(tt).zfill(4)+'_simpledens_'+str(ll).zfill(2)+'layers.dat'
    np.savetxt(ff, fileinput.values, fmt='%f') 
      
  
###### HAS TO BE RUN DIRECTLY AFTER "SIMPLE ALL" TO GET SALINITY IN TRANSITION BETWEEN FYI AND MYI

####################################################

if ee=='70N00W-p2':
  meltseasons = np.array([i for j in (range(0,91),range(577,825),range(1256,1577),range(2036,2293),range(2718,2977)) for i in j])
  freezeseasons = np.array([i for j in (range(91,577),range(825,1256),range(1577,2036),range(2293,2718),range(2977,3285)) for i in j])
elif ee=='75N00W-p2' or ee=='75N00W-p4':
  meltseasons = np.array([i for j in (range(0,57),range(558,867),range(1331,1578),range(2048,2329),range(2771,2998)) for i in j])
  freezeseasons = np.array([i for j in (range(57,558),range(867,1331),range(1578,2048),range(2329,2771),range(2998,3285)) for i in j])
elif ee=='75N180E-p2':
  meltseasons = np.array([i for j in (range(0,92),range(621,941),range(1299,1634),range(2055,2402),range(2792,3145)) for i in j])
  freezeseasons = np.array([i for j in (range(92,621),range(941,1299),range(1634,2055),range(2402,2792),range(3145,3285)) for i in j])
elif ee=='80N00E-p2':
  meltseasons = np.array([i for j in (range(0,79),range(568,878),range(1340,1637),range(2071,2381),range(2816,3144)) for i in j])
  freezeseasons = np.array([i for j in (range(79,568),range(878,1340),range(1637,2071),range(2381,2816),range(3144,3285)) for i in j])
elif ee=='80N90E-p2':
  meltseasons = np.array([i for j in (range(0,92),range(624,856),range(1300,1679),range(2059,2384),range(2823,3161)) for i in j])
  freezeseasons = np.array([i for j in (range(92,624),range(856,1300),range(1679,2059),range(2384,2823),range(3161,3285)) for i in j])
elif ee=='85N180E-p2':
  meltseasons = np.array([i for j in (range(0,98),range(621,870),range(1311,1674),range(2076,2407),range(2806,3144)) for i in j])
  freezeseasons = np.array([i for j in (range(98,621),range(870,1311),range(1674,2076),range(2407,2806),range(3144,3285)) for i in j])
elif ee=='barrow-p2':
  meltseasons = np.array([i for j in (range(0,94),range(609,920),range(1288,1618),range(2031,2362),range(2785,2807)) for i in j])
  freezeseasons = np.array([i for j in (range(94,609),range(920,1288),range(1618,2031),range(2362,2785),range(2807,3285)) for i in j])
elif ee=='NorthPole-p2':
  meltseasons = np.array([i for j in (range(0,105),range(579,893),range(1328,1680),range(2078,2396),range(2807,3158)) for i in j])
  freezeseasons = np.array([i for j in (range(105,579),range(893,1328),range(1680,2078),range(2396,2807),range(3158,3285)) for i in j])
elif ee=='sheba-p2':
  meltseasons = np.array([i for j in (range(0,112),range(615,902),range(1362,1625),range(2090,2373),range(2820,3141)) for i in j])
  freezeseasons = np.array([i for j in (range(112,615),range(902,1362),range(1625,2090),range(2373,2820),range(3141,3285)) for i in j])
elif ee=='NorthPole-p4':
  meltseasons = np.array([i for j in (range(0,105),range(579,840),range(1328,1600),range(2078,2330),range(2807,3070)) for i in j])
  freezeseasons = np.array([i for j in (range(105,579),range(840,1328),range(1600,2078),range(2330,2807),range(3070,3285)) for i in j])
elif ee=='82N50E-p4':
  meltseasons = np.array([i for j in (range(0,131),range(693,820),range(1425,1599),range(2153,2361),range(2885,3060)) for i in j])
  freezeseasons = np.array([i for j in (range(131,693),range(820,1425),range(1599,2153),range(2361,2885),range(3060,3285)) for i in j])
elif ee=='85N50W-p4':
  meltseasons = np.array([i for j in (range(0,121),range(714,875),range(1432,1666),range(2167,2349),range(2887,3118)) for i in j])
  freezeseasons = np.array([i for j in (range(121,714),range(875,1432),range(1666,2167),range(2349,2887),range(3118,3285)) for i in j])
elif ee=='82N120W-p4':
  meltseasons = np.array([i for j in (range(0,132),range(666,852),range(1418,1612),range(2132,2379),range(2868,3143)) for i in j])
  freezeseasons = np.array([i for j in (range(132,666),range(852,1418),range(1612,2132),range(2379,2868),range(3143,3285)) for i in j])
elif ee=='80N160W-p4':
  meltseasons = np.array([i for j in (range(0,126),range(661,847),range(1401,1613),range(2128,2389),range(2866,3147)) for i in j])
  freezeseasons = np.array([i for j in (range(126,661),range(847,1401),range(1613,2128),range(2389,2866),range(3147,3285)) for i in j])
elif ee=='74N170E-p4':
  meltseasons = np.array([i for j in (range(0,128),range(660,857),range(1385,1628),range(2128,2349),range(2868,3036)) for i in j])
  freezeseasons = np.array([i for j in (range(128,660),range(857,1385),range(1628,2128),range(2349,2868),range(3036,3285)) for i in j])
elif ee=='77N39E-p4':
    meltseasons = np.array([i for j in (range(0,141),range(685,889),range(1443,1588),range(2154,2349),range(2887,3089)) for i in j])
    freezeseasons = np.array([i for j in (range(141,685),range(889,1443),range(1588,2154),range(2349,2887),range(3089,3285)) for i in j])



#########################################################
  

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
          ff2 = outputpath+'/timestep'+str(t).zfill(4)+'_simpleall'+saltype+'_'+str(ll).zfill(2)+'layers.dat'
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
        
          ff = outputpath+'/timestep'+str(t).zfill(4)+'_simpleall'+saltype+'_'+str(ll).zfill(2)+'layers.dat'
          np.savetxt(ff, fileinput.values, fmt='%f')  
          
          #### rewrite the simplesal output
          print('rewrite simple sal from ',tmelt,' to ',tfreeze,' now at ',t)
          ff = outputpath+'/timestep'+str(t).zfill(4)+'_complex_'+str(ll).zfill(2)+'layers.dat'
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
              #rho_i6 = icerho2(T_i2,sal_i6,psi_l2)
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

              

        
          ff = outputpath+'/timestep'+str(t).zfill(4)+'_simplesal'+saltype+'_'+str(ll).zfill(2)+'layers.dat'
          np.savetxt(ff, fileinput.values, fmt='%f')  

