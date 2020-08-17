#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 10:01:38 2017

Functions used in run_simplifications

@author: Clara Burgard
"""

import numpy as np
import memls_functions as memls
from scipy import stats
import simplification_functions as sf


def compute_ice_snow_int_temp(sit, snd, tsi):
	"""
    This function computes the temperature at the snow-ice interface
    inspired from Semtner, 1976

    INPUT
    sit : sea-ice thickness in m
    snd : snow thickness in m
    tsi : sea-ice (or snow) surface temperature in K

    OUTPUT
    T_i : temperature at snow-ice interface in K
    """

	k_s = 0.31  # thermal conductivity of snow in W/K/m #from MPIOM paper
	k_i = 2.17  # thermal conductivity of ice in W/K/m #from MPIOM paper
	bottom_temp = tsi * 0 + 273.15 - 1.8
	T_i = ((tsi * (k_s / snd)) + (bottom_temp * (k_i / sit))) / ((k_s / snd) + (k_i / sit))
	return T_i

def icerho2(T,S,lwf): 
	"""
	Define the ice density with liquid water fraction given by SAMSIM. From Notz (2005)
	
	INPUT
	T: temperature in °C or K
	S: salinity in g/kg
	lwf: liquid water (=brine volume) fraction
	
	OUTPUT
	rho_tot: sea-ice density in kg/m**3
	"""

	if T[0] > 100:
		T = T-273.15
	if len(T) == 1:
		#density of seawater
		rho_sw = 1000.3 + 0.78237*S + 2.8008*10**-4*S**2
		rho_tot = rho_sw
	else:
		rho_0 = 916.18 - 0.1403*T
		Sb = memls.Sb(T)
		rho_w = 1000.3 + 0.78237*Sb + 2.8008*(10**-4)*Sb**2
		#rho_w = 998.43 + 0.69722*Sb + 2.5201*(10**-4)*Sb**2
		rho_tot = lwf*rho_w + (1-lwf)*rho_0
	return rho_tot

def icerho(T,S): #von Dirk
	"""
	define the ice density from temperature and salinity. From Notz (2005).
	
	INPUT
	T: temperature in °C or K
	S: salinity in g/kg
	
	OUTPUT
	rho_tot: sea-ice density in kg/m**3
	"""
	if np.any(T>100.):
		T = T-273.15
	#density of pure ice
	rho_0 = 916.18 - 0.1403*T
	Sbr = memls.Sb(T)
	rho_w = 1000.3 + 0.78237*Sbr + 2.8008*(10**-4)*Sbr**2
	lwf=memls.Vb(T,S)
	rho_tot = lwf*rho_w + (1-lwf)*rho_0
	return rho_tot   

def define_fyi_myi(timeseries_depth,tt):
	"""
	define which ice is first-year ice (FYI) and which is multiyear ice (MYI)
	
	INPUT
	timeseries_depth: timeseries of the ice thickness
	tt: current timestep
	
	OUTPUT
	icetype: 'FYI' or 'MYI'
	"""
	
	icedepth = np.array(timeseries_depth)
	if tt > 730:
		iceidx = np.zeros(len(icedepth[tt-730:tt]))
		iceidx[icedepth[tt-730:tt]<0.02] = 1
		iceidx[np.isnan(icedepth[tt-730:tt])] = 1
		#print iceidx
	else:
		iceidx = [1.] 
	#print iceidx
  
	if sum(iceidx) > 0:
		icetype = 'FYI'
		#print icetype
	else:
		icetype = 'MYI'
		#print icetype
    
	return icetype


def sal_approx_fy(norm_z):
	"""
	compute the salinity profile as a function of depth for first-year ice, from Griewank and Notz (2015)
	
	INPUT
	norm_z: normalized depth (0 at the top, 1 at the bottom)
	
	OUTPUT
	sal_fyi: salinity profile in g/kg
	"""	
	a=1.0964
	b=-1.0552
	c=4.41272
	sal_fy = norm_z/(a+b*norm_z)+c
	return sal_fy

 
def sal_approx_my(norm_z):
	"""
	compute the salinity profile as a function of depth for multiyear ice, from Griewank and Notz (2015)
	
	INPUT
	norm_z: normalized depth (0 at the top, 1 at the bottom)
	
	OUTPUT
	sal_myi: salinity profile in g/kg
	"""	
	a=0.17083
	b=0.92762
	c=0.024516
	sal_my = np.zeros(len(norm_z))
	for n in range(len(sal_my)):
		if norm_z[n] > 0.0001:  
			sal_my[n] = norm_z[n]/a + (norm_z[n]/b)**(1/c)
	return sal_my 


def sal_approx_fytomy(norm_z,t):
	"""
	compute the salinity profile as a function of depth in the transition between first-year and multiyear ice, 
	from Griewank and Notz (2015)
	
	INPUT
	norm_z: normalized depth (0 at the top, 1 at the bottom)
	
	OUTPUT
	sal_fytomy: salinity profile in g/kg
	"""	
	sal_fytomy = (1-t)*sal_approx_fy(norm_z) + t*sal_approx_my(norm_z)
	return sal_fytomy 
    

def find_ttransition(icetypes):
	"""
	find the period where the transition between FYI and MYI or between open water and first-year ice occurs
	
	INPUT
	icetypes: list with icetypes
		
	OUTPUT
	ttransition: list with timesteps of the transition
	"""	
	icetype_series = np.array(icetypes)
	ttransition = []  
	for tt in range(1,len(icetype_series)):
		if icetype_series[tt-1] == 'FYI' and icetype_series[tt] == 'MYI':
			ttransition.append(tt)
		if icetype_series[tt-1] == 'OW' and icetype_series[tt] == 'FYI':
			ttransition.append(tt)
	return ttransition


###### find the tmelt and tfreeze for the transition between FYI and MYI
def find_tmelttfreeze(icetypes,kdays,kmonths,timeseries_depth):
	"""
	find the period where the transition between FYI and MYI or between open water and first-year ice occurs
	
	INPUT
	icetypes: list with icetypes
	kdays: how many days of consecutive melting/freezing are needed to be defined as the onset
	kmonths: months to look for the onset before ttrans
	timeseries_depth: timeseries of ice thickness
		
	OUTPUT
	tmelt: approximate timestep where melt starts
	tfreeze: approximate timestep where freezing starts again
	"""	
	ttransition = find_ttransition(icetypes)
	k = kdays #days of melting/freezing for defining that it is the onset
	kk = kmonths*30*2 #months to look for the onset before ttrans
	tmelt = []
	tfreeze = []
	for i in range(len(ttransition)-1,-1,-1):
		if ttransition[i]-kk < 0:
			del ttransition[i]
	for ttrans in ttransition:        
		d_trend = np.zeros(len(range(ttrans-kk,ttrans+kk/2+1)))    
		for i,ts in enumerate(range(ttrans-kk,ttrans+kk/2+1)):
			#compute trend in total thickness over given timespan
			d_trend[i], intercept, r_value, p_value, std_err = stats.linregress(range(k*2),timeseries_depth[ts:ts+k*2])
			#print d_trend[i]

	for i,ts in enumerate(range(ttrans-kk,ttrans+kk/2+1)):
		# print d_trend[i]
		# print ts
		if d_trend[i]<0 and ts==ttrans-kk : #if the trend is negative in the beginning, it could be that the melting season startes even earlier
			print('melt season starts before')
		elif d_trend[i]<0 and ts==ttrans+kk/2+1 : # if the trend is positive in the end, it could be that the freezing season starts even later
			print('melt season ends after')
		#elif ts<(ttrans+kk/2):
			# print 'yep'
		if d_trend[i]>0 and d_trend[i+1]<=0:
			print('Begin of melt season')
			tmelt.append(ts)
		elif d_trend[i]<0 and d_trend[i+1]>=0:
			print('End of melt season')
			tfreeze.append(ts)
	return tmelt, tfreeze




### BARE ICE

def orig_complex_iceonly(tt,orig_T_step, orig_depth_step, orig_sal_step, orig_lwf_step,timeseries_depth):
	"""
	This function writes the SAMSIM results in the form that they can be used as input for MEMLS
	
	INPUT
	tt: timestep
	orig_T_step: SAMSIM temperature profile output in step form
	orig_depth_step: SAMSIM thickness profile output in step form
	orig_sal_step: SAMSIM salinity profile output in step form
	orig_sal_step: SAMSIM liquid water fraction profile output in step form
	timeseries_depth: list to collect the ice thickness over time
	
	OUTPUT
	Tnew: Temperature profile in K
	num: Number of layers 
	W_i: Wetness profile (set to 0)
	rho_i: Density profile in kg/m**3
	pc_i: Correlation length profile in mm
	d_i: Thickness profile in m
	sal_i: Salinity profile in g/kg
	sitype : snow 1 /first year ice 3 /multiyear ice 4
	s_i: snow 0, ice 1
	icetype: 'FYI' or 'MYI' or 'OW'
	psi_l: liquid water fraction
	"""
	
	#####Layer temperature
	T_i0 = np.flipud(orig_T_step[tt,0::2])
	T_i = T_i0[T_i0 != -1.0]
	
	#####Layer number
	num = range(1,len(T_i)+1)  
	#if only one layer and too hot, it is probably ocean 
	if len(num) == 1 and T_i[0]>= -1.0:
		Tnew=np.zeros(len(num))
		Tnew[:] = np.nan
		W_i = np.zeros(len(num))
		W_i[:] = np.nan
		rho_i = np.zeros(len(num))
		rho_i[:] = np.nan
		pc_i = np.zeros(len(num))
		pc_i[:] = np.nan
		d_i = np.zeros(len(num))
		d_i[:] = np.nan
		timeseries_depth.append(0.0)
		sal_i = np.zeros(len(num))
		sal_i[:] = np.nan
		sitype = np.zeros(len(num))
		sitype[:] = np.nan
		s_i = np.zeros(len(num))
		s_i[:] = np.nan
		icetype = np.zeros(len(num))
		icetype = 'OW'
		psi_l = np.zeros(len(num))
		psi_l[:] = np.nan
	else:  
		#####Layer wetness (liquid water fraction)
		W_i = np.zeros(len(num)) 
		 
		#####Layer thickness
		d_i0 = np.flipud(orig_depth_step[tt,1::2])
		d_i0 = d_i0[T_i0 != -1.0]
		d_i = np.zeros(len(num))
		for n in range(len(num)):
			if n==len(num)-1:
				d_i[n] = -d_i0[-1]
			else:  
				d_i[n] = d_i0[n+1] - d_i0[n]
		timeseries_depth.append(-d_i0[0])
		
		#####Layer salinity
		sal_i = np.flipud(orig_sal_step[tt,0::2])
		sal_i = sal_i[T_i0 != -1.0]
		
		#####Liquid water fraction
		psi_l = np.flipud(orig_lwf_step[tt,0::2])
		psi_l = psi_l[T_i0 != -1.0]
		
		#####Layer density  
		rho_i = icerho2(T_i,sal_i,psi_l)
		
		#####Define FYI and MYI
		icetype = define_fyi_myi(timeseries_depth,tt)
		
		#####Correlation length
		pc_i = np.zeros(len(num))
		if icetype == 'FYI':
			for n in range(len(num)):
				if d_i0[n] > -0.20:
					pc_i[n] = 0.35
				else:
					pc_i[n] = 0.25
		elif icetype == 'MYI':
			# more sensitivity studies could be done here as well
		 	pc_i[:] = 1.5

		####Sea ice or snow
		s_i = np.zeros(len(num))
		s_i[:] = 1.0
		sitype = np.zeros(len(num))
		if icetype == 'FYI':
			sitype[:] = 3
		elif icetype == 'MYI':
			sitype[:] = 4

		####avoid that T=0 and S=0
		Tnew = np.ones(len(T_i))*np.nan
		for n, T in enumerate(T_i):
			if T == 0 and sal_i[n] == 0.001:
				if psi_l[n]>0:
					Tnew[n] = - (0.05411*sal_i[n])/(psi_l[n] - 0.001*sal_i[n])
				else:
					Tnew[n] = 0.0
		else:
			Tnew[n] = T

		####transform T in K
		Tnew = Tnew + 273.15

	return Tnew, num, W_i, rho_i, pc_i, d_i, sal_i, sitype, s_i, icetype, psi_l
  
def samelayers(T,d,W,pc,sal,saltype,rho,s,icetype,sitype,thick_sn,thick_ice,Ttop,double_linear):
	"""
	This function simplifies the profiles from SAMSIM as could be inferred from a simple climate model
	output on the same amount of layers as the reference and writes them out to be used with MEMLS
	
	INPUT
	T: reference temperature profile in K
	d: reference thickness profile in m
	W: reference wetness profile (set to 0)
	pc: reference correlation length profile in mm
	sal: reference salinity profile in g/kg
	saltype: "func" if function of depth, "const" if constant at 5 g/kg and 1 g/kg
	rho: reference density profile in kg/m**3
	s: snow 0, ice 1
	icetype: 'FYI' or 'MYI' or 'OW'
	sitype : snow 1 /first year ice 3 /multiyear ice 4
	thick_sn: snow thickness in m
	thick_ice: ice thickness in m
	Ttop: temperature at top of the ice
	double_linear:  if "yes", the linear profile through snow and ice is computed based only on the temperature at 
					the top of the snow
					if "no", we use the ice surface temperature from SAMSIM as a start for the linear profile
	
	
	OUTPUT
	nn: Number of layers 
	d_new: Thickness profile in m
	T_new: Simplified temperature profile in K
	W_new: Simplified wetness profile (set to 0)
	pc_new: Simplified correlation length profile in mm
	sal_new: Simplified salinity profile in g/kg
	rho_new: Simplified density profile in kg/m**3
	s_new: snow 0, ice 1
	sitype_new: snow 1 /first year ice 3 /multiyear ice 4
	icetype_new: 'FYI' or 'MYI' or 'OW'
	"""
	layers = len(d)
	nn = layers
	if layers > 1:
		##### Keep parameters from complex
		d_new = d
		W_new = W
		pc_new = pc
		s_new = s
		icetype_new = icetype
		sitype_new = sitype

		#####Linear temperature profile
	 	#di2 = np.linspace(d[0],np.sum(d[:-1]),layers+5)
		di2 = np.linspace(d[0],np.sum(d),layers+5)
		if double_linear=='yes' and thick_sn > 0:
			new_Ttopice = compute_ice_snow_int_temp(thick_ice,thick_sn,Ttop)
			Ti2 = np.linspace(273.15-1.8,new_Ttopice,layers+5)
		else:
			Ti2 = np.linspace(273.15-1.8,T[-1],layers+5)
			T_new = np.interp(np.cumsum(d),di2,Ti2)

		#####Salinity profiles after functions or constant
		if saltype == 'func':
			norm_z = 1-np.cumsum(d_new/sum(d_new))
			if icetype == 'FYI':
				sal_new=sf.sal_approx_fy(norm_z)
			elif icetype == 'MYI':
				sal_new=sf.sal_approx_my(norm_z)
			#sal_new=flipud(sal_new)
		elif saltype == 'const':
			sal_new = np.zeros(len(d))
			if icetype == 'FYI':
				#sal_new[:]=10. #commented for exp023
		 		sal_new[:]=5.
			elif icetype == 'MYI':
				#sal_new[:]=5.  #commented for exp023
				sal_new[:]=1.

		rho_new = sf.icerho(T_new,sal_new) ######Density only from T and S

	else:
		nn = np.array([1])
		d_new = d
		T_new = T
		W_new = W
		pc_new = pc
		sal_new = sal
		rho_new = rho
		s_new = s
		icetype_new = icetype
		sitype_new = sitype

	return [nn, d_new, T_new, W_new, pc_new, sal_new, rho_new, s_new, sitype_new, icetype_new]
  
def samelayers_fytomy(d,sal,t):
	"""
	This function makes the salinity transition between FYI and MYI
		
	INPUT
	d: thickness profile in m
	sal: simplified salinity profile
	t: current timestep
	
	OUTPUT
	nn: Number of layers 
	sal_new: new simplified salinity profile during transition
	"""

	layers = len(d)
	nn = layers
	if layers > 1:
		norm_z = 1-np.cumsum(d/sum(d))
		sal_new = sal_approx_fytomy(norm_z,t)
	else:
		nn = np.array([1])
		sal_new = sal
	return [nn, sal_new]
  

#difflayers: equidistant layers
def difflayers(layers,T,d,W,pc,sal,saltype,rho,s,icetype,sitype,thick_sn,thick_ice,Ttop,double_linear):
	"""
	This function simplifies the profiles from SAMSIM as could be inferred from a simple climate model
	output on a different amount of layers as the reference and writes them out to be used with MEMLS
	
	INPUT
	layers: number of layers the new profile should be interpolated to
	T: reference temperature profile in K
	d: reference thickness profile in m
	W: reference wetness profile (set to 0)
	pc: reference correlation length profile in mm
	sal: reference salinity profile in g/kg
	saltype: "func" if function of depth, "const" if constant at 5 g/kg and 1 g/kg
	rho: reference density profile in kg/m**3
	s: snow 0, ice 1
	icetype: 'FYI' or 'MYI' or 'OW'
	sitype : snow 1 /first year ice 3 /multiyear ice 4
	thick_sn: snow thickness in m
	thick_ice: ice thickness in m
	Ttop: temperature at top of the ice
	double_linear:  if "yes", the linear profile through snow and ice is computed based only on the temperature at 
					the top of the snow
					if "no", we use the ice surface temperature from SAMSIM as a start for the linear profile
	
	
	OUTPUT
	nn: Number of layers 
	d_new: Thickness profile in m
	T_new: Simplified temperature profile in K
	W_new: Simplified wetness profile (set to 0)
	pc_new: Simplified correlation length profile in mm
	sal_new: Simplified salinity profile in g/kg
	rho_new: Simplified density profile in kg/m**3
	s_new: snow 0, ice 1
	sitype_new: snow 1 /first year ice 3 /multiyear ice 4
	icetype_new: 'FYI' or 'MYI' or 'OW'
	"""
	nn = layers
	num = range(nn)
	if len(d) > 1:

		##### Keep parameters from complex but interpolate to less layers
		d_0 = np.linspace(d[0],np.sum(d),layers)
		d_new = np.zeros(len(num))
		for n in range(len(num)):
			if n==0:
				d_new[n] = d_0[0]
			else:
				d_new[n] = d_0[n] - d_0[n-1]
          
		W_new = np.interp(d_0,np.cumsum(d),W)
		pc_new = np.interp(d_0,np.cumsum(d),pc)  
		s_new = np.interp(d_0,np.cumsum(d),s)
		icetype_new = icetype
		sitype_new = np.interp(d_0,np.cumsum(d),sitype)

		if double_linear=='yes' and thick_sn > 0: #####Linear temperature profile
			new_Ttopice = compute_ice_snow_int_temp(thick_ice,thick_sn,Ttop)
			T_new = np.linspace(273.15-1.8,new_Ttopice,layers)
		else:
			T_new = np.linspace(273.15-1.8,T[-1],layers)

		if saltype == 'func': #####Salinity profiles after functions or constant
			norm_z = 1-np.cumsum(d_new/sum(d_new))
			if icetype == 'FYI':
				sal_new=sal_approx_fy(norm_z)
			elif icetype == 'MYI':
				sal_new=sal_approx_my(norm_z)
				# sal_new=np.flipud(sal_new)
		elif saltype == 'const':
			sal_new = np.zeros(nn)
			if icetype == 'FYI':
				sal_new[:]=10.
			elif icetype == 'MYI':
				sal_new[:]=5.

		######Density only from T and S
		rho_new = icerho(T_new,sal_new)

	else:
		nn = np.array([1])
		d_new = d
		T_new = T
		W_new = W
		pc_new = pc
		sal_new = sal
		rho_new = rho
		s_new = s
		icetype_new = icetype
		sitype_new = sitype

	return [nn, d_new, T_new, W_new, pc_new, sal_new, rho_new, s_new, sitype_new, icetype_new]

def difflayers_fytomy(layers,d,sal,t):
	"""
	This function makes the salinity transition between FYI and MYI on a different amount of layers 
	than the reference
		
	INPUT
	layers: amount of layers of the simplified profiles
	d: thickness profile in m
	sal: simplified salinity profile
	t: current timestep
	
	OUTPUT
	nn: Number of layers 
	sal_new: new simplified salinity profile during transition
	"""

	nn = layers
	if len(d) > 1:
		d_0 = np.linspace(d[0],np.sum(d),layers)
		d_new = np.zeros(nn)
		for n in range(nn):
			if n==0:
				d_new[n] = d_0[0]
			else:
				d_new[n] = d_0[n] - d_0[n-1]
			norm_z = 1-np.cumsum(d_new/sum(d_new))
			sal_new = sal_approx_fytomy(norm_z,t)
	else:
		nn = np.array([1])
		sal_new = sal
	return [nn, sal_new]
