#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 16:56:27 2018

download stuff from ERA-Interim

@author: claraburgard
"""

from ecmwfapi import ECMWFDataServer

    
server = ECMWFDataServer()

outputpath = "/home/mpim/m300411/SatSim/SAMSIM/input/ERA-interim/Clara/"
    
lat_of_int = 75
lon_of_int = 00
#lat_of_int = 90
#lon_of_int = 00

## download forcing for SAMSIM
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2005-07-01/to/2009-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "area": str(lat_of_int)+"/"+str(lon_of_int)+"/"+str(lat_of_int)+"/"+str(lon_of_int),
    "levtype": "sfc",
    "param": "31.128/134.128/144.128/165.128/166.128/167.128/168.128/169.128/175.128/176.128/177.128/228.128",
    "step": "3/6/9/12",
    "stream": "oper",
    "time": "00:00:00/12:00:00",
    "type": "fc",
    "format": "netcdf",
    "target": outputpath+"ERA_interim_"+str(lat_of_int)+"N"+str(lon_of_int)+"E_forecast_2005-2009.nc",
    
})

