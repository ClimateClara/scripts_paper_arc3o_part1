# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 17:44:54 2017

Functions needed for analysis

@author: Clara Burgard
"""

def is_summer(month):
    return (month >= 4) & (month <= 9)

def is_winter(month):
    logic = []
    for j,mm in enumerate(month):
        if mm in [1,2,3,10,11,12]: #[1,2,3,4,10,11,12]
            logic.append(True)
        else:
            logic.append(False)
    return logic