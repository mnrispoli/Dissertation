#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:09:35 2019

@author: matthewrispoli
"""

matica99A = (0.65,0.00,0.00)
matica99B = (0.05,0.53,0.63)
matica99C = (0.44,0.26,0.71)
matica99D = (0.75,0.36,0.13)
matica99E = (0.46,0.56,0.01)
matica99F = (0.66,0.21,0.30)
matica99G = (0.21,0.39,0.80)
matica99H = (0.78,0.52,0.04)
matica99I = (0.52,0.22,0.53)
matica99J = (0.11,0.56,0.42)

matica97A = (0.37,0.51,0.71)
matica97B = (0.88,0.61,0.14)
matica97C = (0.56,0.69,0.19)
matica97D = (0.92,0.39,0.21)
matica97E = (0.53,0.47,0.70)
matica97F = (0.77,0.43,0.10)
matica97G = (0.36,0.62,0.78)
matica97H = (1.00,0.75,0.00)
matica97I = (0.65,0.38,0.61)
matica97J = (0.57,0.59,0.00)

clrWhite = (1,1,1)
clrBlack = (0,0,0)

fontImp = {'family' : 'normal', \
           'size' : 8}
lineSty = {'linewidth' : 1 \
           }

import numpy as np

def blendClr(color1,color2,alphaVal):
    return np.array(color1)*alphaVal + np.array(color2)*(1-alphaVal)


band0Clr=blendClr(matica99B,matica99G,1) 
band1Clr=blendClr(matica99A,matica99F,1)
band2Clr=blendClr(matica99B,matica99G,0.66)
band3Clr=blendClr(matica99A,matica99F,0.66)
band4Clr=blendClr(matica99B,matica99G,0.33)
band5Clr=blendClr(matica99A,matica99F,0.33)
band6Clr=blendClr(matica99B,matica99G,0.15)