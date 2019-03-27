# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:29:05 2019

@author: Matthew
"""

import numpy as np
import matplotlib.pyplot as plt

BHData = np.genfromtxt('BHParameters.csv',delimiter=',')

BHDataR = BHData;
BHDataR[::,3]=BHData[::,3]/BHData[::,2]
BHDataR[::,4]=BHData[::,4]/BHData[::,2]
BHDataR[::,5]=BHData[::,5]/BHData[::,2]
BHDataR[::,6]=BHData[::,6]/BHData[::,2]



fig = plt.figure()
plt.semilogy(BHDataR[0,3::],'-o')
plt.semilogy(BHDataR[1,3::],'-o')
plt.semilogy(BHDataR[3,3::],'-o')
plt.semilogy(BHDataR[5,3::],'-o')
plt.semilogy(BHDataR[7,3::],'-o')
plt.semilogy(BHDataR[9,3::],'-o')