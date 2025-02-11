#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 11:03:31 2020

@author: nicolas
"""

import numpy as np
import matplotlib.pyplot as plt

rho = 125
turns = 15
Nb_ill = 1000
t = np.arange(0, 2*np.pi*turns, step = 2*np.pi*turns/Nb_ill)

x = rho/turns*t/(2*np.pi)*np.cos(t)
y = rho/turns*t/(2*np.pi)*np.sin(t)

plt.plot(x,y,".")
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# Abscisse Curviligne
L = rho/2*(np.log(2*np.pi*turns+np.sqrt((2*np.pi*turns)**2+1))+2*np.pi*turns*np.sqrt((2*np.pi*turns)**2+1))
theta = 0
deltatheta = 0
xc = np.zeros((Nb_ill))
yc = np.zeros((Nb_ill))
for ill in range(Nb_ill):
    theta = theta + deltatheta
    deltatheta = L/Nb_ill*1/(rho*np.sqrt(1+(theta)**2))
    xc[ill]=rho*theta/(2*np.pi*turns)*np.cos(theta)
    yc[ill]=rho*theta/(2*np.pi*turns)*np.sin(theta)

plt.plot(xc,yc,".")
plt.gca().set_aspect('equal', adjustable='box')
plt.show()