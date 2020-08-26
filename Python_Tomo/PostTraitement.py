# -*- coding: utf-8 -*-
"""

@author: Nicolas Verrier
"""
import matplotlib.pyplot as plt
import numpy as np
import FileTools as ft
import Retropropagation as rp
import time
import os

# Path to the parameter file, and parameters reading
DossierAcquis = "/home/nicolas/Acquisitions/ACQUIS_pollen_PN/"
DossierData = f"{DossierAcquis}data/"
ProcessingFolder = f"{DossierData}Reconstruction"
CheminParam = f"{DossierData}Pretraitement/Param.txt"
dimHolo = int(ft.readvalue(CheminParam,'dimHolo'))

# Creation of results folder
GerchbergFolder = f"{DossierData}Gerchberg"
if not os.path.exists(GerchbergFolder):
    os.makedirs(GerchbergFolder)
    
# Paths to the refraction, and absorption of the object
CheminAbsorp = f"{ProcessingFolder}/AbsorptionDiv_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff"
CheminRefrac = f"{ProcessingFolder}/RefractionDiv_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff"
CheminOTF = f"{ProcessingFolder}/OTF_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff"

# Files reading
Absorption = ft.ReadtiffCube(CheminAbsorp)
Refraction = ft.ReadtiffCube(CheminRefrac)
OTF = ft.ReadtiffCube(CheminOTF)
Rec_Object = Refraction + Absorption*1j
plt.imshow(Rec_Object.real[:,:,dimHolo], cmap="gray")
plt.show()
plt.imshow(Rec_Object.imag[:,:,dimHolo], cmap="gray")
plt.show()

# Gerchberg parameters
nbiter = 20    
nmin = 0
nmax = 0.17
kappamin = 0
kappamax = 0

# Gerchberg reconstruction
del Absorption, Refraction
start_time = time.time()
FilteredObj = rp.Gerchberg(Rec_Object,OTF,nmin,nmax,kappamin,kappamax,nbiter)
print("")
print(f"Reconstruction time for {nbiter} iterations: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")
plt.imshow(Rec_Object[:,:,dimHolo].real, cmap="gray")
plt.show()
plt.imshow(Rec_Object[:,:,dimHolo].imag, cmap="gray")
plt.show() 

# Data saving
start_time = time.time()
ft.SAVtiffCube(f"{GerchbergFolder}/RefractionGerch_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff", Rec_Object.real)
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")