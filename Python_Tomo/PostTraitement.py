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
DossierAcquis = "/home/nicolas/Acquisitions/pollen_bouleau_2407220_autre/"
DossierData = f"{DossierAcquis}data/"
ProcessingFolder = f"{DossierData}Reconstruction"
CheminParam = f"{DossierData}Pretraitement/Param.txt"
dimHolo = int(ft.readvalue(CheminParam,'dimHolo'))

# Creation of results folder
GerchbergFolder = f"{DossierData}Gerchberg"
if not os.path.exists(GerchbergFolder):
    os.makedirs(GerchbergFolder)
    
# Paths to the refraction, and absorption of the object
CheminAbsorp = f"{ProcessingFolder}/AbsorptionDiv_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin"
CheminRefrac = f"{ProcessingFolder}/RefractionDiv_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin"
CheminOTF = f"{ProcessingFolder}/OTF_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin"

# Files reading
Absorption = rp.ReadCube(CheminAbsorp, 2*dimHolo, 2*dimHolo, 2*dimHolo, "np.float32")
Refraction = rp.ReadCube(CheminRefrac, 2*dimHolo, 2*dimHolo, 2*dimHolo, "np.float32")
OTF = rp.ReadCube(CheminOTF, 2*dimHolo, 2*dimHolo, 2*dimHolo, "np.int32")
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
print(f"Reconstruction time for {nbiter} iterations: {np.round(time.time() - start_time,decimals=2)} seconds")
plt.imshow(Rec_Object[:,:,dimHolo].real, cmap="gray")
plt.show()
plt.imshow(Rec_Object[:,:,dimHolo].imag, cmap="gray")
plt.show() 

fidRef = open(f"{GerchbergFolder}/RefractionGerch_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
fidAbs = open(f"{GerchbergFolder}/AbsorptionGerch_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
for cpt in range(FilteredObj.shape[2]):
    Refr = Rec_Object[:,:,cpt].real
    Abs = Rec_Object[:,:,cpt].imag
    Refr.tofile(fidRef)
    Abs.tofile(fidAbs) 
fidRef.close()
fidAbs.close()