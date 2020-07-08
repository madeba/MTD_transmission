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
DossierAcquis = "/home/nicolas/Acquisitions/BillesCluster/"
DossierData = f"{DossierAcquis}data/"
ProcessingFolder = f"{DossierData}Reconstruction"
CheminParam = f"{DossierData}Pretraitement/Param.txt"
dimHolo = int(ft.readvalue(CheminParam,'dimHolo'))

# Creation of results folder
GerchbergFolder = f"{DossierData}Gerchberg"
if not os.path.exists(GerchbergFolder):
    os.makedirs(GerchbergFolder)
    
# Paths to the refraction, and absorption of the object
CheminAbsorp = f"{ProcessingFolder}/Absorption_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin"
CheminRefrac = f"{ProcessingFolder}/Refraction_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin"
CheminOTF = f"{ProcessingFolder}/OTF_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin"

# Files reading
Absorption = rp.ReadBornCube(CheminAbsorp, 2*dimHolo, 2*dimHolo, 2*dimHolo, False)
Refraction = rp.ReadBornCube(CheminRefrac, 2*dimHolo, 2*dimHolo, 2*dimHolo,False)
OTF = rp.ReadBornCube(CheminOTF, 2*dimHolo, 2*dimHolo, 2*dimHolo, True)
Rec_Object = Refraction + Absorption*1j
plt.imshow(Refraction[:,:,dimHolo], cmap="gray")
plt.show()
plt.imshow(Absorption[:,:,dimHolo], cmap="gray")
plt.show()

# Gerchberg parameters
nbiter = 100    
nmin = 0.0
nmax = .11
kappamin = 0
kappamax = 0.02

# Gerchberg reconstruction
del Absorption, Refraction
start_time = time.time()
FilteredObj = rp.Gerchberg(Rec_Object,OTF,nmin,nmax,kappamin,kappamax,nbiter)
print(f"Reconstruction time for {nbiter} iterations: {np.round(time.time() - start_time,decimals=2)} seconds")
plt.imshow(FilteredObj[:,:,dimHolo].real, cmap="gray")
plt.show()
plt.imshow(FilteredObj[:,:,dimHolo].imag, cmap="gray")
plt.show() 

fidRef = open(f"{GerchbergFolder}/RefractionGerch_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
fidAbs = open(f"{GerchbergFolder}/AbsorptionGerch_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
for cpt in range(FilteredObj.shape[2]):
    Refr = FilteredObj[:,:,cpt].real
    Abs = FilteredObj[:,:,cpt].imag
    Refr.tofile(fidRef)
    Abs.tofile(fidAbs) 
fidRef.close()
fidAbs.close()