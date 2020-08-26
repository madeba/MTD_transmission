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

# Data folders and config files
DossierAcquis = "/home/nicolas/Acquisitions/ACQUIS_pollen_PN/"
DossierData = f"{DossierAcquis}data/"

# Creating results Folders
ProcessingFolder = f"{DossierData}Reconstruction"
if not os.path.exists(ProcessingFolder):
    os.makedirs(ProcessingFolder)

# Path to the parameter file, and parameters reading
DarkField = False
CheminParam = f"{DossierData}Pretraitement/Param.txt"
REwald = ft.readvalue(CheminParam,'REwald')
nb_angle = int(ft.readvalue(CheminParam,'nb_angle'))
fmaxHolo = int(ft.readvalue(CheminParam,'fmaxHolo'))
dimHolo = int(ft.readvalue(CheminParam,'dimHolo'))
pixTheo = ft.readvalue(CheminParam,'pixTheo')
lambda_0 = 632.8e-9
n_imm = 1.515
UBornPitch = 1/(2*fmaxHolo*pixTheo)
nb_holo = nb_angle

# Paths to the real, and imaginary parts of the field
CheminReUBorn = f"{DossierData}Pretraitement/ReBorn_{dimHolo}.bin"
CheminImUBorn = f"{DossierData}Pretraitement/ImBorn_{dimHolo}.bin"

# Path to the specular coordinates
SpecCoordPath = f"{DossierData}Pretraitement/Centres_{dimHolo}.txt"
fi = rp.Calc_fi(SpecCoordPath, nb_angle, REwald,dimHolo)

# Field files reading
ReUBorn = rp.ReadCube(CheminReUBorn, dimHolo, dimHolo, nb_angle, "np.float64")
ImUBorn = rp.ReadCube(CheminImUBorn, dimHolo, dimHolo, nb_angle, "np.float64")
UBornCplx = ReUBorn + ImUBorn * 1j
del ReUBorn, ImUBorn

start_time = time.time()
f_recon, TFVol, mask_sum = rp.retropropagation(UBornCplx,nb_holo,fi,fmaxHolo,REwald,lambda_0,1.515,2*dimHolo,dimHolo,pixTheo,2*dimHolo,UBornPitch)
print(f"Reconstruction time for a {2*dimHolo}x{2*dimHolo}x{2*dimHolo} volume (3D-FFT included), with {nb_holo} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

Refraction = f_recon.real
Absorption = f_recon.imag
OTF = np.zeros_like(mask_sum)
OTF[mask_sum != 0] = 1

TFVolfilt = np.zeros_like(OTF)
TFVolfilt[TFVol != 0] = 1

plt.imshow(Refraction[:,:,dimHolo], cmap="gray")
plt.show()

# Writting results
start_time = time.time()
ft.SAVtiffCube(f"{ProcessingFolder}/Refraction_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff", Refraction)
ft.SAVtiffCube(f"{ProcessingFolder}/Absorption_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff", Absorption)
# ft.SAVtiffCube(f"{ProcessingFolder}/OTF_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff", OTF)
# ft.SAVtiffCube(f"{ProcessingFolder}/SupRedon_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff", mask_sum)
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")


# Darkfield processing  
if DarkField is True:
    DarkF = (Refraction + Absorption)**2    
    ft.SAVtiffCube(f"{ProcessingFolder}/Darkfield_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.tiff", DarkF.transpose((-1, 0, 1)).astype(np.float32))