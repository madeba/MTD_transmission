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
DossierAcquis = "/home/nicolas/Acquisitions/ACQUIS_pollen_PN18juil/"
DossierData = f"{DossierAcquis}data/"
# DossierAmplitude = 'C:/Users/p1600109/Documents/Recherche/MatlabTomo/Amplitude/'
# DossierPhase = 'C:/Users/p1600109/Documents/Recherche/MatlabTomo/Phase/'
FichierConfig = f"{DossierAcquis}config/config_manip.txt"
CheminMasque = f"{DossierData}Masque.tif"

# Creating results Folders
ProcessingFolder = f"{DossierData}Reconstruction"
if not os.path.exists(ProcessingFolder):
    os.makedirs(ProcessingFolder)

# Path to the parameter file, and parameters reading
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
print(f"Reconstruction time for a {2*dimHolo}x{2*dimHolo}x{2*dimHolo} volume, with {nb_holo} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")

Refraction = f_recon.real
Absorption = f_recon.imag
OTF = np.zeros_like(mask_sum)
OTF[mask_sum != 0] = 1

TFVolfilt = np.zeros_like(OTF)
TFVolfilt[TFVol != 0] = 1

plt.imshow(Refraction[:,:,dimHolo], cmap="gray")
plt.show()
# plt.imshow(OTF[:,:,dimHolo], cmap="gray")
# plt.show()
# plt.imshow(OTF[:,dimHolo,:], cmap="gray")
# plt.show()
# plt.imshow(TFVolfilt[:,dimHolo,:], cmap="gray")
# plt.show()

fidRef = open(f"{ProcessingFolder}/Refraction_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
fidAbs = open(f"{ProcessingFolder}/Absorption_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
fidRedon =  open(f"{ProcessingFolder}/SupRedon_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
fidOTF =  open(f"{ProcessingFolder}/OTF_{2*dimHolo}x{2*dimHolo}x{2*dimHolo}.bin","a")
for cpt in range(Refraction.shape[2]):
    Refr = Refraction[:,:,cpt]
    Abs = Absorption[:,:,cpt]
    Redon = mask_sum[:,:,cpt]
    Support = OTF[:,:,cpt]
    Refr.tofile(fidRef)
    Abs.tofile(fidAbs) 
    Redon.tofile(fidRedon)
    Support.tofile(fidOTF)
fidRef.close()
fidAbs.close()
fidRedon.close()
fidOTF.close()