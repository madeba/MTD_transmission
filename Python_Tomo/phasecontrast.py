#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import time
import os
import matplotlib.pyplot as plt
import numpy as np
import FileTools as ft
import Retropropagation as rp

# Path to the parameter file, and parameters reading
DOSSIERACQUIS = "/home/nicolas/Acquisitions/ACQUIS_pollen_PN/"
DOSSIERDATA = f"{DOSSIERACQUIS}data/"
PROCESSINGFOLDER = f"{DOSSIERDATA}Reconstruction"
CHEMINPARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
DIMHOLO = int(ft.readvalue(CHEMINPARAM, 'dimHolo'))

# Creation of results folder
PHASECONTRAST = f"{DOSSIERDATA}PhaseContrast"
if not os.path.exists(PHASECONTRAST):
    os.makedirs(PHASECONTRAST)

# Paths to the refraction, and absorption of the object
CHEMINABSORP = f"{PROCESSINGFOLDER}/AbsorptionDiv_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff"
CHEMINREFRAC = f"{PROCESSINGFOLDER}/RefractionDiv_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff"

# Files reading
Absorption = ft.ReadtiffCube(CHEMINABSORP)
Refraction = ft.ReadtiffCube(CHEMINREFRAC)
Rec_Object = Refraction + Absorption*1j
plt.imshow(Rec_Object.real[:, :, DIMHOLO], cmap="gray")
plt.show()
plt.imshow(Rec_Object.imag[:, :, DIMHOLO], cmap="gray")
plt.show()

# Darkfield
del Absorption, Refraction
start_time = time.time()
FilteredObj = rp.PhaseContrast(Rec_Object, 10, 500)
print("")
print(f"Reconstruction time: "
      f"{np.round(time.time() - start_time,decimals=2)} seconds")
print("")
plt.imshow(FilteredObj[:, :, DIMHOLO].real, cmap="gray")
plt.show()
plt.imshow(FilteredObj[:, :, DIMHOLO].imag, cmap="gray")
plt.show()

# Data saving
start_time = time.time()
ft.SAVtiffCube(f"{PHASECONTRAST}/PhaseContrast_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
               FilteredObj.imag)
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")
