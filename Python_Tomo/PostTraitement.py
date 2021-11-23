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
import manip
import napari

# Path to the parameter file, and parameters reading
# Data folders and config files
if os.name == 'nt': # Windows
    DOSSIERACQUIS = "C:/Users/p1600109/Documents/Recherche/Acquisitions/Topi/"
else:               # Linux
    DOSSIERACQUIS = "/home/nicolas/Aquisitions/Topi/"
DATA = True # True for data preprocessing, False for white image processing
M = manip.Manip(DOSSIERACQUIS, DATA)
DOSSIERDATA = M.dossier_data
PROCESSINGFOLDER = f"{DOSSIERDATA}Reconstruction"
CHEMINPARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
DIMHOLO = int(ft.readvalue(CHEMINPARAM, 'dimHolo'))
PIXTHEO = float(ft.readvalue(CHEMINPARAM, 'pixTheo'))

# Creation of results folder
GERCHBERGFOLDER = M.dossier_gerchberg
if not os.path.exists(GERCHBERGFOLDER):
    os.makedirs(GERCHBERGFOLDER)

# Paths to the refraction, and absorption of the object
CHEMINABSORP = f"{PROCESSINGFOLDER}/Absorption_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff"
CHEMINREFRAC = f"{PROCESSINGFOLDER}/Refraction_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff"
CHEMINOTF = f"{PROCESSINGFOLDER}/OTF_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff"

# Files reading
Absorption = ft.ReadtiffCube(CHEMINABSORP)
Refraction = ft.ReadtiffCube(CHEMINREFRAC)
OTF = ft.ReadtiffCube(CHEMINOTF)
Rec_Object = Refraction + Absorption*1j
plt.imshow(Rec_Object.real[:, :, DIMHOLO], cmap="gray")
plt.show()
plt.imshow(Rec_Object.imag[:, :, DIMHOLO], cmap="gray")
plt.show()

# Gerchberg parameters
NBITER = 20
NMIN = 0
NMAX = 0.25
KAPPAMIN = 0
KAPPAMAX = 0

# Gerchberg reconstruction
del Absorption, Refraction
start_time = time.time()
FilteredObj = rp.Gerchberg(Rec_Object, OTF, NMIN, NMAX, KAPPAMIN, KAPPAMAX, NBITER)
print("")
print(f"Reconstruction time for {NBITER} iterations: "
      f"{np.round(time.time() - start_time,decimals=2)} seconds")
print("")
plt.imshow(Rec_Object[:, :, DIMHOLO].real, cmap="gray")
plt.show()
plt.imshow(Rec_Object[:, :, DIMHOLO].imag, cmap="gray")
plt.show()

# Data saving
start_time = time.time()
ft.SAVtiffCube(f"{GERCHBERGFOLDER}/RefractionGerch_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
               Rec_Object.real.transpose(1,0,-1), 2*PIXTHEO*1e6)
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")

# Visualization
viewer = napari.view_image(Rec_Object.real.transpose(-1, 1, 0), name='Refraction', colormap='magma', axis_labels=["z", "y", "x"])
