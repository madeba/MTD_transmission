# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import time
import os
# import matplotlib.pyplot as plt
import numpy as np
import FileTools as ft
import Retropropagation as rp
import manip
import napari

# Data folders and config files
DOSSIERACQUIS = "/home/nicolas/Acquisitions/PETIA/PLA_45678/"
DATA = True # True for data preprocessing, False for white image processing
M = manip.Manip(DOSSIERACQUIS, DATA)
if DATA is True:
    DOSSIERDATA = M.dossier_data
else:
    DOSSIERDATA = M.dossier_blanc

# Creating results Folders
PROCESSINGFOLDER = f"{DOSSIERDATA}Reconstruction"
if not os.path.exists(PROCESSINGFOLDER):
    os.makedirs(PROCESSINGFOLDER)

# Path to the parameter file, and parameters reading
DARKFIELD = False
PHASECONTRAST = False
CHEMINPARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
REWALD = float(ft.readvalue(CHEMINPARAM, 'REwald'))
NB_ANGLE = int(ft.readvalue(CHEMINPARAM, 'nb_angle'))
FMAXHOLO = int(ft.readvalue(CHEMINPARAM, 'fmaxHolo'))
DIMHOLO = int(ft.readvalue(CHEMINPARAM, 'dimHolo'))
PIXTHEO = float(ft.readvalue(CHEMINPARAM, 'pixTheo'))
UBornPitch = 1/(2*FMAXHOLO*PIXTHEO)
NB_HOLO = NB_ANGLE

# Paths to the real, and imaginary parts of the field
CHEMIN_RE_UBORN = f"{DOSSIERDATA}Pretraitement/ReBorn_{DIMHOLO}.tiff"
CHEMIN_IM_UBORN = f"{DOSSIERDATA}Pretraitement/ImBorn_{DIMHOLO}.tiff"

# Path to the specular coordinates
SpecCoordPath = f"{DOSSIERDATA}Pretraitement/Centres_{DIMHOLO}.txt"
fi = rp.Calc_fi(SpecCoordPath, NB_ANGLE, DIMHOLO)

# Field files reading
ReUBorn = ft.ReadtiffCube(CHEMIN_RE_UBORN)
ImUBorn = ft.ReadtiffCube(CHEMIN_IM_UBORN)
UBornCplx = ReUBorn + ImUBorn * 1j
del ReUBorn, ImUBorn

start_time = time.time()
f_recon, TFVol, mask_sum = rp.retropropagation(UBornCplx, NB_HOLO, fi, FMAXHOLO,
                                               REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

Refraction = f_recon.real
Absorption = f_recon.imag
OTF = np.zeros_like(mask_sum)
OTF[mask_sum != 0] = 1

TFVolfilt = np.zeros_like(OTF)
TFVolfilt[TFVol != 0] = 1

viewer = napari.view_image(Absorption.transpose(-1, 1, 0), name='Absorption', colormap='magma')
viewer.add_image(Refraction.transpose(-1, 1, 0), name='Refraction',colormap='magma')

# Writting results
start_time = time.time()
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Refraction_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                Refraction, 2*PIXTHEO*1e6)
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Absorption_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                Absorption, 2*PIXTHEO*1e6)
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTF_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                OTF, 1./(Refraction.shape[0]*2*PIXTHEO*1e6))
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")

# Darkfield processing
if DARKFIELD is True:
    DarkF = abs(Refraction + Absorption)**2
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Darkfield_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                    DarkF, 2*PIXTHEO*1e6)
if PHASECONTRAST is True:
    PhaseC = abs(Refraction + Absorption)**2
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Phasecontrast_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                    PhaseC, 2*PIXTHEO*1e6)