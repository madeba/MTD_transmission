# -*- coding: utf-8 -*-
"""
@author: Steve Laroche et Nicolas Verrier
"""

import time
import os
import matplotlib.pyplot as plt
import numpy as np
import FileTools as ft
import Retropropagation as rp
from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
import MultiModalMTD as mmtd
import matplotlib.pyplot as plt
import manip

# Data folders and config files
if os.name == 'nt': # Windows
    DOSSIERACQUIS = "C:/Users/p1600109/Documents/Recherche/Acquisitions/Topi/"
else:               # Linux
    DOSSIERACQUIS = "/home/nicolas/Acquisitions/Topi/"

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

# Choosing method
Method = {0 : "BASE",
          1 : "DARKFIELD",
          2 : "PHASECONTRAST",
          3 : "RHEINBERG",
          4 : "DIC"
          }
MethodUsed = Method[4]

# Path to the parameter file, and parameters reading
CHEMINPARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
REWALD = float(ft.readvalue(CHEMINPARAM, 'REwald'))
NB_ANGLE = int(ft.readvalue(CHEMINPARAM, 'nb_angle'))
FMAXHOLO = int(ft.readvalue(CHEMINPARAM, 'fmaxHolo'))
DIMHOLO = int(ft.readvalue(CHEMINPARAM, 'dimHolo'))
PIXTHEO = float(ft.readvalue(CHEMINPARAM, 'pixTheo'))
UBornPitch = 1/(2*FMAXHOLO*PIXTHEO)
NB_HOLO = NB_ANGLE

# Path to the specular coordinates
SpecCoordPath = f"{DOSSIERDATA}Pretraitement/Centres_{DIMHOLO}.txt"
fi = rp.Calc_fi(SpecCoordPath, NB_ANGLE, DIMHOLO)


# Paths to the real, and imaginary parts of the field
CHEMIN_RE_UBORN = f"{DOSSIERDATA}Pretraitement/ReBorn_{DIMHOLO}.tiff"
CHEMIN_IM_UBORN = f"{DOSSIERDATA}Pretraitement/ImBorn_{DIMHOLO}.tiff"

# Field files reading
ReUBorn = ft.ReadtiffCube(CHEMIN_RE_UBORN)
ImUBorn = ft.ReadtiffCube(CHEMIN_IM_UBORN)
UBornCplx = ReUBorn + ImUBorn * 1j
del ReUBorn, ImUBorn

f_recon, TFVol, mask_sum = rp.retropropagation(UBornCplx, NB_HOLO, fi, FMAXHOLO,
                                               REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)

Refraction = f_recon.real
Absorption = f_recon.imag

if "DIC" == MethodUsed:
    UBornCplxDIC = rp.dic(Refraction + 1j*Absorption, 3, 3, 1, 0.1) 

if "DARKFIELD" == MethodUsed:
    UBornCplxDark = rp.DarkField(Refraction + 1j*Absorption, 5)


if "PHASECONTRAST" == MethodUsed:
    UBornCplxPC = rp.PhaseContrast(Refraction + 1j*Absorption, 5, 250)


if "RHEINBERG" == MethodUsed:
    UBornCplxR, UBornCplxG, UBornCplxB = rp.RheinbergIllumination(Refraction + 1j*Absorption, [30,120], 60, [70,250])

    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/CompoR_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                    np.abs(UBornCplxR)**2, 2*PIXTHEO*1e6)


    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTFR_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
               np.abs(fftshift(fftn(UBornCplxR))), 1./(Refraction.shape[0]*2*PIXTHEO*1e6))

    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/CompoG_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
               np.abs(UBornCplxG)**2, 2*PIXTHEO*1e6)

    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTFG_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
               np.abs(fftshift(fftn(UBornCplxG))), 1./(Refraction.shape[0]*2*PIXTHEO*1e6))


    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/CompoB_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                   np.abs(UBornCplxB)**2, 2*PIXTHEO*1e6)

    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTFB_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
               np.abs(fftshift(fftn(UBornCplxB))), 1./(Refraction.shape[0]*2*PIXTHEO*1e6))


OTF = np.zeros_like(mask_sum)
OTF[mask_sum != 0] = 1

TFVolfilt = np.zeros_like(OTF)
TFVolfilt[TFVol != 0] = 1

# Writting results
start_time = time.time()
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Refraction_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                Refraction, 2*PIXTHEO*1e6)

ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTF_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                OTF, 2*PIXTHEO*1e6)

ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Absorption_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                Absorption, 2*PIXTHEO*1e6)


print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")

# Darkfield processing
if "DARKFIELD" == MethodUsed:
    DarkF = abs(UBornCplxDark)**2
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Darkfield_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                   DarkF, 2*PIXTHEO*1e6)


if "BASE" == MethodUsed:
    Base = np.abs(Refraction + 1j*Absorption)**2
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Brightfield_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                   Base, 2*PIXTHEO*1e6)


if "PHASECONTRAST" == MethodUsed:
    PhaseC = abs(UBornCplxPC)**2
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Phasecontrast_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                    PhaseC, 2*PIXTHEO*1e6)

if "DIC" == MethodUsed:
    DIC = UBornCplxDIC
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/DIC_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                    DIC, 2*PIXTHEO*1e6)