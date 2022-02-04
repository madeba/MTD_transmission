# -*- coding: utf-8 -*-
"""
@author: Steve Laroche et Nicolas Verrier
"""

import time
import os
import numpy as np
import FileTools as ft
import Retropropagation as rp
import manip

# Data folders and config files
if os.name == 'nt': # Windows
    DOSSIERACQUIS = "C:/Users/p1600109/Documents/Recherche/Acquisitions/Topi_pollen_600U/"
else:               # Linux
    DOSSIERACQUIS = "/home/nicolas/Acquisitions/Topi_pollen_600U/"

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
MethodUsed = Method[3]

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

# Rounding tomographic volume dimensions to the next power of 2
pow2 = ft.NextPow2(2*DIMHOLO)
DIMTOMO = 2**pow2


f_recon, TFVol, mask_sum = rp.retropropagation(UBornCplx, NB_HOLO, fi, FMAXHOLO,
                                               REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)

Refraction = f_recon.real
Absorption = f_recon.imag

if "DIC" == MethodUsed:
    print("------------------")
    print("- DIC simulation -")
    print("------------------")
    start_time = time.time()
    UBornCplxDIC = rp.dic(Refraction + 1j*Absorption, 3, 3, 1, 0.1) 
    print(f"DIC processing time: {np.round(time.time() - start_time,decimals=2)} seconds")

if "DARKFIELD" == MethodUsed:
    print("------------------------")
    print("- Darkfield simulation -")
    print("------------------------")
    start_time = time.time()
    UBornCplxDark = rp.DarkField(Refraction + 1j*Absorption, 5)
    print(f"Darkfield processing time: {np.round(time.time() - start_time,decimals=2)} seconds")

if "PHASECONTRAST" == MethodUsed:
    print("------------------------------")
    print("- Phase-Constrast simulation -")
    print("------------------------------")
    start_time = time.time()
    UBornCplxPC = rp.PhaseContrast(Refraction + 1j*Absorption, 5, 250)
    print(f"Phase-Contrast processing time: {np.round(time.time() - start_time,decimals=2)} seconds")

if "RHEINBERG" == MethodUsed:
    print("----------------------------------")
    print("- Rheinberg Constrast simulation -")
    print("----------------------------------")
    start_time = time.time()
    UBornCplxR, UBornCplxG, UBornCplxB = rp.RheinbergIllumination(Refraction + 1j*Absorption, [30,120], 60, [70,250])
    IntensiteRhein = np.zeros((DIMTOMO,DIMTOMO,DIMTOMO,3))
    IntensiteRhein[:,:,:,0] = np.abs(UBornCplxR)**2
    IntensiteRhein[:,:,:,1] = np.abs(UBornCplxG)**2
    IntensiteRhein[:,:,:,2] = np.abs(UBornCplxB)**2
    print(f"Rheinberg Contrast processing time: {np.round(time.time() - start_time,decimals=2)} seconds")
    del UBornCplxR, UBornCplxG, UBornCplxB
    
    start_time = time.time()
    ft.SAVtiffRGBCube(f"{PROCESSINGFOLDER}/IntensiteRheinSpec_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                    IntensiteRhein, 2*PIXTHEO*1e6)
    print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")

# Darkfield processing
if "DARKFIELD" == MethodUsed:
    DarkF = abs(UBornCplxDark)**2
    start_time = time.time()
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Darkfield_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                   DarkF, 2*PIXTHEO*1e6)
    print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")


if "BASE" == MethodUsed:
    Base = np.abs(Refraction + 1j*Absorption)**2
    start_time = time.time()
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Brightfield_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                   Base, 2*PIXTHEO*1e6)
    print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")


if "PHASECONTRAST" == MethodUsed:
    PhaseC = abs(UBornCplxPC)**2
    start_time = time.time()
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/Phasecontrast_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                    PhaseC, 2*PIXTHEO*1e6)
    print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")

if "DIC" == MethodUsed:
    DIC = UBornCplxDIC
    start_time = time.time()
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/DIC_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                    DIC, 2*PIXTHEO*1e6)
    print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")