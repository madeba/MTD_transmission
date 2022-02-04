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

# Filter Radii
GreenRadius = 60
RedRadius = 100

# Rounding tomographic volume dimensions to the next power of 2
pow2 = ft.NextPow2(2*DIMHOLO)
DIMTOMO = 2**pow2


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

UBornCplxGreen = np.zeros_like(UBornCplx)
UBornCplxRed = np.zeros_like(UBornCplx)
UBornCplxBlue = np.zeros_like(UBornCplx)

fiGreen = np.zeros_like(fi)
fiRed = np.zeros_like(fi)
fiBlue = np.zeros_like(fi)

TabFilt=np.zeros((2,NB_HOLO))
cpt = 0
for FiltCount in range(NB_HOLO):
    if fi[1,FiltCount]**2 + fi[0,FiltCount]**2 < GreenRadius**2:
        UBornCplxGreen[:,:,FiltCount] = UBornCplx[:,:,FiltCount]
        fiGreen[:,FiltCount] = fi[:,FiltCount]
        TabFilt[0,FiltCount] = cpt
        TabFilt[1,FiltCount] = 0
    elif fi[1,FiltCount]**2 + fi[0,FiltCount]**2 < RedRadius**2:
        UBornCplxRed[:,:,FiltCount] = UBornCplx[:,:,FiltCount]
        fiRed[:,FiltCount] = fi[:,FiltCount]
        TabFilt[0,FiltCount] = cpt
        TabFilt[1,FiltCount] = 1
    else:
        UBornCplxBlue[:,:,FiltCount] = UBornCplx[:,:,FiltCount]
        fiBlue[:,FiltCount] = fi[:,FiltCount]
        TabFilt[0,FiltCount] = cpt
        TabFilt[1,FiltCount] = 2
    cpt+=1

UBornGreen = UBornCplxGreen[UBornCplxGreen != 0]
UBornGreen = UBornGreen.reshape(DIMHOLO,DIMHOLO,int((UBornGreen.shape[0])/DIMHOLO**2))
# fiG = fiGreen[fiGreen != 0]
fiG = fi[:,TabFilt[1,:] == 0]

UBornRed = UBornCplxRed[UBornCplxRed != 0]
UBornRed = UBornRed.reshape(DIMHOLO,DIMHOLO,int((UBornRed.shape[0])/DIMHOLO**2))
fiR = fi[:,TabFilt[1,:] == 1]

UBornBlue = UBornCplxBlue[UBornCplxBlue != 0]
UBornBlue = UBornBlue.reshape(DIMHOLO,DIMHOLO,int((UBornBlue.shape[0])/DIMHOLO**2))
fiB = fi[:,TabFilt[1,:] == 2]
del UBornCplx, UBornCplxGreen, UBornCplxRed, UBornCplxBlue

# Building composite
IntensiteRhein = np.zeros((DIMTOMO, DIMTOMO, DIMTOMO, 3))
OTFRhein = np.zeros((DIMTOMO, DIMTOMO, DIMTOMO, 3))

start_time = time.time()
print("---------------------------------------")
print("- Reconstruction of the Green Channel -")
print("---------------------------------------")
f_reconG, TFVolG, mask_sum = rp.retropropagation(UBornGreen, UBornGreen.shape[2], fiG, FMAXHOLO,
                                                REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {DIMTOMO}x{DIMTOMO}x{DIMTOMO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

RefractionG = f_reconG.real
AbsorptionG = f_reconG.imag
OTFG = np.zeros_like(mask_sum)
OTFG[mask_sum != 0] = 1
IntensiteG = np.abs(RefractionG+1j*AbsorptionG)**2
IntensiteRhein[:,:,:,1] = IntensiteG
OTFRhein[:,:,:,1] = OTFG
del RefractionG, AbsorptionG, OTFG, IntensiteG

start_time = time.time()
print("-------------------------------------")
print("- Reconstruction of the Red Channel -")
print("-------------------------------------")
f_reconR, TFVolR, mask_sum = rp.retropropagation(UBornRed, UBornRed.shape[2], fiR, FMAXHOLO,
                                                REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {DIMTOMO}x{DIMTOMO}x{DIMTOMO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

RefractionR = f_reconR.real
AbsorptionR = f_reconR.imag
OTFR = np.zeros_like(mask_sum)
OTFR[mask_sum != 0] = 1
IntensiteR = np.abs(RefractionR+1j*AbsorptionR)**2
IntensiteRhein[:,:,:,0] = IntensiteR
OTFRhein[:,:,:,0] = OTFR
del RefractionR, AbsorptionR, OTFR, IntensiteR

start_time = time.time()
print("--------------------------------------")
print("- Reconstruction of the Blue Channel -")
print("--------------------------------------")
f_reconB, TFVolB, mask_sum = rp.retropropagation(UBornBlue, UBornBlue.shape[2], fiB, FMAXHOLO,
                                                REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {DIMTOMO}x{DIMTOMO}x{DIMTOMO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

RefractionB = f_reconB.real
AbsorptionB = f_reconB.imag
OTFB = np.zeros_like(mask_sum)
OTFB[mask_sum != 0] = 1
IntensiteB = np.abs(RefractionB+1j*AbsorptionB)**2
IntensiteRhein[:,:,:,2] = IntensiteB
OTFRhein[:,:,:,2] = OTFB
del RefractionB, AbsorptionB,OTFB, IntensiteB

# Writting results
start_time = time.time()
ft.SAVtiffRGBCube(f"{PROCESSINGFOLDER}/Rheinberg_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                  IntensiteRhein, 2*PIXTHEO*1e6)
ft.SAVtiffRGBCube(f"{PROCESSINGFOLDER}/OTFRheinberg_{DIMTOMO}x{DIMTOMO}x{DIMTOMO}.tiff",
                  OTFRhein, 1./(OTFRhein.shape[0]*2*PIXTHEO*1e6))
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")

viewer = napari.view_image(IntensiteRhein[:,:,:,2].transpose(-1, 1, 0), name='Bleu', colormap='blue')
viewer.add_image(IntensiteRhein[:,:,:,1].transpose(-1, 1, 0), name='Rouge',colormap='red')
viewer.add_image(IntensiteRhein[:,:,:,0].transpose(-1, 1, 0), name='Vert',colormap='green')