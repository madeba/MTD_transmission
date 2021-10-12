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
GreenRadius = 20
RedRadius = 80

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

start_time = time.time()
print("--------------------------------------")
print("- Reconstruction of the Green Channel-")
print("--------------------------------------")
f_reconG, TFVolG, mask_sum = rp.retropropagation(UBornGreen, UBornGreen.shape[2], fiG, FMAXHOLO,
                                                REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

RefractionG = f_reconG.real
AbsorptionG = f_reconG.imag
OTFG = np.zeros_like(mask_sum)
OTFG[mask_sum != 0] = 1
IntensiteG = np.abs(RefractionG+1j*AbsorptionG)**2

# Writting results
start_time = time.time()
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/RefractionG_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 RefractionG, 2*PIXTHEO*1e6)
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/AbsorptionG_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 AbsorptionG, 2*PIXTHEO*1e6)
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/IntensityG_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                IntensiteG, 2*PIXTHEO*1e6)
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTFG_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 OTFG, 1./(RefractionG.shape[0]*2*PIXTHEO*1e6))
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")
del RefractionG, AbsorptionG

start_time = time.time()
print("--------------------------------------")
print("- Reconstruction of the Red Channel-")
print("--------------------------------------")
f_reconR, TFVolR, mask_sum = rp.retropropagation(UBornRed, UBornRed.shape[2], fiR, FMAXHOLO,
                                                REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

RefractionR = f_reconR.real
AbsorptionR = f_reconR.imag
OTFR = np.zeros_like(mask_sum)
OTFR[mask_sum != 0] = 1
IntensiteR = np.abs(RefractionR+1j*AbsorptionR)**2

# Writting results
start_time = time.time()
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/RefractionR_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 RefractionR, 2*PIXTHEO*1e6)
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/AbsorptionR_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 AbsorptionR, 2*PIXTHEO*1e6)
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/IntensityR_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                IntensiteR, 2*PIXTHEO*1e6)
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTFR_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 OTFR, 1./(RefractionG.shape[0]*2*PIXTHEO*1e6))
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")
del RefractionR, AbsorptionR

start_time = time.time()
print("--------------------------------------")
print("- Reconstruction of the Blue Channel-")
print("--------------------------------------")
f_reconB, TFVolB, mask_sum = rp.retropropagation(UBornBlue, UBornBlue.shape[2], fiB, FMAXHOLO,
                                                REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
print(f"Reconstruction time for a {2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO} volume (3D-FFT included), "
      f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
print("")

RefractionB = f_reconB.real
AbsorptionB = f_reconB.imag
OTFB = np.zeros_like(mask_sum)
OTFB[mask_sum != 0] = 1
IntensiteB = np.abs(RefractionB+1j*AbsorptionB)**2

# Writting results
start_time = time.time()
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/RefractionB_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 RefractionB, 2*PIXTHEO*1e6)
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/AbsorptionB_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
#                 AbsorptionB, 2*PIXTHEO*1e6)
ft.SAVtiffCube(f"{PROCESSINGFOLDER}/IntensityB_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                IntensiteB, 2*PIXTHEO*1e6)
# ft.SAVtiffCube(f"{PROCESSINGFOLDER}/OTFB_{2*DIMHOLO}x{2*DIMHOLO}x{2*DIMHOLO}.tiff",
                # OTFB, 1./(RefractionB.shape[0]*2*PIXTHEO*1e6))
print(f"Data saving: {np.round(time.time() - start_time,decimals=2)} seconds")
del RefractionB, AbsorptionB

viewer = napari.view_image(IntensiteB.transpose(-1, 1, 0), name='Bleu', colormap='blue')
viewer.add_image(IntensiteR.transpose(-1, 1, 0), name='Rouge',colormap='red')
viewer.add_image(IntensiteG.transpose(-1, 1, 0), name='Vert',colormap='green')