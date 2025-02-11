# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import FileTools as ft
import manip
import napari

# Path to the parameter file, and parameters reading
DOSSIERACQUIS = "/home/nicolas/Acquisitions/PETIA/PLA_45678/"
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

GerchData = ft.ReadtiffCube(f"{GERCHBERGFOLDER}/indiceGPS.tif")

plt.imshow(GerchData[:, :, 256], cmap="gray")
plt.show()

Rec_ObjectAVG = np.zeros((1,490))
for cpt_avg in range (DIMHOLO):
    Rec_ObjectAVG = Rec_ObjectAVG + GerchData[cpt_avg,10:500,256]
    cpt_avg += 1
Rec_ObjectAVG /= cpt_avg

x = np.linspace(0,489,490)*2*PIXTHEO

plt.plot(x*1e6,Rec_ObjectAVG[0,:])
plt.xlabel(r"$x$ (µm)", fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel(r"$\Delta n$", fontsize=16)
plt.yticks(fontsize=16)
plt.ylim(top=0.015)
plt.savefig("PETIA_45678.svg")
plt.show()

# Visualization
viewer = napari.view_image(GerchData.transpose(-1, 1, 0), name='Refraction', colormap='magma')