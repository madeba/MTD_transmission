# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import time
import os
import numpy as np
import FileTools as ft
import Retropropagation as rp
import manip
import tifffile as tf

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

# Paths to the saved cuts
CHEMINSAV_XY = f"{PROCESSINGFOLDER}/XY_CUT.tiff"
CHEMINSAV_YZ = f"{PROCESSINGFOLDER}/YZ_CUT.tiff"
CHEMINSAVOTF_XY = f"{PROCESSINGFOLDER}/XY_CUT_OTF.tiff"
CHEMINSAVOTF_YZ = f"{PROCESSINGFOLDER}/YZ_CUT_OTF.tiff"
wxy = tf.TiffWriter(CHEMINSAV_XY)
wyz = tf.TiffWriter(CHEMINSAV_YZ)
wxyOTF = tf.TiffWriter(CHEMINSAVOTF_XY)
wyzOTF = tf.TiffWriter(CHEMINSAVOTF_YZ)

# Path to the specular coordinates
SpecCoordPath = f"{DOSSIERDATA}Pretraitement/Centres_{DIMHOLO}.txt"
fi = rp.Calc_fi(SpecCoordPath, NB_ANGLE, DIMHOLO)

# Field files reading
ReUBorn = ft.ReadtiffCube(CHEMIN_RE_UBORN)
ImUBorn = ft.ReadtiffCube(CHEMIN_IM_UBORN)
UBornCplx = ReUBorn + ImUBorn * 1j
del ReUBorn, ImUBorn

# Rounding tomographic volume dimensions to the next power of 2
pow2 = ft.NextPow2(2*DIMHOLO)
DIMTOMO = 2**pow2
for cptview in range(1,NB_ANGLE):
    start_time = time.time()
    f_recon, TFVol, mask_sum = rp.retropropagation(UBornCplx, cptview, fi, FMAXHOLO,
                                                   REWALD, M.LAMBDA, M.NIMM, PIXTHEO, UBornPitch)
    print(f"Reconstruction time for a {DIMTOMO}x{DIMTOMO}x{DIMTOMO} volume (3D-FFT included), "
          f"with {NB_HOLO} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")
    print("")
    OTF = np.zeros_like(mask_sum)
    OTF[mask_sum != 0] = 1
    
    wxy.write(f_recon.real[:,:,int(DIMTOMO/2)-12],contiguous=True)
    wyz.write(f_recon.real[int(DIMTOMO/2)-12,:,:],contiguous=True)
    wxyOTF.write(OTF[:,:,int(DIMTOMO/2)],contiguous=True)
    wyzOTF.write(OTF[int(DIMTOMO/2),:,:],contiguous=True)
    del f_recon, TFVol, mask_sum, OTF
    
wxy.close()
wyz.close()
wxyOTF.close()
wyzOTF.close()