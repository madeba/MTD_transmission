# -*- coding: utf-8 -*-
"""

@author: Nicolas Verrier
"""
import SimuTomo as st
# import matplotlib.pyplot as plt
import time
import numpy as np
import FileTools as ft
from scipy.fftpack import fftn,ifftn,fftshift
import Retropropagation as rp

NA_ill = 1.4
nimm = 1.518
nbead = 1.55
nbangle = 50
kappa = 0.01
Radius = 20
DOSSIERDATA = "/home/nicolas/Simulations/"
CHEMINPARAM = f"{DOSSIERDATA}Param.txt"

# Path to the parameter file, and parameters reading
DARKFIELD = False
PHASECONTRAST = False
REWALD = float(ft.readvalue(CHEMINPARAM, 'REwald'))
NB_ANGLE = int(ft.readvalue(CHEMINPARAM, 'nb_angle'))
FMAXHOLO = int(ft.readvalue(CHEMINPARAM, 'fmaxHolo'))
DIMHOLO = int(ft.readvalue(CHEMINPARAM, 'dimHolo'))
PIXTHEO = float(ft.readvalue(CHEMINPARAM, 'pixTheo'))
UBornPitch = 1/(2*FMAXHOLO*PIXTHEO)
NB_HOLO = NB_ANGLE

# Paths to the real, and imaginary parts of the field
CHEMIN_RE_UBORN = f"{DOSSIERDATA}ReBorn_{DIMHOLO}.tiff"
CHEMIN_IM_UBORN = f"{DOSSIERDATA}ImBorn_{DIMHOLO}.tiff"

# Path to the specular coordinates
SpecCoordPath = f"{DOSSIERDATA}/Centres_{DIMHOLO}.txt"
fi = rp.Calc_fi(SpecCoordPath, NB_ANGLE, DIMHOLO)

total_time = time.time()
start_time = time.time()
if NA_ill/nimm>=1:
    print(f"Simulation impossible: nimm >= {NA_ill}")
else:
    start_time = time.time()
    Bille = st.BeadSimu(Radius, DIMHOLO, nimm, nbead, kappa)
    print(f"Bead simulation time : {np.round(time.time() - start_time, decimals=2)} seconds")

start_time = time.time()
TomoSpectrum = (fftn(Bille))
print(f"3D-FFT computation time : {np.round(time.time() - start_time, decimals=2)} seconds")
fdm, sdzm, dxm, dym = rp.Calc_fd(FMAXHOLO, REWALD)

HOLO = np.zeros((DIMHOLO, DIMHOLO, NB_HOLO), dtype = np.complex128)
start_time = time.time()
for cpt in range(NB_HOLO):
    if cpt % 100 == 0:
        print(f"Hologram = {cpt} out of {NB_HOLO}")
    TF_Holo = np.zeros((DIMHOLO, DIMHOLO), dtype = np.complex128)
    k_inc = np.array([fi[0, cpt], fi[1, cpt]])
    TF_Holo = st.Calc_TF_Holo(TomoSpectrum, TF_Holo, fdm, sdzm, dxm, dym, FMAXHOLO, k_inc, REWALD, 632.8e-9, nimm)
    decal = np.array([-k_inc[0]+DIMHOLO/2, -k_inc[1]+DIMHOLO/2]).astype(int)
    TF_Holo_r = rp.decal_TF_holo(TF_Holo, decal, TF_Holo.shape[0], roll=True)
    HOLO[:, :, cpt] = (ifftn(TF_Holo_r))/(UBornPitch**2*FMAXHOLO**2)
    # plt.imshow(np.abs(HOLO), cmap="gray")
    # plt.show()
print(f"Simulation time for {NB_HOLO} holograms: {np.round(time.time() - start_time, decimals=2)} seconds")
print("")
# ft.SAVtiffCube(CHEMIN_RE_UBORN, HOLO, PIXTHEO*1e6)
ft.SAVtiffCube(CHEMIN_RE_UBORN, HOLO.real, PIXTHEO*1e6)
ft.SAVtiffCube(CHEMIN_IM_UBORN, HOLO.imag, PIXTHEO*1e6)
print(f"Execution time: {np.round(time.time() - total_time, decimals=2)} seconds")
