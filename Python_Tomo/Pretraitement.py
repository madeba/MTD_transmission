# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import time
import os
import matplotlib.pyplot as plt
from scipy.fftpack import fft2, ifft2, ifftshift
from scipy import signal
import numpy as np
import HoloProcessing as holo
import CorrectionAberration as CAber
import MultiModalMTD as mmtd
import imageio as im
import manip

# Data folders and config files
DOSSIERACQUIS = "/home/nicolas/Acquisitions/Topi/"
DATA = True # True for data preprocessing, False for white image processing
M = manip.Manip(DOSSIERACQUIS, DATA)
if DATA is True:
    DOSSIERDATA = M.dossier_data
else:
    DOSSIERDATA = M.dossier_blanc
FICHIERCONFIG = M.fichier_config

# Creating results Folders
PROCESSINGFOLDER = M.dossier_pretraitement
if not os.path.exists(PROCESSINGFOLDER):
    os.makedirs(PROCESSINGFOLDER)

# Acquisition data initialisation
HOLOREF = True
RYTOV = True
DARKFIELD = False
PHASECONTRAST = False
CAMDIM = 1024
Gtot = M.F_TUBE/M.F_OBJ/M.RAPFOC
REwald = CAMDIM*M.PIX/Gtot*M.NIMM/(M.LAMBDA) # Ewald sphere radius (pixel)
fmaxHolo = round(REwald*M.NA/M.NIMM) # Max frequency support (pixel)
dimHolo = int(2*fmaxHolo) # Hologram size
NB_HOLO = M.NB_HOLOTOT # Number of holograms in the sequence
CHEMINSAV_RE = f"{DOSSIERDATA}Pretraitement/ReBorn_{dimHolo}.tiff"
CHEMINSAV_IM = f"{DOSSIERDATA}Pretraitement/ImBorn_{dimHolo}.tiff"
CHEMINSAV_CENTRES = f"{DOSSIERDATA}Pretraitement/Centres_{dimHolo}.tiff"
CHEMINSAV_CENTRESTXT = f"{DOSSIERDATA}Pretraitement/Centres_{dimHolo}.txt"
CHEMINSAV_PARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
fidCentrestxt = open(CHEMINSAV_CENTRESTXT, "a")
fidParams = open(CHEMINSAV_PARAM, "a")
fidParams.write(f"REwald {REwald}\n")

# Initialisation
CPT = 1
CPT_EXIST = 1
Centres = np.zeros((dimHolo, dimHolo))
CentreXShift, CentreYShift = holo.CoordToCoordShift(M.CENTREX, M.CENTREY, CAMDIM, CAMDIM)

# Loading reference image
if HOLOREF is True:
    Refname = f"{DOSSIERDATA}Intensite_ref.pgm"
    Iref = np.sqrt(plt.imread(Refname))

# Tukey Window
TukeyWindow = np.sqrt(np.outer(signal.tukey(CAMDIM, 0.1),
                               signal.tukey(CAMDIM, 0.1)))

# Amplitude and phase correction initialisation
Masque = CAber.InitMasque(M.CHEMINMASQUE, dimHolo)
NBPTOK = CAber.PixInMask(Masque)
DEGREPOLY = 4
NBCOEF = CAber.SizePoly2D(DEGREPOLY)

# Undersampled polynome
Poly_US = np.zeros((NBCOEF, NBPTOK), dtype=np.float64)
Poly_US = CAber.CalcPolyUS_xy(DEGREPOLY, Masque, Poly_US)

# Polynome to be fitted
Poly = np.zeros((NBCOEF, dimHolo*dimHolo), dtype=np.float64)
Poly = CAber.CalcPoly_xy(DEGREPOLY, Masque, Poly)

# File opening for field recording
wreal = im.get_writer(CHEMINSAV_RE)
wimag = im.get_writer(CHEMINSAV_IM)

start_time = time.time()
for hol in range(0, NB_HOLO):
    if CPT % 100 == 0:
        print(f"Hologram = {CPT} out of {NB_HOLO}")
    filename = f"{DOSSIERDATA}i{'%03d' % CPT}.pgm"
    if os.path.isfile(filename):
        if HOLOREF is True:
            Image = (plt.imread(filename))/Iref
        else:
            Image = plt.imread(filename)

        # Hologram spectrum and off-axis filtering
        FImage = fft2(TukeyWindow*Image)
        SpectreFilt = FImage[int(CentreYShift-dimHolo/2):int(CentreYShift+dimHolo/2),
                             int(CentreXShift-dimHolo/2):int(CentreXShift+dimHolo/2)]

        # Specular spot coordinates calculation
        ind = np.unravel_index(np.argmax(np.abs(SpectreFilt), axis=None), SpectreFilt.shape)
        kiy = ind[0]
        kix = ind[1]
        if DARKFIELD is True:
            SpectreFilt = mmtd.darkfield(SpectreFilt, 1, [ind[0]-SpectreFilt.shape[0]/2,
                                                          ind[1]-SpectreFilt.shape[1]/2])
        if PHASECONTRAST is True:
            SpectreFilt = mmtd.phasecontrast(SpectreFilt, 20, 50,
                                             [ind[0]-SpectreFilt.shape[0]/2,
                                              ind[1]-SpectreFilt.shape[1]/2])

        # Coordinate writting
        fidCentrestxt.write(f"{kiy} {kix}\n")

        # Position of the center for further monitoring
        Centres[kiy, kix] = 1

        # Complex Field (UBorn + Ui)
        UBorn = ifft2(ifftshift(SpectreFilt))
        Amp_UBorn = np.abs(UBorn)
        Phase_UBornWrap = np.angle(UBorn)

        # Phase unwrapping
        if RYTOV is True:
            Phase_UBorn = holo.unwrapping(Phase_UBornWrap, M.PIX)
        else:
            Phase_UBorn = Phase_UBornWrap

        # Amplitude correction
        Amp_UBorn = CAber.ampliCorr(Amp_UBorn, Masque, Poly_US, Poly)

        # Phase correction
        if RYTOV is True:
            Phase_UBorn = CAber.aberCorr(Phase_UBorn, Masque, Poly_US, Poly)

        # Field calculation
        if RYTOV is True:
            # Rytov
            Re_UBorn = np.log(np.abs(Amp_UBorn))
            Im_UBorn = Phase_UBorn
            wreal.append_data(np.float32(Re_UBorn))
            wimag.append_data(np.float32(Im_UBorn))
        else:
            # Born
            Re_UBorn = (Amp_UBorn-1)*np.cos(Phase_UBorn)
            Im_UBorn = (Amp_UBorn-1)*np.sin(Phase_UBorn)
            wreal.append_data(np.float32(Re_UBorn))
            wimag.append_data(np.float32(Im_UBorn))
        CPT += 1
        CPT_EXIST += 1
    else:
        CPT += 1
print(f"Pre-Processing time for {CPT_EXIST-1} holograms: "
      f"{np.round(time.time() - start_time,decimals=2)} seconds")

# Center recording and file closing
fidParams.write(f"nb_angle {CPT_EXIST-1}\n")
fidParams.write(f"fmaxHolo {fmaxHolo}\n")
fidParams.write(f"dimHolo {dimHolo}\n")
fidParams.write(f"pixTheo {M.PIX/Gtot}\n")
fidCentres = im.get_writer(CHEMINSAV_CENTRES)
fidCentres.append_data(np.int32(Centres))
fidCentres.close()
fidCentrestxt.close()
fidParams.close()
wreal.close()
wimag.close()
