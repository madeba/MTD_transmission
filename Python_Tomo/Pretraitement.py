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
import FileTools as ft
import MultiModalMTD as mmtd

# Data folders and config files
DOSSIERACQUIS = "/home/nicolas/Acquisitions/ACQUIS_pollen_PN/"
DOSSIERDATA = f"{DOSSIERACQUIS}data/"
FICHIERCONFIG = f"{DOSSIERACQUIS}config/config_manip.txt"

# Creating results Folders
PROCESSINGFOLDER = f"{DOSSIERDATA}Pretraitement"
if not os.path.exists(PROCESSINGFOLDER):
    os.makedirs(PROCESSINGFOLDER)

# Acquisition data initialisation
HOLOREF = True
RYTOV = True
DARKFIELD = False
PHASECONTRAST = True
CAMDIM = 1024
NA = float(ft.readvalue(FICHIERCONFIG, 'NA'))
NIMM = float(ft.readvalue(FICHIERCONFIG, 'N0'))
LAMBDA = float(ft.readvalue(FICHIERCONFIG, 'LAMBDA'))
F_TUBE = float(ft.readvalue(FICHIERCONFIG, 'F_TUBE')) # Tube lens focal length
F_OBJ = float(ft.readvalue(FICHIERCONFIG, 'F_OBJ')) # Microscope objective focal length
PIX = float(ft.readvalue(FICHIERCONFIG, 'TPCAM')) # Physical pixel pitch
RAPFOC = float(ft.readvalue(FICHIERCONFIG, 'RF')) # focal length ratio of the resampling lens dublet
Gtot = F_TUBE/F_OBJ/RAPFOC
REwald = CAMDIM*PIX/Gtot*NIMM/(LAMBDA) # Ewald sphere radius (pixel)
fmaxHolo = round(REwald*NA/NIMM) # Max frequency support (pixel)
dimHolo = int(2*fmaxHolo) # Hologram size
# dimHolo = int(CAMDIM/4) # Hologram size
CHEMINMASQUE = f"{DOSSIERDATA}Mask.tif"
CENTREX = int(ft.readvalue(FICHIERCONFIG, 'CIRCLE_CX')) # Pupil center in Fourier space
CENTREY = int(ft.readvalue(FICHIERCONFIG, 'CIRCLE_CY'))
NB_HOLOTOT = int(ft.readvalue(FICHIERCONFIG, 'NB_HOLO'))
NB_HOLO = NB_HOLOTOT # Number of holograms in the sequence
CHEMINSAV_RE = f"{DOSSIERDATA}Pretraitement/ReBorn_{dimHolo}.bin"
CHEMINSAV_IM = f"{DOSSIERDATA}Pretraitement/ImBorn_{dimHolo}.bin"
CHEMINSAV_CENTRES = f"{DOSSIERDATA}Pretraitement/Centres_{dimHolo}.bin"
CHEMINSAV_CENTRESTXT = f"{DOSSIERDATA}Pretraitement/Centres_{dimHolo}.txt"
CHEMINSAV_PARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
fidCentrestxt = open(CHEMINSAV_CENTRESTXT, "a")
fidParams = open(CHEMINSAV_PARAM, "a")
fidParams.write(f"REwald {REwald}\n")

# Initialisation
CPT = 1
CPT_EXIST = 1
Centres = np.zeros((dimHolo, dimHolo))
CentreXShift, CentreYShift = holo.CoordToCoordShift(CENTREX, CENTREY, CAMDIM, CAMDIM)

# Loading reference image
if HOLOREF is True:
    Refname = f"{DOSSIERDATA}Intensite_ref.pgm"
    Iref = np.sqrt(plt.imread(Refname))

# Tukey Window
TukeyWindow = np.sqrt(np.outer(signal.tukey(CAMDIM, 0.1), signal.tukey(CAMDIM, 0.1)))

# Amplitude and phase correction initialisation
Masque = CAber.InitMasque(CHEMINMASQUE, dimHolo)
NBPTOK = CAber.PixInMask(Masque)
DEGREPOLY = 4
NBCOEF = CAber.SizePoly2D(DEGREPOLY)

# Undersampled polynome
Poly_US = np.zeros((NBCOEF, NBPTOK), dtype=np.float64)
Poly_US = CAber.CalcPolyUS_xy(DEGREPOLY, Masque, Poly_US)

# Polynome to be fitted
Poly = np.zeros((NBCOEF, dimHolo*dimHolo), dtype=np.float64)
Poly = CAber.CalcPoly_xy(DEGREPOLY, Masque, Poly)

start_time = time.time()
for hol in range(0, NB_HOLO):
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
            SpectreFilt = mmtd.phasecontrast(SpectreFilt, 1, 150, [ind[0]-SpectreFilt.shape[0]/2,
                                                          ind[1]-SpectreFilt.shape[1]/2])
            # plt.imshow(np.log10(np.abs(SpectreFilt)))
            # plt.show()

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
            Phase_UBorn = holo.unwrapping(Phase_UBornWrap, PIX)
        else:
            Phase_UBorn = Phase_UBornWrap

        # Amplitude correction
        Amp_UBorn = CAber.ampliCorr(Amp_UBorn, Masque, Poly_US, Poly)

        # Phase correction
        if RYTOV is True:
            Phase_UBorn = CAber.aberCorr(Phase_UBorn, Masque, Poly_US, Poly)

        # File opening for field recording
        fidRe = open(CHEMINSAV_RE, "a")
        fidIm = open(CHEMINSAV_IM, "a")

        # Field calculation
        if RYTOV is True:
            # Rytov
            Re_UBorn = np.log(np.abs(Amp_UBorn))
            Im_UBorn = Phase_UBorn
            Re_UBorn.tofile(fidRe)
            Im_UBorn.tofile(fidIm)
        else:
            # Born
            Re_UBorn = (Amp_UBorn-1)*np.cos(Phase_UBorn)
            Im_UBorn = (Amp_UBorn-1)*np.sin(Phase_UBorn)
            Re_UBorn.tofile(fidRe)
            Im_UBorn.tofile(fidIm)
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
fidParams.write(f"pixTheo {PIX/Gtot}\n")
fidCentres = open(CHEMINSAV_CENTRES, "w")
Centres.tofile(fidCentres)
fidCentres.close()
fidCentrestxt.close()
fidParams.close()
fidRe.close()
fidIm.close()
