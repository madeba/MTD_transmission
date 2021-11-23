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
import tifffile as tf
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
FICHIERCONFIG = M.fichier_config

# Creating results Folders
PROCESSINGFOLDER = M.dossier_pretraitement
if not os.path.exists(PROCESSINGFOLDER):
    os.makedirs(PROCESSINGFOLDER)

# Acquisition data initialisation
HOLOREF = False
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

print("----------------------------------")
print("- Creating/Removing config files -")
print("----------------------------------")
# Former files suppression and creation
if os.path.isfile(CHEMINSAV_CENTRESTXT):
    os.remove(CHEMINSAV_CENTRESTXT)
else:
    print(f"File {CHEMINSAV_CENTRESTXT} does not exist")

print(f"Creating {CHEMINSAV_CENTRESTXT}")
fidCentrestxt = open(CHEMINSAV_CENTRESTXT, "a")

if os.path.isfile(CHEMINSAV_PARAM):
    os.remove(CHEMINSAV_PARAM)
else:
    print(f"File {CHEMINSAV_PARAM} does not exist")
print(f"Creating {CHEMINSAV_PARAM}")
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
print("")
print("------------------------")
print("- Data Mask generation -")
print("------------------------")
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
wreal = tf.TiffWriter(CHEMINSAV_RE)
wimag = tf.TiffWriter(CHEMINSAV_IM)

print("")
print("--------------------------")
print("- Hologram preprocessing -")
print("--------------------------")
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
        Phase_UBorn = holo.unwrapping(Phase_UBornWrap, M.PIX)

        # Amplitude correction
        Amp_UBornC = CAber.ampliCorr(Amp_UBorn, Masque, Poly_US, Poly)
        Amp_UBornC[Amp_UBornC <= -5] = 1
        Amp_UBornC[Amp_UBornC >= 5] = 1

        # Phase correction
        Phase_UBornC = CAber.aberCorr(Phase_UBorn, Masque, Poly_US, Poly)

        # Field calculation
        if RYTOV is True:
            # Rytov
            Re_UBorn = np.log(np.abs(Amp_UBornC))
            Im_UBorn = Phase_UBornC
            wreal.write(np.float32(Re_UBorn), contiguous=True)
            wimag.write(np.float32(Im_UBorn), contiguous=True)
        else:
            # Born
            Re_UBorn = Amp_UBornC*np.cos(Phase_UBornC)-1
            Im_UBorn = Amp_UBornC*np.sin(Phase_UBornC)-1
            wreal.write(np.float32(Re_UBorn), contiguous=True)
            wimag.write(np.float32(Im_UBorn), contiguous=True)
        CPT += 1
        CPT_EXIST += 1
    else:
        CPT += 1
print(f"Pre-Processing time for {CPT_EXIST-1} holograms: "
      f"{np.round(time.time() - start_time,decimals=2)} seconds")
wreal.close()
wimag.close()

# Center recording and file closing
fidParams.write(f"nb_angle {CPT_EXIST-1}\n")
fidParams.write(f"fmaxHolo {fmaxHolo}\n")
fidParams.write(f"dimHolo {dimHolo}\n")
fidParams.write(f"pixTheo {M.PIX/Gtot}\n")
tf.imwrite(CHEMINSAV_CENTRES, np.float32(np.int32(Centres)))
fidCentrestxt.close()
fidParams.close()