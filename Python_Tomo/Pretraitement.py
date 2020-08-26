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
DossierAcquis = "/home/nicolas/Acquisitions/ACQUIS_pollen_PN/"
DossierData = f"{DossierAcquis}data/"
FichierConfig = f"{DossierAcquis}config/config_manip.txt"

# Creating results Folders
ProcessingFolder = f"{DossierData}Pretraitement"
if not os.path.exists(ProcessingFolder):
    os.makedirs(ProcessingFolder)

# Acquisition data initialisation
HoloRef = True
Rytov = True
DarkField = False
CamDim = 1024
NA = ft.readvalue(FichierConfig, 'NA')
nimm = ft.readvalue(FichierConfig, 'N0')
Lambda = ft.readvalue(FichierConfig, 'LAMBDA')
f_tube = ft.readvalue(FichierConfig, 'F_TUBE') # Tube lens focal length
f_obj = ft.readvalue(FichierConfig, 'F_OBJ') # Microscope objective focal length
pix = ft.readvalue(FichierConfig, 'TPCAM') # Physical pixel pitch
RapFoc = ft.readvalue(FichierConfig, 'RF') # focal length ratio of the resampling lens dublet
Gtot = f_tube/f_obj/RapFoc
REwald = CamDim*pix/Gtot*nimm/(Lambda) # Ewald sphere radius (pixel)
fmaxHolo = round(REwald*NA/nimm) # Max frequency support (pixel)
dimHolo = int(2*fmaxHolo) # Hologram size
# dimHolo = int(CamDim/4) # Hologram size
CheminMasque = f"{DossierData}Maske.tif"
CentreX = int(ft.readvalue(FichierConfig, 'CIRCLE_CX')) # Pupil center in Fourier space
CentreY = int(ft.readvalue(FichierConfig, 'CIRCLE_CY'))
nb_holoTot = int(ft.readvalue(FichierConfig, 'NB_HOLO'))
nb_holo = nb_holoTot # Number of holograms in the sequence
# nb_holo = 50 # Number of holograms in the sequence
CheminSAV_Re = f"{DossierData}Pretraitement/ReBorn_{dimHolo}.bin"
CheminSAV_Im = f"{DossierData}Pretraitement/ImBorn_{dimHolo}.bin"
CheminSAV_Centres = f"{DossierData}Pretraitement/Centres_{dimHolo}.bin"
CheminSAV_Centrestxt = f"{DossierData}Pretraitement/Centres_{dimHolo}.txt"
CheminSAV_Param = f"{DossierData}Pretraitement/Param.txt"
fidCentrestxt = open(CheminSAV_Centrestxt, "a")
fidParams = open(CheminSAV_Param, "a")
fidParams.write(f"REwald {REwald}\n")

# Initialisation
cpt = 1
cpt_exist = 1
Centres = np.zeros((dimHolo, dimHolo))
CentreXShift, CentreYShift = holo.CoordToCoordShift(CentreX, CentreY, CamDim, CamDim)

# Loading reference image
if HoloRef is True:
    Refname = f"{DossierData}Intensite_ref.pgm"
    Iref = np.sqrt(plt.imread(Refname))

# Tukey Window
TukeyWindow = np.sqrt(np.outer(signal.tukey(CamDim, 0.1), signal.tukey(CamDim, 0.1)))

# Amplitude and phase correction initialisation
Masque = CAber.InitMasque(CheminMasque, dimHolo)
nbPtOK = CAber.PixInMask(Masque)
degrePoly = 4
nbCoef = CAber.SizePoly2D(degrePoly)

# Undersampled polynome
Poly_US = np.zeros((nbCoef, nbPtOK), dtype=np.float64)
Poly_US = CAber.CalcPolyUS_xy(degrePoly, Masque, Poly_US)

# Polynome to be fitted
Poly = np.zeros((nbCoef, dimHolo*dimHolo), dtype=np.float64)
Poly = CAber.CalcPoly_xy(degrePoly, Masque, Poly)

start_time = time.time()
for hol in range(0, nb_holo):
    filename = f"{DossierData}i{'%03d' % cpt}.pgm"
    if os.path.isfile(filename):
        if HoloRef is True:
            Image = (plt.imread(filename))/Iref
        else:
            Image = plt.imread(filename)

        # Hologram spectrum and off-axis filtering
        FImage = fft2(TukeyWindow*Image)
        SpectreFilt = FImage[int(CentreYShift-dimHolo/2):int(CentreYShift+dimHolo/2), int(CentreXShift-dimHolo/2):int(CentreXShift+dimHolo/2)]

        # Specular spot coordinates calculation
        ind = np.unravel_index(np.argmax(np.abs(SpectreFilt), axis=None), SpectreFilt.shape)
        kiy = ind[0]
        kix = ind[1]
        if DarkField is True:
            SpectreFilt = mmtd.darkfield(SpectreFilt, 1, [ind[0]-SpectreFilt.shape[0]/2, ind[1]-SpectreFilt.shape[1]/2])
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
        if Rytov is True:
            Phase_UBorn = holo.unwrapping(Phase_UBornWrap, pix)
        else:
            Phase_UBorn = Phase_UBornWrap

        # Amplitude correction
        Amp_UBorn = CAber.ampliCorr(Amp_UBorn, Masque, Poly_US, Poly)

        # Phase correction
        if Rytov is True:
            Phase_UBorn = CAber.aberCorr(Phase_UBorn, Masque, Poly_US, Poly)

        # File opening for field recording
        fidRe = open(CheminSAV_Re, "a")
        fidIm = open(CheminSAV_Im, "a")

        # Field calculation
        if Rytov is True:
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
        cpt += 1
        cpt_exist += 1
    else:
        cpt += 1
print(f"Pre-Processing time for {cpt_exist-1} holograms: {np.round(time.time() - start_time,decimals=2)} seconds")

# Center recording and file closing
fidParams.write(f"nb_angle {cpt_exist-1}\n")
fidParams.write(f"fmaxHolo {fmaxHolo}\n")
fidParams.write(f"dimHolo {dimHolo}\n")
fidParams.write(f"pixTheo {pix/Gtot}\n")
fidCentres = open(CheminSAV_Centres, "w")
Centres.tofile(fidCentres)
fidCentres.close()
fidCentrestxt.close()
fidParams.close()
fidRe.close()
fidIm.close()
