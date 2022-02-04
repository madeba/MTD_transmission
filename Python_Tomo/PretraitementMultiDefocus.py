# -*- coding: utf-8 -*-
"""
@author: Steve Laroche et Nicolas Verrier
"""

import time
import os
import matplotlib.pyplot as plt
from scipy.fftpack import fft2, ifft2, ifftshift, fftshift
from scipy import signal
import numpy as np
import FileTools as ft
import HoloProcessing as holo
import CorrectionAberration as CAber
import MultiModalMTD as mmtd
import tifffile as tf
import manip
import math

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
FICHIERCONFIG = M.fichier_config

# Creating results Folders
PROCESSINGFOLDER = M.dossier_pretraitement
if not os.path.exists(PROCESSINGFOLDER):
    os.makedirs(PROCESSINGFOLDER)

# Acquisition data initialisation

HOLOREF = False
RYTOV = False

CAMDIM = 1024
Gtot = M.F_TUBE/M.F_OBJ/M.RAPFOC
REwald = CAMDIM*M.PIX/Gtot*M.NIMM/(M.LAMBDA) # Ewald sphere radius (pixel)
fmaxHolo = round(REwald*M.NA/M.NIMM) # Max frequency support (pixel)
dimHolo = int(2*fmaxHolo) # Hologram size
# NB_HOLO = M.NB_HOLOTOT # Number of holograms in the sequence
NB_HOLO = 50 # Number of holograms in the sequence

CHEMINSAV_RE = f"{DOSSIERDATA}Pretraitement/ReBorn_{dimHolo}.tiff"
CHEMINSAV_RER = f"{DOSSIERDATA}Pretraitement/ReBornR_{dimHolo}.tiff"
CHEMINSAV_REG = f"{DOSSIERDATA}Pretraitement/ReBornG_{dimHolo}.tiff"
CHEMINSAV_REB = f"{DOSSIERDATA}Pretraitement/ReBornB_{dimHolo}.tiff"

CHEMINSAV_IM = f"{DOSSIERDATA}Pretraitement/ImBorn_{dimHolo}.tiff"
CHEMINSAV_IMR = f"{DOSSIERDATA}Pretraitement/ImBornR_{dimHolo}.tiff"
CHEMINSAV_IMG = f"{DOSSIERDATA}Pretraitement/ImBornG_{dimHolo}.tiff"
CHEMINSAV_IMB = f"{DOSSIERDATA}Pretraitement/ImBornB_{dimHolo}.tiff"

CHEMINSAV_CENTRES = f"{DOSSIERDATA}Pretraitement/Centres_{dimHolo}.tiff"
CHEMINSAV_CENTRESTXT = f"{DOSSIERDATA}Pretraitement/Centres_{dimHolo}.txt"
CHEMINSAV_PARAM = f"{DOSSIERDATA}Pretraitement/Param.txt"
fidCentrestxt = open(CHEMINSAV_CENTRESTXT, "a")
fidParams = open(CHEMINSAV_PARAM, "a")
fidParams.write(f"REwald {REwald}\n")

# Initialisation

Centres = np.zeros((dimHolo, dimHolo))
CentreXShift, CentreYShift = holo.CoordToCoordShift(M.CENTREX, M.CENTREY, CAMDIM, CAMDIM)

Method = {0 : "BASE",
          1 : "DARKFIELD",
          2 : "PHASECONTRAST",
          3 : "RHEINBERG",
          4 : "DIC",
          5 : "OBLIQUE"}
MethodUsed = Method[2]

Stack = np.zeros((250,250,30),dtype = float)
StackRGB = np.zeros((250,250,30,3),dtype = float)

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


Illu = np.zeros((250,250),dtype = complex)
Illu[:,fmaxHolo:] = 1

start_time = time.time()
CPTSUM = 0
for Z in range(-1200,1200,80):
    CPT = 1
    CPT_EXIST = 1
    SumHolo = np.zeros((250,250), dtype = float)
    SumHoloR = np.zeros((250,250), dtype = float)
    SumHoloG = np.zeros((250,250), dtype = float)
    SumHoloB = np.zeros((250,250), dtype = float)
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
        
            if "RHEINBERG" == MethodUsed:
                # Spectral filtering
                """CompoR,CompoG,CompoB = mmtd.rheinberg(SpectreFilt, 62, [[62*np.sqrt(2)/2,62*np.sqrt(2)/2]], [[0,0]], [[-62*np.sqrt(2)/2,-62*np.sqrt(2)/2]], [ind[0]-SpectreFilt.shape[0]/2,
                                                                ind[1]-SpectreFilt.shape[1]/2])"""
                
                """CompoR,CompoG,CompoB = mmtd.rheinberg2DTech2(SpectreFilt, 45, [50,0], [0,0], [100,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                                ind[1]-SpectreFilt.shape[1]/2], np.pi/4)"""
            
                """CompoR,CompoG,CompoB = mmtd.rheinberg2DTech3(SpectreFilt, 45, [50,0], [0,0], [50,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                                ind[1]-SpectreFilt.shape[1]/2], np.pi/2, np.pi/4)"""
            
                CompoR,CompoG,CompoB = mmtd.rheinberg2DTech4(SpectreFilt, 45, [50,0], [0,0], [100,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/4)
            
                # Coordinate writting
                fidCentrestxt.write(f"{kiy} {kix}\n")

                # Position of the center for further monitoring
                Centres[kiy, kix] = 1
                
                # Complex Field (UBorn + Ui)
                
                UBornR = CompoR#ifft2(ifftshift(CompoR))
                UBornG = CompoG#ifft2(ifftshift(CompoG))
                UBornB = CompoB#ifft2(ifftshift(CompoB))
                
                # Defocus
                UBornR = holo.reconstruction(UBornR,Z*1e-6,M.LAMBDA,M.PIX)
                UBornG = holo.reconstruction(UBornG,Z*1e-6,M.LAMBDA,M.PIX)
                UBornB = holo.reconstruction(UBornB,Z*1e-6,M.LAMBDA,M.PIX)
                
                Amp_UBornR = np.abs(UBornR)
                Amp_UBornG = np.abs(UBornG)
                Amp_UBornB = np.abs(UBornB)
                Phase_UBornWrapR = np.angle(UBornR)
                Phase_UBornWrapG = np.angle(UBornG)
                Phase_UBornWrapB = np.angle(UBornB)
                
                # Phase unwrapping
                Phase_UBornR = holo.unwrapping(Phase_UBornWrapR, M.PIX)
                Phase_UBornG = holo.unwrapping(Phase_UBornWrapG, M.PIX)
                Phase_UBornB = holo.unwrapping(Phase_UBornWrapB, M.PIX)
                
                # Amplitude correction
                Amp_UBornCR = CAber.ampliCorr(Amp_UBornR, Masque, Poly_US, Poly)
                Amp_UBornCG = CAber.ampliCorr(Amp_UBornG, Masque, Poly_US, Poly)
                Amp_UBornCB = CAber.ampliCorr(Amp_UBornB, Masque, Poly_US, Poly)
                
                Amp_UBornCR[Amp_UBornCR <= -5] = 1
                Amp_UBornCR[Amp_UBornCR >= 5] = 1
                
                Amp_UBornCG[Amp_UBornCG <= -5] = 1
                Amp_UBornCG[Amp_UBornCG >= 5] = 1
                
                Amp_UBornCB[Amp_UBornCB <= -5] = 1
                Amp_UBornCB[Amp_UBornCB >= 5] = 1
                
                # Phase correction
                Phase_UBornCR = CAber.aberCorr(Phase_UBornR, Masque, Poly_US, Poly)
                Phase_UBornCG = CAber.aberCorr(Phase_UBornG, Masque, Poly_US, Poly)
                Phase_UBornCB = CAber.aberCorr(Phase_UBornB, Masque, Poly_US, Poly)
                
                # Field calculation
                if RYTOV is True:
                    # Rytov
                    Re_UBornR = np.log(np.abs(Amp_UBornCR))
                    Re_UBornG = np.log(np.abs(Amp_UBornCG))
                    Re_UBornB = np.log(np.abs(Amp_UBornCB))
                    
                    Im_UBornR = Phase_UBornCR
                    Im_UBornG = Phase_UBornCG
                    Im_UBornB = Phase_UBornCB
                    
                    # Sum of intensities
                    SumHoloR = SumHoloR + np.abs(Re_UBornR+1j*Im_UBornR)**2
                    SumHoloG = SumHoloG + np.abs(Re_UBornG+1j*Im_UBornG)**2
                    SumHoloB = SumHoloB + np.abs(Re_UBornB+1j*Im_UBornB)**2
                    
                else:
                    # Born
                    Re_UBornR = Amp_UBornCR*np.cos(Phase_UBornCR)
                    Re_UBornG = Amp_UBornCG*np.cos(Phase_UBornCG)
                    Re_UBornB = Amp_UBornCB*np.cos(Phase_UBornCB)
                    
                    Im_UBornR = Amp_UBornCR*np.sin(Phase_UBornCR)
                    Im_UBornG = Amp_UBornCG*np.sin(Phase_UBornCG)
                    Im_UBornB = Amp_UBornCB*np.sin(Phase_UBornCB)
                    
                    # Sum of intensities
                    SumHoloR = SumHoloR + np.abs(Re_UBornR+1j*Im_UBornR)**2
                    SumHoloG = SumHoloG + np.abs(Re_UBornG+1j*Im_UBornG)**2
                    SumHoloB = SumHoloB + np.abs(Re_UBornB+1j*Im_UBornB)**2
                
                CPT += 1
                CPT_EXIST += 1
            
            else:
                if "DARKFIELD" == MethodUsed:
                    # Illumination filter
                    FiltDark = mmtd.FiltIlluDarkField(SpectreFilt, 36, [ind[0]-SpectreFilt.shape[0]/2,
                                                                    ind[1]-SpectreFilt.shape[1]/2])
                    # Illumination filtering
                    if FiltDark[kix,kiy]==1:
                        # Spectral filtering
                        SpectreFilt = mmtd.darkfield(SpectreFilt, 55, [ind[0]-SpectreFilt.shape[0]/2,
                                                                  ind[1]-SpectreFilt.shape[1]/2])
                if "PHASECONTRAST" == MethodUsed:
                    # Spectral filtering
                    SpectreFilt = mmtd.phasecontrast(SpectreFilt, 10, 100,
                                                     [ind[0]-SpectreFilt.shape[0]/2,
                                                      ind[1]-SpectreFilt.shape[1]/2])
                if "OBLIQUE" == MethodUsed:
                # Illumination filtering
                    if kiy-fmaxHolo < 0 :
                        SpectreFilt = SpectreFilt*(0+1j*0) 
                        
                # Coordinate writting
                fidCentrestxt.write(f"{kiy} {kix}\n")
                
                # Position of the center for further monitoring
                Centres[kiy, kix] = 1
                
                # Complex Field (UBorn + Ui)
                UBorn = ifft2(ifftshift(SpectreFilt))
                
                # Defocus
                UBorn = holo.reconstruction(UBorn,Z*1e-6,M.LAMBDA,M.PIX)
                
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
                
                if "DIC" == MethodUsed:
                    Amp_UBornC = mmtd.dic((Amp_UBornC + 1j*Phase_UBornC),3,3,0.1)
                
                # Field calculation
                if RYTOV is True:
                    # Rytov
                    Re_UBorn = np.log(np.abs(Amp_UBornC))
                    Im_UBorn = Phase_UBornC
                    
                    # Sum of intensities
                    SumHolo = SumHolo + np.abs(Re_UBorn+1j*Im_UBorn)**2

                else:
                    # Born
                    if "DIC" == MethodUsed:
                        Re_UBorn = Amp_UBornC*np.cos(Phase_UBornC)-1
                        Im_UBorn = Amp_UBornC*np.sin(Phase_UBornC)-1
                        
                        # Sum of intensities
                        SumHolo = SumHolo + np.abs(Re_UBorn+1j*Im_UBorn)
                    else :
                        Re_UBorn = Amp_UBornC*np.cos(Phase_UBornC)
                        Im_UBorn = Amp_UBornC*np.sin(Phase_UBornC)
                        
                        # Sum of intensities
                        SumHolo = SumHolo + np.abs(Re_UBorn+1j*Im_UBorn)**2
                
                CPT += 1
                CPT_EXIST += 1
                
        else:
            CPT += 1
        
    if "RHEINBERG" == MethodUsed:
        StackRGB[:,:,int((Z+1200)/80),0] = SumHoloR
        StackRGB[:,:,int((Z+1200)/80),1] = SumHoloG
        StackRGB[:,:,int((Z+1200)/80),2] = SumHoloB
        
    else :
        Stack[:,:,int((Z+1200)/80)] = SumHolo
    
    CPTSUM += 1
    print(f"Slice {CPTSUM} out of {Stack.shape[2]}")
print(f"Pre-Processing time for {CPTSUM} Slices: "
      f"{np.round(time.time() - start_time,decimals=2)} seconds")


if "RHEINBERG" == MethodUsed:
    ft.SAVtiffRGBCube(f"{PROCESSINGFOLDER}/{MethodUsed}_{2*dimHolo}x{2*dimHolo}x{StackRGB.shape[2]}.tiff",
                        StackRGB, 2*M.PIX*1e6)
else:
    ft.SAVtiffCube(f"{PROCESSINGFOLDER}/{MethodUsed}_{2*dimHolo}x{2*dimHolo}x{Stack.shape[2]}.tiff",
                    Stack, 2*M.PIX*1e6)

# Center recording and file closing
fidParams.write(f"nb_angle {CPT_EXIST-1}\n")
fidParams.write(f"fmaxHolo {fmaxHolo}\n")
fidParams.write(f"dimHolo {dimHolo}\n")
fidParams.write(f"pixTheo {M.PIX/Gtot}\n")
tf.imwrite(CHEMINSAV_CENTRES, np.float32(np.int32(Centres)))
fidCentrestxt.close()
fidParams.close()