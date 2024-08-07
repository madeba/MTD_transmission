# -*- coding: utf-8 -*-
"""
@author: Steve Laroche et Nicolas Verrier
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
NB_HOLO = 50

CHEMINSAV_RE = f"{DOSSIERDATA}Pretraitement/ReBorn_{dimHolo}.tiff"
CHEMINSAV_RER = f"{DOSSIERDATA}Pretraitement/ReBornR_{dimHolo}.tiff"
CHEMINSAV_REG = f"{DOSSIERDATA}Pretraitement/ReBornG_{dimHolo}.tiff"
CHEMINSAV_REB = f"{DOSSIERDATA}Pretraitement/ReBornB_{dimHolo}.tiff"

CHEMINSAV_IM = f"{DOSSIERDATA}Pretraitement/ImBorn_{dimHolo}.tiff"
CHEMINSAV_IMR = f"{DOSSIERDATA}Pretraitement/ImBornR_{dimHolo}.tiff"
CHEMINSAV_IMG = f"{DOSSIERDATA}Pretraitement/ImBornG_{dimHolo}.tiff"
CHEMINSAV_IMB = f"{DOSSIERDATA}Pretraitement/ImBornB_{dimHolo}.tiff"

CHEMINSUM = f"{DOSSIERDATA}Pretraitement/Sum_{dimHolo}.tiff"
CHEMINSUMR = f"{DOSSIERDATA}Pretraitement/SumR_{dimHolo}.tiff"
CHEMINSUMG = f"{DOSSIERDATA}Pretraitement/SumG_{dimHolo}.tiff"
CHEMINSUMB = f"{DOSSIERDATA}Pretraitement/SumB_{dimHolo}.tiff"

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

SumHolo = np.zeros((250,250), dtype = complex)
SumHoloR = np.zeros((250,250), dtype = complex)
SumHoloG = np.zeros((250,250), dtype = complex)
SumHoloB = np.zeros((250,250), dtype = complex)

Method = {0 : "BASE",
          1 : "DARKFIELD",
          2 : "PHASECONTRAST",
          3 : "RHEINBERG",
          4 : "DIC",
          5 : "OBLIQUE"}

MethodUsed = Method[0]

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
wreal = tf.TiffWriter(CHEMINSAV_RE)
wrealR = tf.TiffWriter(CHEMINSAV_RER)
wrealG = tf.TiffWriter(CHEMINSAV_REG)
wrealB = tf.TiffWriter(CHEMINSAV_REB)

wimag = tf.TiffWriter(CHEMINSAV_IM)
wimagR = tf.TiffWriter(CHEMINSAV_IMR)
wimagG = tf.TiffWriter(CHEMINSAV_IMG)
wimagB = tf.TiffWriter(CHEMINSAV_IMB)

wsum = tf.TiffWriter(CHEMINSUM)
wsumR = tf.TiffWriter(CHEMINSUMR)
wsumG = tf.TiffWriter(CHEMINSUMG)
wsumB = tf.TiffWriter(CHEMINSUMB)

start_time = time.time()
for hol in range(0, NB_HOLO):
    if CPT % 100 == 0:
        print(f"Hologram = {CPT} out of {NB_HOLO}")
    filename = f"{DOSSIERDATA}i{'%03d' % hol}.pgm"
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
        R,G,B = mmtd.FiltIlluRheinberg(SpectreFilt, 72,36,108,[ind[0]-SpectreFilt.shape[0]/2,
                                                               ind[1]-SpectreFilt.shape[1]/2])
            
        if "RHEINBERG" == MethodUsed:
            # Spectral filtering
            """CompoR,CompoG,CompoB = mmtd.rheinberg(SpectreFilt, 62, [[-62*np.sqrt(2)/2,62*np.sqrt(2)/2]], [[0,0]], [[62*np.sqrt(2)/2,-62*np.sqrt(2)/2]], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2])"""
            
            """CompoR,CompoG,CompoB = mmtd.rheinberg2DTech2(SpectreFilt, 20, [40,0], [0,0], [60,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/4)"""
            
            """CompoR,CompoG,CompoB = mmtd.rheinberg2DTech3(SpectreFilt, 55, [65,0], [0,0], [65,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/2, np.pi/4)"""
            
            CompoR,CompoG,CompoB = mmtd.rheinberg2DTech4(SpectreFilt, 45, [50,0], [0,0], [100,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/4)
            
            # Case of illumination filtering
            # =========================================================================================================================
            """R,G,B = mmtd.FiltIlluRheinberg(SpectreFilt, 55,30,85,[ind[0]-SpectreFilt.shape[0]/2,
                                                                  ind[1]-SpectreFilt.shape[1]/2])
            
            if R[kix,kiy] == 1:
                CompoR =  mmtd.rheinberg2DTech5(SpectreFilt, 45, [50,0], [0,0], [100,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/4, 'R')
            else :
                CompoR = SpectreFilt*0
                
            if B[kix,kiy] == 1:
                CompoB =  mmtd.rheinberg2DTech5(SpectreFilt, 45, [50,0], [0,0], [100,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/4, 'B')
            else :
                CompoB = SpectreFilt*0
                
            if G[kix,kiy] == 1:
                CompoG =  mmtd.rheinberg2DTech5(SpectreFilt, 45, [50,0], [0,0], [100,0], [ind[0]-SpectreFilt.shape[0]/2,
                                                            ind[1]-SpectreFilt.shape[1]/2], np.pi/4, 'G')
            else :
                CompoG = SpectreFilt*0
                
            """
            # =========================================================================================================================
            
            # Coordinate writting
            fidCentrestxt.write(f"{kiy} {kix}\n")

            # Position of the center for further monitoring
            Centres[kiy, kix] = 1
            
            # Complex Field (UBorn + Ui)
            UBornR = CompoR#ifft2(ifftshift(CompoR))
            UBornG = CompoG#ifft2(ifftshift(CompoG))
            UBornB = CompoB#ifft2(ifftshift(CompoB))
            
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
                
                #Sum of intensities
                SumHoloR = SumHoloR + np.abs(Re_UBornR+1j*Im_UBornR)**2
                SumHoloG = SumHoloG + np.abs(Re_UBornG+1j*Im_UBornG)**2
                SumHoloB = SumHoloB + np.abs(Re_UBornB+1j*Im_UBornB)**2
                
                wrealR.write(np.float32(Re_UBornR), contiguous=True)
                wrealG.write(np.float32(Re_UBornG), contiguous=True)
                wrealB.write(np.float32(Re_UBornB), contiguous=True)
                
                wimagR.write(np.float32(Im_UBornR), contiguous=True)
                wimagG.write(np.float32(Im_UBornG), contiguous=True)
                wimagB.write(np.float32(Im_UBornB), contiguous=True)
                
            else:
                # Born
                Re_UBornR = Amp_UBornCR*np.cos(Phase_UBornCR)
                Re_UBornG = Amp_UBornCG*np.cos(Phase_UBornCG)
                Re_UBornB = Amp_UBornCB*np.cos(Phase_UBornCB)
                
                Im_UBornR = Amp_UBornCR*np.sin(Phase_UBornCR)
                Im_UBornG = Amp_UBornCG*np.sin(Phase_UBornCG)
                Im_UBornB = Amp_UBornCB*np.sin(Phase_UBornCB)
                
                #Sum of intensities
                SumHoloR = SumHoloR + np.abs(Re_UBornR+1j*Im_UBornR)**2
                SumHoloG = SumHoloG + np.abs(Re_UBornG+1j*Im_UBornG)**2
                SumHoloB = SumHoloB + np.abs(Re_UBornB+1j*Im_UBornB)**2
                
                wrealR.write(np.float32(Re_UBornR), contiguous=True)
                wrealG.write(np.float32(Re_UBornG), contiguous=True)
                wrealB.write(np.float32(Re_UBornB), contiguous=True)
                
                wimagR.write(np.float32(Im_UBornR), contiguous=True)
                wimagG.write(np.float32(Im_UBornG), contiguous=True)
                wimagB.write(np.float32(Im_UBornB), contiguous=True)
                
            CPT += 1
            CPT_EXIST += 1
        
        else:
            if "DARKFIELD" == MethodUsed:
                # Illumination filter
                FiltDark = mmtd.FiltIlluDarkField(SpectreFilt, 0, [ind[0]-SpectreFilt.shape[0]/2,
                                                                    ind[1]-SpectreFilt.shape[1]/2])
                # Illumination filtering
                if FiltDark[kix,kiy]==1:
                    # Spectral filtering
                    SpectreFilt = mmtd.darkfield(SpectreFilt, 2, [ind[0]-SpectreFilt.shape[0]/2,
                                                                  ind[1]-SpectreFilt.shape[1]/2])
                else:    
                    SpectreFilt = SpectreFilt*(0+1j*0) 
                    
            if "PHASECONTRAST" == MethodUsed:
                # Illumination filter
                FiltDark = mmtd.FiltIlluDarkField(SpectreFilt, 0, [ind[0]-SpectreFilt.shape[0]/2,
                                                                    ind[1]-SpectreFilt.shape[1]/2])
                # Illumination filtering
                if FiltDark[kix,kiy]==1:
                    # Spectral filtering
                    SpectreFilt = mmtd.phasecontrast(SpectreFilt, 10, 125,[ind[0]-SpectreFilt.shape[0]/2,
                                                                       ind[1]-SpectreFilt.shape[1]/2])
                else:    
                    SpectreFilt = SpectreFilt*(0+1j*0)
                
                
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
                Amp_UBornC = mmtd.dic((Amp_UBornC + 1j*Phase_UBornC),3,3,0.18)
            
            # Field calculation
            if RYTOV is True:
                # Rytov
                Re_UBorn = np.log(np.abs(Amp_UBornC))
                Im_UBorn = Phase_UBornC
                
                # Sum of intensities
                SumHolo = SumHolo + np.abs(Re_UBorn+1j*Im_UBorn)**2
                
                wreal.write(np.float32(Re_UBorn), contiguous=True)
                wimag.write(np.float32(Im_UBorn), contiguous=True)
            else:
                # Born
                if "DIC" == MethodUsed:
                    Re_UBorn = Amp_UBornC*np.cos(Phase_UBornC)-1
                    Im_UBorn = Amp_UBornC*np.sin(Phase_UBornC)-1
                    
                    # Sum of intensities
                    SumHolo = SumHolo + (Re_UBorn+1j*Im_UBorn)
                else :
                    Re_UBorn = Amp_UBornC*np.cos(Phase_UBornC)
                    Im_UBorn = Amp_UBornC*np.sin(Phase_UBornC)
                    
                    #Sum of intensities
                    SumHolo = SumHolo + np.abs(Re_UBorn+1j*Im_UBorn)**2
                wreal.write(np.float32(Re_UBorn), contiguous=True)
                wimag.write(np.float32(Im_UBorn), contiguous=True)
            CPT += 1
            CPT_EXIST += 1
    else:
        CPT += 1
print(f"Pre-Processing time for {CPT_EXIST-1} holograms: "
      f"{np.round(time.time() - start_time,decimals=2)} seconds")

wsum.write(np.float32(np.abs(SumHolo)))
wsumR.write(np.float32(SumHoloR))
wsumG.write(np.float32(SumHoloG))
wsumB.write(np.float32(SumHoloB))
wreal.close()
wrealR.close()
wrealG.close()
wrealB.close()

wimag.close()
wimagR.close()
wimagG.close()
wimagB.close()
wsum.close()
wsumR.close()
wsumG.close()
wsumB.close()

# Center recording and file closing
fidParams.write(f"nb_angle {CPT_EXIST-1}\n")
fidParams.write(f"fmaxHolo {fmaxHolo}\n")
fidParams.write(f"dimHolo {dimHolo}\n")
fidParams.write(f"pixTheo {M.PIX/Gtot}\n")
tf.imwrite(CHEMINSAV_CENTRES, np.float32(np.int32(Centres)))
fidCentrestxt.close()
fidParams.close()