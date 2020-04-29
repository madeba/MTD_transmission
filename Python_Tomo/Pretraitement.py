# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy.fft as nfft
import numpy as np
import HoloProcessing as holo
import CorrectionAberration as CAber
import os
import FileTools as ft
import time
# from libtiff import TIFFfile, TIFFimage

# Dossier de Données et fichier de configuration
DossierData = '../PollenAziz/'
# DossierData = '../Test/'
# nbFichier = len(os.listdir(DossierData))
# print("Nombre de fichiers =", nbFichier)
DossierAmplitude = 'C:/Users/p1600109/Documents/Recherche/MatlabTomo/Amplitude/'
DossierPhase = 'C:/Users/p1600109/Documents/Recherche/MatlabTomo/Phase/'
FichierConfig = f"{DossierData}config_manip.txt" # Prévoir lecture des paramètres depuis le fichier texte
CheminMasque = 'Masque.tif'


# Données de l'acquisition
Rytov = True
CamDim = 1024
NA = ft.readvalue(FichierConfig,'NA')
nimm = ft.readvalue(FichierConfig,'N0')
Lambda = ft.readvalue(FichierConfig,'LAMBDA')
f_tube = ft.readvalue(FichierConfig,'F_TUBE') # focale de la lentille de tube
f_obj = ft.readvalue(FichierConfig,'F_OBJ') # focale de l'objectif
pix = ft.readvalue(FichierConfig,'TPCAM') # taille des pixels du capteur
RapFoc = ft.readvalue(FichierConfig,'RF') # rapport des focales du doublet de rééchantillonnage
Gtot = f_tube/f_obj/0.7
REwald = CamDim*pix/Gtot*nimm/(Lambda) # rayon support de fréquence accessible en pixel
fmaxHolo = round(REwald*NA/nimm) # support max de fréquence
dimHolo = int(2*fmaxHolo) # dimension de l'hologramme
CentreX = int(ft.readvalue(FichierConfig,'CIRCLE_CX')) # position du centre de la pupille dans l'espace de Fourier
CentreY = int(ft.readvalue(FichierConfig,'CIRCLE_CY'))
nb_holoTot = int(ft.readvalue(FichierConfig,'NB_HOLO'))
nb_holo = 1 # nombre d'hologrammes à traiter
CheminSAV_Re = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/ReBorn_{dimHolo}.bin"
CheminSAV_Im = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/ImBorn_{dimHolo}.bin"
CheminSAV_Centres = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/Centres_{dimHolo}.bin"

# Traitement de la séquence
cpt = 1
cpt_exist = 1
Centres = np.zeros((dimHolo,dimHolo))

# Initialisation correction amplitude et phase
Masque = CAber.InitMasque(CheminMasque,dimHolo)
nbPtOK = CAber.PixInMask(Masque)
degrePoly = 4
nbCoef = CAber.SizePoly2D(degrePoly)

# Calcul du polynome sous-échantillonné        
Poly_US = np.zeros((nbCoef,nbPtOK),dtype=np.float64)        
Poly_US = CAber.CalcPolyUS_xy(degrePoly,Masque,Poly_US)

# Calcul du polynome total
Poly = np.zeros((nbCoef,dimHolo*dimHolo),dtype=np.float64)          
Poly = CAber.CalcPoly_xy(degrePoly,Masque,Poly)

start_time = time.time()
for hol in range(0,nb_holo):
    filename = f"{DossierData}i{'%03d' % cpt}.pgm"
    if os.path.isfile(filename):
        Image = plt.imread(filename)
        
        # Spectre hologramme
        FImage = nfft.fftshift(nfft.fft2(Image))
        nfy,nfx = np.shape(FImage)
        
        # Découpe hors-axe
        Spectre = holo.filtrage(FImage,CentreX,CentreY,fmaxHolo)
        SpectreFilt=Spectre[int((nfx-dimHolo)/2):int((nfx-dimHolo)/2+dimHolo),int((nfy-dimHolo)/2):int((nfy-dimHolo)/2+dimHolo)]
        
        # Coordonnées du speculaire
        ind = np.unravel_index(np.argmax(SpectreFilt, axis=None), SpectreFilt.shape)
        kiy = ind[0]
        kix = ind[1]
        Centres[kiy,kix] = 1
        
        # Champ complexe (UBorn + Ui)
        UBorn = nfft.ifft2(nfft.ifftshift(SpectreFilt))
        Amp_UBorn = np.abs(UBorn)
        Phase_UBornWrap = np.angle(UBorn)
        
        # Depliement de la phase
        Phase_UBorn = holo.unwrapping(Phase_UBornWrap, pix)
            
        # Correction de l'amplitude
        Amp_UBornCorr = CAber.ampliCorr(Amp_UBorn,Masque,Poly_US,Poly)
        
        # Correction de la phase
        Phase_UBornCorr = CAber.aberCorr(Phase_UBorn,Masque,Poly_US,Poly)
        
        # # Enregistrement des résultats
        # CheminAmp = f"{DossierAmplitude}AmpUBorn_{'%03d' % cpt}.tiff"
        # CheminPh = f"{DossierPhase}PhaseUBorn_{'%03d' % cpt}.tiff"
        
        # # Enregistrement de l'amplitude et la phase dépliée (avant correction pour le test)
        # plt.imsave(CheminAmp,Amp_UBornCorr,cmap=plt.cm.gray)
        # plt.imsave(CheminPh,Phase_UBornCorr,cmap=plt.cm.gray)     
        
        fidRe = open(CheminSAV_Re,"a")
        fidIm = open(CheminSAV_Im,"a")
        # Calcul du Champ
        if Rytov is True:
            # print("Rytov")
            Re_UBorn = np.log(np.abs(Amp_UBornCorr))
            Im_UBorn = Phase_UBornCorr
            Re_UBorn.tofile(fidRe)
            Im_UBorn.tofile(fidIm)
        else:
            # print("Born")
            Re_UBorn = (Amp_UBornCorr-1)*np.cos(Phase_UBornCorr)
            Im_UBorn = (Amp_UBornCorr-1)*np.sin(Phase_UBornCorr)
            Re_UBorn.tofile(fidRe)
            Im_UBorn.tofile(fidIm)
        cpt += 1
        cpt_exist += 1
    else:
        cpt += 1

# Sauvegarde des centres   
fidCentres = open(CheminSAV_Centres,"w")
Centres.tofile(fidCentres)
fidCentres.close()
fidRe.close()
fidIm.close()
print(f"Temps d'éxecution : {np.round(time.time() - start_time,decimals=2)} seconds")

