import matplotlib.pyplot as plt
import numpy.fft as nfft
# import pyfftw
import numpy as np
import HoloProcessing as holo
import CorrectionAberration as CAber
import os
import time
import cv2

# Dossier de Données et fichier de configuration
DossierData = '../PollenAziz/'
DossierAmplitude = 'C:/Users/p1600109/Documents/Recherche/MatlabTomo/Amplitude/'
DossierPhase = 'C:/Users/p1600109/Documents/Recherche/MatlabTomo/Phase/'
FichierConfig = DossierData + 'config_manip.txt' # Prévoir lecture des paramètres depuis le fichier texte
CheminMasque = 'Masque.tif'

# Données de l'acquisition (Bientôt lues depuis le fichier de données)
NA = 1.4
nimm = 1.515
Lambda = 632.8e-9
f_tube = 15e-2 # focale de la lentille de tube
f_obj = 1.8e-3 # focale de l'objectif
pix = 5.5e-6 # taille des pixels du capteur
RapFoc = 0.7 # rapport des focales du doublet de rééchantillonnage
SupFreq = 128 # rayon support de fréquence accessible en pixel
dimHolo = 2*SupFreq # dimension de l'hologramme
fmaxHolo = round(SupFreq*NA/nimm) # support max de fréquence
CentreX = 773 # position du centre de la pupille dans l'espace de Fourier
CentreY = 210
nb_holo = 1 # nombre d'hologrammes

# Détermination de la taille de l'image
# I = plt.imread(DossierData + 'i001.pgm');
I = cv2.imread(DossierData + 'i001.pgm',0);
# del I;

# Traitement de la séquence
cpt = 1
cpt_exist = 1
Centres = np.zeros((dimHolo,dimHolo))
OTFTomo = np.zeros((2*dimHolo,2*dimHolo,2*dimHolo))
SupRedon = np.zeros((2*dimHolo,2*dimHolo,2*dimHolo))
# SpectreFilt = pyfftw.empty_aligned((dimHolo, dimHolo), dtype='complex128')
# UBorn = pyfftw.empty_aligned((dimHolo, dimHolo), dtype='complex128')
# Image = pyfftw.empty_aligned((Nx, Ny), dtype='complex128')
# FImage = pyfftw.empty_aligned((Nx, Ny), dtype='complex128')
Support2D = np.ones((dimHolo,dimHolo))

start_time = time.time()
for hol in range(0,nb_holo):
    filename = DossierData + 'i' + str('%03d' % cpt ) + '.pgm'
    if os.path.isfile(filename):
        Image = plt.imread(filename)
        # Spectre hologramme
        FImage = nfft.fftshift(nfft.fft2(Image))
        # fft_obj = pyfftw.FFTW(Image,FImage,axes=(0,1),direction='FFTW_FORWARD',threads=4)
        # fft_obj.execute()
        nfy,nfx = np.shape(FImage)
        Spectre = holo.filtrage(FImage,CentreX,CentreY,fmaxHolo)
        SpectreFilt=Spectre[int((nfx-2*SupFreq)/2):int((nfx-2*SupFreq)/2+2*SupFreq),int((nfy-2*SupFreq)/2):int((nfy-2*SupFreq)/2+2*SupFreq)]
        # Coordonnées du speculaire
        ind = np.unravel_index(np.argmax(SpectreFilt, axis=None), SpectreFilt.shape)
        kiy = ind[0]
        kix = ind[1]
        Centres[kiy,kix] = 1
        # Champ complexe (UBorn + Ui)
        UBorn = nfft.ifft2(nfft.ifftshift(SpectreFilt))
        # fft_obj2 = pyfftw.FFTW(SpectreFilt,UBorn,axes=(0,1),direction='FFTW_BACKWARD',threads=4)
        # fft_obj2.execute()
        Amp_UBorn = np.abs(UBorn)
        Phase_UBornWrap = np.angle(UBorn)
        # Depliement de la phase
        Phase_UBorn = holo.unwrapping(Phase_UBornWrap, pix)
        # Correction amplitude et phase
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
        
        Amp_UBornCorr = CAber.ampliCorr(Amp_UBorn,Masque,Poly_US,Poly)
        plt.imshow(Amp_UBornCorr, cmap=plt.cm.gray)
        plt.colorbar()
        plt.show()
        
        CheminAmp = DossierAmplitude + 'AmpUBorn_' + str('%03d' % cpt ) + '.tiff'
        CheminPh = DossierPhase + 'PhaseUBorn_' + str('%03d' % cpt ) + '.tiff'
        # Enregistrement de l'amplitude et la phase dépliée (avant correction pour le test)
        # plt.imsave(CheminAmp,Amp_UBorn,cmap=plt.cm.gray)
        # plt.imsave(CheminPh,Phase_UBorn,cmap=plt.cm.gray)
        # Amp_UBornNorm = np.zeros((dimHolo,dimHolo))
        # cv2.normalize(Amp_UBorn, Amp_UBornNorm, 1.0, 0.0, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
        # cv2.imwrite(CheminAmp,Amp_UBornNorm)
        # cv2.imwrite(CheminPh,Phase_UBorn)      
        
        cpt = cpt + 1
        cpt_exist = cpt_exist + 1
    else:
        cpt = cpt + 1
print("--- %s seconds ---" % np.round(time.time() - start_time))
# plt.imshow(Centres, cmap=plt.cm.gray)
# plt.colorbar()
# plt.show()
# holo.SAVbin(np.float32(Centres),'Centres',str(dimHolo))
