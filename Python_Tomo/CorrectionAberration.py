# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os
import cv2
import numba


# Initialisation du masque de filtrage
# @numba.jit(nopython=True)
def InitMasque(Chemin,dimMasque):
    """
    Initialization of the mask for aberration correction. If the path is not valid, a blank mask of dimension dimMasque * dimMasque is generated.
    Parameters
    ----------
    Chemin : str
        Path + filename + extension to the binary phase mask.
    dimMasque : int
        width (assumed to equal heigth) of the mask to be applied in pixel.
    Returns
    -------
    Masque : uint8
        Binary mask applied to both amplitude and phase images.
    """
    if os.path.isfile(Chemin):
        Masque = np.uint8(np.array(plt.imread(Chemin)))
        print("Ouverture du masque binaire")
    else:
        Masque = np.uint8(255 * np.ones((dimMasque,dimMasque)))
        print("Pas de masque trouvé ... chargement masque unité")
    return Masque
    
# Détermination de la taille du polynome
# @numba.jit(nopython=True)
def SizePoly2D(deg):
    """
    Estimation of the size of the computed polynom. Further used in estimation of the column number in PolynomToFit function
    Parameters
    ----------
    deg : int
        Degree of the polynom to be fitted.
    Returns
    -------
    size : int
        Size of the polynom (equivalent to the number of columns in PolynomToFit function).
    """
    j = size = 0
    for i in range(0,deg-1):
        while i+j <=deg:
            size += 1
            j += 1
        j = 0
    return size

# Comptage des pixels dans le masque
# @numba.jit(nopython=True)
def PixInMask(Masque):
    """
    Estimation of the number of pixels within the mask
    Parameters
    ----------
    Masque : uint8
        Mask generated with the function \" InitMasque \".
    Returns
    -------
    nbPtRand : int
        Number of pixels used for the background estimation. The background is undersampled by a factor \" step \".
    Masque : uint8
        Modified mask.    
    """
    step = 5
    count = 0
    dimMasqueX, dimMasqueY = np.shape(Masque)
    for i in range(0,dimMasqueY):
        for j in range(0,dimMasqueX):
            if Masque[j,i]>45:
                count += 1
    nbPtRand = 0
    for i in range(0,dimMasqueY,step):
        for j in range(0,dimMasqueX,step):
            if Masque[j,i]>45:
                Masque[j,i]=220
                nbPtRand += 1
    return nbPtRand

# Génération du polynome à ajuster pour tout x,y en dehors du masque (sous échantillonné)
@numba.jit(nopython=True)
def CalcPolyUS_xy(degre_poly, Masque, polynomeUS_to_fit):
    """
    Generation of the undersampled polynome to fit for estimation of the background image
    Parameters
    ----------
    degre_poly : int
        Degree of the polynome to be fitted.
    Masque : uint8
        Mask applied to raw image for background estimation.
    polynomeUS_to_fit : float64
        Undersampled polynome to be calculated.
    Returns
    -------
    polynomeUS_to_fit : float64
        Undersampled polynome.
    """
    nbCoef,nbPtPolyUS = np.shape(polynomeUS_to_fit)
    dimMasqueX, dimMasqueY = np.shape(Masque) 
    if nbPtPolyUS > nbCoef:
        num_coef = Coord1D = 0  
        for y in range(0,dimMasqueY):
            for x in range(0,dimMasqueX):
                num_coef=0
                if Masque[y,x] == 220:                    
                    for powX in range(0,degre_poly-1):
                        powY = 0
                        while powX + powY <= degre_poly:                            
                            polynomeUS_to_fit[num_coef,Coord1D] = pow(x,powX)*pow(y,powY)
                            num_coef += 1
                            powY += 1                             
                    Coord1D += 1
    return polynomeUS_to_fit

# Génération du polynome à ajuster pour tout x,y en dehors du masque
@numba.jit(nopython=True)
def CalcPoly_xy(degre_poly, Masque, polynome_to_fit):
    """
    Generation of the polynome to fit for estimation of the background image
    Parameters
    ----------
    degre_poly : int
        Degree of the polynome to be fitted.
    Masque : uint8
        Mask applied to raw image for background estimation.
    polynome_to_fit : float64
        Polynome to be calculated.
    Returns
    -------
    polynome_to_fit : float64
        Calculated polynome.
    """
    nbCoef,nbPtPoly = np.shape(polynome_to_fit)   
    dimMasqueX, dimMasqueY = np.shape(Masque)
    if nbPtPoly > nbCoef:
        num_coef = Coord1D = 0
        for y in range(0,dimMasqueY):
            for x in range(0,dimMasqueX):
                num_coef=0
                for powX in range(0,degre_poly-1):
                    powY = 0
                    while powX + powY <= degre_poly:
                        polynome_to_fit[num_coef,Coord1D] = pow(x,powX)*pow(y,powY)
                        num_coef += 1
                        powY += 1
                Coord1D += 1
    return polynome_to_fit

# Calcul du fond
@numba.jit(nopython=True)
def compuBackgr(coefficients,polynome_to_fit,PolyBackgr):
    """
    Estimation of the hologram background
    Parameters
    ----------
    coefficients : float64
        Coefficients of the fitted polynome.
    polynome_to_fit : float64
        Calculated polynome.
    PolyBackgr : float64
        Background polynome.
    Returns
    -------
    None.
    """
    UsCoord1D = 0
    nbCoef, NbPtImg = np.shape(polynome_to_fit)
    result_mltpn=np.zeros((nbCoef,1),dtype=np.float64)
    BackgrRows, BackgrCols = np.shape(PolyBackgr)
    #PolySize = nbCoef
    for y in range(0,BackgrCols):
        for x in range(0,BackgrRows):
            result_mltpn[:,0]=coefficients[:,0] * polynome_to_fit[:,UsCoord1D]
            PolyBackgr[y,x] =np.sum(result_mltpn)
            UsCoord1D += 1
    # UsCoord1D = 0
    # nbCoef, NbPtImg = np.shape(polynome_to_fit)
    # BackgrRows, BackgrCols = np.shape(PolyBackgr)
    # PolySize = nbCoef
    # for y in range(0,BackgrRows):
    #     for x in range(0,BackgrCols):    
    #         sum = 0.
    #         for num_coef in range(0,PolySize):             
    #             sum += coefficients[num_coef]*polynome_to_fit[num_coef,UsCoord1D]       
    #         PolyBackgr[y,x] = sum
    #         UsCoord1D += 1
      
# Calcul des coefficients du polynome
# @numba.jit(nopython=True)
def compuCoefPoly(ImageBrut,Masque,coef_polynomial,polynomeUS_to_fit):
    """
    Least-Square computation of the coefficient of the fitted polynome
    Parameters
    ----------
    ImageBrut : uint8
        Source image.
    Masque : uint8
        Binary mask applied to the data.
    coef_polynomial : float64
        Coefficients of the fit (copied from coef therein).
    polynomeUS_to_fit : float64
        Polynome to be fitted (undersampled version).
    Returns
    -------
    None.
    """
    ImgRows,ImgCols = np.shape(ImageBrut)
    nbCoef,nbPtUs = np.shape(polynomeUS_to_fit)
    undersampled_background = np.zeros((nbPtUs,1),dtype=np.float64)
    x_poly,y_poly = np.shape(undersampled_background)
    if nbPtUs>nbCoef:
        UsCoord = 0
        for y in range(0,ImgRows):
            for x in range(0,ImgCols):
                if Masque[y,x] == 220:
                    undersampled_background[UsCoord] = np.float64(ImageBrut[y,x])
                    UsCoord += 1
        tabCoef = np.zeros((nbCoef,1),dtype=np.float64)                
        cv2.solve(np.transpose(polynomeUS_to_fit), undersampled_background, tabCoef, cv2.DECOMP_NORMAL)
        np.copyto(coef_polynomial,tabCoef)
    else:
        tabCoef = np.zeros((nbCoef,1),dtype=np.float64)

# Correction de la phase du fond
# @numba.jit(nopython=True)
def aberCorr(Image,Masque,polynomeUS_to_fit,polynome_to_fit):
    """
    Correction of the phase aberrations
    Parameters
    ----------
    Image : uint8
        Phase image.
    Masque : uint8
        Background mask.
    polynomeUS_to_fit : float64
        Polynome to be fitted (under-sampled for computation speeding-up).
    polynome_to_fit : float64
        Polynome to be fitted.
    Returns
    -------
    resultat_final : float64
        Corrected phase image.
    """
    nbCoef,nbPtUs = np.shape(polynomeUS_to_fit)
    ImgRows,ImgCols = np.shape(Image)
    coefsolve = np.zeros((nbCoef,1),dtype=np.float64)
    compuCoefPoly(Image,Masque,coefsolve,polynomeUS_to_fit)
    resultatpolyBG = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    resultat_final = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    compuBackgr(coefsolve,polynome_to_fit,resultatpolyBG)
    resultat_final=Image-resultatpolyBG
    return resultat_final

# Normalisation de l'amplitude du fond
# @numba.jit(nopython=True)
def ampliCorr(Image,Masque,polynomeUS_to_fit,polynome_to_fit):
    """
    Amplitude normalization
    Parameters
    ----------
    Image : uint8
        Amplitude image.
    Masque : uint8
        Background mask.
    polynomeUS_to_fit : float64
        Polynome to be fitted (under-sampled for computation speeding-up).
    polynome_to_fit : float64
        Polynome to be fitted.
    Returns
    -------
    resultat_final : float64
        Normalized amplitude image.
    """
    ImgRows,ImgCols = np.shape(Image)
    nbCoef,nbPtUs = np.shape(polynomeUS_to_fit)
    coefsolve = np.zeros((nbCoef,1),dtype=np.float64)
    compuCoefPoly(Image,Masque,coefsolve,polynomeUS_to_fit)
    resultatpoly = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    resultat = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    compuBackgr(coefsolve,polynome_to_fit,resultatpoly)
    resultat=Image/(resultatpoly+np.finfo(np.float32).eps)
    return resultat
