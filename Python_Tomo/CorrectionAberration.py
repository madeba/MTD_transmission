import numpy as np
import matplotlib.pyplot as plt
import os
import cv2


# Initialisation du masque de filtrage
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
        # Ouverture si le masque existe
        Masque = np.uint8(np.array(plt.imread(Chemin)))
        print("Ouverture du masque binaire")
    else:
        # Masque à 255 s'il n'est pas donné
        Masque = np.uint8(255 * np.ones((dimMasque,dimMasque)))
        print("Pas de masque trouvé ... chargement masque unité")
    return Masque
    
# Détermination de la taille du polynome
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
    for i in range(0,dimMasqueY-1):
        for j in range(0,dimMasqueX-1):
            if Masque[j,i]>45:
                count += 1
    nbPtRand = 0
    for i in range(0,dimMasqueY-1,step):
        for j in range(0,dimMasqueX-1,step):
            if Masque[j,i]>45:
                Masque[j,i]=220
                nbPtRand += 1
    return nbPtRand

# Génération du polynome à ajuster pour tout x,y en dehors du masque (sous échantillonné)
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
    # nbPtMask,Masque = PixInMask(Masque)
    nbRows,nbCols = np.shape(polynomeUS_to_fit)
    dimMasqueX, dimMasqueY = np.shape(Masque)
   
    if nbRows > 9:
        num_coef = Coord1D = 0
        # print("degre_poly=",degre_poly)  
        for y in range(0,dimMasqueY-1):
            for x in range(0,dimMasqueX-1):
                num_coef=0
                if Masque[y,x] == 220:                    
                    for powX in range(0,degre_poly-1):
                        powY = 0
                        while powX + powY <= degre_poly:                            
                            polynomeUS_to_fit[num_coef,Coord1D] = pow(x,powX)*pow(y,powY)
                            num_coef += 1
                            powY += 1                             
                    Coord1D += 1
                #print("coord",Coord1D)
    return polynomeUS_to_fit

# Génération du polynome à ajuster pour tout x,y en dehors du masque
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
    # nbPtMask,Masque = PixInMask(Masque)
    nbRows,nbCols = np.shape(polynome_to_fit)
    dimMasqueX, dimMasqueY = np.shape(Masque)
    if nbRows > 9:
        num_coef = Coord1D = 0
        for y in range(0,dimMasqueY-1):
            for x in range(0,dimMasqueX-1):
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
    coefX,coefY=np.shape(coefficients)
    
    UsCoord1D = 0
    nbCoef, NbPtImg = np.shape(polynome_to_fit)
    #print("PolyRows, PolyCols=",PolyRows,PolyCols)
    BackgrRows, BackgrCols = np.shape(PolyBackgr)
    PolySize = nbCoef
    for y in range(0,BackgrRows-1):
        for x in range(0,BackgrCols-1):
            sum = 0.
            for num_coef in range(0,PolySize-1):             
                sum += coefficients[num_coef] * polynome_to_fit[num_coef,UsCoord1D]       
            PolyBackgr[y,x] = sum
            UsCoord1D += 1
    # plt.imshow(PolyBackgr, cmap=plt.cm.gray)
    # plt.colorbar()
    # plt.show()
      
# Calcul des coefficients du polynome
def compuCoefPoly(ImageBrut,Masque,coef_polynomial,polynomeUS_to_fit,methode):
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
    methode : bool
        Choice of the method. If methode is true, lsq inversion. If methode is false, direct matrix inversion.

    Returns
    -------
    None.

    """
    ImgRows,ImgCols = np.shape(ImageBrut)
    nbCoef,nbPtUs = np.shape(polynomeUS_to_fit)
    # print("nbPtUs dimXpolyus==",nbCoef)
    # print("nbCoef dimYpolyUS==",nbPtUs)
    Bt = np.zeros((nbCoef,nbPtUs),dtype=np.float64)
    undersampled_background = np.zeros((nbPtUs,1),dtype=np.float64)
    x_poly,y_poly = np.shape(undersampled_background)
    # print("xusbg==",x_poly)
    # print("yusbg==",y_poly)
    if nbCoef>9:
        UsCoord = 0
        for y in range(0,ImgRows):
            for x in range(0,ImgCols):
                if Masque[y,x] == 220:
                    undersampled_background[UsCoord] = np.float64(ImageBrut[y,x])
                    UsCoord += 1
                    # print("USCoord",UsCoord)
        tabCoef = np.zeros((nbCoef,1),dtype=np.float64)
        
        tabCoefX,tabCoefY=np.shape(tabCoef)
        # print("tabCoefX=",tabCoefX)
        # print("tabCoefY=",tabCoefY)
        # print("type(tabCoef)=",type(tabCoef))
        # print("bckgrd",undersampled_background)
        
        D = np.zeros((nbPtUs,nbCoef),dtype=np.float64)
        invD = np.zeros((nbCoef,nbPtUs),dtype=np.float64)
        if methode:
            # print("solve==")
            cv2.solve(np.transpose(polynomeUS_to_fit), (undersampled_background), (tabCoef), cv2.DECOMP_SVD)
            # for cpt in range(0,nbCoef-1):
            #     print(tabCoef[cpt])
        
        else:
            cv2.transpose(polynomeUS_to_fit,Bt)
            D = Bt * polynomeUS_to_fit
            cv2.invert(D,invD)
            tabCoef = (invD * Bt) * undersampled_background
        # coef_polynomial=tabCoef
        # print("tabCoef=",tabCoef)
        np.copyto(coef_polynomial,tabCoef)
    else:
        tabCoef = np.zeros((nbCoef,1),dtype=np.float64)

# Correction de la phase du fond
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
    # print("typecoefsolve=",type(coefsolve))
    compuCoefPoly(Image,Masque,coefsolve,polynomeUS_to_fit,True)
    resultatpolyBG = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    resultat_final = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    compuBackgr(coefsolve,polynome_to_fit,resultatpolyBG)
    # plt.imshow(resultatpolyBG, cmap=plt.cm.gray)
    # plt.colorbar()
    # plt.show()
    resultat_final=Image-resultatpolyBG
    return resultat_final

# Normalisation de l'amplitude du fond
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
    compuCoefPoly(Image,Masque,coefsolve,polynomeUS_to_fit,True)
    resultatpoly = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    resultat = np.zeros((ImgRows,ImgCols),dtype=np.float64)
    compuBackgr(coefsolve,polynome_to_fit,resultatpoly)
    resultat=Image/(resultatpoly+np.finfo(float).eps)
    return resultat
