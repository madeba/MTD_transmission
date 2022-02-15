# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""
import numpy as np
from scipy.fftpack import fft2, ifft2

def darkfield(Spectrum, FiltRadius, SpecCoord):
    """
    Emulation of darkfield measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    FiltRadius : int
        Radius of the filter.
    SpecCoord : int32
        Coordinates of the specular spot.

    Returns
    -------
    Spectrum : complex128
        Processed spectrum.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    Spectrum[(kx-SpecCoord[1])**2+(ky-SpecCoord[0])**2 < FiltRadius**2] = 0
    return Spectrum

def FiltIlluDarkField(Spectrum, Radius, SpecCoord):
    
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    
    Filt = np.ones((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    
    Filt[(kx)**2+(ky)**2 < Radius**2] = 0   
    
    return Filt

def phasecontrast(Spectrum, FiltRadius1, FiltRadius2, SpecCoord):
    """
    Emulation of phase contrast measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    FiltRadius1 : int
        Radius of the inner filter.
    FiltRadius2 : int
        Radius of the outer filter.
    SpecCoord : int32
        Coordinates of the specular spot.

    Returns
    -------
    Spectrum : complex128
        Processed spectrum.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    Mask = np.ones((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    Mask[(kx-SpecCoord[1])**2+(ky-SpecCoord[0])**2 < FiltRadius2**2] = 1 * np.exp(1j*np.pi/2)
    Mask[(kx-SpecCoord[1])**2+(ky-SpecCoord[0])**2 < FiltRadius1**2] = 0.7*(1+ 1*1j)
    Spectrum = Spectrum * Mask
    return Spectrum

def rheinberg(Spectrum, Radius, CentreR, CentreG, CentreB, SpecCoord):
    """
    Emulation of Rheinberg illumination measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    Radius : int
        Radius of the Red, Green and Blue Filter.
    CentreR : int
        Center of the Red filter.
    CentreG : int
        Center of the Green filter.
    CentreB : int
        Center of the Blue filter.
    SpecCoord : int32
        Coordinates of the specular spot.

    Returns
    -------
    R, G, B : complex128
        Processed spectrum.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    for i in range (len(CentreR)):
        FiltR[(kx+CentreR[i][0]-SpecCoord[1])**2+(ky+CentreR[i][1]-SpecCoord[0])**2 < Radius**2] = 1
    for j in range (len(CentreG)):
        FiltG[(kx+CentreG[j][0]-SpecCoord[1])**2+(ky+CentreG[j][1]-SpecCoord[0])**2 < Radius**2] = 3
    for k in range (len(CentreB)):
        FiltB[(kx+CentreB[k][0]-SpecCoord[1])**2+(ky+CentreB[k][1]-SpecCoord[0])**2 < Radius**2] = 1
    
    R = (FiltR * Spectrum)
    G = (FiltG * Spectrum)
    B = (FiltB * Spectrum)
    
    return R, G, B

def rheinberg2DTech2(Spectrum, Radius, CentreR, CentreG, CentreB, SpecCoord, Angle):
    """
     An other emulation of Rheinberg illumination measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    Radius : int
        Radius of the Red, Green and Blue Filter.
    CentreR : int
        Center of the Red filter.
    CentreG : int
        Center of the Green filter.
    CentreB : int
        Center of the Blue filter.
    SpecCoord : int32
        Coordinates of the specular spot.
    Angle : float
        Angle for calculate the position of filter

    Returns
    -------
    R, G, B : complex128
        Processed spectrum.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    
    FiltG[(kx+CentreG[0]-SpecCoord[1])**2+(ky+CentreG[1]-SpecCoord[0])**2 < Radius**2] = 1
    
    
    for i in range (8):
        FiltR[(kx+(CentreR[0]*(np.cos(Angle*i)))-SpecCoord[1])**2+(ky+(CentreR[0]*(np.sin(Angle*i)))-SpecCoord[0])**2 < Radius**2] = 1
        FiltB[(kx+(CentreB[0]*(np.cos(Angle*i)))-SpecCoord[1])**2+(ky+(CentreB[0]*(np.sin(Angle*i)))-SpecCoord[0])**2 < Radius**2] = 1
        
    R = (FiltR * Spectrum)
    G = (FiltG * Spectrum)
    B = (FiltB * Spectrum)

    return R, G, B

def rheinberg2DTech3(Spectrum, Radius, CentreR, CentreG, CentreB, SpecCoord, AngleR, AngleB):
    """
     An other emulation of Rheinberg illumination measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    Radius : int
        Radius of the Red, Green and Blue Filter.
    CentreR : int
        Center of the Red filter.
    CentreG : int
        Center of the Green filter.
    CentreB : int
        Center of the Blue filter.
    SpecCoord : int32
        Coordinates of the specular spot.
    AngleR : float
        Angle for calculate the position of Red filter
    AngleB : float
        Angle for calculate the position of Blue filter

    Returns
    -------
    R, G, B : complex128
        Intensity image

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    ListFiltR = [0]*8
    ListFiltB = [0]*8
    
    FiltG[(kx+CentreG[0]-SpecCoord[1])**2+(ky+CentreG[1]-SpecCoord[0])**2 < 50**2] = 1
    
    
    for i in range(4):
        FiltR[(kx+(CentreR[0]*(np.cos(AngleR*i)))-SpecCoord[1])**2+(ky+(CentreR[0]*(np.sin(AngleR*i)))-SpecCoord[0])**2 < Radius**2] = 1
        PartFiltR = np.copy(FiltR)
        ListFiltR[i] = PartFiltR
        FiltR[:,:] = 0
        
    for j in [1,3,5,7]:
        FiltB[(kx+(CentreB[0]*(np.cos(AngleB*j)))-SpecCoord[1])**2+(ky+(CentreB[0]*(np.sin(AngleB*j)))-SpecCoord[0])**2 < Radius**2] = 1
        PartFiltB = np.copy(FiltB)
        ListFiltB[i] = PartFiltB
        FiltB[:,:] = 0
        
    ListFiltR = [np.abs(ifft2(n * Spectrum))**2 for n in ListFiltR]
    ListFiltB = [np.abs(ifft2(n * Spectrum))**2 for n in ListFiltB]
    
    R = sum(ListFiltR)
    B = sum(ListFiltB) 
        
    # R = (FiltR * Spectrum)
    G =  np.abs(ifft2(FiltG * Spectrum))**2
    # B = (FiltB * Spectrum)

    return R, G, B

def rheinberg2DTech4(Spectrum, Radius, CentreR, CentreG, CentreB, SpecCoord, Angle):
    """
     An other emulation of Rheinberg illumination measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    Radius : int
        Radius of the Red, Green and Blue Filter.
    CentreR : int
        Center of the Red filter.
    CentreG : int
        Center of the Green filter.
    CentreB : int
        Center of the Blue filter.
    SpecCoord : int32
        Coordinates of the specular spot.
    Angle : float
        Angle for calculate the position of filter

    Returns
    -------
    R, G, B : complex128
        Intensity image.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    ListFiltR = [0]*8
    ListFiltB = [0]*8
    
    for i in range(8):
        FiltR[(kx+(CentreR[0]*(np.cos(Angle*i)))-SpecCoord[1])**2+(ky+(CentreR[0]*(np.sin(Angle*i)))-SpecCoord[0])**2 < Radius**2] = 1
        PartFiltR = np.copy(FiltR)
        ListFiltR[i] = PartFiltR
        FiltR[:,:] = 0
        
        FiltB[(kx+(CentreB[0]*(np.cos(Angle*i)))-SpecCoord[1])**2+(ky+(CentreB[0]*(np.sin(Angle*i)))-SpecCoord[0])**2 < Radius**2] = 1
        PartFiltB = np.copy(FiltB)
        ListFiltB[i] = PartFiltB
        FiltB[:,:] = 0
        
    FiltG[(kx+CentreG[0]-SpecCoord[1])**2+(ky+CentreG[1]-SpecCoord[0])**2 < 50**2] = 1
        
    ListFiltR = [np.abs(ifft2(n * Spectrum))**2 for n in ListFiltR]
    ListFiltB = [np.abs(ifft2(n * Spectrum))**2 for n in ListFiltB]
        
    R = sum(ListFiltR)
    B = sum(ListFiltB)    
    G = np.abs(ifft2(FiltG * Spectrum))**2
    
    
    
    return R, G, B

def FiltIlluRheinberg(Spectrum, RadiusR,RadiusG,RadiusB, SpecCoord):
    
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    
    
    FiltR[(kx)**2+(ky)**2 < RadiusR**2] = 1    
    FiltG[(kx)**2+(ky)**2 < RadiusG**2] = 1
    FiltB[(kx)**2+(ky)**2 < RadiusB**2] = 1
    
    FiltB = FiltB - FiltR
    FiltR = FiltR - FiltG
    
    
    return FiltR, FiltG, FiltB
    
def rheinberg2DTech5(Spectrum, Radius, CentreR, CentreG, CentreB, SpecCoord, Angle, Compo):
    """
     An other emulation of Rheinberg illumination measurement from holographic acquisitions

    Parameters
    ----------
    Spectrum : complex128
        Spectrum of the field to be processed.
    Radius : int
        Radius of the Red, Green and Blue Filter.
    CentreR : int
        Center of the Red filter.
    CentreG : int
        Center of the Green filter.
    CentreB : int
        Center of the Blue filter.
    SpecCoord : int32
        Coordinates of the specular spot.
    Angle : float
        Angle for calculate the position of filter
    Compo : String
        Caracter for choose which component (RGB)

    Returns
    -------
    SpectrumFilt : complex128
        Processed spectrum.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                         np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)))
    
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    ListFiltR = [0]*8
    ListFiltB = [0]*8
    
    for i in range(8):
        FiltR[(kx+(CentreR[0]*(np.cos(Angle*i)))-SpecCoord[1])**2+(ky+(CentreR[0]*(np.sin(Angle*i)))-SpecCoord[0])**2 < Radius**2] = 1
        PartFiltR = np.copy(FiltR)
        ListFiltR[i] = PartFiltR
        FiltR[:,:] = 0
        
        FiltB[(kx+(CentreB[0]*(np.cos(Angle*i)))-SpecCoord[1])**2+(ky+(CentreB[0]*(np.sin(Angle*i)))-SpecCoord[0])**2 < Radius**2] = 1
        PartFiltB = np.copy(FiltB)
        ListFiltB[i] = PartFiltB
        FiltB[:,:] = 0
        
    FiltG[(kx+CentreG[0]-SpecCoord[1])**2+(ky+CentreG[1]-SpecCoord[0])**2 < Radius**2] = 1
        
    ListFiltR = [np.abs(ifft2(n * Spectrum))**2 for n in ListFiltR]
    ListFiltB = [np.abs(ifft2(n * Spectrum))**2 for n in ListFiltB]
        
    R = sum(ListFiltR)
    B = sum(ListFiltB)    
    G = np.abs(ifft2(FiltG * Spectrum))**2
    
    if Compo == 'R':
        SpectrumFilt = fft2(R)
    if Compo == 'G':
        SpectrumFilt = fft2(G)
    if Compo == 'B':
        SpectrumFilt = fft2(B)
    
    return SpectrumFilt


    
def dic(ComplexField, DeltaX, DeltaY, Phi):
    """
    Emulation of DIC measurement from holographic acquisitions

    Parameters
    ----------
    ComplexField : complex128
        Field to be processed.
    DeltaX : int
        Field shift on X.
    DeltaY : int
        Field shift on Y.
    Phi : float
        phase shift angle
    SpecCoord : int32
        Coordinates of the specular spot.

    Returns
    -------
    Spectrum : complex128
        Processed spectrum.

    """
    ComplexFieldShift = np.roll(ComplexField, DeltaX, axis=1)
    ComplexFieldShift = np.roll(ComplexFieldShift, DeltaY, axis=0)
    IDic = np.abs(ComplexField - (ComplexFieldShift*np.exp(1j*Phi)))**2
    return IDic


