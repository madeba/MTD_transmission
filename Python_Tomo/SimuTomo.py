# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:31:30 2020

@author: Nicolas Verrier
"""
import numpy as np

def OTF_Cap(dimHolo, NA_ill, nimm, SpecCoord):
    """
    OTF simulation for a single cap of sphere

    Parameters
    ----------
    dimHolo : int
        Lateral dimension of.
    NA_ill : float
        Numerical aperture of the collection objective.
    nimm : float
        Immersion medium.
    SpecCoord : (int,int)
        (y,x) Coordinates of the specular spot.

    Returns
    -------
    OTF_Simu : int32
        Optical Transfert Function.    

    """
    REwald = dimHolo/2
    OTF_Simu = np.zeros((2*dimHolo, 2*dimHolo, 2*dimHolo))
    
    # Maximal frequency   
    Fmax = np.round(REwald*NA_ill/nimm)
    
    # Definition of diffracted wavevectors
    kdx = np.arange(-Fmax,Fmax)
    kdy = np.arange(-Fmax,Fmax)
    KDX, KDY = np.meshgrid(kdx, kdy)
    
    # Illumination wavevectors
    KIX = SpecCoord[1]
    KIY = SpecCoord[0]
    KIZ = np.round(np.sqrt(REwald**2-KIX**2-KIY**2))
    
    index = np.where(REwald**2-KDX**2-KDY**2 >= 0)
    KDZ = np.zeros(KDX.shape)
    KDZ[index] = np.round(np.sqrt(REwald**2-KDX[index]**2-KDY[index]**2))
        
    KOX = KDX - KIX
    KOY = KDY - KIY
    KOZ = KDZ - KIZ

    OTF_Simu[np.int32(KOY[index]+dimHolo),np.int32(KOX[index]+dimHolo),np.int32(KOZ[index]+dimHolo)]=1
    return OTF_Simu
    

def OTF_Flower(dimHolo, NA_ill, nimm, nbangle):
    """
    OTF simulation with a flower illumination swweping

    Parameters
    ----------
    dimHolo : int
        Lateral dimension of.
    NA_ill : float
        Numerical aperture of the collection objective.
    nimm : float
        Immersion medium.
    nbangle : int
        Number of illumination angles.

    Returns
    -------
    OTF_Simu : int32
        Optical Transfert Function.

    """
    REwald = dimHolo/2
    SpecCoord = np.zeros((dimHolo, dimHolo))
    OTF_Simu = np.zeros((2*dimHolo, 2*dimHolo, 2*dimHolo))
    
    # MaxAngle = np.arcsin(NA_ill/nimm)
    Fmax = np.round(REwald*NA_ill/nimm)
    
    # Flower parametrization
    n = 4
    d = 1
    NbPetals = n/d
    theta = np.linspace(0,2*np.pi,nbangle)
    # Diffracted frequencies
    kdx = np.arange(-Fmax,Fmax)
    kdy = np.arange(-Fmax,Fmax)
    KDX, KDY, THETA = np.meshgrid(kdx, kdy, theta)
    illum_number = np.linspace(1,len(theta),num=len(theta))

    KIX = np.round(Fmax*np.cos(NbPetals*THETA)*np.cos(THETA))
    KIY = np.round(Fmax*np.cos(NbPetals*THETA)*np.sin(THETA))
    KIZ = np.round(np.sqrt(REwald**2-KIX**2-KIY**2))
    
    SpecCoord[np.int32(KIY[0,0,:]+dimHolo/2+1),np.int32(KIX[0,0,:]+dimHolo/2+1)] = illum_number
    index = np.where(REwald**2-KDX**2-KDY**2 >= 0)
    KDZ = np.zeros(KDX.shape)
    KDZ[index] = np.round(np.sqrt(REwald**2-KDX[index]**2-KDY[index]**2))
        
    KOX = KDX - KIX
    KOY = KDY - KIY
    KOZ = KDZ - KIZ

    OTF_Simu[np.int32(KOY[index]+dimHolo),np.int32(KOX[index]+dimHolo),np.int32(KOZ[index]+dimHolo)]=1
    return OTF_Simu

def BeadSimu(Radius,dimHolo,nimm,nbead,kappa):
    """
    Simulation of a spherical bead for testing purposes

    Parameters
    ----------
    Radius : int
        Pixel radius of the bead.
    dimHolo : int
        Lateral dimension of the hologram.
    nimm : float
        Refractive index of the immersion medium.
    nbead : float
        Refractive index of the bead
    kappa : float
        Absorption of the bead

    Returns
    -------
    Bead : Complex128
        Simulated bead.

    """
    Refraction = np.zeros((2*dimHolo,2*dimHolo,2*dimHolo))
    Absorption = np.zeros((2*dimHolo,2*dimHolo,2*dimHolo))
    [X,Y,Z] = np.meshgrid(np.arange(0,Refraction.shape[1]),np.arange(0,Refraction.shape[0]),np.arange(0,Refraction.shape[2]))
    Center = np.array([Refraction.shape[1]/2,Refraction.shape[0]/2,Refraction.shape[2]/2])
    Refraction = (X - Center[1])**2 + (Y - Center[0])**2 + (Z - Center[2])**2 <= Radius**2
    Absorption = (X - Center[1]+10)**2 + (Y - Center[0]+10)**2 + (Z - Center[2])**2 <= (0.5*Radius)**2
    Absorption = kappa * Absorption
    Refraction = Refraction * (nbead-nimm) + nimm
    Bead = Refraction + 1j* Absorption
    return Bead

def Calc_TF_Holo(TF_vol3D, TF_Holo, fd_m, sdz_m, dx_m, dy_m, Nmax, k_inc, R_Ewald, lambda_v, n0):
    """
    Projection of the cap of sphere on a 2D plane

    Parameters
    ----------
    TF_vol3D : complex128
        3D spectrum of the object.
    TF_Holo : complex128
        2D spectrum of the object.
    fd_m : int
        Masked coordinates of the diffracted vectors.
    sdz_m : float
        Normalization of the diffracted vectors.
    dx_m : int
        Abcissa of the filtered coordinate.
    dy_m : int
        Ordonnae of the filtered coordinate.
    Nmax : int
        Maximal spatial frequency.
    k_inc : int32
        Array with the y, x coordinates of the specular spot.
    REwald : float64
        Radius of the Ewald sphere.
    lambda_v : float64
        Wavelength of light in vacuum.
    n0 : float64
        Refractive index of the immersion medium.

    Returns
    -------
    TF_Holo : complex128
        2D spectrum of the object.

    """
    fi = np.array([k_inc[0], k_inc[1], np.round(np.sqrt(R_Ewald**2 - k_inc[0]**2 - k_inc[1]**2))])
    
    kv = 2*np.pi/lambda_v
    k0 = kv * n0

    cteInd2Pot = -2*kv**2*n0
    ctePot2UBorn = -1j*np.pi/k0
    cteNormalisation = -1/(2*np.pi)

    fobj_m = (fd_m-np.round(fi[:, np.newaxis])/2).astype(int)
    
    TF_Holo[dx_m+Nmax,dy_m+Nmax] = (cteInd2Pot * ctePot2UBorn * cteNormalisation)/sdz_m * TF_vol3D[fobj_m[0,:], fobj_m[1,:], fobj_m[2,:]]
    return TF_Holo