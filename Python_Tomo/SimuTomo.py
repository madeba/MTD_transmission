# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:31:30 2020

@author: Nicolas Verrier
"""
import numpy as np

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

def BeadSimu(Radius,dimHolo,nimm,nbead):
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

    Returns
    -------
    Bead : float64
        Simulated bead.

    """
    Bead = np.zeros((2*dimHolo,2*dimHolo,2*dimHolo))
    [X,Y,Z] = np.meshgrid(np.arange(0,Bead.shape[1]),np.arange(0,Bead.shape[0]),np.arange(0,Bead.shape[2]))
    Center = np.array([Bead.shape[1]/2,Bead.shape[0]/2,Bead.shape[2]/2])
    Bead = (X - Center[1])**2 + (Y - Center[0])**2 + (Z - Center[2])**2 <= Radius**2
    Bead = Bead * (nbead-nimm) + nimm
    return Bead