# -*- coding: utf-8 -*-
"""

@author: Nicolas Verrier
"""
import SimuTomo as st
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy.fftpack import fftn,ifftn,fftshift
import Retropropagation as rp

dimHolo = 210
NA_ill = 1.4
nimm = 1.518
nbead = 1.55
nbangle = 600
Radius = 30

start_time = time.time()
if NA_ill/nimm>=1:
    print(f"Simulation impossible: nimm >= NA_ill")
else:
    OTF = st.OTF_Flower(dimHolo, NA_ill, nimm, nbangle)
    print(f"OTF simulation time for {nbangle} angles: {np.round(time.time() - start_time,decimals=2)} seconds")
    
    # plt.imshow(SphericalCoord,cmap = "gray")
    # plt.show()
    # plt.imshow(OTF[:,:,dimHolo],cmap = "gray")
    # plt.show()
    # plt.imshow(OTF[:,dimHolo,:],cmap = "gray")
    # plt.show()
    
    start_time = time.time()
    Bille = st.BeadSimu(Radius,dimHolo,nimm,nbead)
    print(f"Bead simulation time : {np.round(time.time() - start_time,decimals=2)} seconds")
    plt.imshow(Bille[:,:,dimHolo],cmap = "gray")
    plt.show()
    
    TomoSpectrum = fftn(Bille)*fftshift(OTF)
    Bille = ifftn(TomoSpectrum)
    plt.imshow(Bille[:,dimHolo,:].real,cmap = "gray")
    plt.show()
    plt.imshow(Bille[:,dimHolo,:].imag,cmap = "gray")
    plt.show()
    
    # Gerchberg parameters
    nbiter = 10    
    nmin = 1.518
    nmax = 1.55
    kappamin = 0
    kappamax = 0
    
    # Gerchberg reconstruction
    start_time = time.time()
    FilteredObj = rp.Gerchberg(Bille,OTF,nmin,nmax,kappamin,kappamax,nbiter)
    print(f"Reconstruction time for {nbiter} iterations: {np.round(time.time() - start_time,decimals=2)} seconds")
    plt.imshow(FilteredObj[:,dimHolo,:].real, cmap="gray")
    plt.show()
    plt.imshow(FilteredObj[:,dimHolo,:].imag, cmap="gray")
    plt.show() 