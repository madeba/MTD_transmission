# -*- coding: utf-8 -*-
import numpy.fft as nfft
import numpy as np

"""
Collection of functions helping for hologram processing.
"""

# Filtrage et recentrage de l'information
def filtrage(spectrum, dx, dy, radius):
    """
    
    Bandpass filtering of the hologram spectrum. A circular mask is applied to the
    Fourier space data. Filtered spectrum is finally shifted of the quantities dx and dy
    Parameters
    ----------
    spectrum : complex128
        Spectrum to be filtered.
    dx : int
        Shift in pixels of the filtered spectrum in \"x\" direction in pixel.
    dy : int
        Shift in pixels of the filtered spectrum in \"y\" direction in pixel.
    radius : int
        Radius of the bandpass filter in pixel.
    Returns
    -------
    filtSpectrum : complex128
        Filtered and recentered spectrum.
    """
    ny, nx = np.shape(spectrum)
    xmin = 0
    xmax = nx - 1
    ymin = 0
    ymax = ny - 1

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    xpix2d, ypix2d = np.meshgrid(x, y)  # Mesh en pixel
    shiftx = int(nx/2-dx)
    shifty = int(ny/2-dy)

    mask = np.ones((ny, nx)) - (np.hypot(xpix2d - dx, ypix2d - dy) > radius)
    filtSpectrum = spectrum * mask
    filtSpectrum = np.roll(filtSpectrum, shiftx, axis=1)
    filtSpectrum = np.roll(filtSpectrum, shifty, axis=0)
    return filtSpectrum

# Pour l'holographie : méthodes de reconstruction
def reconstruction(I, z, Lambda, pix):
    """
    
    Implementation of the Fresnel transform for hologram reconstruction. Depending on the object to sensor distance, reconstruction kernel is estimated in image or Fourier
    space.
    Parameters
    ----------
    I : float64
    Hologram to be reconstructed.
    z : float
        Object to sensor distance.
    Lambda : float
        Illumination wavelength.
    pix : float
        Sensor pixel size. If working with microscope objective, effective pixel size
        (accounting for total magnification) has to be considered.
    Returns
    -------
    Irest : complex128
        Reconstructed hologram. Amplitude of the wavefront can be estimated via np.abs(Irecons), while wrapped phase is obtained via np.angle(Irecons).
    """
    nx, ny = np.shape(I)  # Taille de l'image
    xmin_m = -nx / 2
    xmax_m = nx / 2 - 1
    ymin_m = -ny / 2
    ymax_m = ny / 2 - 1
    x_m = np.linspace(xmin_m, xmax_m, nx) * pix
    y_m = np.linspace(ymin_m, ymax_m, ny) * pix
    X_m, Y_m = np.meshgrid(x_m, y_m)  # Mesh metrique

    if np.abs(z) <= nx*pix**2/Lambda:
        pasu = 1 / (nx * pix)
        pasv = 1 / (ny * pix)
        kx = np.linspace(xmin_m, xmax_m, nx) * pasu
        ky = np.linspace(ymin_m, ymax_m, ny) * pasv

        KX, KY = np.meshgrid(kx, ky)  # Mesh frequentiel

        H = np.exp(-1j * z * np.pi * Lambda * (KX ** 2 + KY ** 2))
        C11 = nfft.fft2(I) * H
        Irecons = nfft.ifft2(C11)
        return Irecons
    else:
        phaseQ = np.exp(1j*np.pi/(Lambda*z)*(X_m**2+Y_m**2))
        C22 = nfft.fft2(I)*nfft.fft2(phaseQ)/nx
        D22 = (nfft.ifft2(C22))
        F22 = D22/1j
        Irecons=(F22)
        return Irecons

# Dépliement de Phase Volkov
def unwrapping(PhiW, pix, approx=True):
    """
    
    Phase unwrapping algorithm as proposed in Vyacheslav V. Volkov and Yimei Zhu, \"Deterministic phase unwrapping in the presence of noise,\" Opt. Lett. 28, 2156-2158 (2003).
    Parameters
    ----------
    PhiW : float64
        Wrapped phase image.
    pix : float64
        Sensor pixel size. If working with microscope objective, effective pixel size
        (accounting for total magnification) has to be considered.
    approx : bool, optional
        Choice of the unwrapping method. If approx is true, PhiUW is calculating using Eq. (4) to approximate phase gradients.
    Returns
    -------
    PhiUW : float64
        Unwrapped phase image.
    """
    nx, ny = np.shape(PhiW)  # Taille de l'image
    xmin_m = -nx / 2
    xmax_m = nx / 2 - 1
    ymin_m = -ny / 2
    ymax_m = ny / 2 - 1

    # Maillage spatial
    x_m = np.linspace(xmin_m, xmax_m, nx) * pix
    y_m = np.linspace(ymin_m, ymax_m, ny) * pix
    X_m, Y_m = np.meshgrid(x_m, y_m)  # Mesh metrique

    # Maillage frequentiel
    pasu = 1 / (nx * pix)
    pasv = 1 / (ny * pix)
    kx = np.linspace(xmin_m, xmax_m, nx) * pasu
    ky = np.linspace(ymin_m, ymax_m, ny) * pasv

    KX, KY = np.meshgrid(kx, ky)  # Mesh frequentiel
    NormK = nfft.ifftshift(KX**2)+nfft.ifftshift(KY**2)
    NormK[0, 0] = 1/3*(NormK[0, 1]+NormK[1, 0]+NormK[1,1])

    # Gestion des sauts de phase
    Zphi = np.exp(1j*PhiW)
    Gradz_x = 2*1j*np.pi*nfft.ifft2(nfft.fft2(Zphi)*nfft.ifftshift(KX))
    Gradz_y = 2*1j*np.pi*nfft.ifft2(nfft.fft2(Zphi)*nfft.ifftshift(KY))
    if approx is True:
        Eq_x = nfft.fft2(np.conj(1j*Zphi)*Gradz_x)*nfft.ifftshift(KX)/NormK
        Eq_y = nfft.fft2(np.conj(1j*Zphi)*Gradz_y)*nfft.ifftshift(KY)/NormK
        PhiUW = np.real(1/(2*1j*np.pi)*nfft.ifft2(Eq_x+Eq_y))
    else:
        Gradpw_x = 2*1j*np.pi*nfft.ifft2(nfft.fft2(PhiW)*nfft.ifftshift(KX))
        Gradpw_y = 2*1j*np.pi*nfft.ifft2(nfft.fft2(PhiW)*nfft.ifftshift(KY))
        Gradk_x = 1/(2*np.pi)*(np.real(Gradz_x/(1j*Zphi))-Gradpw_x)
        Gradk_y = 1/(2*np.pi)*(np.real(Gradz_y/(1j*Zphi))-Gradpw_y)
        q_x = nfft.fft2(Gradk_x)*nfft.ifftshift(KX)/NormK
        q_y = nfft.fft2(Gradk_y)*nfft.ifftshift(KY)/NormK
        kphi = np.real(1/(2*1j*np.pi)*nfft.ifft2(q_x+q_y))
        PhiUW = PhiW+2*np.pi*kphi
    return PhiUW
