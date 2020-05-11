# -*- coding: utf-8 -*-
# import numpy.fft as nfft
import scipy.fftpack as sfft
import numpy as np
import numba

"""
Collection of functions helping for hologram processing.
"""

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
    xpix2d, ypix2d = np.meshgrid(x, y)
    shiftx = int(nx/2-dx)
    shifty = int(ny/2-dy)

    mask = np.ones((ny, nx)) - (np.hypot(xpix2d - dx, ypix2d - dy) > radius)
    filtSpectrum = spectrum * mask
    filtSpectrum = np.roll(filtSpectrum, shiftx, axis=1)
    filtSpectrum = np.roll(filtSpectrum, shifty, axis=0)
    return filtSpectrum

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
    nx, ny = np.shape(I)
    xmin_m = -nx / 2
    xmax_m = nx / 2 - 1
    ymin_m = -ny / 2
    ymax_m = ny / 2 - 1
    x_m = np.linspace(xmin_m, xmax_m, nx) * pix
    y_m = np.linspace(ymin_m, ymax_m, ny) * pix
    X_m, Y_m = np.meshgrid(x_m, y_m)

    if np.abs(z) <= nx*pix**2/Lambda:
        pasu = 1 / (nx * pix)
        pasv = 1 / (ny * pix)
        kx = np.linspace(xmin_m, xmax_m, nx) * pasu
        ky = np.linspace(ymin_m, ymax_m, ny) * pasv

        KX, KY = np.meshgrid(kx, ky)

        H = np.exp(-1j * z * np.pi * Lambda * (KX ** 2 + KY ** 2))
        C11 = sfft.fft2(I) * H
        Irecons = sfft.ifft2(C11)
        return Irecons
    else:
        phaseQ = np.exp(1j*np.pi/(Lambda*z)*(X_m**2+Y_m**2))
        C22 = sfft.fft2(I)*sfft.fft2(phaseQ)/nx
        D22 = (sfft.ifft2(C22))
        F22 = D22/1j
        Irecons=(F22)
        return Irecons

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
    nx, ny = np.shape(PhiW)
    xmin_m = -nx / 2
    xmax_m = nx / 2 - 1
    ymin_m = -ny / 2
    ymax_m = ny / 2 - 1

    # Spatial mesh
    x_m = np.linspace(xmin_m, xmax_m, nx) * pix
    y_m = np.linspace(ymin_m, ymax_m, ny) * pix
    X_m, Y_m = np.meshgrid(x_m, y_m)

    # Frequency mesh
    pasu = 1 / (nx * pix)
    pasv = 1 / (ny * pix)
    kx = np.linspace(xmin_m, xmax_m, nx) * pasu
    ky = np.linspace(ymin_m, ymax_m, ny) * pasv

    KX, KY = np.meshgrid(kx, ky)
    KX = sfft.ifftshift(KX)
    KY = sfft.ifftshift(KY)
    NormK = KX**2+KY**2
    NormK[0, 0] = 1/3*(NormK[0, 1]+NormK[1, 0]+NormK[1,1])

    # Phase jump management
    Zphi = np.exp(1j*PhiW)
    Gradz_x = 2*1j*np.pi*sfft.ifft2(sfft.fft2(Zphi)*KX)
    Gradz_y = 2*1j*np.pi*sfft.ifft2(sfft.fft2(Zphi)*KY)
    if approx is True:
        Eq_x = sfft.fft2(np.conj(1j*Zphi)*Gradz_x)*KX/NormK
        Eq_y = sfft.fft2(np.conj(1j*Zphi)*Gradz_y)*KY/NormK
        PhiUW = np.real(1/(2*1j*np.pi)*sfft.ifft2(Eq_x+Eq_y))
    else:
        Gradpw_x = 2*1j*np.pi*sfft.ifft2(sfft.fft2(PhiW)*KX)
        Gradpw_y = 2*1j*np.pi*sfft.ifft2(sfft.fft2(PhiW)*KY)
        Gradk_x = 1/(2*np.pi)*(np.real(Gradz_x/(1j*Zphi))-Gradpw_x)
        Gradk_y = 1/(2*np.pi)*(np.real(Gradz_y/(1j*Zphi))-Gradpw_y)
        q_x = sfft.fft2(Gradk_x)*KX/NormK
        q_y = sfft.fft2(Gradk_y)*KY/NormK
        kphi = np.real(1/(2*1j*np.pi)*sfft.ifft2(q_x+q_y))
        PhiUW = PhiW+2*np.pi*kphi
    return PhiUW

@numba.jit(nopython=True)
def CoordToCoordShift(CoordX,CoordY,dimImgX,dimImgY):
    """
    Convert position of the filter center to fftshifted ones. Allows to get rid of the fftshift operation in hologram filtering

    Parameters
    ----------
    CoordX : int
        X coordinate of the filter center.
    CoordY : int
        Y coordinate of the filter center.
    dimImgX : int
        Dimension of the spectrum along kx axis.
    dimImgY : int
        Dimension of the spectrum along ky axis.

    Returns
    -------
    CoordShiftX : int
        Shifted X coordinate of the filter center.
    CoordShiftY : int
        Shifted Y coordinate of the filter center.

    """
    if CoordX-dimImgX/2>0:
        CoordShiftX = CoordX-dimImgX/2
    else:
        CoordShiftX = dimImgX/2+CoordX
    if CoordY-dimImgY/2>0:
        CoordShiftY = CoordY-dimImgY/2
    else:
        CoordShiftY = dimImgY/2+CoordY
    return int(CoordShiftX),int(CoordShiftY)