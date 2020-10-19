# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""
import numpy as np

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
    Mask = np.zeros((Spectrum.shape[1], Spectrum.shape[0]), dtype=complex)
    Mask[(kx-SpecCoord[1])**2+(ky-SpecCoord[0])**2 < FiltRadius2**2] = 0.7 * np.exp(-1j*np.pi/2)
    Mask[(kx-SpecCoord[1])**2+(ky-SpecCoord[0])**2 < FiltRadius1**2] = 1 + 1j
    Spectrum = Spectrum * Mask
    return Spectrum
