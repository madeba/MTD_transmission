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
        Processed filter.

    """
    kx, ky = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2),int(Spectrum.shape[1]/2)), np.arange(-int(Spectrum.shape[0]/2),int(Spectrum.shape[0]/2)))
    Spectrum[(kx-SpecCoord[1])**2+(ky-SpecCoord[0])**2<FiltRadius**2] = 0
    return Spectrum