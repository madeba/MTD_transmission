# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier et Steve Laroche
"""

from linecache import getline
import time
from scipy import signal
from scipy.fftpack import fftn, ifftn, fftshift, ifftshift
import numpy as np
from FileTools import NextPow2

def ReadCube(Chemin, dimX, dimY, nb_img, datatype):
    """
    Opening raw values of UBorn files and transforming it as a data cube with (dimX,dimY,nb_img) dimensions

    Parameters
    ----------
    Chemin : str
        Path to the file containing raw data.
    dimX : int
        X dimension of the datacube.
    dimY : int
        Y dimension of the datacube.
    nb_img : int
        Number of images in the datacube.
    datatype : str
        Data type.

    Returns
    -------
    DataCube : Raw data rearranged as a (y,x,z)- datacube.

    """
    with open(Chemin, 'r') as fid:
        DataCube = np.fromfile(fid, eval(datatype))
    DataCube = DataCube.reshape((nb_img, dimX, dimY)).transpose(1, -1, 0)
    return DataCube

def Calc_fi(Chemin, nb_angle, dimHolo):
    """
    Extraction of center coordinates

    Parameters
    ----------
    Chemin : str
        Path to the file containing center coordinates.
    nb_angle : int
        Number of points to be processed.
    dimHolo : int
        Size of the hologram

    Returns
    -------
    k_inc : int
        Array with the y, x coordinates of the specular spot.

    """
    CoordSpec = np.zeros((nb_angle, 2), dtype=int)
    for line in range(1, nb_angle+1):
        Ligne = getline(Chemin, line).split()
        CoordSpec[line-1, :] = [int(Ligne[0]), int(Ligne[1])]
    kix = CoordSpec[:, 1]-dimHolo/2
    kiy = CoordSpec[:, 0]-dimHolo/2
    k_inc = np.array([kiy, kix])
    return k_inc

def Calc_fd(Nmax, REwald):
    """
    Estimation of the diffracted vector coordinates

    Parameters
    ----------
    Nmax : int
        Maximal accessible frequency.
    REwald : float64
        Radius of the Ewald sphere.

    Returns
    -------
    fd_m : int
        Masked coordinates of the diffracted vectors.
    sdz_m : float
        Normalization of the diffracted vectors.
    dx_m : int
        Abcissa of the filtered coordinate.
    dy_m : int
        Ordonnae of the filtered coordinate.

    """
    dy, dx = np.meshgrid(np.arange(-Nmax, Nmax), np.arange(-Nmax, Nmax))
    Mask = dx**2+dy**2 < Nmax**2
    dxm = dx[Mask]
    # print(f"{dxm}")
    dym = dy[Mask]
    dzm = np.round(np.sqrt(REwald**2 - dxm**2 - dym**2))
    sdzm = np.sqrt(REwald**2 - dxm**2 - dym**2) / REwald
    fdm = np.array([dym, dxm, dzm])
    return fdm, sdzm, dxm, dym

def decal_TF_holo(TF_holo, decal, P_holo, roll=True):
    """
    Shifting the holograms according to illumination frequency

    Parameters
    ----------
    TF_holo : complex128
        2D Fourier transform of the hologram.
    decal : int
        Amount of shift in pixels.
    P_holo : int
        Lateral dimension of the holograms.
    roll : bool
        Using np.roll to perform registration. Default is False

    Returns
    -------
    TF_holo_shift : complex128
        Shifted 2D Fourier transform.

    """
    TF_holo_shift = np.zeros_like(TF_holo)
    decal[0] = decal[0]%P_holo
    decal[1] = decal[1]%P_holo
    if decal[0] < 0:
        decal[0] += P_holo
    if decal[1] < 0:
        decal[1] += P_holo
    if roll is True:
            TF_holo_shift = np.roll(TF_holo, decal[1], axis=1)
            TF_holo_shift = np.roll(TF_holo_shift, decal[0], axis=0)
    else:
        TF_holo_shift[decal[0]:, decal[1]:] = TF_holo[0:P_holo-decal[0], 0:P_holo-decal[1]]
        TF_holo_shift[:decal[0], :decal[1]] = TF_holo[P_holo-decal[0]:, P_holo-decal[1]:]
        TF_holo_shift[decal[0]:, :decal[1]] = TF_holo[0:P_holo-decal[0], P_holo-decal[1]:]
        TF_holo_shift[:decal[0], decal[1]:] = TF_holo[P_holo-decal[0]:, 0:P_holo-decal[1]]
    return TF_holo_shift

def calc_tf_calotte(TF_vol3D, mask_calotte, TF_holo, fd_m, sdz_m, dx_m, dy_m, Nmax, k_inc, R_Ewald, lambda_v, n0):
    """
    Calculation of the cap of sphere

    Parameters
    ----------
    TF_vol3D : complex128
        3D Fourier transform of the object.
    mask_calotte : int32
        3D Mask corresponding to the OTF.
    TF_holo : complex128
        2D Fourier transform of the hologram.
    fd_m : int
        Masked coordinates of the diffracted vectors.
    sdz_m : float
        Normalization of the diffracted vectors.
    dx_m : int
        Abcissa of the filtered coordinate.
    dy_m : int
        Ordonnae of the filtered coordinate.
    Nmax : int
        Maximal accessible frequency.
    k_inc : int
        Coordinates of the illumination vector.
    R_Ewald : float
        Ewald sphere radius.
    lambda_v : float
        Wavelength in vacuum of the illumination source.
    n0 : float
        Refractive index of the immersion medium.

    Returns
    -------
    TF_vol3D : complex128
        3D Fourier transform of the object.
    mask_calotte : int32
        3D Mask corresponding to the OTF.

    """
    # print(f"{Nmax}")
    fi = np.array([k_inc[0], k_inc[1], np.round(np.sqrt(R_Ewald**2 - k_inc[0]**2 - k_inc[1]**2))])
    # print(f"{fi}")

    kv = 2*np.pi/lambda_v
    k0 = kv * n0

    cteUBorn2Pot = 1j*k0/np.pi
    cteNormalisation = -2*np.pi
    # cteInd2Pot = -1/(2*kv**2*n0)
    # ctePot2UBorn = -1j*k0/np.pi
    # cteNormalisation = -(2*np.pi)

    fobj_m = (fd_m-(fi[:, np.newaxis])).astype(int)
    TF_vol3D[fobj_m[0, :], fobj_m[1, :], fobj_m[2, :]] += TF_holo[dx_m+Nmax, dy_m+Nmax]*sdz_m*cteUBorn2Pot*cteNormalisation
    mask_calotte[fobj_m[0, :], fobj_m[1, :], fobj_m[2, :]] += 1

    return TF_vol3D, mask_calotte

def retropropagation(holo_pile, nb_holo, SpecCoord, Nmax, R_Ewald, lambda_v, n0, Tp_Tomo, Delta_f_Uborn):
    """
    Retropropagation algorithm

    Parameters
    ----------
    holo_pile : complex128
        Holograms to be processed. The datacube has dimHolo x dimHolo x nb_angle voxels
    SpecCoord : int
        Coordinates of the specular spot.
    Nmax : int
        Maximal accessible frequency.
    R_Ewald : float
        Ewald sphere radius.
    lambda_v : float
        Wavelength in vaccum of the illumination source.
    n0 : float
        Refractive index of the immersion medium.
    Tp_Tomo : float
        Pixel pitch for the tomographic experiment (accounting for magnification).
    Delta_f_Uborn : float
        Frequency pitch of tomography.

    Returns
    -------
    f_recon : complex128
        Reconstructed 3D object.
    TF_vol : complex128
        3D Fourier transform of the object masked according to experimental parameters.
    mask_sum : int32
        3D OTF

    """
    if nb_holo > holo_pile.shape[2]:
        nb_holo = holo_pile.shape[2]
    
    pow2 = NextPow2(2*holo_pile.shape[0])
    
    kv = 2*np.pi/lambda_v
    ctePot2Ind = -1/(2*kv**2*n0)

    # TF_vol = np.zeros(shape=(2*holo_pile.shape[0],
    #                          2*holo_pile.shape[0],
    #                          2*holo_pile.shape[0]), dtype=complex)
    # mask_sum = np.zeros(shape=(2*holo_pile.shape[0],
    #                            2*holo_pile.shape[0],
    #                            2*holo_pile.shape[0]), dtype=np.int32)
    TF_vol = np.zeros(shape=(2**pow2, 2**pow2, 2**pow2), dtype=complex)
    mask_sum = np.zeros(shape=(2**pow2, 2**pow2, 2**pow2), dtype=np.int32)    
    fd_m, sdz_m, dx_m, dy_m = Calc_fd(Nmax, R_Ewald)

    # Tukey Window
    TukeyWindow = np.sqrt(np.outer(signal.tukey(holo_pile.shape[0], 0.1),
                                   signal.tukey(holo_pile.shape[0], 0.1)))
    print("")
    print("-------------------------")
    print("- Mapping Fourier space -")
    print("-------------------------")
    for cpt in range(nb_holo):
        if cpt % 100 == 0:
            print(f"Hologram = {cpt} out of {nb_holo}")
        k_inc = np.array([SpecCoord[1, cpt], SpecCoord[0, cpt]])
        TF_holo_shift_r = fftshift(fftn(TukeyWindow*holo_pile[:, :, cpt]))/(Delta_f_Uborn**2*Nmax**2)
        decal = np.array([k_inc[1], k_inc[0]]).astype(int)
        TF_holo_r = decal_TF_holo(TF_holo_shift_r, decal, holo_pile.shape[0], roll=True)
        TF_vol3D, mask_sum = calc_tf_calotte(TF_vol, mask_sum, TF_holo_r, fd_m, sdz_m,
                                           dx_m, dy_m, Nmax, k_inc, R_Ewald, lambda_v, n0)
    TF_vol3D[mask_sum != 0] = TF_vol3D[mask_sum != 0] / mask_sum[mask_sum != 0]
    TF_vol3D[mask_sum == 0] = 0
    print("")
    print("----------")
    print("- FFT 3D -")
    print("----------")
    start_time = time.time()
    f_recon = fftshift(ifftn(TF_vol3D))/(Tp_Tomo**3)
    print("")
    print(f"3D-FFT calculation: {np.round(time.time() - start_time,decimals=2)} seconds")
    print("")
    mask_sum = fftshift(mask_sum)
    TF_vol3D = fftshift(TF_vol3D)
    f_recon = ctePot2Ind*fftshift(f_recon, [0, 1])
    return f_recon, TF_vol, mask_sum

def Gerchberg(ReconsObjOrig, OTF, Delta_nmin, Delta_nmax, Kappa_min, Kappa_max, nbiter):
    """
    Iterative Gerchberg processing

    Parameters
    ----------
    ReconsObjOrig : complex128
        3D distribution of the complex refractive index of the object.
    OTF : int32
        Frequency support for regularization.
    Delta_nmin : float
        Minimal value of the real part of the refractive index.
    Delta_nmax : float
        Maximal value of the real part of the refractive index.
    Kappa_min : flaot
        Minimal value of the imaginary part of the refractive index (absorption).
    Kappa_max : float
        Maximal value of the imaginary part of the refractive index (absorption).
    nbiter : int
        Number of iterations.

    Returns
    -------
    ReconsObjEst : complex128
        Estimated reconstruction.

    """
    OTF = fftshift(OTF)
    ReconsFieldOrig = fftn(ReconsObjOrig)
    ReconsObjEst = ReconsObjOrig
    for cpt in range(nbiter):
        start_time = time.time()
        print(f"Gerchberg iteration {cpt}")
        Refr_index = ReconsObjEst.real
        Absorb_val = ReconsObjEst.imag

        Refr_index[Refr_index < Delta_nmin] = Delta_nmin
        Refr_index[Refr_index > Delta_nmax] = Delta_nmax
        Absorb_val[Absorb_val < Kappa_min] = Kappa_min
        Absorb_val[Absorb_val > Kappa_max] = Kappa_max

        ReconsObjEst = Refr_index + Absorb_val*1j
        ReconsFieldEst = fftn(ReconsObjEst)

        if cpt < nbiter/2:
            index = np.nonzero(OTF[0: int(OTF.shape[0]/2-1), :, :])
        if cpt >= nbiter/2:
            index = np.nonzero(OTF)

        ReconsFieldEst[index] = ReconsFieldOrig[index]
        ReconsObjEst = ifftn(ReconsFieldEst)
        print(f"Iteration duration: {np.round(time.time() - start_time,decimals=2)} seconds")
    return ReconsObjEst

def DarkField(ComplexField, CutOff):
    """
    Processing darkfield images from tomographic acquisitions

    Parameters
    ----------
    ComplexField : complex128
        Reconstructed complex refractive index distribution.
    CutOff : int
        Cut-off frequency of the darkfield filter (assumed to be a circulat mask).

    Returns
    -------
    Field : complex128
        Darkfield filtered refractive index distribution.

    """
    Spectrum = fftshift(fftn(ComplexField))
    kx, ky, kz = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                             np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)),
                             np.arange(-int(Spectrum.shape[2]/2), int(Spectrum.shape[2]/2)))
    Spectrum[kx**2 + ky**2 + kz**2 < CutOff**2] = 0
    Field = ifftn(ifftshift(Spectrum))
    return Field

def PhaseContrast(ComplexField, CutOff1, CutOff2):
    """
    Processing phase-contrast images from tomographic acquisitions

    Parameters
    ----------
    ComplexField : complex128
        Reconstructed complex refractive index distribution.
    CutOff : int
        Cut-off frequency of the phase-contrast filter (assumed to be a circulat mask).

    Returns
    -------
    Field : complex128
        Phase-contrast filtered refractive index distribution.

    """
    Spectrum = fftshift(fftn(ComplexField))
    kx, ky, kz = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                             np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)),
                             np.arange(-int(Spectrum.shape[2]/2), int(Spectrum.shape[2]/2)))
    Mask = np.zeros((Spectrum.shape[1], Spectrum.shape[0], Spectrum.shape[2]), dtype=complex)
    Mask[kx**2 + ky**2 + kz**2 < CutOff2**2] = 0.1*(1 + 1*1j)
    Mask[kx**2 + ky**2 + kz**2 < CutOff1**2] = np.exp(-1j*np.pi/2)
    Spectrum = Spectrum * Mask
    Field = ifftn(ifftshift(Spectrum))
    return Field

def RheinbergIllumination(ComplexField, CutR, CutG, CutB):
    """
    Processing Rheinberg illumination images from tomographic acquisitions

    Parameters
    ----------
    ComplexField : complex128
        Reconstructed complex refractive index distribution.
    CutR : int
        Cut-off frequency of the Red filter (assumed to be a circulat mask).
    CutG : int
        Cut-off frequency of the Green filter (assumed to be a circulat mask).
    CutB : int
        Cut-off frequency of the Blue filter (assumed to be a circulat mask).

    Returns
    -------
    Field : complex128
        Phase-contrast filtered refractive index distribution.

    """
    Spectrum = fftshift(fftn(ComplexField))
    kx, ky, kz = np.meshgrid(np.arange(-int(Spectrum.shape[1]/2), int(Spectrum.shape[1]/2)),
                             np.arange(-int(Spectrum.shape[0]/2), int(Spectrum.shape[0]/2)),
                             np.arange(-int(Spectrum.shape[2]/2), int(Spectrum.shape[2]/2)))
    FiltR = np.zeros((Spectrum.shape[1], Spectrum.shape[0], Spectrum.shape[2]), dtype=complex)
    FiltG = np.zeros((Spectrum.shape[1], Spectrum.shape[0], Spectrum.shape[2]), dtype=complex)
    FiltB = np.zeros((Spectrum.shape[1], Spectrum.shape[0], Spectrum.shape[2]), dtype=complex)

    FiltR[kx**2 + ky**2 + kz**2 < CutR[1]**2] = 1
    FiltR[kx**2 + ky**2 + kz**2 < CutR[0]**2] = 0

    FiltG[kx**2 + ky**2 + kz**2 < CutG**2] = 1

    FiltB[kx**2 + ky**2 + kz**2 < CutB[1]**2] = 1
    FiltB[kx**2 + ky**2 + kz**2 < CutB[0]**2] = 0

    R = Spectrum * FiltR
    G = Spectrum * FiltG
    B = Spectrum * FiltB

    FieldR = ifftn(ifftshift(R))
    FieldG = ifftn(ifftshift(G))
    FieldB = ifftn(ifftshift(B))

    return FieldR, FieldG, FieldB

def dic(ComplexField, DeltaX, DeltaY, DeltaZ, Phi):
    """
    Processing DIC images from tomographic acquisitions

    Parameters
    ----------
    ComplexField : complex128
        Reconstructed complex refractive index distribution.
    DeltaX : int
        Hologram shearing in x-direction
    DeltaY : int
        Hologram shearing in y-direction
    DeltaZ : int
        Hologram shearing in z-direction        
    Phi : float64
        Hologram phase-shift.

    Returns
    -------
    Field : complex128
        Phase-contrast filtered refractive index distribution.

    """
    ComplexFieldShift = np.roll(ComplexField, DeltaX, axis=1)
    ComplexFieldShift = np.roll(ComplexFieldShift, DeltaY, axis=0)
    ComplexFieldShift = np.roll(ComplexFieldShift, DeltaZ, axis=2)    
    IDic = ComplexField - (ComplexFieldShift*np.exp(1j*Phi))
    return IDic