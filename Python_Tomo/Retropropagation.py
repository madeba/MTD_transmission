# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

from linecache import getline
from scipy.fftpack import fftn,ifftn,fftshift
import numpy as np
import time

def ReadBornCube(Chemin,dimX,dimY,nb_img,isint):
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
    isint : bool
        Data type. If True, data to be loaded are int32, if False, data to be loaded are float64

    Returns
    -------
    DataCube : int32 or float64 depending on isint value
        Raw data rearranged as a datacube.

    """
    if isint == False:
        with open(Chemin,'r') as fid:
            DataCube = np.fromfile(fid, np.float64)
        DataCube = DataCube.reshape((nb_img,dimX,dimY)).transpose(1,2,0)
    else:
        with open(Chemin,'r') as fid:
            DataCube = np.fromfile(fid, np.int32)
        DataCube = DataCube.reshape((nb_img,dimX,dimY)).transpose(1,2,0)        
    return DataCube

def Calc_fi(Chemin,nb_angle, REwald, dimHolo):
    """
    Extraction of center coordinates

    Parameters
    ----------
    Chemin : str
        Path to the file containing center coordinates.
    nb_angle : int
        Number of points to be processed.
    REwald : float64
        Radius of the Ewald sphere
    dimHolo : int
        Size of the hologram

    Returns
    -------
    k_inc : int
        Array with the y, x coordinates of the specular spot.

    """
    CoordSpec = np.zeros((nb_angle,2),dtype=int)
    for line in range (1,nb_angle+1):
        Ligne = getline(Chemin,line).split()
        CoordSpec[line-1,:] = [int(Ligne[0]),int(Ligne[1])]
    kix = CoordSpec[:,1]-dimHolo/2
    kiy = CoordSpec[:,0]-dimHolo/2
    k_inc = np.array([kiy, kix])
    return k_inc

def Calc_fd(Nmax,REwald):
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
    dy, dx, _ = np.meshgrid(np.arange(-Nmax,Nmax), np.arange(-Nmax,Nmax), np.arange(-Nmax,Nmax))
    Mask = dx**2+dy**2 < Nmax**2
    dxm = dx[Mask]
    dym = dy[Mask]
    dzm = np.round(np.sqrt(REwald**2 - dxm**2 - dym**2))
    sdzm = np.sqrt(REwald**2 - dxm**2 - dym**2) / REwald
    fdm =  np.array([dym ,dxm, dzm])
    return fdm, sdzm, dxm, dym

def decal_TF_holo(TF_holo,decal,P_holo):
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

    Returns
    -------
    TF_holo_shift : complex128
        Shifted 2D Fourier transform.

    """
    TF_holo_shift = np.zeros_like(TF_holo)
    decal[0] = decal[0]%P_holo
    decal[1] = decal[1]%P_holo
    if decal[0] < 0 : decal[0] +=P_holo
    if decal[1] < 0 : decal[1] +=P_holo
    TF_holo_shift[decal[0]:,decal[1]:] = TF_holo[0:P_holo-decal[0], 0:P_holo-decal[1]]  
    TF_holo_shift[:decal[0],:decal[1]] = TF_holo[P_holo-decal[0]:, P_holo-decal[1]:]
    TF_holo_shift[decal[0]:,:decal[1]] = TF_holo[0:P_holo-decal[0], P_holo-decal[1]:]
    TF_holo_shift[:decal[0],decal[1]:] = TF_holo[P_holo-decal[0]:, 0:P_holo-decal[1]]
    return TF_holo_shift

def calc_tf_calotte(TF_vol3D,mask_calotte,TF_holo,fd_m, sdz_m, dx_m, dy_m,P,P_holo,Nmax,k_inc,R_Ewald,lambda_v,n0):
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
    P : int
        Width of the tomographic volume.
    P_holo : int
        Width of the hologram.
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
    fi = np.array([k_inc[0], k_inc[1], np.round(np.sqrt(R_Ewald**2 - k_inc[0]**2 - k_inc[1]**2))])
    
    kv = 2*np.pi/lambda_v
    k0 =  kv * n0
    
    cteInd2Pot =np.complex(-2*kv**2*n0,0)
    ctePot2UBorn  =  np.complex(0,-np.pi/k0)
    cteNormalisation = np.complex(-1/(2*np.pi),0)

    fobj_m = (fd_m-(fi[:,np.newaxis])).astype(int)
    TF_vol3D[fobj_m[0,:],fobj_m[1,:],fobj_m[2,:]] += TF_holo[dx_m+Nmax,dy_m+Nmax]*sdz_m/(cteInd2Pot*ctePot2UBorn*cteNormalisation)
    mask_calotte[fobj_m[0,:],fobj_m[1,:],fobj_m[2,:]] += 1
    
    return TF_vol3D,mask_calotte

def retropropagation(holo_pile,nb_holo,SpecCoord,Nmax,R_Ewald,lambda_v,n0,P,P_holo,Tp_Tomo,dim_Uborn,Delta_f_Uborn):
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
    P : int
        Width of the tomographic volume.
    P_holo : int
        Width of the hologram.
    Tp_Tomo : float
        Pixel pitch for the tomographic experiment (accounting for magnification).
    dim_Uborn : int
        Width of the hologram.
    Delta_f_Uborn : float
        Frequency pitch of tomography.

    Returns
    -------
    f_recon : complex128
        Reconstructed 3D object.
    TF_vol : complex128
        3D Fourier transform of the object masked according to experimental parameters.
    mask_sum : int32
        3D OTF.
        
    """
    if nb_holo>holo_pile.shape[2]:
        nb_holo = holo_pile.shape[2]

    TF_vol = np.zeros(shape=(P,P,P),dtype=complex)
    mask_sum = np.zeros(shape=(P,P,P),dtype=int)
    fd_m, sdz_m, dx_m, dy_m = Calc_fd(Nmax,R_Ewald)
        
    for i in range(nb_holo):
        k_inc = np.array([SpecCoord[1,i], SpecCoord[0,i]])
        TF_holo_shift_r = fftshift(fftn(holo_pile[:,:,i]))/(Delta_f_Uborn**2*Nmax**2)
        decal = np.array([k_inc[1], k_inc[0]]).astype(int)
        TF_holo_r = decal_TF_holo(TF_holo_shift_r,decal,P_holo)
        TF_vol,mask_sum = calc_tf_calotte(TF_vol,mask_sum,TF_holo_r, fd_m, sdz_m, dx_m, dy_m, P, P_holo, Nmax, k_inc, R_Ewald, lambda_v, n0)
    TF_vol[mask_sum !=0] = TF_vol[mask_sum!=0] / mask_sum[mask_sum!=0]
    f_recon = fftshift(ifftn(TF_vol))/(Tp_Tomo**3)
    mask_sum = fftshift(mask_sum)
    f_recon = fftshift(f_recon,[0,1])
    return f_recon, TF_vol, mask_sum

def Gerchberg(ReconsObjOrig,OTF,Delta_nmin,Delta_nmax,Kappa_min,Kappa_max,nbiter):
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
    ReconsObjEst=ReconsObjOrig
    for cpt in range(nbiter):
        start_time = time.time()
        print(f"Gerchberg iteration {cpt}")
        Refr_index = ReconsObjEst.real
        Absorb_val = ReconsObjEst.imag
        
        # Refr_index[np.unravel_index(Refr_index<Delta_nmin, Refr_index.shape)]=Delta_nmin
        # Refr_index[np.unravel_index(Refr_index>Delta_nmax, Refr_index.shape)]=Delta_nmax
        # Absorb_val[np.unravel_index(Absorb_val<Kappa_min, Absorb_val.shape)]=Kappa_min
        # Absorb_val[np.unravel_index(Absorb_val>Kappa_max, Absorb_val.shape)]=Kappa_max
        
        Refr_index[Refr_index<Delta_nmin]=Delta_nmin
        Refr_index[Refr_index>Delta_nmax]=Delta_nmax
        Absorb_val[Absorb_val<Kappa_min]=Kappa_min
        Absorb_val[Absorb_val>Kappa_max]=Kappa_max
        
        ReconsObjEst = Refr_index + Absorb_val*1j
        ReconsFieldEst = fftn(ReconsObjEst)
        
        if cpt<nbiter/2:
            # index = np.unravel_index(OTF[0:int(OTF.shape[0]/2-1),:,:], OTF.shape)
            index = np.nonzero(OTF[0:int(OTF.shape[0]/2-1),:,:])
            
        if cpt>=nbiter/2:
            # index = np.unravel_index(OTF, OTF.shape)
            index = np.nonzero(OTF)
        ReconsFieldEst[index] = ReconsFieldOrig[index]
        ReconsObjEst = ifftn(ReconsFieldEst)
        print(f"Iteration duration: {np.round(time.time() - start_time,decimals=2)} seconds")
    return ReconsObjEst