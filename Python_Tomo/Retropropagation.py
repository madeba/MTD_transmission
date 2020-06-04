# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

from linecache import getline
from scipy.fftpack import fftn,ifftn,fftshift
import numpy as np

def ReadBornCube(Chemin,dimX,dimY,nb_img):
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

    Returns
    -------
    DataCube : float64
        Raw data rearranged as a datacube.

    """
    with open(Chemin,'r') as fid:
        DataCube = np.fromfile(fid, np.float64)
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

def calc_tf_calotte(TF_holo,fd_m, sdz_m, dx_m, dy_m,P,P_holo,Nmax,k_inc,R_Ewald,lambda_v,n0):
    """
    Calculation of the cap of sphere

    Parameters
    ----------
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
    TF_vol3D = np.zeros(shape=(P,P,P),dtype=complex)
    mask_calotte = np.zeros_like(TF_vol3D,dtype=bool)
    fi = np.array([k_inc[0], k_inc[1], (np.sqrt(R_Ewald**2 - k_inc[0]**2 - k_inc[1]**2))])
    
    kv = 2*np.pi/lambda_v
    k0 =  kv * n0
    
    cteInd2Pot =np.complex(-2*kv**2*n0,0)
    ctePot2UBorn  =  np.complex(0,-np.pi/k0)
    cteNormalisation = np.complex(-1/(2*np.pi),0)

    fobj_m = (np.round(fd_m)-np.round(fi[:,np.newaxis])).astype(int)
    # print(f"Vecteur Objet de dimension : {fobj_m.shape}")

    TF_vol3D[fobj_m[0,:],fobj_m[1,:],fobj_m[2,:]] = TF_holo[dx_m+Nmax,dy_m+Nmax] * sdz_m
    
    mask_calotte[fobj_m[0,:],fobj_m[1,:],fobj_m[2,:]] = 1
    TF_vol3D /= (cteInd2Pot*ctePot2UBorn*cteNormalisation)

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
        3D Fourier transform of the object masked accordinf to experimental parameters.
    mask_sum : int32
        3D OTF.
        
    """
    if nb_holo>holo_pile.shape[2]:
        nb_holo = holo_pile.shape[2]

    TF_vol = np.zeros(shape=(P,P,P),dtype=complex)
    mask_sum = np.zeros(shape=(P,P,P),dtype=int) # changer pour OTF_filling
    fd_m, sdz_m, dx_m, dy_m = Calc_fd(Nmax,R_Ewald)
        
    for i in range(nb_holo):
        k_inc = np.array([SpecCoord[1,i], SpecCoord[0,i]])
        TF_holo_shift_r = fftshift(fftn(holo_pile[:,:,i]))/(Delta_f_Uborn**2*Nmax**2) # Nmax !
        decal = np.round(np.array([-k_inc[1]+dim_Uborn/2, -k_inc[0]+dim_Uborn/2])).astype(int)
        TF_holo_r = decal_TF_holo(TF_holo_shift_r,decal,P_holo)
        TF_calotte,mask_calotte = calc_tf_calotte(TF_holo_r, fd_m, sdz_m, dx_m, dy_m, P, P_holo, Nmax, k_inc, R_Ewald, lambda_v, n0)
        TF_vol += TF_calotte
        mask_sum += mask_calotte
    TF_vol[mask_sum !=0] = TF_vol[mask_sum!=0] / mask_sum[mask_sum!=0]
    f_recon = fftshift(ifftn(TF_vol))/(Tp_Tomo**3)
    return f_recon, TF_vol, mask_sum