# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import matplotlib.pyplot as plt
import numpy as np
import FileTools as ft
import Retropropagation as rp
import time
# from PIL import Image
# from scipy.fftpack import fftn,fftshift

# Path to the parameter file, and parameters reading
CheminParam = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/Param.txt"
REwald = ft.readvalue(CheminParam,'REwald')
nb_angle = int(ft.readvalue(CheminParam,'nb_angle'))
fmaxHolo = int(ft.readvalue(CheminParam,'fmaxHolo'))
dimHolo = int(ft.readvalue(CheminParam,'dimHolo'))
pixTheo = ft.readvalue(CheminParam,'pixTheo')
lambda_0 = 632.8e-9
n_imm = 1.515
UBornPitch = 1/(2*fmaxHolo*pixTheo)

# Paths to the real, and imaginary parts of the field
CheminReUBorn = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/ReBorn_{dimHolo}.bin"
CheminImUBorn = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/ImBorn_{dimHolo}.bin"

# Path to the specular coordinates
SpecCoordPath = f"C:/Users/p1600109/Documents/Recherche/MatlabTomo/Centres_{dimHolo}.txt"
fi = rp.Calc_fi(SpecCoordPath, nb_angle, REwald, dimHolo)

# Estimation of the diffracted vector coordinates
fd, sdz, dxm, dym = rp.Calc_fd(fmaxHolo,REwald)

# Field files reading
ReUBorn = rp.ReadBornCube(CheminReUBorn, dimHolo, dimHolo, nb_angle)
ImUBorn = rp.ReadBornCube(CheminImUBorn, dimHolo, dimHolo, nb_angle)
UBornCplx = ReUBorn + ImUBorn * 1j

# FU = fftshift(fftn(UBornCplx[:,:,1]))
# plt.imshow(np.log(abs(FU)),cmap="gray")
# plt.show()


start_time = time.time()
f_recon, TF_vol, mask_sum = rp.retropropagation(UBornCplx,fi,fmaxHolo,REwald,lambda_0,1.515,2*dimHolo,dimHolo,pixTheo,2*dimHolo,UBornPitch)
print(f"Reconstruction time for a {2*dimHolo}x{2*dimHolo}x{2*dimHolo} volume: {np.round(time.time() - start_time,decimals=2)} seconds")
Refraction = f_recon.real
Absorption = f_recon.imag

fidRef = open("Refraction_420.bin","a")
fidAbs = open("Absorption_420.bin","a")
for cpt in range(Refraction.shape[2]):
    Refr = Refraction[:,:,cpt]
    Abs = Absorption[:,:,cpt]
    Refr.tofile(fidRef)
    Abs.tofile(fidAbs)    
fidRef.close()
fidAbs.close()
# plt.imshow(Absorption[:,:,219],cmap="gray")
# plt.show()
# plt.imshow(Refraction[:,:,219],cmap="gray")
# plt.show()
# plt.imshow(mask_sum[:,:,0],cmap="gray")
# plt.show()
# ft.SAVbin(Refraction, "Refraction", str(2*dimHolo))
# # imAbs = Image.fromarray(Absorption)
# # imRefr = Image.fromarray(Refraction)
# # imAbs.save("Absorption.tiff")
# # imRefr.save("Refraction.tiff")