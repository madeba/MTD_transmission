# MTD_transmission
Code for tomographic diffractive microscopy : off axis holography, aberration correction, phase unwrapping, 3D reconstruction...

from Projet_tomo directory, launch
1. sh COMPIL.sh
to compile all the programs
2. sh INSTALL.sh 
to install them in the system path (root/sudo permission required)


C++ Dependencies : opencv4, fftw3, libtiff,   libboost-thread, libboost-chrono, , libboost-thread, libboost-system, wxWidgets (Graphical user interface). For GPU code, the API Arrayfire and a backend (CUDA or OpenCL), must be installed. `TomoManip` is used to cvontrol the setupo and will need libusb-1.0, exodriver labjack (https://labjack.com/support/software/installers/exodriver), pleora sdk5. An GUI can be used to control the set up, the parameters and the reconstruction process. wx-Widget must be in,stalled to compile it. 

The process essentially relies on 2 binaries : 
- `pretraitement_CPU` (or GPU) to extract the phase and amplitude of  an off-axis hologram, correct aberrations, unwrap the phase, normalise holograms....
- `Reconstruction` to calculate the final complex 3D image, saved in Tiff format.

These binaries use 3 config files :
1. `gui_tomo.conf`, usually placed in ~/.config. This is the first file read, needed to obtain the path to data or saved results.
2. `config_manip`, which gives the set-up parameters (magnification, wavelentgh, pixel size, etc.) usually placed in the data directory
3. `recon.txt`, which controls all paramters for the reconstruction

