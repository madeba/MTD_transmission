# MTD_transmission
Code for tomographic diffractive microscopy : off axis holography, aberration correction, phase unwrapping, 3D reconstruction...

From the Projet_tomo directory, run:

1. sh COMPIL.sh
to compile all the programs
2. sh INSTALL.sh 
to install them in the system path (root/sudo privileges required)


# C++ Dependencies

The following libraries are needed tu compile the code : 
OpenCV 4, FFTW3, libtiff, libboost-thread, libboost-chrono, libboost-system, and wxWidgets (graphical user interface).

For GPU code, the API Arrayfire and a backend (CUDA or OpenCL), must be installed. 

`TomoManip` is used to control the setup and requires libusb-1.0, the LabJack Exodriver (https://labjack.com/support/software/installers/exodriver), and the pleora SDK5. A GUI can be used to control the setup, the parameters, and the reconstruction process. wxWidgets must be installed to compile it.

The process essentially relies on 2 binaries : 
- `pretraitement_CPU` (or GPU) to extract the phase and amplitude of  an off-axis hologram, correct aberrations, unwrap the phase, normalise holograms....
- `Reconstruction`  to calculate the final complex 3D image, saved in TIFF format.

These binaries use 3 config files :
1. `gui_tomo.conf`, usually placed in ~/.config. This is the first file read and it's used to obtain the paths to the data and saved results.
2. `config_manip`, which gives the set-up parameters (magnification, wavelentgh, pixel size, etc.) usually placed in the data directory.
3. `recon.txt`, which controls all the reconstruction parameters.

