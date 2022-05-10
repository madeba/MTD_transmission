#!/bin/bash



DIR_TOMO=`pwd`
cd ${DIR_TOMO}/Tomo_prog/Gui_Tomo/
make clean
make
cd ${DIR_TOMO}/Tomo_prog/manipTomo5_pipeline/
make clean
make
cd ${DIR_TOMO}/Tomo_prog/pretraitement_hors_axe/
make clean
make
cd ${DIR_TOMO}/Tomo_prog/pretraitement_hors_axe_GPU/
make clean
make
cd ${DIR_TOMO}/Tomo_prog/Reconstruction
make clean
make
cd ${DIR_TOMO}/Tomo_prog/Tomo_GPS
make clean
make
cd ${DIR_TOMO}/utilTomo/show_fourier_PP
make clean
make




