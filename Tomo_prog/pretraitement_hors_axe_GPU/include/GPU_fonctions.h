#ifndef __GPU_FONCTIONS__
#define __GPU_FONCTIONS__

#include <vector>
#include <complex>
#include "vecteur.h"
#include <arrayfire.h>
#include "struct.h"
#include "projet.h"

void display(af::array img);
af::array tukey2D(size_t dimx,size_t dimy, float alpha);
af::array crop(af::array img_orig,Var2D dim_crop, Var2D coord_coin);
void chargeImageGPU(af::array &img, std::string nomImg, Var2D coin);
void chargeImageGPU_via_OCV(af::array &img_array, std::string nomImg, Var2D coin);
void holo2TF_UBorn_GPU(af::array  &holo1,af::array &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, af::array  const &tukeyHolo, af::array const &tukeyBorn);
;
#endif
