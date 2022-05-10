#ifndef __DEROUL_VOLKOV3_GPU__
#define __DEROUL_VOLKOV3_GPU__
//#include <fftw3.h>
//#include <vector>
//#include <complex>
//#include "vecteur.h"
//#include "FFTW_init.h"
//#include "src/vecteur.h"


void deroul_volkov3_GPU2(af::array const &phase_enroul, af::array  &phase_deroul,af::array  const &kvect_shiftX,af::array const &kvect_shiftY);
void gradient_fft3_GPU2(af::array const &entree, af::array &gradx, af::array &grady, af::array const &kvect_shiftX,af::array const  &kvect_shiftY) ;
void integ_grad3_GPU2(af::array const &gradx, af::array const& grady, af::array &sortie,af::array const &kvect_shiftX, af::array const &kvect_shiftY);
af::array init_kvect_shiftX_GPU(int dim2DHA);
af::array init_kvect_shiftY_GPU(int dim2DHA);
#endif
