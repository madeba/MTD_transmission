#ifndef __DEROUL_VOLKOV3__
#define __DEROUL_VOLKOV3__
#include <fftw3.h>
#include <vector>
#include <complex>
#include "vecteur.h"
#include "FFTW_init.h"
//#include "src/vecteur.h"

//fonction FFTW_init

//void deroul_volkov2(std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul,FFTW_init paramC2rHA);
//std::vector<vecteur> init_kvect_shift(Var2D dim2DHA); //déplacé dans fft_fonctions
//std::vector<double> init_kvect_mod2Shift(std::vector<vecteur> &kvect_shift);//déplacé dans fft_fonctions
void deroul_volkov3(std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul, std::vector<vecteur> &kvect_shift,  FFTW_init &param_c2c);

void gradient_fft3(std::vector<double> const &entree, std::vector<std::complex<double>> &gradx,std::vector<std::complex<double>> &grady, std::vector<vecteur>  &kvect_shift, FFTW_init &paramC2rHA);

void gradient_fft3(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady,std::vector<vecteur>  &kvect_shift,FFTW_init &paramC2rHA);
void integ_grad3(std::vector<double> const &gradx, std::vector<double> const&grady, std::vector<std::complex<double>> &sortie,std::vector<vecteur> &kvect_shift, FFTW_init &paramC2rHA);


/*
void gradient2(std::vector<double> src,std::vector<double> &grad, std::string direction, unsigned int dim);
void gradient_central2(std::vector<double> src, std::vector<double> &grad,std::string direction,unsigned short int dim);
*/
#endif
