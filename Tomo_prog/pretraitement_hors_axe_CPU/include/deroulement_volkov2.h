#ifndef __DEROUL_VOLKOV2__
#define __DEROUL_VOLKOV2__
#include <fftw3.h>
#include <vector>
#include <complex>
#include "FFTW_init.h"
//#include "src/vecteur.h"

//fonction FFTW_init

//void deroul_volkov2(std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul,FFTW_init paramC2rHA);

void deroul_volkov2(std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul , FFTW_init param_c2c, FFTW_init param_r2c_Uborn);
void gradient_fft2(std::vector<double> const &entree, std::vector<std::complex<double>> &gradx,std::vector<std::complex<double>> &grady, FFTW_init &paramC2rHA,FFTW_init &param_r2c_Uborn);

void gradient_fft2(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady,FFTW_init &paramC2rHA);
void integ_grad2(std::vector<double> const &gradx, std::vector<double> const&grady, std::vector<std::complex<double>> &sortie,FFTW_init &paramC2rHA,FFTW_init &param_r2c_Uborn);


/*
void gradient2(std::vector<double> src,std::vector<double> &grad, std::string direction, unsigned int dim);
void gradient_central2(std::vector<double> src, std::vector<double> &grad,std::string direction,unsigned short int dim);
*/
#endif
