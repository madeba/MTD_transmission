#ifndef __DEROUL_VOLKOV4__
#define __DEROUL_VOLKOV4__
#include <fftw3.h>
#include <vector>
#include <complex>
#include "vecteur.h"
#include "FFTW_init.h"
//#include "src/vecteur.h"

//fonction FFTW_init

//void deroul_volkov2(std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul,FFTW_init paramC2rHA);
std::vector<vecteur> init_kvect_shift(Var2D dim2DHA);
std::vector<double> init_kvect_mod2Shift(std::vector<vecteur> &kvect_shift);
//void deroul_volkov4_AS(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul, std::vector<vecteur> &kvect_shift,  FFTW_init &param_c2c);
//void deroul_volkov4_gradient_central(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul);
//void deroul_volkov4_total_sym(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul,std::vector<vecteur>  &kvect_shift,  FFTW_init &param_c2c);
//void deroul_volkov4_total_sym_mix(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul);
//void deroul_volkov4_total_sym_AS(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul);
//void deroul_volkov4_total_sym_paire(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul);
void deroul_volkov4_total_sym_paire(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul,std::vector <vecteur> double_kvect_shift,FFTW_init &param_c2c_double);
void gradient_fft4(std::vector<double>  &entree, std::vector<std::complex<double>> &gradx,std::vector<std::complex<double>> &grady, std::vector<vecteur>  &kvect_shift, FFTW_init &paramC2rHA);

void gradient_fft4(std::vector<std::complex<double>>  &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady,std::vector<vecteur>  &kvect_shift,FFTW_init &paramC2rHA);
void integ_grad4(std::vector<double> const &gradx, std::vector<double> const&grady, std::vector<std::complex<double>> &sortie,std::vector<vecteur> &kvect_shift, FFTW_init &paramC2rHA);

std::vector<double>  SymetriseY(std::vector<double> const &monImg, std::vector<double> &monImgSymetricY);
std::vector<double> SymetriseX(std::vector<double> const &monImg, std::vector<double> &monImgSymetricX);
void Symetrise_mirror(std::vector<double> const &monImg, std::vector<double> &monImgSymetric);
//std::vector<double> ASymetrise_OddX(std::vector<double> const &monImg, std::vector<double> &monImgSymetricX);
//std::vector<double> ASymetrise_EvenY(std::vector<double> const &monImg, std::vector<double> &monImgSymetricY);
std::vector<double> cut_quad4(std::vector<double> const &monImg4Quad);
std::vector<double> cut_quad4(std::vector<std::complex<double>> const &monImg4Quad);
void gradient_central(const std::vector<double> src, std::vector<double> &grad,std::string direction);

void gradient_central(const std::vector<std::complex<double>> src, std::vector<std::complex<double>> &grad,std::string direction);

void gradient_back(const std::vector<double> src, std::vector<double> &grad,std::string direction);

void gradient_back(const std::vector<std::complex<double>> src, std::vector<std::complex<double>> &grad,std::string direction);
/*
void gradient2(std::vector<double> src,std::vector<double> &grad, std::string direction, unsigned int dim);
void gradient_central2(std::vector<double> src, std::vector<double> &grad,std::string direction,unsigned short int dim);
*/
#endif
