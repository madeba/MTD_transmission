#ifndef __DEROUL_VOLKOV4__
#define __DEROUL_VOLKOV4__
#include <fftw3.h>
#include <vector>
#include <complex>
#include "vecteur.h"
#include "FFTW_init.h"
#include "gradient.h"
//#include "src/vecteur.h"

//fonction FFTW_init

//void deroul_volkov2(std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul,FFTW_init paramC2rHA);
std::vector<vecteur> init_kvect_shift(Var2D dim2DHA);
std::vector<double> init_kvect_mod2Shift(std::vector<vecteur> &kvect_shift);
void deroul_volkov4_total_sym_paire(std::vector<double>  &phase_enroul,std::vector<double> &phase_deroul,std::vector <vecteur> double_kvect_shift,FFTW_init &param_c2c_double);
//std::vector<double> deroul_volkov5_sym_paire_gradu(std::vector<std::complex<double>> const & UBorn, std::vector <vecteur> double_kvect_shift,FFTW_init &param_c2c_double);
std::vector<double> deroul_volkov5_sym_paire_gradu(std::vector<std::complex<double>>  & UBorn,std::vector<vecteur> kvect_shift, std::vector <vecteur> double_kvect_shift,FFTW_init &param_c2c,FFTW_init &param_c2c_double);

std::vector<double> deroul_volkov5_total_sym_paire_gradu(std::vector<std::complex<double>> & UBorn, std::vector<vecteur> double_kvect_shift,FFTW_init &param_c2c_double, double alpha_damp);
void gradient_fft4(std::vector<double>  &entree, std::vector<std::complex<double>> &gradx,std::vector<std::complex<double>> &grady, std::vector<vecteur>  &kvect_shift, FFTW_init &paramC2rHA);

void gradient_fft4(std::vector<std::complex<double>>  &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady,std::vector<vecteur>  &kvect_shift,FFTW_init &paramC2rHA);

std::vector<double>  SymetriseY(std::vector<double> const &monImg, std::vector<double> &monImgSymetricY);
std::vector<double> SymetriseX(std::vector<double> const &monImg, std::vector<double> &monImgSymetricX);
void Symetrise_mirror(std::vector<double> const &monImg, std::vector<double> &monImgSymetric);





#endif
