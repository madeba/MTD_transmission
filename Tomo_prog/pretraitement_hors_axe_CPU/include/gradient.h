#ifndef __GRADIENT__
#define __GRADIENT__
#include <fftw3.h>
#include <vector>
#include "vecteur.h"
#include <complex>
#include "FFTW_init.h"
///gradient par fft, entrée réelle
///calculate gradient by fft, real input
void gradient_fft4(std::vector<double> &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady, std::vector<vecteur> &kvect_shift, FFTW_init &param_c2c);
///gradient par fft, entrée complexe
///calculate gradient by fft, complex input (overloaded function)
void gradient_fft4(std::vector<std::complex<double>> &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady,std::vector<vecteur>  &kvect_shift, FFTW_init &param_c2c);

////--------------------------gradient par différence dans l'espace image------------------------
void gradient_central(const std::vector<double> src, std::vector<double> &grad,std::string direction);
void gradient_central(const std::vector<std::complex<double>> src, std::vector<std::complex<double>> &grad,std::string direction);
void gradient_back(const std::vector<double> src, std::vector<double> &grad,std::string direction);
void gradient_back(const std::vector<std::complex<double>> src, std::vector<std::complex<double>> &grad,std::string direction);
#endif

