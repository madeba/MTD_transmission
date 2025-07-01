#include <fftw3.h>
#include <vector>
#include <complex>
//#include "src/vecteur.h"


void gradient_fft(fftw_complex *in,fftw_complex *out, std::vector<std::complex<double>> const &entree, std::vector<std::complex<double>> &gradx, std::vector<std::complex<double>> &grady,fftw_plan p_forward, fftw_plan p_backward);
void gradient_fft(fftw_complex *in, fftw_complex *out, std::vector<double> const &entree, std::vector<std::complex<double>> &gradx,std::vector<std::complex<double>> &grady, fftw_plan p_forward, fftw_plan p_backward);
void gradient(std::vector<double> src,std::vector<double> &grad, std::string direction, unsigned int dim);
void gradient_central(std::vector<double> src, std::vector<double> &grad,std::string direction,unsigned short int dim);
void deroul_volkov(fftw_complex *in, fftw_complex *out,std::vector<double> const &phase_enroul,std::vector<double> &phase_deroul , fftw_plan p_forward, fftw_plan p_backward);
void integ_grad(fftw_complex *in, fftw_complex *out, std::vector<double> const &gradx, std::vector<double> const&grady, std::vector<std::complex<double>> &sortie,fftw_plan p_forward, fftw_plan p_backward);
