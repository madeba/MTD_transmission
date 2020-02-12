#ifndef __FFT_FONCTIONS__
#define __FFT_FONCTIONS__

#include <vector>
#include <complex>
#include <fftw3.h>

std::vector<std::complex<double> >  fftshift2D(std::vector<std::complex<double> > &entree);
std::vector<double>  fftshift2D(std::vector<double> &entree);
//std::vector<vecteur>  fftshift2D(std::vector<vecteur> &entree);

void TF2D_vec(fftw_complex *in,fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);

void TF2Dcplx_vec(fftw_complex *in, fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2Dcplx_vec_INV(fftw_complex *in,fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p);

void TF2Dcplx_vec(fftw_complex *in,fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2Dcplx_vec_INV(fftw_complex *in, fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);

void TF2Dcplx_vec_INPLACE(fftw_complex *in_out, std::vector<double> entree,  std::vector< std::complex<double> > &sortie, fftw_plan p);

std::vector<double> tukey2D(int dimx,int dimy, float alpha);
#endif
