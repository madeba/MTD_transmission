#ifndef __FFT_FONCTIONS__
#define __FFT_FONCTIONS__

#include <vector>
#include <complex>
#include <fftw3.h>
#include "src/vecteur.h"
#include "FFT_encaps.h"

std::vector<std::complex<double> >  fftshift2D(std::vector<std::complex<double> > &entree);
std::vector<double>  fftshift2D(std::vector<double> &entree);
std::vector<vecteur>  fftshift2D(std::vector<vecteur> &entree);

void TF2D_vec(fftw_complex *in,fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2D_r2c_symetric(std::vector<double> const &entree, std::vector<std::complex<double> > &sortie, FFT_encaps const &tf2D_Re, double delta_x);

void TF2Dcplx_vec(fftw_complex *in,fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2Dcplx_vec_INV(fftw_complex *in, fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
///Surcharge pour utiliser la classe encapsulant la FFTW
void TF2Dcplx_vec(std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie,FFT_encaps &tf2D);
//surcharge pour utiliser double
void TF2Dcplx_vec(fftw_complex *in, fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2Dcplx_vec_INV(fftw_complex *in,fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p);


//test avec l'option FFTW INPLACE
void TF2Dcplx_vec_INPLACE(fftw_complex *in_out, std::vector<double> entree,  std::vector< std::complex<double> > &sortie, fftw_plan p);

std::vector<double> tukey2D(int dimx,int dimy, float alpha);
#endif
