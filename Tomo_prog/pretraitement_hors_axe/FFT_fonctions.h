#ifndef __FFT_FONCTIONS__
#define __FFT_FONCTIONS__

#include <vector>
#include <complex>

#include <fftw3.h>
#include "src/vecteur.h"
#include "FFTW_init.h"


std::vector<std::complex<double> >  fftshift2D(std::vector<std::complex<double> > const &entree);
std::vector<double>  fftshift2D(std::vector<double> const &entree);
void fftshift2D(std::vector<double> const &entree, std::vector<double> &sortie);
void  fftshift2D(std::vector<std::complex<double>> const &entree,std::vector<std::complex<double>>  &sortie);
std::vector<vecteur>  fftshift2D(std::vector<vecteur> &entree);

///fftw c2r complex input
void TF2D_vec(fftw_complex *in,fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);

void TF2Dcplx_vec(fftw_complex *in, fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2Dcplx_vec_INV(fftw_complex *in,fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p);

///fftw c2r  real  input

void TF2Dcplx_vec(std::vector<double> const &entree, std::vector<std::complex<double>> &sortie, FFTW_init &tf2D_Holo_c2r);
void TF2Dcplx_vec(fftw_complex *in,fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);
void TF2Dcplx_vec_INV(fftw_complex *in, fftw_complex *out, std::vector<double> entree, std::vector<std::complex<double> > &sortie, fftw_plan p);

void TF2Dcplx_vec_INPLACE(fftw_complex *in_out, std::vector<double> entree,  std::vector< std::complex<double> > &sortie, fftw_plan p);

std::vector<double> tukey2D(int dimx,int dimy, float alpha);


///----------------------fonctions avec FFTW_init---------------
void TF2D(std::vector<double> const & entree, std::vector<std::complex<double> > &sortie, FFTW_init const &tf2D, double delta_x);
void TF2D_r2c_symetric(std::vector<double> const &entree, std::vector<std::complex<double> > &sortie, FFTW_init  &tf2D_Re);
void TF2D_r2c(std::vector<double> const &entree, std::vector<std::complex<double> > &sortie, FFTW_init const &tf2D_Re, double delta_x);

void TF2Dcplx(std::vector<std::complex<double> > const & entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D, double delta_x);

void TF2Dcplx_INV(std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D, double Delta_f);
#endif
