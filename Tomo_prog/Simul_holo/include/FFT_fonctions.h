#ifndef __FFT_FONCTIONS__
#define __FFT_FONCTIONS__

#include <vector>
#include <complex>
#include <fftw3.h>
#include "manip.h"
#include "FFTW_init.h"
#include "vecteur.h"
#include "struct.h"
//std::vector<std::complex<double> >  fftshift2D(std::vector<std::complex<double> > &entree);
std::vector<vecteur> init_kvect_shift(Var2D dim2DHA);
///--------------------------overload for fftshift----------------------------------------------
std::vector<double>  fftshift2D(std::vector<double> &entree);
std::vector<std::complex<double>> fftshift2D(std::vector<std::complex<double>> &entree);
std::vector<std::complex<double>> fftshift3D(std::vector<std::complex<double>> &entree);
std::vector<vecteur>  fftshift2D(std::vector<vecteur> &entree);

///-----------------------2D FFT and exact FT----------------------------------------------
void TF2Dcplx(std::vector<double> const &entree, std::vector<std::complex<double>> &sortie, FFTW_init &tf2D_c2r);//fft real input
void TF2Dcplx(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double>> &sortie,FFTW_init &param_c2c);//fft complex input
void TF2D(std::vector<double>  entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D, double delta_x);

void TF2Dcplx(std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D, double delta_x);//exact FT, complex input

void TF2Dcplx_INV(std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D, double Delta_f);
void TF2Dcplx_INV(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D);


void TF2Dcplx_INPLACE(fftw_complex *in_out, std::vector<double> entree,  std::vector< std::complex<double> > sortie, fftw_plan p, double delta_x);

///-----------------------3D FFT----------------------------------------------
void TF3Dcplx_inplace(fftw_complex *in_out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p3d, double delta_x);
void TF3Dcplx(fftw_complex *in, fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p3d, double delta_x);
void TF3Dcplx_INV(fftw_complex *in, fftw_complex *out, std::vector<std::complex<double> > entree, std::vector<std::complex<double> > &sortie, fftw_plan p3d, double delta_f);


std::vector<double> tukey2D(int dimx,int dimy, float alpha);
#endif
