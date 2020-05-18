#ifndef __FFT_FONCTIONS__
#define __FFT_FONCTIONS__

#include <vector>
#include <complex>

#include <fftw3.h>

#include "FFTW_init.h"


std::vector<std::complex<double> >  fftshift2D(std::vector<std::complex<double> > const &entree);
std::vector<double>  fftshift2D(std::vector<double> const &entree);
//void fftshift2D(std::vector<double> const &entree, std::vector<double> &sortie);
//void  fftshift2D(std::vector<std::complex<double>> const &entree,std::vector<std::complex<double>>  &sortie);
///surcharge utilisant std:copy
std::vector<double>  fftshift2D2(std::vector<double> const &entree);
std::vector<std::complex<double>>  fftshift2D2(std::vector<std::complex<double>> const &entree);
std::vector<std::complex<double>> fftshift3D(std::vector<std::complex<double>> const &entree);

std::vector<double> tukey2D(int dimx,int dimy, float alpha);


///----------------------fonctions avec FFTW_init---------------


//------------------------------------------------------------------------------------------

void FFT2Dcplx(std::vector<double> const &entree, std::vector<std::complex<double>> &sortie, FFTW_init &tf2D_c2c);//calcul fft, real input
void TF2D(std::vector<double> const & entree, std::vector<std::complex<double> > &sortie, FFTW_init const &tf2D_c2c, double delta_x);//calcul TF, real input
void TF2Dcplx(std::vector<std::complex<double>> const & entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D_c2c, double delta_x);//calcul TF, cplx input
void FFT2Dcplx(std::vector<std::complex<double>> const & entree, std::vector<std::complex<double> > &sortie,FFTW_init &param_c2c);///->renommer fft?

void FFT2Dcplx_INV(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D_c2c);//calcul fft inv
void TF2Dcplx_INV(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double> > &sortie, FFTW_init &tf2D_c2c, double Delta_f);//calcul TF inv
void FFT3Dcplx_Inplace_INV(std::vector<std::complex<double>> const & entree, std::vector<std::complex<double> > &sortie, FFTW_init &param_Inplace_c2c);
void FFT3Dcplx_Inplace(std::vector<std::complex<double>> const & entree, std::vector<std::complex<double> > &sortie, FFTW_init &param_Inplace_c2c);
#endif
