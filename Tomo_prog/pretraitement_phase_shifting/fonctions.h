#ifndef __FONCTIONS__
#define __FONCTIONS__

/* --------------------------------------------------------------------------- */
// incl
/* --------------------------------------------------------------------------- */
#include "struct.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <tiffio.h>
#include <cstdlib>
#include <strstream>
//#include <Magick++.h>
#include <fftw3.h>
//#include <core/core_c.h>
//#include <highgui/highgui_c.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "projet.h"
//#include "vecteur.h"
#include <vector>
#include <complex>
#include "manip.h"
void antigaussienne(std::vector<std::complex<double>> &tab, int sigma, float A, int Exy);
void calcPhase_mpi_pi_atan2(std::vector<std::complex<double>> obj, std::vector<double> &phaseMod2pi);
void circshift2DCplx(std::vector<std::complex<double>> entree, std::vector<std::complex<double>> &result, Var2D dim,Var2D decal);
void decal2DCplxGen(std::vector<std::complex<double>> entree,std::vector<std::complex<double>> &result, Var2D dim,Var2D decal);
void calc_Uborn(std::vector<std::complex<double>> TF_UBorn,std::vector<std::complex<double>> &UBorn,Var2D dim2DHA,Var2D PosSpec,fftw_complex *in,fftw_complex *out,fftw_plan p);
int coordSpec(std::vector<std::complex<double>> TF_UBorn, std::vector<double> &TF_champMod,Var2D NMAX);
//void holo2TF_UBorn(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA,
                 //  size_t NumAngle, std::vector<double> masque);std::
void holo2TF_UBorn_PS(std::vector<std::complex <double>> holo1, std::vector<std::complex<double>> &TF_UBornTot, size_t NumAngle, std::vector<double> tukey_holo, fftw_complex *in,fftw_complex *out,fftw_plan p_forward_holo,manip m1);

void holo2TF_UBorn_old(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, std::vector<double> tukey_holo);
void holo2TF_UBorn(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle,
                   std::vector<double> tukey_holo, fftw_complex *in, fftw_complex *out,fftw_plan p_forward_holo);
void holo2TF_UBorn_INPLACE(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle,
                   std::vector<double> tukey_holo, fftw_complex *in_out,fftw_plan p_forward_holo);

void SAV2(std::vector<double> v, std::string chemin, enum PRECISION2 precision, char options[]);
void SAV2(double *var_sav, int NbPix2D, std::string chemin, enum PRECISION2 precision, char options[]);
void SAVCplx(std::vector<std::complex<double> > v, std::string partie, std::string chemin, enum PRECISION2 precision, char options[]);
void SAV_Tiff2D(std::vector<double> var_sav,  std::string chemin, double taille_pixel);
void SAV_Tiff2DCplx(std::vector<std::complex<double>> var_sav, std::string partie, std::string chemin, double taille_pixel);
void SAV3D_Tiff(std::vector<std::complex <double>> var_sav, std::string partie, std::string chemin, double taille_pixel);
void charger_image2D_OCV(std::vector<double> &imgTab, std::string imgFile, Var2D coin,Var2D taille);

void coupeCplx(std::vector<std::complex<double>> src, std::vector<std::complex<double>> &dest, Var2D coin);
void coupeCplx2Stack(std::vector<std::complex<double>> src, std::vector<std::complex<double>> &dest, Var2D dim_dest, Var2D coin, size_t NumAngle);

float extract_val(std::string token,  std::string chemin_fic);
std::string extract_string(std::string token,  std::string chemin_fic);



#endif
