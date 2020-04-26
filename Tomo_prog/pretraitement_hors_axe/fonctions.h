#ifndef __FONCTIONS__
#define __FONCTIONS__

/* --------------------------------------------------------------------------- */
// includes...
/* --------------------------------------------------------------------------- */
#include "struct.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <strstream>
//#include <Magick++.h>
#include <fftw3.h>
#include "FFTW_init.h"
//#include <core/core_c.h>
//#include <highgui/highgui_c.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "projet.h"

#include "Point2D.h"
//#include "vecteur.h"
#include <vector>
#include <complex>
#include <tiffio.h>
//*calculate the wrapped phase from -pi to pi*/
void calcPhase_mpi_pi_atan2(std::vector<std::complex<double>> const& cplx_field, std::vector<double> &phaseMod2pi);
void circshift2DCplx(std::vector<std::complex<double>> const &entree, std::vector<std::complex<double>> &result, Var2D dim,Var2D decal);
void decal2DCplxGen(std::vector<std::complex<double>> const &entree,std::vector<std::complex<double>> &result, Var2D dim,Var2D decal);
void calc_Uborn(std::vector<std::complex<double>> const &TF_UBorn,std::vector<std::complex<double>> &UBorn,Var2D dim2DHA,Var2D PosSpec,fftw_complex *in,fftw_complex *out,fftw_plan p);
void calc_Uborn2(std::vector<std::complex<double>> const &TF_UBorn,std::vector<std::complex<double>> &UBorn,Var2D dim2DHA,Var2D PosSpec,FFTW_init &param_c2c);
int coordSpec(std::vector<std::complex<double>> const &TF_UBorn, std::vector<double> &TF_champMod,Var2D NMAX);
//void holo2TF_UBorn(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA,
                 //  size_t NumAngle, std::vector<double> masque);
//void holo2TF_UBorn_old(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, std::vector<double> tukey_holo);
void holo2TF_UBorn(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle,
                   std::vector<double> tukey_holo, fftw_complex *in, fftw_complex *out,fftw_plan p_forward_holo);

void holo2TF_UBorn(std::vector<double> &holo1, std::vector<std::complex<double>> &TF_UBornTot, Var2D dimROI, Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, std::vector<double> const &masqueTukeyHolo, FFTW_init &tf2D_Holo_c2r);//surcharge FFTW_init pour c2r
//void holo2TF_UBorn_INPLACE(std::vector<double> holo1, std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle,
                  // std::vector<double> tukey_holo, fftw_complex *in_out,fftw_plan p_forward_holo);

//void holo2TF_UBorn2(std::vector<double> const holo1,std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, std::vector<double> const &masqueTukeyHolo,
//                      FFTW_init const &FFTparam_fftw2DHolo);
///extract complex field form off-axis hologram
void holo2TF_UBorn2(std::vector<double>  &holo1,std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, std::vector<double> const &tukeyHolo,FFTW_init  &FFTparam_fftw2DHolo);
//void holo2TF_UBorn2(std::vector<double>  &holo1,std::vector<complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, std::vector<double> const &tukeyHolo,FFTW_init const &FFTparam_fftw2DHolo);
void SAV2(std::vector<double>  &v, std::string chemin, enum PRECISION2 precision, char options[]);
void SAV2(double *var_sav, int NbPix2D, std::string chemin, enum PRECISION2 precision, char options[]);
void SAVCplx(std::vector<std::complex<double>> const &v, std::string partie, std::string chemin, enum PRECISION2 precision, char options[]);
void charger_image2D_OCV(std::vector<double> &imgTab, std::string imgFile, Var2D coin,Var2D taille);

void coupeCplx(std::vector<std::complex<double>> const &src, Var2D dim_src, std::vector<std::complex<double>> &dest, Var2D dim_dest, Var2D coin,size_t NumAngle);

float extract_val(std::string token,  std::string chemin_fic);
std::string extract_string(std::string token,  std::string chemin_fic);
bool is_readable( const std::string & file )  ;
void SAV3D_Tiff(std::vector<std::complex <double>> const &var_sav, std::string partie, std::string chemin, double taille_pixel);
void SAV_Tiff2D(std::vector<double> const &var_sav, std::string chemin, double taille_pixel);
void SAV_Tiff2D(std::vector<std::complex<double>> const &var_sav, std::string partie, std::string chemin, double taille_pixel);



#endif
