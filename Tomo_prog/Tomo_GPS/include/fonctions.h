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
//#include "projet.h"
#include "IO_fonctions.h"
#include "Point2D.h"
//#include "vecteur.h"
#include <vector>
#include <complex>
#include <tiffio.h>
//*calculate the wrapped phase from -pi to pi*/


Var2D coord_to_coordShift(Var2D coord2D, Var2D dimROI);
void calcPhase_mpi_pi_atan2(std::vector<std::complex<double>> const& cplx_field, std::vector<double> &phaseMod2pi);
double calc_quad_err(std::vector <std::complex <double> > &indice_actual_iter, std::vector< std::complex <double> > &indice_prev_iter);
void  decal2DCplxGen2(std::vector<std::complex<double> > const &entree,std::vector <std::complex<double> >  &result, Var2D &decalGen);
std::vector<double> initRef(std::string chemin_ref,Var2D coin, Var2D dimROI);
void calc_Uborn2(std::vector<std::complex<double>> const &TF_UBorn,std::vector<std::complex<double>> &UBorn,Var2D dim2DHA,Var2D PosSpec,FFTW_init &param_c2c);
int coordSpec(std::vector<std::complex<double>> const &TF_UBorn, std::vector<double> &TF_champMod,Var2D NMAX);
///extract complex field form off-axis hologram
void holo2TF_UBorn2(std::vector<double>  &holo1,std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, std::vector<double> const &tukeyHolo,FFTW_init  &FFTparam_fftw2DHolo);
//////extract complex field form off-axis hologram, but shifted (zero frequency on the top left). Warning you have to pass a shifted coin_HA !
void holo2TF_UBorn2_shift(std::vector<double>  &holo1,std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA_shift, size_t NbAngleOk, std::vector<double> const &tukeyHolo,FFTW_init  &FFTparam_fftw2DHolo);
void holo2TF_UBorn2_shift_r2c(std::vector<double>  &holo1,std::vector<std::complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA_shift, size_t NbAngleOk, std::vector<double> const &tukeyHolo,FFTW_init  &param_fftw2D_r2c_Holo);
void coupeCplx(std::vector<std::complex<double>> const &src, Var2D dim_src, std::vector<std::complex<double>> &dest, Var2D dim_dest, Var2D coin,size_t NumAngle);
void coupe2D_RefI_to3D(std::vector<std::complex<double>> const &src, std::vector<std::complex<double>> &dest, Var2D dim_dest, unsigned short int numAngle);
void coupe2D_I_to_H3D(std::vector<std::complex<double>> const &src2D, std::vector<std::complex<double>> &dest3D,Var2D dim_dest2D, unsigned short int numAngle);
///-----------entree/sortie+parseur---------------------------------
/*float extract_val(std::string token,  std::string chemin_fic);
std::string extract_string(std::string token,  std::string chemin_fic);
bool is_readable( const std::string & file )  ;

void charger_image2D_OCV(std::vector<double> &imgTab, std::string imgFile, Var2D coin,Var2D taille);

void SAV2(std::vector<double>  &v, std::string chemin, enum PRECISION2 precision, char options[]);
void SAV2(double *var_sav, int NbPix2D, std::string chemin, enum PRECISION2 precision, char options[]);
void SAVCplx(std::vector<std::complex<double>> const &v, std::string partie, std::string chemin, enum PRECISION2 precision, char options[]);

void SAV3D_Tiff(std::vector<std::complex <double>> const &var_sav, std::string partie, std::string chemin, double taille_pixel);
void SAV_Tiff2D(std::vector<double> const &var_sav, std::string chemin, double taille_pixel);
void SAV_Tiff2D(std::vector<std::complex<double>> const &var_sav, std::string partie, std::string chemin, double taille_pixel);

*/


#endif
