#ifndef __IO_FONCTIONS__
#define __IO_FONCTIONS__

/* --------------------------------------------------------------------------- */
// includes...
/* --------------------------------------------------------------------------- */
//#include "struct.h"
//#include <time.h>
//#include <math.h>
//#include <iostream>
//#include <cstdlib>
//#include <strstream>
//#include <Magick++.h>
//#include <fftw3.h>
#include "FFTW_init.h"
//#include <core/core_c.h>
//#include <highgui/highgui_c.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

//#include <cstring>
//#include <fstream>
//#include <sstream>
//#include <assert.h>

//#include "Point2D.h"
//#include "vecteur.h"
#include <vector>
#include <complex>
#include <tiffio.h>
#include "projet.h"
#include <arrayfire.h>

///-----------entree/sortie+parseur---------------------------------
void deleteCplxField(std::string chemin_result, std::string dimImg);//delete complex field file
float extract_val(std::string token,  std::string chemin_fic);
std::string extract_string(std::string token,  std::string chemin_fic);
bool is_readable( const std::string & file);
void charger_image2D_OCV(std::vector<double> &imgTab, std::string imgFile, Var2D coin,Var2D taille);

void SAV2(std::vector<double>  &v, std::string chemin, enum PRECISION2 precision, char options[]);
void SAV2(cv::Mat &imgOpenCV, std::string chemin, enum PRECISION2 precision, char options[]);
void SAV2(double *var_sav, int NbPix2D, std::string chemin, enum PRECISION2 precision, char options[]);

void SAVCplx(std::vector<std::complex<double>> const &v, std::string partie, std::string chemin, enum PRECISION2 precision, char options[]);

void SAV3D_Tiff(std::vector<std::complex <double>> const &var_sav, std::string partie, std::string chemin, double taille_pixel);
void SAV_Tiff2D(std::vector<double> const &var_sav, std::string chemin, double taille_pixel);
void SAV_Tiff2D(std::vector<std::complex<double>> const &var_sav, std::string partie, std::string chemin, double taille_pixel);
void SAV_Tiff2D(af::array  af_img, std::string chemin, double taille_pixel);
void SAV2(cv::Mat const &imgCrop, std::string chemin, enum PRECISION2 precision, char options[]);
#endif
