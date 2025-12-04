#ifndef __FONCTIONS__
#define __FONCTIONS__

/* --------------------------------------------------------------------------- */
// incl
/* --------------------------------------------------------------------------- */
#include "struct.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <strstream>
#include <pthread.h>
#include <omp.h>
//#include <Magick++.h>
#include <fftw3.h>
//#include <core/core_c.h>
//#include <highgui/highgui_c.h>

#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/objdetect/objdetect.hpp>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "projet.h"
#include "vecteur.h"
#include <vector>
#include <complex>
#include "Point3D.h"
#include "OTF.h"
#include "manip.h"
//void SAV_Tiff2D2(vector<double> v, std::string chemin, unsigned int dim);

void delete_file(std::string file_path);
std::vector<double> wrap_phase(std::vector<double> &unwrapped_phase);
void interp_lin3D(std::vector <std::complex<double>> &volume_interp_3D);
Point2D maj_fleur(Point2D Vin, float rho, int nbHolo, double *theta, manip m1);
void calcHolo(Point2D spec,std::vector<std::complex<double>> const& TF_vol3D,std::vector<std::complex<double>> &TF_hologramme,manip const &m1);
void Conv_Ewald(Point2D spec,std::vector<std::complex<double>> TF_vol3D,std::vector<std::complex<double>> &TF_conv3D, manip m1);
void calcPhase_mpi_pi_atan2(vector<complex<double>> obj, vector<double> &phaseMod2pi);///calcul phase -PI-PI
void SAV2(std::vector<double> v, std::string chemin, enum PRECISION precision, char options[]);
void SAVCplx(std::vector<complex<double> > v, std::string partie, std::string chemin, enum PRECISION precision, char options[]);
void SAV3D_Tiff(vector <double> var_sav, string chemin, double taille_pixel);
void SAV3D_Tiff(vector<complex <double>> var_sav, string partie,string chemin, double taille_pixel);
void SAV2D_Tiff(std::vector<complex<double>> var_sav, string partie, string chemin,double taille_pixel);
void SAV_Tiff2D(std::vector<double> var_sav, string chemin, double taille_pixel);

void decal2DCplxGen(vector<complex<double>> &entree, vector<complex<double>> &result, Var2D decal);
string type2str(int type);
void gradient(std::vector<double> src,std::vector<double> &grad, string direction, unsigned int dim);
void gradient_central(std::vector<double> src, std::vector<double> &grad,string direction,unsigned short int dim);
void mat2vector(cv::Mat src,vector<double> &dst);
cv::Mat vector2mat(std::vector<double> &src, unsigned short dim);
void sobel_filtre2D(cv::Mat src, std::vector<double> &grad, string direction);
void lire_bin_vector(string chemin, vector<double> &dst, unsigned short int dim,  unsigned int NbPix);
vector<complex<double> > fftshift2D(vector<complex<double> > &entree);
void genere_bille(std::vector <std::complex<double>> &vol_bille, Point3D centre, size_t rayon,std::complex<double> delta_indice,size_t dim_espace);

void genere_barre(std::vector<std::complex<double>> &vol_obj,  Point3D coordMin, Point3D coordMax, std::complex<double> delta_indice, manip m1);
float extract_val(string token,  string chemin_fic);
float extract_val(string token,  string chemin_fic, double defaut);
string extract_string(std::string token,  std::string chemin_fic);
#endif
