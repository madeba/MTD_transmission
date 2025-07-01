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
//#include <Magick++.h>
#include <fftw3.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <tiffio.h>
#include <vector>
#include <complex>
//using namespace Magick;
//using namespace std;

#include "projet.h"
#include "manip.h"
//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
typedef struct{
clock_t init, fin;
float   total;
}temps;

/* --------------------------------------------------------------------------- */
// Prototypes
/* --------------------------------------------------------------------------- */
std::vector<std::complex<double>> extractSliceZ(std::vector<std::complex<double>> &index3D,std::string axis, size_t z_height);
float extract_val(std::string token,  std::string chemin_fic);
float extract_val(std::string token,  std::string chemin_fic, double defaut);
void interp_lin3D(std::vector <std::complex <double> > &volume_interp_3D);
//void divCplx(nbCplx* imageA,nbCplx* imageB,nbCplx* resultat, int NbPix);
//void multiplierCplx(nbCplx* image1,nbCplx* image2,nbCplx* resultat, int NbPix);
//void decal2DGen(double* entree, double* result, Var2D dim,Var2D decal);
//void decal2DCplxGen(vector<complex<double>> &entree, vector<complex<double>> &result, Var2D dim,Var2D decal);
Var2D corr_crois(double* obj2D_A, double* obj2D_B, Var2D dim);

void methodeCarre(int NbPixROI2d, double *holo1,  double *holo2,  double *holo3,  double *holo4);
void changeDim2D(double* tab, double* tabFinal, Var2D dimInit,Var2D dimFin);
void changeDim2DCplx(nbCplx *tab, nbCplx* tabFinal, Var2D dimInit, Var2D dimFin);
//double *circshift(double *entree,int dimx,int dimy,int decal_x,int decal_y);
double bruit(int attenuation);
void calcPhase2pi(nbCplx* obj, Var2D taille,double* phaseMod2pi);
//void calc_Uborn(nbCplx *TF_UBorn,nbCplx *UBorn,Var2D dim2DHA,Var2D PosSpec);
void Chrono(temps *t, std::string message);
void coupeCplx(nbCplx *src, nbCplx *dest, Var2D dim_src, Var2D dim_dest, Var2D coin);
void circshift2(double* entree, double* result, Var2D dim,Var2D decal);
void circshift3(double* entree, double* result, Var2D dim,Var2D decal);
void circshift2DCplx(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
void circshift3DCplx(std::vector <std::complex<double>> const &volume3D, std::vector <std::complex<double>> &volume3D_shift, Var3D dimFinal3D, Var3D decal3D);

int chargeBin(float *objet, std::string chemin,  int NbPix);

void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
int CreerZoneFresnel(double *FresnelRe,double * FresnelIm, Var2D dim, Var2D centre, float d, float lambda);
std::string extract_string(std::string const &token,  std::string chemin_fic);
void genere_rectang3D(nbCplx *objet,Var3D posI_Coin,Var3D dimRect, Var3D dim);
void genere_rectang2D(double *objet,Var2D posI_Coin,Var2D dimRect,Var2D dim);
void genere_OTF_RB_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon);
void genere_OTF_T_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon);
void genere_OTF_RH_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon);

void lire_bin(std::string chemin, double resultat[], short int precision, const size_t NbParam);
int get_bin_file_size(std::string chemin);
void Import3D_Tiff(std::vector<double> &imgTiff, std::string chemin, double taille_pixel);
void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY);
void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);
void ecrire_rapport(int NXMAX,float rayon,float Rf, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,std::string chemin,int nb_proj,float n1,float NA,float Tp, int G);

void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY);
//void concatener(char chaine1[],char chaine2[],char resultat[]);

void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY);
void multiplier_masque2Cplx(nbCplx image[], double masque[], int t_image, int t_mask, Var2D posCentre);
void multiplier_masqueCplx2(nbCplx *image, nbCplx *masque, int t_image, int t_mask, Var2D CentreI);
void prepare_wisdom2D(Var2D dim, const char *chemin);
void prepare_wisdom3D(Var3D dim, char *chemin);

void SAV2(double *var_sav, int NbPix2D, std::string chemin, enum PRECISION2 precision, char options[]);
//void SAVCplx(std::vector<std::complex<double> > var_sav, std::string partie, std::string chemin, enum PRECISION2 precision, char options[]);

void SAV3D_Tiff(std::vector<std::complex <double>> var_sav, std::string partie, std::string chemin, double taille_pixel);
void SAV3D_Tiff(std::vector<std::complex <double>> var_sav, Var3D const dim, std::string partie, std::string chemin, double taille_pixel);
void SAV3D_Tiff(std::vector<double> var_sav, std::string chemin, double taille_pixel);
void SAV3D_Tiff_Optimized(const std::vector<std::complex<double>> &var_sav, const std::string &partie, const std::string &chemin, double taille_pixel);
void SAV_Tiff2D(double *var_sav, std::string chemin, const size_t dim);

void SAV_Tiff2D(std::vector<double> var_sav, std::string chemin, double taille_pixel);
void SAV_Tiff2DCplx(std::vector<std::complex<double>> var_sav, std::string partie, std::string chemin, double taille_pixel);

double max(double* entree, int dim);
void phase2pi(nbCplx* obj, Var2D taille,double* WrappedImage);
//int retroPropag_Born(nbCplx *TF3D_PotObj, nbCplx *TF_Uborn_norm, double * sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon);
//int retroPropag_Born(vector <complex <double>> &TF3D_PotObj, vector <complex <double>> &TF_Uborn_norm, double * sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon);
void retroPropag_Born(std::vector <std::complex<double>> &TF3D_PotObj, std::vector<std::complex<double>> const &TF_Uborn_norm, std::vector<double> &sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon, manip m1);
void retroPropagSA(int deltaZ, nbCplx *fft_shift_norm, nbCplx *planObjet, Var3D decal, Var2D NMAX, double rayon);
void decalCoupeCplx(nbCplx *fft, nbCplx *fft_tmp, Var2D NMAX,Var2D dimCCD);
void Plan_ds_VolCplx(nbCplx *Vol3D, nbCplx *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di);
void Plan_ds_Vol(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di);
void Vol_ds_Plan(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di);
int coordMaxMod2D(nbCplx *entree, size_t const  tailleTab);
void InitTabCplx(nbCplx *z,int taille);
//void decal2DCplxGen(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal);
void decal2DCplxGen(std::vector<std::complex<double>> const &entree,std:: vector<std::complex<double>> &result, Var2D dim,Var2D decal);

void recale(nbCplx* obj,nbCplx* objDecal,nbCplx *objRecal, Var3D dimVol);

//void recal_obj(nbCplx *FFTobj, nbCplx *FFTobjDecal,nbCplx *objRecal, Var3D dimVol);
#endif

