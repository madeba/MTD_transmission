#ifndef __FONCTIONS__
#define __FONCTIONS__

/* --------------------------------------------------------------------------- */
// incl
/* --------------------------------------------------------------------------- */
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <strstream>
#include <Magick++.h>
#include <fftw3.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>


using namespace Magick;
using namespace std;

#include "projet.h"

/* --------------------------------------------------------------------------- */
// Prototypes
/* --------------------------------------------------------------------------- */

void changeDim2D(double* tab, double* tabFinal, Var2D dimInit,Var2D dimFin);
//double *circshift(double *entree,int dimx,int dimy,int decal_x,int decal_y);

void circshift2(double* entree, double* result, Var2D dim,Var2D decal);
void circshift3(double* entree, double* result, Var2D dim,Var2D decal);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);

//void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille);
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y);
void TF2D_INV(double entree_reelle[],double entree_imag[], double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y);
double *tukey2D(int dimx,int dimy, float alpha);
void charger_image2D(unsigned char* phasei, string chemin,Var2D coin,Var2D taille);

void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY);
void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);
void ecrire_rapport(int NXMAX,float rayon,float Rf, float K, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,string chemin,int nb_proj,float n1,float NA,float Tp, int G);
//void fftw_plan_with_nthreads(int nthreads);
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY);
//void concatener(char chaine1[],char chaine2[],char resultat[]);
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy);
void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY);

void crop_window(double* img_src, size_t src_dimx, size_t src_dimy, double* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y);
void SAV(double *var_sav, int taille, char* chemin, enum PRECISION,char options[]);

#endif
