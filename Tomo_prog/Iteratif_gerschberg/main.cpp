#include <iostream>
#include <fstream>
#include <sstream>
#define cimg_use_tiff
#define cimg_use_fftw3
#include <math.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include "CImg.h"
#include "fonctions.h"
#define NITER 20

using namespace cimg_library;
using namespace std;


CImg<float>indice;
CImg<float>abso;
CImg<float>masque_support;

int main(int argc, char *argv[])
{
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
    string chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    string fic_cfg_recon=repertoire_config+"/recon.txt";
    string indice_filename_str=chemin_result+"/indice.tif";
    string absorpt_filename_str=chemin_result+"/absorption.tif";
    string OTF_filename_str=chemin_result+"/OTF3D.tif";
    string resultat_indice_file_str=chemin_result+"/indice_GSP.tif";
    string resultat_kappa_file_str=chemin_result+"/kappa_GSP.tif";
    string resultat_OTF_GPS_str=chemin_result+"/OTF_GPS.tif";
    const float deltaNmax=extract_val("DELTA_NMAX",fic_cfg_recon);
    const float deltaNmin=extract_val("DELTA_NMIN",fic_cfg_recon);
    const float kappaMax=extract_val("KAPPA_MAX",fic_cfg_recon);
    const float kappaMin=extract_val("KAPPA_MIN",fic_cfg_recon);
    const int niterations=extract_val("NB_ITER_GPS",fic_cfg_recon);
    cout<<"fichier recon="<<fic_cfg_recon<<endl;

//  const char* indicefilename = cimg_option("-i",chemin_result.c_str()+"/indice.tif","Input index image file");
 // const char* absofilename = cimg_option("-i","/home/mat/tomo_test/absorption.tif","Input absorption image file");
 // const char* supportfilename = cimg_option("-sup","/home/mat/tomo_test/OTF3D.tif","Input support file");
 // const char* outputfilename   = cimg_option("-o","/home/mat/tomo_test/resultatindice.tif","Output image file");
    const char* outputfilename_indice   = resultat_indice_file_str.c_str();
    const char* outputfilename_kappa   = resultat_kappa_file_str.c_str();
    const char* outputfilename_OTF_GPS =resultat_OTF_GPS_str.c_str();

    const double sigma   = cimg_option("-s",1.0,"Standard variation of the gaussian smoothing");
    // const  int niterations  = cimg_option("-n",50,"Number of iterations");
    const bool hidden    = cimg_option("-hidden",false,0);      // This is a hidden option

    cout<<"Tiff:"<<TIFFLIB_VERSION<<endl;


   // indice.load_tiff(indicefilename);
    CImg<float> indice(indice_filename_str.c_str());
    abso.load_tiff(absorpt_filename_str.c_str());
    masque_support.load_tiff(OTF_filename_str.c_str());

    //indice.load_tiff("/media/bruno/Donnees/OSIRIS_NP-clean/cleaned/indice_deg_0.tif");
   // abso.load_tiff("/media/bruno/Donnees/OSIRIS_NP/absorption_deg_0.bin.tiff");
   // masque_support.load_tiff("/media/bruno/Donnees/OSIRIS_NP/sup_redon_C_deg_0.bin.tiff");


    if (!hidden) indice.display();

//   masque_support.shift(masque_support.width()/2,masque_support.height()/2,masque_support.depth()/2,0,2);

  if (!hidden) masque_support.display();
 //  indice.shift(indice.width()/2,indice.height()/2,indice.depth()/2,0,2);
   // CImgList<float> Image(indice*1.0,abso.fill(0));  //Estimation de l'objet initialisée avec le champ mesuré
    CImgList<float> Image(indice*1.0,abso*1.0);  //Estimation de l'objet initialisée avec le champ mesuré
    CImgList<float> F_mesure(Image);            //Champ mesuré

    F_mesure.FFT(false);

cout<<"GS-Papoulis.."<<endl;

    for(int i=1; i<niterations; i++)
    {
        cout<<"Iteration "<< i <<endl;

        cimg_forXYZ(Image[0],x,y,z)
        {
            if (Image[0](x,y,z)<deltaNmin) Image[0](x,y,z)=deltaNmin; // contrainte sur la variation d'indice: le fond est fixé à zero le delta est positif
            if (Image[0](x,y,z)>deltaNmax) Image[0](x,y,z)=deltaNmax;

            if (Image[1](x,y,z)<kappaMin) Image[1](x,y,z)=kappaMin; // contrainte sur le coefficient extinction
            if (Image[1](x,y,z)>kappaMax) Image[1](x,y,z)=kappaMax;
        }

        Image.FFT(false);

        if (i<niterations/2) cimg_forXYZ(Image[0],x,y,z)
        {
            if ((masque_support(x,y,z)>0)&&(x<masque_support.width()/2)) // contrainte du support : tout ce qui est dans le support est fixé à ce qui a été mesuré
            {                                                              // pour la première moitié des itérations, on réalise la contrainte uniquement sur la moitié du support (on évite la convergence vers la solution conjugée)
                Image[0](x,y,z)=F_mesure[0](x,y,z);
                Image[1](x,y,z)=F_mesure[1](x,y,z);
            }
        }

        if (i>=niterations/2) cimg_forXYZ(Image[0],x,y,z) // contrainte du support : tout ce qui est dans le support est fixé à ce qui a été mesuré
        {
            if (masque_support(x,y,z)>0)
            {
                Image[0](x,y,z)=F_mesure[0](x,y,z);
                Image[1](x,y,z)=F_mesure[1](x,y,z);
            }
        }


        Image.FFT(true);
      //  Image[1].fill(0);
        if ((i % 5) == 0) Image[0].save_tiff(outputfilename_indice);
    }
    Image[0].save_tiff(outputfilename_indice);
    Image[1].save_tiff(outputfilename_kappa);

    if (!hidden) Image.display();


    cout<<"OK";

    return 0;
}
