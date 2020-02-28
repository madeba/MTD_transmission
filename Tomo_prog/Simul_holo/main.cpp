#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>  //setprecision
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>
#include "fonctions.h"
#include <complex>
#include "include/Point3D.h"
#include <fftw3.h>
#include "include/OTF.h"
#define pi M_PI
#include "Point3D.h"
#include "FFT_fonctions.h" //fonctions fftw
#include "FFT_encaps.h" //gestion init fftw
#include "manip.h" //gestion manip

using namespace std;
using namespace cv;

/** @function main */
int main( int argc, char** argv )
{
    string home=getenv("HOME");
    string dir_sav=home+"/tomo_test/";
    cout<<"dirsav="<<dir_sav<<endl;
///--------------- Chargement ou création de l'objet (bille, spectre etc.)------------
    int dimROI=512;
    manip m1(dimROI);//Initialiser la manip avec la taille du champ holographique (acquisition) en pixel
    Point3D dim3D(m1.dim_final,m1.dim_final,m1.dim_final);
    Point2D dim2D(m1.dim_Uborn,m1.dim_Uborn,m1.dim_Uborn);
    unsigned int nbPix3D=pow(m1.dim_final,3),nbPix2D=pow(m1.dim_Uborn,2);
///variables 3D
    vector<complex<double>> vol_bille(nbPix3D),  TF_bille(nbPix3D), OTF_shift(nbPix3D);
    vector<complex<double>> test(nbPix3D);
   // vector<complex<double>> spectre_conv(nbPix3D);
    vector<complex<double>> obj_conv(nbPix3D);
    vector<complex<double>> TF_convEwald(nbPix3D);
///variables 2D
    vector<complex<double>> TF_holo(nbPix2D);
    vector<complex<double>> TF_holo_shift(nbPix2D);
    vector<complex<double>> holo(nbPix2D);
    vector<double> centres(nbPix2D);
    vector<double> tabPosSpec(m1.nbHolo*2);  ///nbhoo*2 coordonnées. stockage des speculaires pour exportation vers reconstruction
///Génération de l'objet (bille polystyrène, n=1.5983, absorption nulle)
    double Rboule_metrique=2.75*pow(10,-6);///rayon bille en m
    int rayon_boule_pix=round(Rboule_metrique/m1.Tp_Tomo);///rayon bille en pixel
    Point3D centre_boule(dim3D.x/2,dim3D.x/2,dim3D.x/2,dim3D.x);//bille centrée dans l'image
    double indice=1.45,kappa=0.001;//indice + coef d'extinction
    complex<double> nObj= {indice,kappa},n0= {m1.n0,0.0},Delta_n=nObj-n0;///init propriété bille
    genere_bille(vol_bille,centre_boule, rayon_boule_pix,nObj-n0,n0-n0,dim3D.x);
    SAV3D_Tiff((vol_bille),"Re",dir_sav+"bille_Re.tif",m1.Tp_Tomo);

///--------------- Données physiques (en µm)----------------------------

    double phase_au_centre=Delta_n.real()*2*Rboule_metrique*2*pi/m1.lambda_v;
    cout<<setprecision(3);
    cout<<"*-------------Données physiques objet---------------*"<<endl;
    cout<<"| Dim            | "<<m1.dim_Uborn<<" pixel                         |"<<endl;
    cout<<"| Rayon boule    | "<<rayon_boule_pix*m1.Tp_Tomo*pow(10,6)<<" µm                           |"<<endl;
    cout<<"| Delta_n        | ("<<Delta_n.real()<<","<<Delta_n.imag()<<")                        |"<<endl;
    cout<<"| phase au centre| "<<phase_au_centre<<" rad                          |"<<endl;
    cout<<"*---------------------------------------------------*"<<endl;





///Calcul spectre objet
    FFT_encaps tf3D(dim3D),tf2D(dim2D); ///init fftw pour spectre 3D et 2D
    TF3Dcplx(tf3D.in,tf3D.out,fftshift3D(vol_bille),TF_bille,tf3D.p_forward_OUT,m1.Tp_Tomo);
    SAV3D_Tiff(fftshift3D(TF_bille),"Re",dir_sav+"TF_bille_Re.tif", m1.Delta_f_tomo*pow(10,-6));

///Balayage : onde incidente k_inc
    double phi_i=0,theta_i=0;//angle polaire et azimuthal.

    Point2D f_inc(round(m1.R_EwaldPix*cos(phi_i)*sin(theta_i)),round(m1.R_EwaldPix*sin(phi_i)*sin(theta_i)),dim3D.x);//coordonnées ki
    Point2D spec(round(m1.R_EwaldPix*cos(phi_i)*sin(theta_i)),round(m1.R_EwaldPix*sin(phi_i)*sin(theta_i)),m1.dim_Uborn);//corrdonnée spec2D identique ki lateral mais definit en 2D pourt cpt2D
///boucle de balayage sur les 400 hologrammes
   for(int holo_numero=1;holo_numero<m1.nbHolo+1;holo_numero++)
   {
      spec=maj_fleur(spec,m1.NXMAX,m1.nbHolo,&theta_i,m1);//rosace
      tabPosSpec[holo_numero]=(double)spec.x;
      tabPosSpec[holo_numero+m1.nbHolo]=(double)spec.y;
     //cout<<"spec.x="<<spec.x<<"| spec.y="<<spec.y<<endl;
     int cpt2D=round(spec.y)*m1.dim_Uborn+spec.x;
    // cout<<cpt2D<<endl;
     if(holo_numero<200)
      centres[spec.coordI().cpt2D()]=holo_numero;
///Calcul des TF2D et des hologrammes à partir du spectre 3D de l'objet.
    calcHolo(spec,fftshift3D(TF_bille),TF_holo,m1);//extraction des sphères d'Ewald+placage 2D
   // SAV2D_Tiff(TF_holo,"Im",dir_sav+"Tf_holo_Im.tif",m1.Tp_Uborn);
      Var2D decal={-f_inc.x-m1.dim_Uborn/2,f_inc.y-m1.dim_Uborn/2};

  decal2DCplxGen(TF_holo,TF_holo_shift,decal);
    //TF2Dcplx_INV(fftshift2D(TF_holo_shift),holo,tf2D,m1.Delta_f_Uborn);
     TF2Dcplx_INV(TF_holo_shift,holo,tf2D,m1.Delta_f_Uborn);
   // SAV2D_Tiff(TF_holo_shift,"Mod",dir_sav+"Tf_holo_Mod_shift.tif",m1.Tp_Uborn);
    //SAV2D_Tiff(fftshift2D(holo),"Im",dir_sav+"holo_Im.tif",m1.Tp_Uborn);
    SAVCplx(fftshift2D(holo),"Im",dir_sav+"holo_Im.raw",t_double,"a+b");
   }
 SAV_Tiff2D(centres,dir_sav+"centres.tif",m1.Tp_holo);
///calcul objet convolué
    cout<<"nb_holo="<<m1.nbHolo<<endl;
  /*  OTF mon_OTF(dimROI, m1);
    mon_OTF.bFleur();

    OTF_shift=fftshift3D(mon_OTF.Valeur);
    SAV3D_Tiff(fftshift3D(OTF_shift),"Re",dir_sav+"OTF_shift_Re.tif",m1.Tp_Tomo);
    for(int cpt=0;cpt<pow(m1.dim_final,3);cpt++)
        test[cpt]=OTF_shift[cpt]*TF_bille[cpt];

    TF3Dcplx_INV(tf3D.in,tf3D.out,test,obj_conv,tf3D.p_forward_OUT,m1.Delta_f_tomo);
    SAV3D_Tiff(fftshift3D(obj_conv),"Re",dir_sav+"obj_conf_Re.tif",m1.Tp_Tomo);
*/
    SAV2(tabPosSpec,dir_sav+"/tab_posSpec.raw",t_double,"wb");
    return 0;
}
