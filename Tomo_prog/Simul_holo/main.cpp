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
#include "FFTW_init.h" //gestion init fftw
#include "manip.h" //gestion manip

using namespace std;
using namespace cv;

/** @function main */
int main( int argc, char** argv )
{
    string home=getenv("HOME");
///--------------- Chargement ou création de l'objet (bille, spectre etc.)------------
    int dimROI=1024;
    manip m1(dimROI);//Initialiser la manip avec la taille du champ holographique (acquisition) en pixel
    Point3D dim3D(m1.dim_final,m1.dim_final,m1.dim_final);
    Point2D dim2D((double)m1.dim_Uborn,(double)m1.dim_Uborn,round(m1.dim_Uborn));
    unsigned int nbPix3D=pow(m1.dim_final,3),nbPix2D=pow(m1.dim_Uborn,2);
    ///variables 3D : objet+spectre
    vector<complex<double>> vol_bille(nbPix3D),  TF_bille(nbPix3D);
    vector<complex<double>> obj_conv(nbPix3D), SpectreObjConv(nbPix3D);

    ///variables 2D : champ complexe+spectre
    vector<complex<double>> TF_holo(nbPix2D), TF_holo_shift(nbPix2D);
    vector<complex<double>> holo(nbPix2D), holo_centre(nbPix2D);
    ///exportation des centres, pour contrôle
    vector<double> centres(nbPix2D);


    ///Génération de l'objet (bille polystyrène, n=1.5983, absorption=?)
    double Rboule_metrique=2.75*pow(10,-6);///rayon bille en m
    int rayon_boule_pix=round(Rboule_metrique/m1.Tp_Tomo);///rayon bille en pixel
    Point3D centre_boule(dim3D.x/2,dim3D.x/2,dim3D.x/2,dim3D.x);//bille centrée dans l'image
    double indice=1.45,kappa=0.105;//indice + coef d'extinction
    complex<double> nObj= {indice,kappa},n0= {m1.n0,0.0},Delta_n=nObj-n0;///init propriété bille
    genere_bille(vol_bille,centre_boule, rayon_boule_pix,nObj-n0,n0-n0,dim3D.x);
    SAV3D_Tiff((vol_bille),"Re",m1.chemin_result+"bille_Re.tif",m1.Tp_Tomo);
    SAV3D_Tiff((vol_bille),"im",m1.chemin_result+"bille_im.tif",m1.Tp_Tomo);

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
    FFTW_init tf3D(dim3D),tf2D(dim2D); ///init fftw pour spectre 3D et 2D
    TF3Dcplx(tf3D.in,tf3D.out,fftshift3D(vol_bille),TF_bille,tf3D.p_forward_OUT,m1.Tp_Tomo);
    vector<complex<double>>().swap(vol_bille);
    TF_bille=fftshift3D(TF_bille);
    //SAV3D_Tiff(TF_bille,"Re",m1.chemin_result+"TF_bille_Re.tif", m1.Delta_f_tomo*pow(10,-6));

    Point2D spec_H(0.0,0.0,m1.dim_Uborn);//coordonnées spec2D dans une image de dimension dimUBorn
    const int nbAngle=m1.nbHolo;
    cout<<"nbAngle="<<nbAngle;
    ///Generate OTF and 2D center (specular beam).
    //vector<Var2D> CoordSpec(m1.nbHolo);//table of specular coordinates
    vector<Point2D> CoordSpec_H(nbAngle, spec_H);
    short unsigned int const nbAxes=4;//nombre de branche de la fleur

    cout<<"------------------------Calcul OTF et objet convolué"<<endl;
    OTF mon_OTF(m1);//init OTF
    mon_OTF.bFleur(CoordSpec_H, nbAxes);//retrieve tabular of specular beam from class OTF & create 3D OTF;
    //CoordSpec2=mon_OTF.bFleur(nbAxes);//retrieve tabular of specular beam from class OTF & create 3D OTF;
    SAV3D_Tiff(mon_OTF.Valeur,"Re",m1.chemin_result+"OTF_simule_Re.tif",m1.Tp_Tomo);

    ///calcul objet convolué
    for(int cpt=0; cpt<pow(m1.dim_final,3); cpt++){
     //SpectreObjConv[cpt]=mon_OTF.Valeur[cpt]*TF_bille[cpt];
     SpectreObjConv[cpt].real(mon_OTF.Valeur[cpt].real()*TF_bille[cpt].real());
     SpectreObjConv[cpt].imag(mon_OTF.Valeur[cpt].imag()*TF_bille[cpt].imag());
    }
    vector<complex<double>>().swap(mon_OTF.Valeur);
    TF3Dcplx_INV(tf3D.in,tf3D.out,fftshift3D(SpectreObjConv),obj_conv,tf3D.p_forward_OUT,m1.Delta_f_tomo);
    vector<complex<double>>().swap(SpectreObjConv);
    SAV3D_Tiff(fftshift3D(obj_conv),"Re",m1.chemin_result+"obj_conv_Re.tif",m1.Tp_Tomo);
    SAV3D_Tiff(fftshift3D(obj_conv),"Im",m1.chemin_result+"obj_conv_Im.tif",m1.Tp_Tomo);
    vector<complex<double>>().swap(obj_conv);


    cout<<"extraction hologramme"<<endl;
    ///extract complex field from 3D spectrum
    for(int holo_numero=0; holo_numero<m1.nbHolo; holo_numero++){

        spec_H.x=CoordSpec_H[holo_numero].x;//the old code use spec, so we convert OTF.centre to spec
        spec_H.y=CoordSpec_H[holo_numero].y;

        int cpt2D=round(spec_H.y)*m1.dim_Uborn+round(spec_H.x);
        centres[spec_H.coordI().cpt2D()]=holo_numero;//used save centres in a image file, for quick visualisation

        ///Calcul des TF2D  des hologrammes à partir du spectre 3D de l'objet.
        calcHolo(spec_H,TF_bille,TF_holo,m1);//extract Ewald spheres + projection on a 2D plane
        // SAV2D_Tiff(TF_holo,"Im",dir_sav+"Tf_holo_Im.tif",m1.Tp_Uborn);
        SAVCplx(TF_holo,"Re",m1.chemin_result+"TF_holo_H_Re_208x208x600.raw",t_float,"a+b");
        Var2D decal2centreI= {-spec_H.x+m1.dim_Uborn/2,spec_H.y+m1.dim_Uborn/2};///centering all the spectrum in computer axes
       // SAVCplx(TF_holo,"Re",m1.chemin_result+"TF_holo_I_Re_208x208x600.raw",t_float,"a+b");
        decal2DCplxGen(TF_holo,TF_holo_shift,decal2centreI);
        ///calcul des hologrammes après recalage (hologrammes centrés)//calculation of hologramms, after shifting
        TF2Dcplx_INV(TF_holo_shift,holo,tf2D,m1.Delta_f_Uborn);
        SAVCplx(fftshift2D(holo),"Im",m1.chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw",t_double,"a+b");
        SAVCplx(fftshift2D(holo),"Re",m1.chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw",t_double,"a+b");
    }
    //sauver les centres sous forme d'image, pour contrôle
    SAV_Tiff2D(centres,m1.chemin_result+"centres_simul.tif",m1.Tp_holo);
    //repasser les spéculaires en coordonnées informatique (attendu par tomo_reconstruction)
    vector<double> tabPosSpec_I(m1.nbHolo*2);  ///nbhoo*2 coordonnées. stockage des speculaires pour exportation vers reconstruction
    for(int holo_numero=0;holo_numero<m1.nbHolo;holo_numero++){
        tabPosSpec_I[holo_numero]=CoordSpec_H[holo_numero].x+m1.dim_Uborn/2; //save tab_posSpec to a format readable by Tomo_reconstruction (computer )
        tabPosSpec_I[holo_numero+m1.nbHolo]=-CoordSpec_H[holo_numero].y+m1.dim_Uborn/2;
    }
    ///donnée utiles à la reconstruction
    SAV2(tabPosSpec_I,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");
    vector<double> param{m1.NXMAX,m1.nbHolo,m1.R_EwaldPix,dimROI,m1.Tp_holo};//devenu inutile avec fichier de config
    SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");
    return 0;
}
