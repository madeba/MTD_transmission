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

    // vector<complex<double>> spectre_conv(nbPix3D);
    vector<complex<double>> obj_conv(nbPix3D);
    vector<complex<double>> SpectreObjConv(nbPix3D);
    //vector<complex<double>> TF_convEwald(nbPix3D);
///variables 2D : champ complexe
    vector<complex<double>> TF_holo(nbPix2D);
    vector<complex<double>> TF_holo_shift(nbPix2D);
    vector<complex<double>> holo(nbPix2D);
    vector<complex<double>> holo_centre(nbPix2D);
    vector<double> centres(nbPix2D);
    vector<double> tabPosSpec(m1.nbHolo*2);  ///nbhoo*2 coordonnées. stockage des speculaires pour exportation vers reconstruction

///Génération de l'objet (bille polystyrène, n=1.5983, absorption=?)
    double Rboule_metrique=2.75*pow(10,-6);///rayon bille en m
    int rayon_boule_pix=round(Rboule_metrique/m1.Tp_Tomo);///rayon bille en pixel
    Point3D centre_boule(dim3D.x/2,dim3D.x/2,dim3D.x/2,dim3D.x);//bille centrée dans l'image
    double indice=1.45,kappa=0.001;//indice + coef d'extinction
    complex<double> nObj= {indice,kappa},n0= {m1.n0,0.0},Delta_n=nObj-n0;///init propriété bille
    genere_bille(vol_bille,centre_boule, rayon_boule_pix,nObj-n0,n0-n0,dim3D.x);
    //  SAV3D_Tiff((vol_bille),"Re",m1.chemin_result+"bille_Re.tif",m1.Tp_Tomo);

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
//   SAV3D_Tiff(TF_bille,"Re",m1.chemin_result+"TF_bille_Re.tif", m1.Delta_f_tomo*pow(10,-6));

    Point2D spec(0.0,0.0,m1.dim_Uborn);//corrdonnée spec2D identique ki lateral mais definit
    //en 2D pour cpt2D
    ///Generate OTF and 2D center (specular beam).
    //vector<Var2D> CoordSpec(m1.nbHolo);//table of specular coordinates

    vector<Point2D> CoordSpec2(m1.nbHolo, spec);

    cout<<"creation point 2D m1.nbHolo="<<m1.nbHolo<<endl;
    cout<<"------------------------Calcul OTF et objet convolué"<<endl;
    OTF mon_OTF(m1);//init OTF
    //mon_OTF.bFleur(CoordSpec2);//retrieve tabular of specular beam from class OTF & create 3D OTF;
    mon_OTF.bFleur(CoordSpec2);//retrieve tabular of specular beam from class OTF & create 3D OTF;
   //for(int cpt=0;cpt<m1.nbHolo;cpt++)
     //   cout<<"main, cpt="<<cpt<<" : CoorSpec=("<<CoordSpec2[cpt].x<<","<<CoordSpec2[cpt].y<<")"<<endl;
    SAV3D_Tiff(mon_OTF.Valeur,"Re",m1.chemin_result+"OTF_Re.tif",m1.Tp_Tomo);

     ///calcul objet convolué
    //OTF_shift=fftshift3D(mon_OTF.Valeur);

    for(int cpt=0; cpt<pow(m1.dim_final,3); cpt++)
        SpectreObjConv[cpt]=mon_OTF.Valeur[cpt]*TF_bille[cpt];
   vector<complex<double>>().swap(mon_OTF.Valeur);
    TF3Dcplx_INV(tf3D.in,tf3D.out,fftshift3D(SpectreObjConv),obj_conv,tf3D.p_forward_OUT,m1.Delta_f_tomo);
   vector<complex<double>>().swap(SpectreObjConv);
    SAV3D_Tiff(fftshift3D(obj_conv),"Re",m1.chemin_result+"obj_conv_Re.tif",m1.Tp_Tomo);
   vector<complex<double>>().swap(obj_conv);
    ///extract complex field from 3D spectrum
    cout<<"extraction hologramme"<<endl;
    for(int holo_numero=0; holo_numero<m1.nbHolo; holo_numero++){
        cout<<"holo_numero="<<holo_numero<<endl;
        spec.x=CoordSpec2[holo_numero].x;//the old code use spec, so we convert OTF.centre to spec
        spec.y=CoordSpec2[holo_numero].y;


       cout<<"main spec=("<<spec.x<<","<<spec.y<<")"<<endl;
        tabPosSpec[holo_numero]=(double)CoordSpec2[holo_numero].x; //save tab_posSpec to a format readable by Tomo_reconstruction
        tabPosSpec[holo_numero+m1.nbHolo]=(double)spec.y;
        int cpt2D=round(spec.y)*m1.dim_Uborn+round(spec.x);
        centres[spec.coordI().cpt2D()]=holo_numero;//used save centres in a image file, for quick visualisation

        ///Calcul des TF2D  des hologrammes à partir du spectre 3D de l'objet.
        calcHolo(spec,TF_bille,TF_holo,m1);//extract Ewald spheres + projection on a 2D plane
        // SAV2D_Tiff(TF_holo,"Im",dir_sav+"Tf_holo_Im.tif",m1.Tp_Uborn);
        Var2D decal= {-spec.x-m1.dim_Uborn/2,spec.y-m1.dim_Uborn/2};

        decal2DCplxGen(TF_holo,TF_holo_shift,decal);
        //TF2Dcplx_INV(fftshift2D(TF_holo_shift),holo,tf2D,m1.Delta_f_Uborn);
        ///calcul des hologrammes après recalage (hologrammes centrés)//calculation of hologramms, after shifting
        TF2Dcplx_INV(TF_holo_shift,holo,tf2D,m1.Delta_f_Uborn);
        SAVCplx(fftshift2D(holo),"Im",m1.chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw",t_double,"a+b");
        SAVCplx(fftshift2D(holo),"Re",m1.chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw",t_double,"a+b");
    }

///boucle de balayage sur les 400 hologrammes : ancienne version avec mise à jour iterative de la rosace
    /*   for(int holo_numero=1; holo_numero<m1.nbHolo+1; holo_numero++)
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
           calcHolo(spec,TF_bille,TF_holo,m1);//extraction des sphères d'Ewald+placage 2D

           Var2D decal= {-f_inc.x-m1.dim_Uborn/2,f_inc.y-m1.dim_Uborn/2};

           decal2DCplxGen(TF_holo,TF_holo_shift,decal);
            SAVCplx(TF_holo,"Im",m1.chemin_result+"Tf_holo_Im_122x122x200x32"+".raw",t_float,"a+b");
            SAVCplx(TF_holo_shift,"Im",m1.chemin_result+"Tf_holo_shift_Im_122x122x200x32"+".raw",t_float,"a+b");
          // TF2Dcplx_INV(fftshift2D(TF_holo_shift),holo,tf2D,m1.Delta_f_Uborn);
           TF2Dcplx_INV(TF_holo_shift,holo,tf2D,m1.Delta_f_Uborn);
           // SAV2D_Tiff(TF_holo_shift,"Mod",dir_sav+"Tf_holo_Mod_shift.tif",m1.Tp_Uborn);
           //SAV2D_Tiff(fftshift2D(holo),"Im",dir_sav+"holo_Im.tif",m1.Tp_Uborn);

           SAVCplx(fftshift2D(holo),"Im",m1.chemin_result+"holo_Im_122x122x200x32"+".raw",t_float,"a+b");
       }*/
    SAV_Tiff2D(centres,m1.chemin_result+"centres.tif",m1.Tp_holo);



    SAV2(tabPosSpec,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");
    vector<double> param{m1.NXMAX,m1.nbHolo,m1.R_EwaldPix,dimROI,m1.Tp_holo};
    SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");
    return 0;
}
