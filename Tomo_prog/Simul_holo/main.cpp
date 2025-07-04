//######Generate a 3D complex object (bead or box), calculate its spectrum and extract hologram thanks to Ewald sphere.######//
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>  //setprecision
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include <stdlib.h>
#include <stdio.h>
#include "fonctions.h"
#include "tiff_functions.h"
#include <complex>
#include "include/Point3D.h"
#include <fftw3.h>
#include "include/OTF.h"
#define pi M_PI
#include "Point3D.h"
#include "FFT_fonctions.h" //fonctions fftw
#include "FFTW_init.h" //gestion init fftw
#include "manip.h" //gestion manip
#include "champCplx_functions.h"
using namespace std;
using namespace cv;

/** @function main */
int main( int argc, char** argv )
{
    string home=getenv("HOME");
///--------------- Chargement ou création de l'objet (bille, spectre etc.)------------
    int dimROI=2048;
    manip m1(dimROI);//Initialiser la manip avec la taille du champ holographique (acquisition) en pixel. Attention, cela effacera la valeur du fichier de config !
    Point3D dim3D(m1.dim_final,m1.dim_final,m1.dim_final);
    Var3D dim={m1.dim_final,m1.dim_final,m1.dim_final};
    Point2D dim2D((double)m1.dim_Uborn,(double)m1.dim_Uborn,round(m1.dim_Uborn));
    unsigned int nbPix3D=pow(m1.dim_final,3),nbPix2D=pow(m1.dim_Uborn,2);
    ///variables 3D : objet+spectre
    vector<complex<double>> vol_obj(nbPix3D),  TF_obj(nbPix3D);
    vector<complex<double>> obj_conv(nbPix3D), SpectreObjConv(nbPix3D);

    ///variables 2D : champ complexe+spectre
    vector<complex<double>> TF_holo(nbPix2D), TF_holo_shift(nbPix2D);
    vector<complex<double>> holo(nbPix2D), holo_centre(nbPix2D);
    ///exportation des centres, pour contrôle
    vector<double> centres(nbPix2D);


    ///Génération de l'objet (bille polystyrène, n=1.5983, absorption=?)
    double Rboule_metrique=5*pow(10,-6);///rayon bille en m
    int rayon_boule_pix=round(Rboule_metrique/m1.Tp_Tomo);///rayon bille en pixel
    Point3D centre_boule(dim3D.x/2,dim3D.x/2,dim3D.x/2,dim3D.x);//bille centrée dans l'image
    double indice=1.395,kappa=0.000;//indice + coef d'extinction
    complex<double> nObj= {indice,kappa},n0= {m1.n0,0.0},Delta_n=nObj-n0;///init propriété bille

   // genere_bille(vol_obj,centre_boule, rayon_boule_pix,nObj-n0,dim3D.x);
   ///first abscissa, we add a width (x_width) and a repetition distance (delta_x)
   /* double x0=-2.5*pow(10,-6),delta_x=10*pow(10,-6),x_width=5*pow(10,-6);
    double x1=x0+x_width;
    double y0=-2.5*pow(10,-6),delta_y=0*pow(10,-6),y_width= 5*pow(10,-6);
    double y1=y0+y_width;
    double z_width=5*pow(10,-6);*/

    double x0=-15*pow(10,-6),delta_x=10*pow(10,-6),x_width=30*pow(10,-6);
    double x1=x0+x_width;
    double y0=-15*pow(10,-6),delta_y=0*pow(10,-6),y_width= 30*pow(10,-6);
    double y1=y0+y_width;
    double z_width=5*pow(10,-6);

    Point3D coordMin(x0,y0,-z_width/2,dim3D.x),coordMax(x1,y1,z_width/2,dim3D.x);

   genere_barre(vol_obj,coordMin,coordMax,nObj-n0, m1);

  /*  coordMin.set_coord3D(x1+delta_x,y1+delta_y,-z_width/2);
    coordMax.set_coord3D(x1+delta_x+x_width,y1+delta_y+y_width,z_width/2);
    genere_barre(vol_obj,coordMin,coordMax,nObj-n0, m1);*/

  /*  coordMin.set_coord3D(x0+2*(delta_x+x_width),-15*pow(10,-6),-z_width/2);
    coordMax.set_coord3D(x0+2*(delta_x+x_width)+x_width,15*pow(10,-6),z_width/2);
    genere_barre(vol_obj,coordMin,coordMax,nObj-n0, m1);*/



  // coordMin.set_coord3D(-25*pow(10,-6),-10*pow(10,-6),2.31*pow(10,-6));
      //  coordMax.set_coord3D(25*pow(10,-6),10*pow(10,-6),2.61*pow(10,-6));
  //  genere_barre(vol_bille,coordMin,coordMax,-0.02, m1);

   // coordMin.set_coord3D(-25*pow(10,-6),-10*pow(10,-6),3.31*pow(10,-6));
   // coordMax.set_coord3D(25*pow(10,-6),10*pow(10,-6),3.61*pow(10,-6));
  /*  genere_barre(vol_bille,coordMin,coordMax,-0.06, m1);*/


    SAV3D_Tiff(vol_obj,"Re",m1.chemin_result+"/obj_Re.tif",m1.Tp_Tomo);
   // SAV3D_Tiff((vol_bille),"im",m1.chemin_result+"/bille_im.tif",m1.Tp_Tomo);
    double energy = computeEnergy(vol_obj);
    std::cout << "Energy = " << energy/vol_obj.size() << std::endl;
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
    TF3Dcplx(tf3D.in,tf3D.out,fftshift3D(vol_obj),TF_obj,tf3D.p_forward_OUT,m1.Tp_Tomo);
    vector<complex<double>>().swap(vol_obj);
    TF_obj=fftshift3D(TF_obj);
    SAV3D_Tiff(TF_obj,"Re",m1.chemin_result+"/TF_obj_Re.tif", m1.Delta_f_tomo*pow(10,-6));

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
    //mon_OTF.bMultiCercleUNI(10);
   // mon_OTF.bFermat(nbAngle);
    //CoordSpec2=mon_OTF.bFleur(nbAxes);//retrieve tabular of specular beam from class OTF & create 3D OTF;

   // interp_lin3D(mon_OTF.Valeur);
    //SAV3D_Tiff(mon_OTF.Valeur,"Re",m1.chemin_result+"/OTF_simule_Re.tif",m1.Tp_Tomo);
    cout<<"Tp_tomo"<<m1.Tp_Tomo<<endl;
    write3D_Tiff(mon_OTF.Valeur,dim, "Re",m1.chemin_result+"/OTF_simule_Re.tif",m1.Tp_Tomo,"OTF partie reelle");



    ///calcul objet convolué
    for(int cpt=0; cpt<pow(m1.dim_final,3); cpt++){
     //SpectreObjConv[cpt]=mon_OTF.Valeur[cpt]*TF_bille[cpt];
     SpectreObjConv[cpt].real(mon_OTF.Valeur[cpt].real()*TF_obj[cpt].real());
     SpectreObjConv[cpt].imag(mon_OTF.Valeur[cpt].imag()*TF_obj[cpt].imag());
    }

    vector<complex<double>>().swap(mon_OTF.Valeur);
    TF3Dcplx_INV(tf3D.in,tf3D.out,fftshift3D(SpectreObjConv),obj_conv,tf3D.p_forward_OUT,m1.Delta_f_tomo);
    vector<complex<double>>().swap(SpectreObjConv);
    string description="partie réeelle convoluée,carre 30 µm, Δz=3.5, Δn=0.12\n";
    double energy_conv=computeEnergy(obj_conv);
    cout<<"Energie avant convolution="<<energy/nbPix3D<<endl;
    cout<<"Energie après convolution="<<energy_conv/obj_conv.size()<<endl;
    cout<<"ratio energie après/avant="<<energy_conv/energy<<endl;
    write3D_Tiff(fftshift3D(obj_conv),dim, "Re",m1.chemin_result+"/obj_conv_Re.tif",m1.Tp_Tomo,description.c_str());
    write3D_Tiff(fftshift3D(obj_conv),dim, "Im",m1.chemin_result+"/obj_conv_Im.tif",m1.Tp_Tomo,description.c_str());

    vector<complex<double>>().swap(obj_conv);

    vector<double> wrappedPhase(nbPix2D);
    vector<double> unwrappedPhase(nbPix2D);
    vector<double> amplitudeRytov(nbPix2D);
    vector<double> amplitudeBorn(nbPix2D);

    cout<<"extraction hologramme"<<endl;
    ///extract complex field from 3D spectrum

    for(int holo_numero=0; holo_numero<m1.nbHolo; holo_numero++){

        spec_H.x=CoordSpec_H[holo_numero].x;//the old code use spec, so we convert OTF.centre to spec
        spec_H.y=CoordSpec_H[holo_numero].y;

        int cpt2D=round(spec_H.y)*m1.dim_Uborn+round(spec_H.x);
        centres[spec_H.coordI().cpt2D()]=holo_numero;//used save centres in a image file, for quick visualisation

        ///Calcul des TF2D  des hologrammes à partir du spectre 3D de l'objet.
        calcHolo(spec_H,TF_obj,TF_holo,m1);//extract Ewald spheres + projection on a 2D plane
        // SAV2D_Tiff(TF_holo,"Im",dir_sav+"Tf_holo_Im.tif",m1.Tp_Uborn);
        //SAVCplx(TF_holo,"Re",m1.chemin_result+"TF_holo_H_Re_208x208x600.raw",t_float,"a+b");
        Var2D decal2centreI= {-spec_H.x+m1.dim_Uborn/2,spec_H.y+m1.dim_Uborn/2};///centering all the spectrum in "computer axes"
       // SAVCplx(TF_holo,"Re",m1.chemin_result+"TF_holo_I_Re_208x208x600.raw",t_float,"a+b");
        decal2DCplxGen(TF_holo,TF_holo_shift,decal2centreI);
        ///calcul des hologrammes après recalage (hologrammes centrés)//calculation of hologramms, after shifting
        TF2Dcplx_INV(TF_holo_shift,holo,tf2D,m1.Delta_f_Uborn);
        SAVCplx(fftshift2D(holo),"Im",m1.chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw",t_double,"a+b");
        SAVCplx(fftshift2D(holo),"Re",m1.chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw",t_double,"a+b");

        //calcPhase_mpi_pi_atan2(fftshift2D(holo),phase);
        for(int cpt=0;cpt<nbPix2D;cpt++){
            amplitudeRytov[cpt]=exp(sqrt(holo[cpt].imag()*holo[cpt].imag()+holo[cpt].real()*holo[cpt].real()));
            amplitudeBorn[cpt]=(sqrt(holo[cpt].imag()*holo[cpt].imag()+holo[cpt].real()*holo[cpt].real()));
            unwrappedPhase[cpt]=holo[cpt].imag();
            //phaseWrap[cpt]=wrapPhase(holo[cpt].imag());
        }
        //wrappedPhase=wrap_phase(unwrappedPhase);

        SAV2(fftshift2D(wrappedPhase),m1.chemin_result+"/wrapped_phase_rytov"+m1.dimImg+".raw",t_double,"a+b");

        SAV2(fftshift2D(unwrappedPhase),m1.chemin_result+"/phase_rytov"+m1.dimImg+".raw",t_double,"a+b");
       // SAV2(fftshift2D(amplitudeRytov),m1.chemin_result+"/amplitude_rytov"+m1.dimImg+".raw",t_double,"a+b");
        //SAV2(fftshift2D(amplitudeBorn),m1.chemin_result+"/amplitude_born"+m1.dimImg+".raw",t_double,"a+b");
    }

    //sauver les centres sous forme d'image, pour contrôle
    SAV_Tiff2D(centres,m1.chemin_result+"centres_simul.tif",m1.Tp_holo);
    //repasser les spéculaires en coordonnées informatique (car attendu par tomo_reconstruction)
    vector<double> tabPosSpec_I(m1.nbHolo*2);  ///nbHolo*2 coordonnées. stockage des speculaires pour exportation vers reconstruction
    for(int holo_numero=0;holo_numero<m1.nbHolo;holo_numero++){
        tabPosSpec_I[holo_numero]=CoordSpec_H[holo_numero].x+m1.dim_Uborn/2; //save tab_posSpec_X to a format readable by Tomo_reconstruction (computer )
        tabPosSpec_I[holo_numero+m1.nbHolo]=-CoordSpec_H[holo_numero].y+m1.dim_Uborn/2; //tab_posSpec_Y
    }
    cout<<"Résultats écrits dans "<<m1.chemin_result<<endl;
    ///donnée utiles à la reconstruction
    SAV2(tabPosSpec_I,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");
    vector<double> param{m1.NXMAX,m1.nbHolo,m1.R_EwaldPix,dimROI,m1.Tp_holo};//devenu inutile avec fichier de config
    SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");

    return 0;
}
