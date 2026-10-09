//######Generate a 3D complex object (bead or box), calculate its spectrum and extract hologram thanks to Ewald sphere.######//
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>  //setprecision
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include <cstdlib>
#include <cstdio>
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
#include "bruit.h"
#include "deroulement_volkov4.h"
#include "Correction_aberration.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace std;
using namespace cv;

/** @function main */
int main( int argc, char** argv )
{     int dimROI=512;
    string home=getenv("HOME");
///--------------- Chargement ou création de l'objet (bille, spectre etc.)------------
    manip m1(dimROI);//Initialiser la manip avec la taille du champ holographique (acquisition) en pixel. Attention, cela effacera la valeur du fichier de config !
    if(m1.dim_final<2*m1.dim_Uborn)
    {std::cerr << "Erreur : dim_final="<<m1.dim_final<<" < 2dimUborn="<<2*m1.dim_Uborn<< "\n";
        return EXIT_FAILURE; //
    }
    cout<<"chemin_acquis="<<m1.chemin_acquis<<endl;
    const int nbAngle=m1.nbHolo;
    cout<<"nbAngle="<<nbAngle;

    Point3D dim3D(m1.dim_final,m1.dim_final,m1.dim_final);
    Var3D dim={m1.dim_final,m1.dim_final,m1.dim_final};
    Var3D dimStack={m1.dim_Uborn,m1.dim_Uborn,m1.nbHolo};//dimension of complex field stack.
    Point2D dim2D((double)m1.dim_Uborn,(double)m1.dim_Uborn,round(m1.dim_Uborn));
    unsigned int nbPix3D=pow(m1.dim_final,3),nbPix2D=pow(m1.dim_Uborn,2);
    ///variables 3D : objet+spectre
    vector<complex<double>> vol_obj(nbPix3D),  TF_obj(nbPix3D);
    vector<complex<double>> obj_conv(nbPix3D), SpectreObjConv(nbPix3D);
    vector<double> unwrappedPhase(nbPix2D);
    ///variables 2D : champ complexe+spectre
    vector<complex<double>> TF_holo(nbPix2D), TF_holo_shift(nbPix2D);
    vector<complex<double>> holo(nbPix2D),  holo_centre(nbPix2D);
    ///exportation des centres, pour contrôle
    vector<double> centres(nbPix2D);

    ///Génération de l'objet (bille polystyrène, n=1.5983, absorption=?)
    double Rboule_metrique=1.95*pow(10,-6);///rayon bille en m
    int rayon_boule_pix=round(Rboule_metrique/m1.Tp_Tomo);///rayon bille en pixel
    Point3D centre_boule(dim3D.x/2,dim3D.x/2,dim3D.x/2,dim3D.x);//bille centrée dans l'image
    double indice=1.41,kappa=0.000;//indice + coef d'extinction
    complex<double> nObj= {indice,kappa},nM= {m1.nM,0.0},Delta_n=nObj-nM;///init propriété bille

    genere_bille(vol_obj,centre_boule, rayon_boule_pix,nObj-nM,dim3D.x);
   ///first abscissa, we add a width (x_width) and a repetition distance (delta_x)
   /* double x0=-2.5*pow(10,-6),delta_x=10*pow(10,-6),x_width=5*pow(10,-6);
    double x1=x0+x_width;
    double y0=-2.5*pow(10,-6),delta_y=0*pow(10,-6),y_width= 5*pow(10,-6);
    double y1=y0+y_width;
    double z_width=5*pow(10,-6);*/

    double x0=-2.5*pow(10,-6),delta_x=10*pow(10,-6),x_width=5*pow(10,-6);
    double x1=x0+x_width;
    double y0=-2.5*pow(10,-6),delta_y=0*pow(10,-6),y_width= 5*pow(10,-6);
    double y1=y0+y_width;
    double z_width=5*pow(10,-6);

    Point3D coordMin(x0,y0,-z_width/2,dim3D.x),coordMax(x1,y1,z_width/2,dim3D.x);

  // genere_barre(vol_obj,coordMin,coordMax,nObj-n0, m1);

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

    ///Calcul spectre objetabfeyn
    FFTW_init tf3D(dim3D),tf2D(dim2D); ///init fftw pour spectre 3D et 2D
    TF3Dcplx(tf3D.in,tf3D.out,fftshift3D(vol_obj),TF_obj,tf3D.p_forward_OUT,m1.Tp_Tomo);
    vector<complex<double>>().swap(vol_obj);
    TF_obj=fftshift3D(TF_obj);
    SAV3D_Tiff(TF_obj,"Re",m1.chemin_result+"/TF_obj_Re.tif", m1.Delta_f_tomo*pow(10,-6));
    SAV3D_Tiff(TF_obj,"Im",m1.chemin_result+"/TF_obj_Im.tif", m1.Delta_f_tomo*pow(10,-6));
    Point2D spec_H(0.0,0.0,m1.dim_Uborn);//coordonnées spec2D dans une image de dimension dimUBorn

    ///Generate OTF and 2D center (specular beam).
    //vector<Var2D> CoordSpec(m1.nbHolo);//table of specular coordinates
    vector<Point2D> CoordSpec_H(nbAngle, spec_H);
   // short unsigned int const nbAxes=4;//nombre de branche de la fleur

    cout<<"------------------------Calcul OTF et objet convolué"<<endl;
    OTF mon_OTF(m1);//init OTF
    //mon_OTF.bFleur(CoordSpec_H, nbAxes);//retrieve tabular of specular beam from class OTF & create 3D OTF;
    mon_OTF.scan_uniform3D(CoordSpec_H,0.95);
    //mon_OTF.bSpiral();
   /*for(int cpt=0;cpt<m1.nbHolo;cpt++)
    {  cout<<"---------------------"<<endl;
        cout<<"CoordSpec_H.x["<<cpt<<"]="<<CoordSpec_H[cpt].x<<endl;
        cout<<"CoordSpec_H.y["<<cpt<<"]="<<CoordSpec_H[cpt].y<<endl;
    }*/
    if(m1.b_no_absorption==true){
            cout<<"NO ABSORPTION : symetrisation of OTF"<<endl;
            mon_OTF.symetrize_xoy();
    }
    //mon_OTF.bMultiCercleUNI(10);
   // mon_OTF.bFermat(nbAngle);
    //CoordSpec2=mon_OTF.bFleur(nbAxes);//retrieve tabular of specular beam from class OTF & create 3D OTF;

   // interp_lin3D(mon_OTF.Valeur);
    //SAV3D_Tiff(mon_OTF.Valeur,"Re",m1.chemin_result+"/OTF_simule_Re.tif",m1.Tp_Tomo);
    cout<<"Tp_tomo"<<m1.Tp_Tomo<<endl;
    write3D_Tiff(mon_OTF.Valeur,dim, "Re",m1.chemin_result+"/OTF_simule_Re.tif",m1.Tp_Tomo,"OTF partie reelle");
    write3D_Tiff(mon_OTF.Valeur,dim, "Im",m1.chemin_result+"/OTF_simule_Im.tif",m1.Tp_Tomo,"OTF partie imag");



    ///calcul objet convolué
    for(int cpt=0; cpt<pow(m1.dim_final,3); cpt++){
     //SpectreObjConv[cpt]=mon_OTF.Valeur[cpt]*TF_bille[cpt];
     SpectreObjConv[cpt].real(mon_OTF.Valeur[cpt].real()*TF_obj[cpt].real());
     SpectreObjConv[cpt].imag(mon_OTF.Valeur[cpt].imag()*TF_obj[cpt].imag());
    }
     write3D_Tiff(SpectreObjConv,dim, "Re",m1.chemin_result+"/spectre_conv_simule_Re.tif",m1.Tp_Tomo,"spectre convolué partie reelle");
     write3D_Tiff(SpectreObjConv,dim, "Im",m1.chemin_result+"/spectre_conv_simule_Im.tif",m1.Tp_Tomo,"spectre convolué partie Imag");
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
   // vector<double> phase(nbPix2D);
    vector<double> wrappedPhase(nbPix2D);
    vector<double> amplitude(nbPix2D);


    cout<<"extraction hologramme"<<endl;

    ///---------------------------------------correction aberration : INIT---------------
        ///initialisation
        string Chemin_mask=m1.chemin_acquis+"/Image_mask.pgm";
        cout<<"Chemin_mask"<<Chemin_mask<<endl;
        vector<double>  ampli_ref(nbPix2D);
        Mat src=Mat(1, ampli_ref.size(), CV_64F, ampli_ref.data());
        Var2D dim2DUborn{m1.dim_Uborn,m1.dim_Uborn};
        Mat  mask_aber=init_mask_aber(Chemin_mask,m1.chemin_acquis,dim2DUborn);
        if(fs::exists(Chemin_mask)){
        string info="An aberration Mask has been used--";
        cout<<"chemin_mask==============="<<Chemin_mask<<endl;
        }
        size_t NbPtOk=countM(mask_aber),  degre_poly=4, nbCoef = sizePoly2D(degre_poly);//Nb coef poly
        Mat polynomeUs_to_fit(Size(nbCoef,NbPtOk), CV_64F);///(undersampled) Polynome to fit= function to fit (We use a polynome). we have to generate a table containing polynome_to_fit=[1,x,x^2,xy,y^2] for each coordinate (x,y)
        Mat polynome_to_fit(Size(nbCoef,m1.dim_Uborn*m1.dim_Uborn), CV_64F);

        string str_degre_poly="degré poly aberration="+to_string(degre_poly);
        initCorrAber(Chemin_mask, mask_aber, degre_poly,dim2DUborn,polynome_to_fit,polynomeUs_to_fit);
///loop on hologramms
   for(size_t holo_numero=0; holo_numero<m1.nbHolo; holo_numero++){

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
        holo=fftshift2D(holo);
        SAVCplx(holo,"Im",m1.chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw",t_double,"a+b");
        SAVCplx(holo,"Re",m1.chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw",t_double,"a+b");
       // write3D_Tiff(amplitude,dimStack,m1.chemin_result+"/amplitude.tif",m1.Tp_Tomo,description.c_str());
        ///test du déroulement et eventuellement, d'un décalage de phase/amplitude
     //  double maxAmplitude = *max_element(amplitude.begin(), amplitude.end());
       // calcPhase_mpi_pi_atan2(fftshift2D(holo),maxAmplitude,wrappedPhase);///wrapped phase
       double deltaPhi=0.5;
      /* for(size_t cpt=0;cpt<holo.size();++cpt)
       {    //Re'=Re*cos(dPhi)-Im*sin(dPhi);
           holo[cpt].real(holo[cpt].real()*cos(deltaPhi)-holo[cpt].imag()*sin(deltaPhi));
           //Im'=Im*cos(dPhi)+Re*sin(dPhi)
           holo[cpt].imag(holo[cpt].imag()*cos(deltaPhi)+holo[cpt].real()*sin(deltaPhi));
       }*/
     /*  double epsilon=0.2;
         for(size_t cpt=0;cpt<holo.size();++cpt)
       {    //Re'=Re*cos(dPhi)-Im*sin(dPhi);
           holo[cpt].real(holo[cpt].real()+epsilon*cos(deltaPhi));
           //Im'=Im*cos(dPhi)+Re*sin(dPhi)
           holo[cpt].imag(holo[cpt].imag()+epsilon*sin(deltaPhi));
       }*/


        calcPhase_mpi_pi_atan2(holo,wrappedPhase);///wrapped phase


        ///---phase unwrapping-----
        ///kvect_shift for symetrized image
        //vector<vecteur>  double_kvect_shift(4*m1.dim_Uborn);


        vector<vecteur> double_kvect_shift=init_kvect_shift({2*m1.dim_Uborn,2*m1.dim_Uborn});
       // cout<<"double_kvect_shift size="<<double_kvect_shift.size()<<endl;
        ///---------------------phase unwrapping--------------------

        vector<double> phase_2Pi_vec_double(4*nbPix2D);//variable created only to init param_fftw2D_r2c_HA_double
        FFTW_init param_fftw2D_c2r_HA_double(phase_2Pi_vec_double,m1.nbThreads);//OUTPLACE
       // cout<<"fft size (dimension)="<<sqrt(param_fftw2D_c2r_HA_double.getFFTSize())<<endl;
        //cout<<"dimension image"<<sqrt(4*nbPix2D)<<endl;
        deroul_volkov4_total_sym_paire(wrappedPhase,unwrappedPhase, double_kvect_shift, param_fftw2D_c2r_HA_double);
        ///-------------------
        SAV2(wrappedPhase,m1.chemin_result+"/wrappedPhase"+m1.dimImg+".raw",t_double,"a+b");
        SAV2(unwrappedPhase,m1.chemin_result+"/unwrappedPhase"+m1.dimImg+".raw",t_double,"a+b");
        SAV2(fftshift2D(amplitude),m1.chemin_result+"/amplitude"+m1.dimImg+".raw",t_double,"a+b");


        ///correction en elle même--------------
        src=Mat(m1.dim_Uborn,m1.dim_Uborn,CV_64F, unwrappedPhase.data());
        Mat Phase_corr(aberCorr2(src, mask_aber,polynomeUs_to_fit,polynome_to_fit));
        ///-------------------------------------------------------------------------------------------
        // write3D_Tiff(amplitude,dim, "Re",m1.chemin_result+"/OTF_simule_Re.tif",m1.Tp_Tomo,"OTF partie reelle");
        bruit bruitPhase,bruitAmplitude;
        bruitPhase.randomDoubleUnit();
        bruitAmplitude.valeur=1;

        for(int cpt=0;cpt<nbPix2D;cpt++){
        amplitude[cpt]=abs(holo[cpt]);
        }
        vector<double> PhaseFinal(nbPix2D);
        for(size_t y=0; y<m1.dim_Uborn; y++){ // reconstruire l'onde complexe/Recalculate the complex field
            for(size_t x=0; x<m1.dim_Uborn; x++){
            size_t cpt=x+y*m1.dim_Uborn;
            PhaseFinal[cpt]=Phase_corr.at<double>(y,x) ;//copie opencV->Tableau
            //UBornAmpFinal[cpt]=UBornAmp_corr.at<double>(y,x);
            }
        }
         SAV2(PhaseFinal,m1.chemin_result+"/PhaseFinal"+m1.dimImg+".raw",t_double,"a+b");
        ///calculate Rytov complex field and save it for 3D reconstruction
        vector<complex<double>> chpRytov(nbPix2D);
        vector<complex<double>> chpBorn(nbPix2D);
      /*  for(int cpt=0;cpt<nbPix2D;cpt++){
            chpRytov[cpt].real(log(amplitude[cpt]));
           // chpRytov[cpt].imag(unwrappedPhase[cpt]);
            chpRytov[cpt].imag(PhaseFinal[cpt]);
        }*/

      //  SAVCplx(chpRytov,"Im",m1.chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw",t_double,"a+b");
     //   SAVCplx(chpRytov,"Re",m1.chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw",t_double,"a+b");
        //wrappedPhase=wrap_phase(unwrappedPhase);

       // SAV2(fftshift2D(wrappedPhase),m1.chemin_result+"/wrapped_phase_rytov"+m1.dimImg+".raw",t_double,"a+b");
       // SAV2(fftshift2D(unwrappedPhase),m1.chemin_result+"/phase_rytov"+m1.dimImg+".raw",t_double,"a+b");
        //SAV2(fftshift2D(amplitudeRytov),m1.chemin_result+"/amplitude_rytov"+m1.dimImg+".raw",t_double,"a+b");
        //SAV2(fftshift2D(amplitudeBorn),m1.chemin_result+"/amplitude_born"+m1.dimImg+".raw",t_double,"a+b");
    }

    //sauver les centres sous forme d'image, pour contrôle
    SAV_Tiff2D(centres,m1.chemin_result+"centres_simul.tif",m1.Tp_holo);
    //repasser les spéculaires en coordonnées informatique (car attendu par tomo_reconstruction)
    vector<double> tabPosSpec_I(m1.nbHolo*2);  ///nbHolo*2 coordonnées. stockage des speculaires pour exportation vers reconstruction
    bruit randomizeSpecular;
    randomizeSpecular.valeur=0;
    for(int holo_numero=0;holo_numero<m1.nbHolo;holo_numero++){
        //randomizeSpecular.randomDoubleUnit();
       /// cout<<round(randomizeSpecular.valeur*5)<<endl;
        tabPosSpec_I[holo_numero]=CoordSpec_H[holo_numero].x+m1.dim_Uborn/2+round(randomizeSpecular.valeur*5); //save tab_posSpec_X to a format readable by Tomo_reconstruction (computer )
       /// randomizeSpecular.randomDoubleUnit();
        ///cout<<round(randomizeSpecular.valeur*5)<<endl;
        tabPosSpec_I[holo_numero+m1.nbHolo]=-CoordSpec_H[holo_numero].y+m1.dim_Uborn/2+round(randomizeSpecular.valeur*5); //tab_posSpec_Y
    }
    cout<<"Résultats écrits dans "<<m1.chemin_result<<endl;
    ///donnée utiles à la reconstruction
    SAV2(tabPosSpec_I,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");
    vector<double> param{m1.NXMAX,m1.nbHolo,m1.R_EwaldPix,dimROI,m1.Tp_holo};
    SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");

    return 0;
}
