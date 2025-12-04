#include <iostream>
//#include <Magick++.h>

#include <vector>
#include "struct.h"
#include <complex>
#include "manip.h"
#include <string>
#include "projet.h"
#include "fonctions.h"
//#include "IO_fonctions.h"
//#include "FFTW_init.h"
#include "FFT_fonctions.h"
#include "deroulement_volkov4.h"
#include "deroulement_volkov3.h"
//#include "deroulement_volkov2.h"
//#include "deroulement_volkov.h"
#include "deroulement_herraez.h"
#include "Correction_aberration.h"
#include <chrono>
#include <opencv2/core/utility.hpp>
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"
#include <filesystem>
#include <cmath>
namespace fs = std::filesystem;
//#include <cv.h>
using namespace std;

int main(int argc,char *argv[]){
    //std::cout << "Nombre d'arguments : " << argc << "\n";
    string etat_polar;
    string gui_config_name;
    bool b_polar=false;//test if data are from polarisation tomography
   // cout<<"b_polar="<<b_polar<<endl;
    for (int i = 0; i < argc; ++i) {
        //std::cout << "Argument " << i << " : " << argv[i] << "\n";

        std::string arg = argv[i];
        if(arg == "--help"){
                cout<<"-polar : mode polar"<<endl;
                cout<<"-etat : indiquez un des 4 états :"<<endl;
                cout<<"/LC/polar0/"<<endl;
                cout<<"/LC/polar90/"<<endl;
                cout<<"/RC/polar0/"<<endl;
                cout<<"/RC/polar90/"<<endl;

        return 0;
        }
        else if(arg=="-polar"){
                cout<<"Mode polarisation activé, reconstruction des 4 états : RC0,RC90,LC0,LC90"<<endl;
                b_polar=true;
                gui_config_name="gui_tomo_polar.conf";
        }else if(arg == "-etat" && i + 1 < argc) {
        etat_polar = argv[++i];
        cout<<"etat polar="<<etat_polar <<endl;
        gui_config_name="gui_tomo_polar.conf";
        }
        else gui_config_name="gui_tomo.conf";
    }

    manip m1(gui_config_name); //créer un objet manip, avec fichier de config soit TDM classique, soit polar

    string chemin_result=m1.chemin_result,chemin_acquis;
    if(b_polar==true) {
            cout<<"b_polar="<<b_polar<<endl;
            chemin_acquis=m1.chemin_racine+"demosaic"+etat_polar;
            chemin_result=m1.chemin_racine+"demosaic"+etat_polar+"/UBorn/";
            cout<<"modification du chemin acquisitions : "<<chemin_acquis<<endl;
            cout<<"modification du chemin resultats : "<<chemin_result<<endl;
    }
    else chemin_acquis=m1.chemin_acquis;
    string Chemin_mask=chemin_acquis+"/Image_mask.pgm";
    if(b_polar==true)//if polar set up, pick up the root folder for aberration mask
    {
        Chemin_mask=m1.chemin_racine+"Image_mask.pgm";
    }

    cout<<"dans main"<<chemin_acquis<<endl;
    string str_sav_param_path=chemin_result+"/SAV_param_manip.txt";

    cout<<"camdimROI="<<m1.CamDimROI<<endl;
    Var2D const dimROI= {m1.CamDimROI,m1.CamDimROI}, coin= {0,0};
    Point2D const dimHolo(m1.CamDimROI,m1.CamDimROI,m1.CamDimROI);
    size_t const NbPixROI2d=dimROI.x*dimROI.y;
    //tableaux hologramme et réference en 1024x1024
    vector<double> static holo1(NbPixROI2d), intensite_ref(NbPixROI2d), ampli_ref(NbPixROI2d);
    char charAngle[4+1];
    ///-----------Init FFTW Holo---------------
    size_t nb_thread_fftw=m1.nbThreads;
    int fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(nb_thread_fftw);
    FFTW_init param_fftw2D_r2c_Holo(holo1,"r2c",nb_thread_fftw);/// /!\ init r2c->surcharge avec image2D réelle en entrée!
    FFTW_init param_fftw2D_c2r_Holo(holo1,nb_thread_fftw);

        ///--------------Init FFTW Holo-------------------------------------------------

        ///prepare fftw plan+tableaux-----------------
        fftw_plan p_forward_holo;
        fftw_complex *in_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//
        fftw_complex *out_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//
        p_forward_holo=fftw_plan_dft_2d(dimROI.x, dimROI.y, in_holo, out_holo, FFTW_FORWARD,FFTW_MEASURE);

    vector<double> static const masqueTukeyHolo(tukey2D(dimROI.x,dimROI.y,0.05));
    //masqueTukeyHolo=tukey2D(dimROI.x,dimROI.y,0.05);
    ///--------------------Init Référence Amplitude pour correction plan central-------------------------------------------

   // ampli_ref=initRef(chemin_acquis+"/Intensite_ref.pgm",coin,dimROI);
    //ampli_ref=initRef("/opt/Acquis/phytolites/acquis2/ACQUIS/Intensite_ref.pgm",coin,dimROI);

    ///------Init variable champ complexe après découpe hors axe-----------------------------
    const size_t  NbPixUBorn=4*m1.NXMAX*m1.NXMAX, NbAngle=m1.NbAngle;//dimensions
    cout<<"NbAngle="<<NbAngle<<endl<<"m1.NXMAX="<<m1.NXMAX<<endl;//nombre d'hologrammes
    size_t NbAngleOk=0;//Nbangle reellement utilisé
    Var2D const dim2DHA= {(size_t)2*m1.NXMAX,(size_t)2*m1.NXMAX},coinHA= {m1.circle_cx-m1.NXMAX,m1.circle_cy-m1.NXMAX},coinHA_shift={m1.fPortShift.x-m1.NXMAX,m1.fPortShift.y-m1.NXMAX};
    Var2D posSpec= {0,0},decal2DHA= {m1.NXMAX,m1.NXMAX};;
    cout<<"dimension decoupe dimDHA="<<dim2DHA.x<<endl;
    vector<complex<double>> TF_UBornTot(NbPixUBorn*NbAngle);///variable stockant les N champs complexes decoupés depuis la zone 1024 (pour utiliser wisdom en 1024)
    vector<double>  TF_champMod(NbPixUBorn), centre(NbPixUBorn);///centre pour controler balayage

auto start_decoupeHA = std::chrono::system_clock::now();///démarrage chrono Hors-axe
///Charger les acqusitions
//#pragma omp parallel for reduction(cpt)
FILE* test_existence;//tester l'existence des fichiers
unsigned short int cptAngle=0;
//#pragma omp parallel forTF2D_r2c
//cout<<"hello"<<endl;
///loop on the holograms
for(cptAngle=0; cptAngle<NbAngle; cptAngle++){
        if((cptAngle-100*(cptAngle/100))==0)    cout<<cptAngle<<endl;
        sprintf(charAngle,"%03i",cptAngle);
        string nomFichierHolo=chemin_acquis+"/i"+charAngle+".pgm";
        test_existence = fopen(nomFichierHolo.c_str(), "rb");
        // cout<<nomFichierHolo<<endl;
        if(test_existence!=NULL)
        {
            fclose(test_existence);
            //  cout<<"cahrger..."<<endl;
            //    auto start_holo2tf = chrono::steady_clock::now();
            charger_image2D_OCV(holo1,nomFichierHolo, coin, dimROI);
            //  auto end_holo2tf = chrono::steady_clock::now();
            //auto elapsed = end_holo2tf - start_holo2tf;
            //std::cout <<"Temps pour chargement= "<< chrono::duration_cast<chrono::microseconds>(elapsed).count()<< '\n';
            for(size_t cpt=0; cpt<NbPixROI2d; cpt++)
            {
                //holo1[cpt]=holo1[cpt]/ampli_ref[cpt];
                holo1[cpt]=holo1[cpt];
      }
     // SAV2(holo1,chemin_result+"holo1_divise_1024x1024.raw",t_float,"a+b");

    //auto start_holo2tf = chrono::steady_clock::now();
    //off-axis extraction
    holo2TF_UBorn( holo1, TF_UBornTot,dimROI, dim2DHA, coinHA, NbAngleOk,masqueTukeyHolo, in_holo,out_holo, p_forward_holo);
       //     auto end_holo2tf = chrono::steady_clock::now();
    //auto elapsed = end_holo2tf - start_holo2tf;
    //std::cout <<"Temps pour holo2tf= "<< chrono::duration_cast<chrono::microseconds>(elapsed).count()<< '\n';
     ///nouvelles versions
    //holo2TF_UBorn2(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA,NbAngleOk, masqueTukeyHolo,param_fftw2D_r2c_Holo);///calcul TF holo+ découpe dans TF symétrisée, repère humain
    //holo2TF_UBorn2_shift(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA_shift,NbAngleOk, masqueTukeyHolo,param_fftw2D_r2c_Holo);///calcul TF hologrammes +  découper de dimROI à 2NXMAX dans TF symétrisée, repère informatique
    //holo2TF_UBorn2_shift_r2c(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA_shift,NbAngleOk, masqueTukeyHolo,param_fftw2D_r2c_Holo);///calcul TF hologrammes +  découper de dimROI à 2NXMAX dans TF NON symétrique, repère info
      NbAngleOk++;
    }
    else cout<<"fichier "<<cptAngle<<" inexistant\n";
}

auto end_decoupeHA = std::chrono::system_clock::now();
auto elapsed = end_decoupeHA - start_decoupeHA;
std::cout <<"Temps pour FFT holo+découpe Spectre= "<< elapsed.count()/(pow(10,9)) << '\n';
 //SAVCplx(TF_UBornTot,"Re",chemin_result+"/TF_Uborn_Tot_GPU_250x250x599x64_orig.raw",t_float,"wb");
m1.dimImg=to_string(dim2DHA.x)+"x"+to_string(dim2DHA.y)+"x"+to_string(NbAngleOk);
cout<<"m1.dimImg="<<m1.dimImg<<endl;
deleteCplxField(chemin_result, m1.dimImg);
///initialiser les variables de champ, phase, amplitude etc.

vector<complex<double>> TF_UBorn(NbPixUBorn),  UBorn(NbPixUBorn);
vector<double> phase_2Pi_vec(NbPixUBorn),  UnwrappedPhase(NbPixUBorn),PhaseFinal(NbPixUBorn);
double *UnwrappedPhase_herraez=new double[NbPixUBorn];

float alpha=0.1;//coeff pour le masque de tuckey
vector<double> masqueTukeyHA(tukey2D(dim2DHA.x,dim2DHA.y,alpha));
FFTW_init param_fftw2D_c2r_HA(TF_UBorn,m1.nbThreads);//OUTPLACE

//FFTW_init param_fftw2D_c2r_HA(dim2DHA,1,m1.nbThreads);//INPLACE
//FFTW_init param_fftw2D_r2c_HA(phase_2Pi_vec,"r2c",m1.nbThreads);

///---------------------use symetry for phase unwrapping--------------------
vector<double> phase_2Pi_vec_double(4*NbPixUBorn);//variable created only to init param_fftw2D_r2c_HA_double
FFTW_init param_fftw2D_c2r_HA_double(phase_2Pi_vec_double,m1.nbThreads);//OUTPLACE
///variable pour correction aberration
cout<<"chemin_mask==============="<<Chemin_mask<<endl;
Mat src=Mat(1, ampli_ref.size(), CV_64F, ampli_ref.data()), mask_aber=init_mask_aber(Chemin_mask,chemin_acquis,dim2DHA);


if(fs::exists(Chemin_mask)){
        string info="un masque d'aberration a été utilisé----------------------------------------------------------------";
sav_param2D(info,str_sav_param_path);
}
//imshow("Image",mask_aber);
//waitKey(0);
size_t NbPtOk=countM(mask_aber),  degre_poly=4, nbCoef = sizePoly2D(degre_poly);//Nb coef poly
Mat polynomeUs_to_fit(Size(nbCoef,NbPtOk), CV_64F);///(undersampled) Polynome to fit= function to fit (We use a polynome). we have to generate a table containing polynome_to_fit=[1,x,x^2,xy,y^2] for each coordinate (x,y)
Mat polynome_to_fit(Size(nbCoef,dim2DHA.x*dim2DHA.y), CV_64F);

string str_degre_poly="degré poly aberration="+to_string(degre_poly);
sav_param2D(str_degre_poly,str_sav_param_path);

initCorrAber(Chemin_mask, mask_aber, degre_poly,dim2DHA,polynome_to_fit,polynomeUs_to_fit);
cout<<"COUCOU"<<endl;

cout<<"\n#########################Calcul champs cplx 2D Uborn/Rytov + eventuelle Correction aberrations#############################"<<endl;
cout<<"NbPixUborn="<<NbPixUBorn<<endl;
vector<complex<double>> UBornFinal(NbPixUBorn), UBornFinalDecal(NbPixUBorn), TF_UBorn_norm(NbPixUBorn);
cout<<"COUCOU3"<<endl;
vector<double> UBornAmpFinal(NbPixUBorn),  UBornAmp(NbPixUBorn);
cout<<"COUCOU4"<<endl;
vector<double> tabPosSpec(NbAngleOk*2);  ///stockage des spéculaires pour exportation vers reconstruction
vector<vecteur>  kvect_shift(init_kvect_shift(dim2DHA));///init opérateur differentiation kvect
vector<double> kvect_mod2Shift(init_kvect_mod2Shift(kvect_shift));
///kvect_shift for symetrized image
vector<vecteur>  double_kvect_shift(4*dim2DHA.x);
double_kvect_shift=init_kvect_shift({2*dim2DHA.x,2*dim2DHA.y});


cout<<"COUCOU4"<<endl;
auto start_part2= std::chrono::system_clock::now();
//#pragma omp parallel for

cout<<"avant boucle"<<endl;
for(size_t cpt_angle=0; cpt_angle<NbAngleOk; cpt_angle++){ //boucle sur tous les angles : correction aberrations
  TF_UBorn.assign(TF_UBornTot.begin()+NbPixUBorn*cpt_angle, TF_UBornTot.begin()+NbPixUBorn*(cpt_angle+1));  ///Récupérer la TF2D dans la pile de spectre2D//get back 2D spectrum in the stack (dimxdimxNumber_holograms)
 // SAVCplx(TF_UBorn,"Re",chemin_result+"/TF_Uborn_iterateur_Re_220x220x599x32.raw",t_float,"a+b");
 // SAVCplx(TF_UBorn,"Im",chemin_result+"/TF_Uborn_iterateur_Im_220x220x599x32.raw",t_float,"a+b");
  //Recherche de la valeur maximum du module dans ref non centré-----------------------------------------
  size_t cpt_max=coordSpec(TF_UBorn, TF_champMod,decal2DHA);
  double  max_part_reel = TF_UBorn[cpt_max].real(),///sauvegarde de la valeur cplx du spéculaire/save complex value of the specular beam
  max_part_imag = TF_UBorn[cpt_max].imag(),
  max_module = sqrt(TF_UBorn[cpt_max].imag()*TF_UBorn[cpt_max].imag()+TF_UBorn[cpt_max].real()*TF_UBorn[cpt_max].real());
  const int kxmi=cpt_max%(2*m1.NXMAX), kymi=cpt_max/(2*m1.NXMAX);
  posSpec= {kxmi,kymi}; ///coord informatique speculaire

  ///calculate phi and theta (angles of illumination)
  Var2D posSpecH={kxmi-m1.NXMAX,kymi-m1.NXMAX};
  float kiz=sqrt(pow(m1.rayon,2)-pow(posSpecH.x,2)-pow(posSpecH.y,2));
  float theta=acos(kiz/m1.rayon)*180/M_PI;
  // cout<<"theta="<<theta<<endl;
  float phi=atan2(static_cast<double>(posSpecH.y),static_cast<double>(posSpecH.x));
  // cout<<"(kixs,kiy,kiz)=("<<posSpecH.x<<","<<posSpecH.y<<","<<kiz<<")"<<endl;
 // cout<<"num angle="<<cpt_angle<<", phi="<<phi*180/M_PI<<endl;
  tabPosSpec[cpt_angle]=(double)posSpec.x;
  tabPosSpec[cpt_angle+NbAngleOk]=(double)posSpec.y;
  centre[kxmi*2*m1.NXMAX+kymi]=cpt_angle;
///calculer la phase, dérouler, corriger les aberrations
  if(m1.b_CorrAber==true){
    calc_Uborn2(TF_UBorn,UBorn,dim2DHA,posSpec,param_fftw2D_c2r_HA);
  //  SAVCplx(UBorn,"Re",chemin_result+"/ampli_UBorn_debut_extract.raw",t_float,"a+b");
    ///----------Calcul phase + deroulement--------------------------------
    calcPhase_mpi_pi_atan2(UBorn,phase_2Pi_vec); ///fonction atan2
    //SAV2(phase_2Pi_vec,chemin_result+"/phasePI_atan2.raw",t_float,"a+b");
    //phase2pi(UBorn, dim2DHA,phase2Pi);//asin
    if(m1.b_Deroul==true){
      if(m1.b_volkov==0){
        phaseUnwrapping_Mat(dim2DHA, phase_2Pi_vec, UnwrappedPhase_herraez);
        for(size_t cpt=0;cpt<NbPixUBorn;cpt++)
          UnwrappedPhase[cpt]=UnwrappedPhase_herraez[cpt];///plutôt passer pointeur ?
        }
        else{
           // auto start_volkov = std::chrono::system_clock::now();
           //deroul_volkov3(phase_2Pi_vec,UnwrappedPhase, kvect_shift, param_fftw2D_c2r_HA);
           ///symetrized function to avoid FFT artifact (simple mirror symetry before unwrapping)
           deroul_volkov4_total_sym_paire(phase_2Pi_vec,UnwrappedPhase, double_kvect_shift, param_fftw2D_c2r_HA_double);
             // auto end_volkov = std::chrono::system_clock::now();
          //  auto elapsed_volkov = end_volkov - start_volkov;
           // std::cout <<"Temps pour deroul volkov2= "<< elapsed_volkov.count()/(pow(10,9)) << '\n';
        }
    }
    else UnwrappedPhase=phase_2Pi_vec;
     // SAV2(UnwrappedPhase,chemin_result+"/phase_deroul_volkov_avant_corr_aber_vvvvolkov3.raw",t_float,"a+b");
      //-------------Correction aberration phase-------------------------------
      src=Mat(dim2DHA.x,dim2DHA.y,CV_64F, UnwrappedPhase.data());
      // auto start_calcAber = std::chrono::system_clock::now();
      Mat Phase_corr(aberCorr2(src, mask_aber,polynomeUs_to_fit,polynome_to_fit));
        //  auto end_calcAber = std::chrono::system_clock::now();
        //auto elapsed = end_calcAber - start_calcAber;
        //std::cout <<"Temps pour FFT holo+découpe Spectre= "<< elapsed.count()/(pow(10,9)) << '\n';

            // SAV2((double*)Phase_corr.data,Phase_corr.rows*Phase_corr.cols,"/home/mat/tmp/phaseCorr_main.raw",t_float,"a+b");
            ///---------------Correction amplitude----------------------------------------
            for(size_t cpt=0; cpt<(NbPixUBorn); cpt++)
                UBornAmp[cpt]=abs(UBorn[cpt]);
            //SAV2(UBornAmp,chemin_result+"/UBornAmp.raw",t_float,"a+b");
            Mat srcAmp=Mat(dim2DHA.x, dim2DHA.y, CV_64F, UBornAmp.data());///image source
            Mat UBornAmp_corr(ampliCorr2(srcAmp, polynomeUs_to_fit, polynome_to_fit, mask_aber));///résultat amplitude
    ///Fin Correction amplitude----------------------------------------*/
int flag=0;
    for(size_t y=0; y<dim2DHA.y; y++){ // reconstruire l'onde complexe/Recalculate the complex field
      for(size_t x=0; x<dim2DHA.x; x++){
        size_t cpt=x+y*dim2DHA.x;
        PhaseFinal[cpt]=Phase_corr.at<double>(y,x) ;//copie opencV->Tableau
        UBornAmpFinal[cpt]=UBornAmp_corr.at<double>(y,x);
        if(UBornAmpFinal[cpt]>10) UBornAmpFinal[cpt]=1;//Amplitude should always be around 1 after correction
        if(UBornAmpFinal[cpt]<-10)  UBornAmpFinal[cpt]=1;
       // if(UBornAmpFinal[cpt]<-1 && flag==0) { cout<<"Holo numero"<<cpt_angle<<endl;flag=1;}
        if(m1.b_Born==true){ //UBORN=U_tot-U_inc=u_tot_norm-1
         // UBornFinal[cpt].real( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])- 1 )*cos(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);///correction amplitude
         // UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]) )*sin(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);
          UBornFinal[cpt].real( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*cos(PhaseFinal[cpt]) -1);//*masqueTukeyHolo[cpt]);///correction amplitude
         // UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*sin(PhaseFinal[cpt]) -1);//*masqueTukeyHolo[cpt]);
          UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*sin(PhaseFinal[cpt]));
        }
        else{ //RYTOV URytov = log a_t/a_i (=log a_t après correction AmpliCorr)
          UBornFinal[cpt].real(log(sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])));///racine(U^2) pour éliminer les éventuelles amplitudes négatives liées au bruit
          UBornFinal[cpt].imag(PhaseFinal[cpt]);
        }
      }
    }
    //cout<<"hello01"<<endl;
    //SAV2(PhaseFinal,chemin_result+"/phase_finale.raw",t_float,"a+b");
    //SAV2(UBornAmpFinal,chemin_result+"/UBornAmpFinal.raw",t_float,"a+b");
    ///Recalculer la TF décalée pour le programme principal.
    Var2D recal= {kxmi,kymi};
    decal2DCplxGen2(UBornFinal,UBornFinalDecal, decal2DHA);
    //TF2Dcplx_vec(in_HA,out_HA,UBornFinalDecal,TF_UBorn_norm,p_forward_HA);
    TF2Dcplx(UBornFinalDecal,TF_UBorn_norm,param_fftw2D_c2r_HA);
    SAVCplx(UBornFinal,"Re", chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw", t_double, "a+b");
    SAVCplx(UBornFinal,"Im", chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw", t_double, "a+b");
  }
  else{ ///sauvegarde onde avec aberration
   if(m1.b_Born==0)
    cout<<"Attention : Rytov nécessite un déroulement"<<endl;
    //cout<<"normalise spec"<<endl;
    for(size_t cpt=0; cpt<(4*m1.NXMAX*m1.NXMAX); cpt++){ //correction phase à l'ordre zéro et normalisatoin amplitude par ampli_spec=ampli_inc*ampli_ref
      TF_UBorn_norm[cpt].real((TF_UBorn[cpt].real()*max_part_reel+TF_UBorn[cpt].imag()*max_part_imag)/max_module);
      TF_UBorn_norm[cpt].imag((TF_UBorn[cpt].imag()*max_part_reel-TF_UBorn[cpt].real()*max_part_imag)/max_module);
    }
    //SAVCplx(TF_UBorn,"Re", chemin_result+"/TF_Uborn_Re.raw", t_double, "a+b");
    calc_Uborn2(TF_UBorn_norm,UBorn,dim2DHA,posSpec,param_fftw2D_c2r_HA);

    SAVCplx(UBorn,"Re", chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw", t_double, "a+b");
    SAVCplx(UBorn,"Im", chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw", t_double, "a+b");
  }
}//fin de boucle for sur tous les angles
//SAV2(UBornAmpFinal,chemin_result+"/UBornAmpFinal.raw",t_float,"a+b");
auto end_part2= std::chrono::system_clock::now();
auto elapsed_part2 = end_part2 - start_part2;
std::cout <<"Temps pour part2= "<< elapsed_part2.count()/(pow(10,9)) << '\n';
delete[] UnwrappedPhase_herraez;
 //cout<<"hello02"<<endl;
///------------SAUVER LES PARAMETRES UTILES À RECONSTRUCTION ou au contrôle//Save useful paramaters for tomo_reconstruction
cout<<"NXMAX sauve="<<m1.NXMAX<<endl;
vector<double> param{m1.NXMAX,NbAngleOk,m1.rayon,dimROI.x,m1.tailleTheoPixelHolo};
SAV2(param,chemin_result+"/parametres.raw", t_double, "wb");
/*for(int cpt=0;cpt<NbAngle;cpt++)
{
    cout<<"tabPosSpec.y="<<tabPosSpec[cpt]<<","<<"tabPosSpec.y="<<tabPosSpec[cpt+NbAngle]<<endl;
}*/
SAV2(tabPosSpec,chemin_result+"/tab_posSpec.raw",t_double,"wb");///sav speculaire en coord informatique [0->2NXMAX-1]
SAV_Tiff2D(centre,chemin_result+"/centres.tif",m1.NA/m1.NXMAX); //exportation des spéculaires en "unité NA"
cout<<"Fin prétraitement"<<endl;
cout<<"Results in "<<chemin_result<<endl;

return 0;
}
