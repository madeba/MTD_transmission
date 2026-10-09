#include <iostream>
#include "zernike.h"
#include <vector>
#include "struct.h"
#include <complex>
#include "manip.h"
#include "math_functions.h"
#include <string>
#include "projet.h"
#include "fonctions.h"
//#include "IO_fonctions.h"
//#include "FFTW_init.h"
#include "FFT_fonctions.h"
#include "deroulement_volkov4.h"
#include "deroulement_herraez.h"
#include "Correction_aberration2.h"
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

    string etat_polar;
    string gui_config_name;
    bool b_polar=false;//test if data are from polarisation tomography
    // cout<<"b_polar="<<b_polar<<endl;
    for (int i = 0; i < argc; ++i)
    {
        //std::cout << "Argument " << i << " : " << argv[i] << "\n";

        std::string arg = argv[i];
        if(arg == "--help")
        {
            cout<<"-polar : polar mode"<<endl;
            cout<<"-etat : indiquez un des 4 états :"<<endl;
            cout<<"/LC/polar0/"<<endl;
            cout<<"/LC/polar90/"<<endl;
            cout<<"/RC/polar0/"<<endl;
            cout<<"/RC/polar90/"<<endl;

            return 0;
        }
        else if(arg=="-polar")
        {
            cout<<"Polarisation mode  activated, reconstruction  of the 4 states : RC0,RC90,LC0,LC90"<<endl;
            b_polar=true;
            gui_config_name="gui_tomo_polar.conf";
        }
        else if(arg == "-etat" && i + 1 < argc)
        {
            etat_polar = argv[++i];
            cout<<"etat polar="<<etat_polar <<endl;
            gui_config_name="gui_tomo_polar.conf";
        }
        else gui_config_name="gui_tomo.conf";
    }

    manip m1(gui_config_name); //Class containing all useful information about the set up (read in config_files.txt) and reconstruction (read in recon.txt) parameters
    string chemin_result=m1.chemin_result,chemin_acquis;
    ///-----------------init Polarisation calculation if needed-------------------------------
    if(b_polar==true)
    {
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
    string str_sav_param_path=chemin_result+"/SAV_param_manip.txt";

    cout<<"camdimROI="<<m1.CamDimROI<<endl;
    Var2D const dimROI= {m1.CamDimROI,m1.CamDimROI}, coin= {0,0};
    Point2D const dimHolo(m1.CamDimROI,m1.CamDimROI,m1.CamDimROI);
    size_t const NbPixROI2d=dimROI.x*dimROI.y;
    //tableaux hologramme et réference en 1024x1024
    vector<double> static holo1(NbPixROI2d), intensite_ref(NbPixROI2d);
    char charAngle[4+1];
    ///-----------Init FFTW Holo---------------
    size_t nb_thread_fftw=m1.nbThreads;
    int fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(nb_thread_fftw);

    FFTW_init param_fftw2D_r2c_Holo(holo1,"r2c",nb_thread_fftw);/// /!\ overload of init fftw for real to complex (r2c) calculation

    vector<double> static const masqueTukeyHolo(tukey2D(dimROI.x,dimROI.y,0.1));//Tukey mask for fft apodisation. Alpha=0.1 control the apodisation width

    ///------Init variable champ complexe used after off axis extraction-----------------------------
    const size_t  NbPixUBorn=2*m1.NXMAX*2*m1.NXMAX, nbAngle=m1.NbAngle;//dimensions
    cout<<"nbAngle="<<nbAngle<<endl<<"m1.NXMAX="<<m1.NXMAX<<endl;//number of holograms
    size_t nbAngleOk=0;//Nbangle reellement utilisé
    Var2D const dim2DHA= {(size_t)2*m1.NXMAX,(size_t)2*m1.NXMAX},coinHA={m1.circle_cx-m1.NXMAX,m1.circle_cy-m1.NXMAX},coinHA_shift= {m1.fPortShift.x-m1.NXMAX,m1.fPortShift.y-m1.NXMAX};
    Var2D posSpec= {0,0},decal2DHA= {m1.NXMAX,m1.NXMAX};
    cout<<"off axis="<<dim2DHA.x<<endl;

    vector<double>  TF_champMod(NbPixUBorn), centre(NbPixUBorn), visibility2Dimg(NbPixUBorn);///Uborn modulus, specular coordinate, fringes visibility

    ///--------------------Init Reference Amplitude (amplitude only correction, useless with a blank acquisition)--------------------------------------
    vector<double> ampli_ref(dim2DHA.x*dim2DHA.y, 1.0);  // valeur neutre par défaut
    ampli_ref=initRef(chemin_acquis+"/Intensite_ref.pgm", coin, dimROI, dim2DHA);

    auto start_decoupeHA = std::chrono::system_clock::now();///démarrage chrono Hors-axe

    //#pragma omp parallel for reduction(cpt)
    FILE* test_existence;//test for fiule existance
    unsigned short int cptAngle=0;
    //#pragma omp parallel forTF2D_r2c
    float alpha=0.1;//coefficient to adjust tukey mask width
    vector<double> masqueTukeyHA(tukey2D(dim2DHA.x,dim2DHA.y,alpha));
///---------------------------------------Check holograms visibility-------------------------------
    vector<double> visibilityRaw(m1.NbAngle);//visiblity in a raw table
    visibilityRaw=checkVisibility(holo1,  m1, nbAngleOk,dim2DHA,coinHA, masqueTukeyHolo, param_fftw2D_r2c_Holo);

   /* string path_to_visibility=m1.chemin_acquis+"visibilityRaw.raw";
    cout<<"path to visibility"<<path_to_visibility<<endl;
    if(is_readable(path_to_visibility)){
       cout<<"visibility table already calculated"<<endl;
       visibilityRaw=lire_bin(path_to_visibility,64,nbAngle);
    }
    else{
        cout<<"Checking fringes visibility"<<endl;
        visibilityRaw=calcVisibility(holo1,m1, dim2DHA,coinHA,masqueTukeyHolo,param_fftw2D_r2c_Holo);
        vector<double> calcVisibility(vector<double>  &holo1, manip const &m1, Var2D dim2DHA,Var2D coinHA, vector<double> const &tukeyHolo, FFTW_init  &param_fftw2DHolo);
        SAV2(visibilityRaw,chemin_acquis+"/visibilityRaw.raw",t_double,"wb");
    }
    double seuilcontraste=m1.minVisibility;
    for(int cptHolo=0;cptHolo<nbAngle;cptHolo++){
      if(visibilityRaw[cptHolo]>seuilcontraste) nbAngleOk++;
    }
    cout<<"nbAngleOk="<<nbAngleOk<<endl;*/
///----------------------------------------Off axis extraction------------------------------
    vector<complex<double>> TF_UBornTot(NbPixUBorn*nbAngle);///stack of wrapped complex fields
    for(cptAngle=0; cptAngle<nbAngle; cptAngle++)
    {
        if((cptAngle-100*(cptAngle/100))==0)    cout<<cptAngle<<endl;
        sprintf(charAngle,"%03i",cptAngle);
        string nomFichierHolo=chemin_acquis+"/i"+charAngle+".pgm";
        test_existence = fopen(nomFichierHolo.c_str(), "rb");
        if(test_existence==NULL)
        {
            continue;
        }
        fclose(test_existence);
        charger_image2D_OCV_UNI(holo1,nomFichierHolo, coin, dimROI);
        for(size_t cpt=0; cpt<NbPixROI2d; cpt++)
        {
            holo1[cpt]=holo1[cpt]*masqueTukeyHolo[cpt];
        }

        holo2TF_UBornTukeyHA_r2c(holo1, TF_UBornTot,dimROI, dim2DHA, coinHA, cptAngle,masqueTukeyHA, param_fftw2D_r2c_Holo);
    }
///----------------------------------------END Off axis extraction------------------------------

    auto end_decoupeHA = std::chrono::system_clock::now();
    auto elapsed = end_decoupeHA - start_decoupeHA;
    std::cout <<"Temps pour FFT holo+découpe Spectre= "<< elapsed.count()/(pow(10,9)) << '\n';
    //SAVCplx(TF_UBornTot,"Re",chemin_result+"/TF_Uborn_Tot_GPU_250x250x599x64_orig.raw",t_float,"wb");


    m1.dimImg=to_string(dim2DHA.x)+"x"+to_string(dim2DHA.y)+"x"+to_string(nbAngleOk);//string for off axis field dimension
    cout<<"m1.dimImg="<<m1.dimImg<<endl;
    deleteCplxField(chemin_result, m1.dimImg);

///---------------------------------------------init phase, filed, amplitude variables etc.--------

    vector<complex<double>> TF_UBorn(NbPixUBorn),  UBorn(NbPixUBorn);
    vector<complex<double>> UBorn_unwrapped(NbPixUBorn), TF_UBorn_unwrapped(NbPixUBorn);//Ubornunwrapped to control spectrum after unwrapping. useless for reconstruciton.
    vector<double> UBornAmpCorr(NbPixUBorn);//Ubornunwrapped to control spectrum after unwrapping. useless for reconstruciton.
    vector<double> phase_2Pi_vec(NbPixUBorn),  UnwrappedPhase(NbPixUBorn),PhaseFinal(NbPixUBorn);
    double *UnwrappedPhase_herraez=new double[NbPixUBorn];

    FFTW_init param_fftw2D_c2r_HA(TF_UBorn,m1.nbThreads);//OUTPLACE

    //FFTW_init param_fftw2D_c2r_HA(dim2DHA,1,m1.nbThreads);//INPLACE
    //FFTW_init param_fftw2D_r2c_HA(phase_2Pi_vec,"r2c",m1.nbThreads);

///---------------------phase unwrapping + aberration correction---------------------------------------------------------------------------
    vector<double> phase_2Pi_vec_double(4*NbPixUBorn);//variable created only to init param_fftw2D_r2c_HA_double
    FFTW_init param_fftw2D_c2r_HA_double(phase_2Pi_vec_double,m1.nbThreads);
    ///variable pour correction aberration
    cout<<"chemin_mask==============="<<Chemin_mask<<endl;
    Mat src=Mat(1, ampli_ref.size(), CV_64F, ampli_ref.data()), mask_aber=init_mask_aber(Chemin_mask,chemin_acquis,dim2DHA);

    if(fs::exists(Chemin_mask)) //mask used to exclude object from aberration correction
    {
        string info="An aberration Mask has been used----------------------------------------------------------------";
        sav_param2D(info,str_sav_param_path);
    }

    size_t NbPtOk=countM(mask_aber),  degre_poly=5, nbCoef = sizePoly2D(degre_poly);//Nb coef poly
    Mat polynomeUs_to_fit(Size(nbCoef,NbPtOk), CV_64F);///(undersampled) Polynome to fit= function to fit (We use a polynome). we have to generate a table containing polynome_to_fit=[1,x,x^2,xy,y^2] for each coordinate (x,y)
    Mat polynome_to_fit(Size(nbCoef,dim2DHA.x*dim2DHA.y), CV_64F);

    string str_degre_poly="degré poly aberration="+to_string(degre_poly);
    sav_param2D(str_degre_poly,str_sav_param_path);

    initCorrAber(Chemin_mask, mask_aber, degre_poly,dim2DHA,polynome_to_fit,polynomeUs_to_fit);

///------------------------------------------Init variables used to calculate unwrapped complex field------------------------------------------------------------------
    cout<<"\n#########################Calculate complex 2D field (BNorn or Rytov) + aberrations correction#############################"<<endl;
    cout<<"NbPixUborn="<<NbPixUBorn<<endl;
    vector<complex<double>> UBornFinal(NbPixUBorn), UBornFinalDecal(NbPixUBorn), TF_UBorn_norm(NbPixUBorn);

    vector<double> UBornAmpFinal(NbPixUBorn),  UBornAmp(NbPixUBorn);

    vector<double> tabPosSpec(nbAngleOk*2);  ///table of speucalr coordinate (need for reconstruction binary)
    vector<vecteur>  kvect_shift(init_kvect_shift(dim2DHA));///init  differentiation operator kvect, used to calculate gradient with fft
    vector<double> kvect_mod2Shift(init_kvect_mod2Shift(kvect_shift));
    vector<vecteur> double_kvect_shift=init_kvect_shift({2*dim2DHA.x,2*dim2DHA.y}); ///------kvect_shift for symetrized image
    Var2D doubled_dimHA= {dim2DHA.x*2,dim2DHA.y*2};
    FFTW_init fftw_c2c_doubleHA(doubled_dimHA,m1.nbThreads);///init fftw for symetrized complex image
    auto start_part2= std::chrono::system_clock::now();
    double alpha_damp=3e-4;///damping factor for grad U/U regularisation in pahse unwrapping

///------------------------------------------Loop on wrapped complex field------------------------------------------------------------------
    //#pragma omp parallel for
    for(size_t cpt_angle=0; cpt_angle<nbAngleOk; cpt_angle++)  //boucle sur tous les angles : correction aberrations
    {
        ///Récupérer la TF2D dans la pile de spectre2D//get back 2D spectrum in the stack (dim x dim x Number_holograms)
        TF_UBorn.assign(TF_UBornTot.begin()+NbPixUBorn*cpt_angle, TF_UBornTot.begin()+NbPixUBorn*(cpt_angle+1));
        //SAVCplx(TF_UBorn,"Re",chemin_result+"/TF_Uborn_iterateur_Re_208x208x500x32.raw",t_float,"a+b");

        //Recherche de la valeur maximum du module dans ref non centré-----------------------------------------
        size_t cpt_max=coordSpec(TF_UBorn, TF_champMod,decal2DHA);
        double  max_part_reel = TF_UBorn[cpt_max].real(),///sauvegarde de la valeur cplx du spéculaire/save complex value of the specular beam
                max_part_imag = TF_UBorn[cpt_max].imag(),
                max_module = sqrt(TF_UBorn[cpt_max].imag()*TF_UBorn[cpt_max].imag()+TF_UBorn[cpt_max].real()*TF_UBorn[cpt_max].real());
        const int kxmi=cpt_max%(2*m1.NXMAX), kymi=cpt_max/(2*m1.NXMAX);
        posSpec= {kxmi,kymi}; ///coord informatique speculaire

        visibility2Dimg[kxmi*2*m1.NXMAX+kymi]=visibilityRaw[cpt_angle];
        ///calculate phi and theta (angles of illumination)
        Var2D posSpecH= {kxmi-m1.NXMAX,kymi-m1.NXMAX};
        /* float kiz=sqrt(pow(m1.rayon,2)-pow(posSpecH.x,2)-pow(posSpecH.y,2));
         float theta_bis=acos(kiz/m1.rayon)*180/M_PI;
         float phi_bis=atan2(static_cast<double>(posSpecH.y),static_cast<double>(posSpecH.x));
          cout<<"(kix,kiy,kiz)=("<<posSpecH.x<<","<<posSpecH.y<<","<<kiz<<")"<<endl;
          cout<<"num angle="<<cpt_angle<<", phi_bis="<<phi_bis*180/M_PI<<",theta_bis="<<theta_bis<<endl;*/
        ///calculate angle theta, phi with vecteur class
        vecteur kmi(posSpecH.x,posSpecH.y,round(sqrt(m1.rayon*m1.rayon-posSpecH.x*posSpecH.x-posSpecH.y*posSpecH.y)));//incident vector (pixel)
        kmi.setNorm(m1.rayon);
        kmi.calc_angle();
        double theta=kmi.theta;
        double phi=kmi.phi;

        tabPosSpec[cpt_angle]=(double)posSpec.x;
        tabPosSpec[cpt_angle+nbAngleOk]=(double)posSpec.y;
        centre[kxmi*2*m1.NXMAX+kymi]=cpt_angle;
        ///calculate phase and correct aberration on the illumination beam
        if(m1.b_CorrAber==true)
        {
            calc_Uborn2(TF_UBorn,UBorn,dim2DHA,posSpec,param_fftw2D_c2r_HA);
            //vector<complex<double>> TF_UBorn_I(dim2DHA.x*dim2DHA.y);
            //TF_UBorn_I=calc_Uborn2exportTF(TF_UBorn,UBorn,dim2DHA,posSpec,param_fftw2D_c2r_HA);
            if(m1.b_ampliRef==1)
            {
                correctAmpliRef(UBorn,ampli_ref);
            }
            // SAVCplx(UBorn,"Im",chemin_result+"/UBorn_Im_debut_extract.raw",t_float,"a+b");
            ///----------Calculate phase + unwrapping---------------------------------------
            calcPhase_mpi_pi_atan2(UBorn,phase_2Pi_vec); ///fonction atan2 for Herraez unwrapping
            //SAV2(phase_2Pi_vec,chemin_result+"/phasePI_atan2.raw",t_float,"a+b");
            //phase2pi(UBorn, dim2DHA,phase2Pi);//
            if(m1.b_Deroul==true)
            {
                if(m1.b_volkov==0)
                {
                    phaseUnwrapping_Mat(dim2DHA, phase_2Pi_vec, UnwrappedPhase_herraez);
                    for(size_t cpt=0; cpt<NbPixUBorn; cpt++)
                        UnwrappedPhase[cpt]=UnwrappedPhase_herraez[cpt];///plutôt passer pointeur ?
                }
                else
                {
                    //deroul_volkov4_total_sym_paire(phase_2Pi_vec,UnwrappedPhase, double_kvect_shift, param_fftw2D_c2r_HA_double);
                    UnwrappedPhase=deroul_volkov6_total_sym_paire_gradu(UBorn,double_kvect_shift,param_fftw2D_c2r_HA_double,alpha_damp);//"exact" unwrapping with damped division of (grad U)/U
                    /*  for(int cpt=0;cpt<UBorn.size();cpt++){
                      UBorn_unwrapped[cpt].real(abs(UBorn[cpt]) *cos(UnwrappedPhase[cpt]));//*masqueTukeyHolo[cpt]);
                      UBorn_unwrapped[cpt].imag(abs(UBorn[cpt]) *sin(UnwrappedPhase[cpt]));//*masqueTukeyHolo[cpt]);/
                      }*/
                    // SAVCplx(UBorn_unwrapped,"Im",chemin_result+"/UBorn_Im_unwrraped_208x208x500.raw",t_float,"a+b");
                    // TF2Dcplx(UBorn_unwrapped,TF_UBorn_unwrapped,param_fftw2D_c2r_HA);
                    //SAVCplx(fftshift2D(TF_UBorn_unwrapped),"Im",chemin_result+"/TF_UBorn_Im_unwrapped_208x208x500.raw",t_float,"a+b");
                }
            }
            else UnwrappedPhase=phase_2Pi_vec;
            // SAV2(UnwrappedPhase,chemin_result+"/phase_deroul_volkov_avant_corr_aber_vvvvolkov3.raw",t_float,"a+b");
            //---------------------polynomial correction-------------------------------
            src=Mat(dim2DHA.x,dim2DHA.y,CV_64F, UnwrappedPhase.data());
            // auto start_calcAber = std::chrono::system_clock::now();
            Mat Phase_corr(aberCorr2(src, mask_aber,polynomeUs_to_fit,polynome_to_fit));

            //  auto end_calcAber = std::chrono::system_clock::now();
            //auto elapsed = end_calcAber - start_calcAber;
            //std::cout <<"Temps pour FFT holo+découpe Spectre= "<< elapsed.count()/(pow(10,9)) << '\n';
            // SAV2((double*)Phase_corr.data,Phase_corr.rows*Phase_corr.cols,"/home/mat/tmp/phaseCorr_main.raw",t_float,"a+b");

            ///--------------- Amplitude Normalisation-----------------------------------------------------------------------------

            for(size_t cpt=0; cpt<(NbPixUBorn); cpt++)
                UBornAmp[cpt]=abs(UBorn[cpt]);
            SAV2(UBornAmp,chemin_result+"/UBornAmp.raw",t_float,"a+b");
            Mat srcAmp=Mat(dim2DHA.x, dim2DHA.y, CV_64F, UBornAmp.data());///image source
            // Mat UBornAmp_corr(ampliCorr2(srcAmp, polynomeUs_to_fit, polynome_to_fit, mask_aber));///amplitude result
            double gamma=0.01; //epsilon controlant la division=gamma*max(amplitude)
            Mat UBornAmp_corr(ampliCorr3(srcAmp, polynomeUs_to_fit, polynome_to_fit, mask_aber,gamma));///résultat amplitude avec gestion de la division+inversion robuste

            ///------------------End amplitude normalisation-----------------------------------------------------------------------------
            ///export spectrum to control
            for(size_t y=0; y<dim2DHA.y; y++)  // reconstruire l'onde complexe/Recalculate the complex field
            {
                for(size_t x=0; x<dim2DHA.x; x++)
                {
                    size_t cpt=x+y*dim2DHA.x;
                    PhaseFinal[cpt]=Phase_corr.at<double>(y,x) ;//copie opencV->Tableau
                    UBornAmpCorr[cpt]=UBornAmp_corr.at<double>(y,x) ;//copie opencV->Tableau

                    UBorn_unwrapped[cpt].real( UBornAmpCorr[cpt]*cos(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);///correction amplitude
                    UBorn_unwrapped[cpt].imag( UBornAmpCorr[cpt]*sin(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);///correction amplitude
                }
            }
            TF2Dcplx(UBorn_unwrapped,TF_UBorn_unwrapped,param_fftw2D_c2r_HA);
            //SAVCplx(fftshift2D(TF_UBorn_unwrapped),"Im",chemin_result+"/TF_UBorn_Im_AMPcorr_208x208x500.raw",t_float,"a+b");
            int flag=0;
            for(size_t y=0; y<dim2DHA.y; y++)  // reconstruire l'onde complexe/Recalculate the complex field
            {
                for(size_t x=0; x<dim2DHA.x; x++)
                {
                    size_t cpt=x+y*dim2DHA.x;
                    PhaseFinal[cpt]=Phase_corr.at<double>(y,x) ;//copie opencV->Tableau
                    UBornAmpFinal[cpt]=UBornAmp_corr.at<double>(y,x);
                    if(UBornAmpFinal[cpt]>10) UBornAmpFinal[cpt]=1;//Amplitude should always be around 1 after correction since Amplitude illumination>>Amplitude ohject
                    if(UBornAmpFinal[cpt]<-10)  UBornAmpFinal[cpt]=1;
                    // if(UBornAmpFinal[cpt]<-1 && flag==0) { cout<<"Holo numero"<<cpt_angle<<endl;flag=1;}
                    if(m1.b_Born==true)  //UBORN=U_tot-U_inc=u_tot_norm-1
                    {
                        // UBornFinal[cpt].real( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])- 1 )*cos(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);///correction amplitude
                        // UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]) )*sin(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);
                        UBornFinal[cpt].real( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*cos(PhaseFinal[cpt]) -1);//*masqueTukeyHolo[cpt]);///correction amplitude
                        // UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*sin(PhaseFinal[cpt]) -1);//*masqueTukeyHolo[cpt]);
                        UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*sin(PhaseFinal[cpt]));
                    }
                    else  //RYTOV URytov = log a_t/a_i (=log a_t après correction AmpliCorr)
                    {
                        UBornFinal[cpt].real(log(sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])));
                        UBornFinal[cpt].imag(PhaseFinal[cpt]);
                    }
                }
            }
            //SAV2(PhaseFinal,chemin_result+"/phase_finale.raw",t_float,"a+b");
            ///Recalculer la TF décalée pour le programme principal.
            Var2D recal= {kxmi,kymi};
            // decal2DCplxGen2(UBornFinal,UBornFinalDecal, decal2DHA);
            //TF2Dcplx_vec(in_HA,out_HA,UBornFinalDecal,TF_UBorn_norm,p_forward_HA);
            //  TF2Dcplx(UBornFinalDecal,TF_UBorn_norm,param_fftw2D_c2r_HA);
            //  SAVCplx(TF_UBorn_norm,"Re", chemin_result+"/TF_UBorn_norm_Re"+m1.dimImg+".raw", t_double, "a+b");
            SAVCplx(UBornFinal,"Re", chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw", t_double, "a+b");
            SAVCplx(UBornFinal,"Im", chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw", t_double, "a+b");
        }
        else  ///sauvegarde onde avec aberration
        {
            cout<<"normalisation avec spec"<<endl;
            for(size_t cpt=0; cpt<(4*m1.NXMAX*m1.NXMAX); cpt++)  //correction phase à l'ordre zéro et normalisatoin amplitude par ampli_spec=ampli_inc*ampli_ref
            {
                TF_UBorn_norm[cpt].real((TF_UBorn[cpt].real()*max_part_reel+TF_UBorn[cpt].imag()*max_part_imag)/max_module);
                TF_UBorn_norm[cpt].imag((TF_UBorn[cpt].imag()*max_part_reel-TF_UBorn[cpt].real()*max_part_imag)/max_module);
                calc_Uborn2(TF_UBorn_norm,UBorn,dim2DHA,posSpec,param_fftw2D_c2r_HA);
            }
            if(m1.b_Born==1)
            {
                //SAVCplx(TF_UBorn,"Re", chemin_result+"/TF_Uborn_Re.raw", t_double, "a+b");
                SAVCplx(UBorn,"Re", chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw", t_double, "a+b");
                SAVCplx(UBorn,"Im", chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw", t_double, "a+b");
            }
            else
            {
            cout<<"Rytov approximation need a phase unwrapping"<<endl;
            }
        }
    }//fin de boucle for sur tous les angles

//SAV2(UBornAmpFinal,chemin_result+"/UBornAmpFinal.raw",t_float,"a+b");
    auto end_part2= std::chrono::system_clock::now();
    auto elapsed_part2 = end_part2 - start_part2;
    std::cout <<"Temps pour part2= "<< elapsed_part2.count()/(pow(10,9)) << '\n';
    delete[] UnwrappedPhase_herraez;

///-------------Save useful paramaters for tomo_reconstruction and log-------------------------
    cout<<"NXMAX sauve="<<m1.NXMAX<<endl;
    vector<double> param{m1.NXMAX,nbAngleOk,m1.rayon,dimROI.x,m1.tailleTheoPixelHolo};
    SAV2(param,chemin_result+"/parametres.raw", t_double, "wb");
    SAV2(tabPosSpec,chemin_result+"/tab_posSpec.raw",t_double,"wb");///speculaire coordinate with computer coordinate [0->2NXMAX-1]
    SAV_Tiff2D(centre,chemin_result+"/centres.tif",m1.NA/m1.NXMAX); //export speucalr position in a  2D image
    SAV_Tiff2D(visibility2Dimg,chemin_result+"/visibility.tif",m1.NA/m1.NXMAX);
    cout<<"End preprocessing, just run Reconstruction now"<<endl;
    cout<<"Results in "<<chemin_result<<endl;
    return 0;
}
