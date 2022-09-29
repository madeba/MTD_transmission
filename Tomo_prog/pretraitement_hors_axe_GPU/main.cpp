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
#include "deroulement_volkov3_GPU.h"
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
#include "GPU_fonctions.h"
#define GPU
//#include <cv.h>
using namespace std;

int main(int argc, char** argv){
    manip m1; //créer un objet manip
    string chemin_result=m1.chemin_result, chemin_acquis=m1.chemin_acquis, Chemin_mask=m1.chemin_acquis+"/Image_mask.pgm";
    //dimensions hologrammes
    cout<<"camdimROI"<<m1.CamDimROI;
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
    FFTW_init param_fftw2D_r2c_Holo(holo1,"r2c",nb_thread_fftw);/// /!\ init r2c->surcharge avec image2D en entrée!
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
    ampli_ref=initRef(chemin_acquis+"/Intensite_ref.pgm",coin,dimROI);
    //ampli_ref=initRef("/opt/Acquis/phytolites/acquis2/ACQUIS/Intensite_ref.pgm",coin,dimROI);

    ///------Init variable champ complexe après découpe hors axe-----------------------------
    const size_t  NbPixUBorn=4*m1.NXMAX*m1.NXMAX, NbAngle=m1.NbAngle;//dimensions
    cout<<"NbAngle="<<NbAngle<<endl<<"m1.NXMAX="<<m1.NXMAX<<endl;//nombre d'hologrammes
    size_t NbAngleOk=0;//Nbangle reellement utilisé
    Var2D const dim2DHA= {(size_t)2*m1.NXMAX,(size_t)2*m1.NXMAX},coinHA= {m1.circle_cx-m1.NXMAX,m1.circle_cy-m1.NXMAX},coinHA_shift={m1.fPortShift.x-m1.NXMAX,m1.fPortShift.y-m1.NXMAX};
    Var2D posSpec= {0,0},decal2DHA= {m1.NXMAX,m1.NXMAX};;
    cout<<"dimension decoupe dimDHA="<<dim2DHA.x<<endl;
    vector<complex<double>> TF_UBornTot(NbPixUBorn*NbAngle);///variable stockant les N champs complexes decoupés depuis la zone 1024 (pour utilsier wisdom en 1024)
    vector<double>  TF_champMod(NbPixUBorn), centre(NbPixUBorn);///centre pour controler balayage

    float alpha=0.025;//coeff pour le masque de tuckey
vector<double> masqueTukeyHA(tukey2D(dim2DHA.x,dim2DHA.y,alpha));

af::array masqueTukey_GPU=tukey2D_GPU(dimROI.x,dimROI.y,alpha);
af::array TukeyBorn_GPU=tukey2D_GPU(dim2DHA.x,dim2DHA.y,alpha);
af::array TukeyBorn_cplx_GPU=af::complex(TukeyBorn_GPU,TukeyBorn_GPU);
//display(masqueTukey_GPU);
af::array TF_UBornTot_GPU(dim2DHA.x,dim2DHA.y,NbAngle,c64);

//display(real(TukeyBorn_cplx_GPU));
FILE* test_existence;//tester l'existence des fichiers
unsigned short int cptAngle=0;
        auto start_decoupeHA = std::chrono::system_clock::now();///démarrage chrono Hors-axe

   try {
        int device = argc > 1 ? atoi(argv[1]) : 0;
        af::setDevice(device);
        af::info();
        af::array holoGPU(dimHolo.x,dimHolo.y,f64);

        for(cptAngle=0; cptAngle<NbAngle; cptAngle++)
        {
            if((cptAngle-100*(cptAngle/100))==0)    cout<<cptAngle<<endl;
            sprintf(charAngle,"%03i",cptAngle);
            string nomFichierHolo=m1.chemin_acquis+"/i"+charAngle+".pgm";
            test_existence = fopen(nomFichierHolo.c_str(), "rb");
            // cout<<nomFichierHolo<<endl;
            if(test_existence!=NULL)
            {
                fclose(test_existence);

               // chargeImageGPU(holoGPU,nomFichierHolo,coin);
               #ifdef GPU
                chargeImageGPU_via_OCV(holoGPU,nomFichierHolo,coin);

               // holoGPU=transpose(holoGPU);///to compare with non GPU code using opencv
                 // af::timer start1=af::timer::start();
                holo2TF_UBorn_GPU(holoGPU, TF_UBornTot_GPU,dimROI, dim2DHA, coinHA, NbAngleOk,masqueTukey_GPU,TukeyBorn_cplx_GPU);
                  //  printf("elapsed seconds: %g\n", af::timer::stop(start1));
                  #endif // GPU
                NbAngleOk++;
            }
            else cout<<"fichier "<<cptAngle<<" inexistant\n";
        }

            ///------------------------------------------

            //string tampon="/UBornfinal_Im.raw";
        m1.dimImg=to_string(dim2DHA.x)+"x"+to_string(dim2DHA.y)+"x"+to_string(NbAngleOk);
        deleteCplxField(chemin_result, m1.dimImg);
        TF_UBornTot_GPU.host(TF_UBornTot.data());

      //  SAV3D_Tiff(TF_UBornTot,"Re",chemin_result+"/TF_Uborn_Tot_GPU_250x250x599x64.raw",1);
      // SAVCplx(TF_UBornTot,"Re",chemin_result+"/TF_Uborn_Tot_RE_GPU_250x250x599x64.raw",t_float,"wb");
       // SAVCplx(TF_UBornTot,"Im",chemin_result+"/TF_Uborn_Tot_Im_GPU_250x250x599x64.raw",t_float,"wb");
      //  auto end_decoupeHA = std::chrono::system_clock::now();
      //  auto elapsed = end_decoupeHA - start_decoupeHA;
       // std::cout <<"Temps pour FFT holo+découpe Spectre= "<< elapsed.count()/(pow(10,9)) << '\n';
   }
      catch (af::exception& e)
    {
        fprintf(stderr, "%s\n", e.what());
    }


///initialiser les variables de champ, phase, amplitude etc.

vector<complex<double>> TF_UBorn(NbPixUBorn),  UBorn(NbPixUBorn);
vector<double> phase_2Pi_vec(NbPixUBorn),  UnwrappedPhase(NbPixUBorn),PhaseFinal(NbPixUBorn);
double *UnwrappedPhase_herraez=new double[NbPixUBorn];


auto start_init_var = std::chrono::system_clock::now();///démarrage chrono Hors-axe

FFTW_init param_fftw2D_c2r_HA(TF_UBorn,m1.nbThreads);//OUTPLACE
//FFTW_init param_fftw2D_c2r_HA(dim2DHA,1,m1.nbThreads);//INPLACE
//FFTW_init param_fftw2D_r2c_HA(phase_2Pi_vec,"r2c",m1.nbThreads);

auto end_init_var = std::chrono::system_clock::now();
        auto elapsed = end_init_var - start_init_var;
        std::cout <<"Temps pour init var= "<< elapsed.count()/(pow(10,9)) << '\n';


///variable pour correction aberration
Mat src=Mat(1, ampli_ref.size(), CV_64F, ampli_ref.data()), mask_aber=init_mask_aber(Chemin_mask,chemin_acquis,dim2DHA);


size_t NbPtOk=countM(mask_aber),  degre_poly=3, nbCoef = sizePoly2D(degre_poly);//Nb coef poly
Mat polynomeUs_to_fit(Size(nbCoef,NbPtOk), CV_64F);///(undersampled) Polynome to fit= function to fit (We use a polynome). we have to generate a table containing polynome_to_fit=[1,x,x^2,xy,y^2] for each coordinate (x,y)
Mat polynome_to_fit(Size(nbCoef,dim2DHA.x*dim2DHA.y), CV_64F);

initCorrAber(Chemin_mask, mask_aber, degre_poly,dim2DHA,polynome_to_fit,polynomeUs_to_fit);

cout<<"\n#########################Calcul champs cplx 2D Uborn/Rytov + eventuelle Correction aberrations#############################"<<endl;

vector<complex<double>> UBornFinal(NbPixUBorn), UBornFinalDecal(NbPixUBorn), TF_UBorn_norm(NbPixUBorn);
vector<double> UBornAmpFinal(NbPixUBorn),  UBornAmp(NbPixUBorn);
vector<double> tabPosSpec(NbAngleOk*2);  ///stockage des spéculaires pour exportation vers reconstruction

///init opérateur différentiation kvect
vector<vecteur>  kvect_shift(init_kvect_shift(dim2DHA));
/*vector<double> kvect_shift_x(NbPixUBorn);
for(int cpt=0;cpt<NbPixUBorn;cpt++)
{
    kvect_shift_x[cpt]=kvect_shift[cpt].x;
}*/
vector<double> kvect_mod2Shift(init_kvect_mod2Shift(kvect_shift));


//SAV_Tiff2D(kvect_mod2Shift,"/home/mat/tomo_test/kvectmod2Shift.tif",1);
//SAV_Tiff2D(kvect_shift_x,"/home/mat/tomo_test/kvect_shift_X.tif",1);

auto start_part2= std::chrono::system_clock::now();
                    // display(phase_2Pi_vec_GPU);
                    af::array UnwrappedPhase_GPU(dim2DHA.x,dim2DHA.y,f64);//->a mettre en retour de deroul_volkov?
                    af::array kvect_shiftX_GPU=init_kvect_shiftX_GPU(dim2DHA.x);
                    af::array kvect_shiftY_GPU=transpose(kvect_shiftX_GPU);
                    af::array phase_2Pi_vec_GPU(dim2DHA.x,dim2DHA.y,f64);

for(size_t cpt_angle=0; cpt_angle<NbAngleOk; cpt_angle++){ //boucle sur tous les angles : correction aberrations
        TF_UBorn.assign(TF_UBornTot.begin()+NbPixUBorn*cpt_angle, TF_UBornTot.begin()+NbPixUBorn*(cpt_angle+1));  ///Récupérer la TF2D dans la pile de spectre2D
        // SAVCplx(TF_UBorn,"Re",chemin_result+"/TF_Uborn_iterateur_Re_220x220x599x32.raw",t_float,"a+b");
        //SAVCplx(TF_UBorn,"Im",chemin_result+"/TF_Uborn_iterateur_Im_250x250x599x32.raw",t_float,"a+b");
        //Recherche de la valeur maximum du module dans ref non centré-----------------------------------------
        size_t cpt_max=coordSpec(TF_UBorn, TF_champMod,decal2DHA);
        double  max_part_reel = TF_UBorn[cpt_max].real(),///sauvegarde de la valeur cplx des  spéculaires
        max_part_imag = TF_UBorn[cpt_max].imag(),
        max_module = sqrt(TF_UBorn[cpt_max].imag()*TF_UBorn[cpt_max].imag()+TF_UBorn[cpt_max].real()*TF_UBorn[cpt_max].real());
        const int kxmi=cpt_max%(2*m1.NXMAX), kymi=cpt_max/(2*m1.NXMAX);
        posSpec= {kxmi,kymi}; ///coord informatique speculaire
        tabPosSpec[cpt_angle]=(double)posSpec.x;
        tabPosSpec[cpt_angle+NbAngleOk]=(double)posSpec.y;
        centre[kxmi*2*m1.NXMAX+kymi]=cpt_angle;

        if(m1.b_CorrAber==true)
        {
            calc_Uborn2(TF_UBorn,UBorn,dim2DHA,posSpec,param_fftw2D_c2r_HA);
            //  SAVCplx(UBorn,"Re",chemin_result+"/ampli_UBorn_debut_extract.raw",t_float,"a+b");
            ///----------Calcul phase + deroulement--------------------------------
            calcPhase_mpi_pi_atan2(UBorn,phase_2Pi_vec); ///fonction atan2
            //  SAV2(phase_2Pi_vec,chemin_result+"/phasePI_atan2.raw",t_float,"a+b");
            //phase2pi(UBorn, dim2DHA,phase2Pi);//asin
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
                   /// auto start_volkov = std::chrono::system_clock::now();
                    // deroul_volkov2(phase_2Pi_vec,UnwrappedPhase, param_fftw2D_c2r_HA,param_fftw2D_r2c_HA);    //  auto end_calcAber = std::chrono::system_clock::now();
                    //deroul_volkov3(phase_2Pi_vec,UnwrappedPhase, kvect_shift, param_fftw2D_c2r_HA);
                    phase_2Pi_vec_GPU=af::array(dim2DHA.x,dim2DHA.y,phase_2Pi_vec.data());
                    //  SAV_Tiff2D(kvect_shiftX_GPU,"/home/mat/tomo_test/kvect_shift_X_GPU.tif",1);
                    deroul_volkov3_GPU2(phase_2Pi_vec_GPU,UnwrappedPhase_GPU,kvect_shiftX_GPU,kvect_shiftY_GPU);
                    //    SAV_Tiff2D(UnwrappedPhase_GPU,"/home/mat/tomo_test/Unwrapped_GPU.tif",1);
                    UnwrappedPhase_GPU.host(UnwrappedPhase.data());
                    //  SAV2(UnwrappedPhase,"/home/mat/UnWrapped_GPU.bin",t_float,"a+b");
                    //deroul_volkov4_AS(phase_2Pi_vec,UnwrappedPhase, kvect_shift, param_fftw2D_c2r_HA);//antisymétrie avant intégration
                    //deroul_volkov4(phase_2Pi_vec,UnwrappedPhase, kvect_shift, param_fftw2D_c2r_HA);
                  ///  auto end_volkov = std::chrono::system_clock::now();
                  ///  auto elapsed_volkov = end_volkov - start_volkov;
                   /// std::cout <<"Temps pour deroul volkov2= "<< elapsed_volkov.count()/(pow(10,9)) << '\n';
                }
            }
    else UnwrappedPhase=phase_2Pi_vec;
            //      SAV2(UnwrappedPhase,chemin_result+"/phase_deroul_volkov_avant_corr_aber_vvvvolkov3.raw",t_float,"a+b");

            //auto start_aber = std::chrono::system_clock::now();
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

            /*  auto end_aber = std::chrono::system_clock::now();
                       auto elapsed_aber = end_aber - start_aber;
                        std::cout <<"Temps pour corrAber= "<< elapsed_aber.count()/(pow(10,9)) << '\n';*/
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
         // UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]) - 1)*sin(PhaseFinal[cpt]) );//*masqueTukeyHolo[cpt]);

          UBornFinal[cpt].real( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*cos(PhaseFinal[cpt]) -1);//*masqueTukeyHolo[cpt]);///correction amplitude
          UBornFinal[cpt].imag( (sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*sin(PhaseFinal[cpt]) -1);//*masqueTukeyHolo[cpt]);
        }
        else{ //RYTOV URytov = log a_t/a_i (=log a_t après correction AmpliCorr)
          UBornFinal[cpt].real(log(sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])));///racine(U^2) pour éliminer les éventuelles amplitudes négatives
          UBornFinal[cpt].imag(PhaseFinal[cpt]);
        }
      }
    }
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
   if(m1.b_Born==0){
    cout<<"Attention : Rytov nécessite un déroulement"<<endl;}
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
SAV2(UBornAmpFinal,chemin_result+"/UBornAmpFinal.raw",t_float,"a+b");
auto end_part2= std::chrono::system_clock::now();
auto elapsed_part2 = end_part2 - start_part2;
std::cout <<"Temps pour part2= "<< elapsed_part2.count()/(pow(10,9)) << '\n';
delete[] UnwrappedPhase_herraez;

///------------SAUVER LES PARAMETRES UTILES À RECONSTRUCTION ou au contrôle//Save useful paramaters for tomo_reconstruction
cout<<"NXMAX sauve="<<m1.NXMAX<<endl;
vector<double> param{m1.NXMAX,NbAngleOk,m1.rayon,dimROI.x,m1.tailleTheoPixelHolo};
SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");
/*for(int cpt=0;cpt<NbAngle;cpt++)
{
    cout<<"tabPosSpec.y="<<tabPosSpec[cpt]<<","<<"tabPosSpec.y="<<tabPosSpec[cpt+NbAngle]<<endl;
}*/
SAV2(tabPosSpec,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");///sav speculaire en coord informatique [0->2NXMAX-1]
SAV_Tiff2D(centre,m1.chemin_result+"/centres.tif",m1.NA_obj/m1.NXMAX); //exportation des spéculaires en "unité NA"
cout<<"Fin prétraitement"<<endl;
return 0;
}
