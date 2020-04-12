#include <iostream>
//#include <Magick++.h>

#include <vector>
#include "struct.h"
#include <complex>
#include "manip.h"
#include <string>
#include "projet.h"
#include "fonctions.h"
#include "FFT_fonctions.h"
#include "deroulement_volkov.h"
#include "deroulement_herraez.h"
#include "Correction_aberration.h"
#include <chrono>
#include <cv.h>
using namespace std;

int main()
{
    manip m1; //créer un objet manip
   // m1.init();//initialiser les constantes de la manip
    string chemin_result=m1.chemin_result;
    string chemin_acquis=m1.chemin_acquis;
    string Chemin_mask=m1.chemin_acquis+"Image_mask.pgm";



    Var2D dimROI= {1024,1024}, coin= {0,0};
    size_t NbPixROI2d=dimROI.x*dimROI.y;
    vector<double> holo1(NbPixROI2d);
    vector<double> intensite_ref(NbPixROI2d);
    vector<double> ampli_ref(NbPixROI2d);
    char charAngle[4+1];

    ///--------------Init FFTW Holo-------------------------------------------------
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        fftw_plan_with_nthreads(4);
        ///prepare fftw plan+tableaux-----------------
        fftw_plan p_forward_holo;
        //fftw_complex *in_out_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//in=out pour transformation "inplace".
        fftw_complex *in_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//
        fftw_complex *out_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//
        p_forward_holo=fftw_plan_dft_2d(dimROI.x, dimROI.y, in_holo, out_holo, FFTW_FORWARD,FFTW_ESTIMATE);
        vector<double> masqueTukeyHolo(NbPixROI2d);
        masqueTukeyHolo=tukey2D(dimROI.x,dimROI.y,0.05);

    ///--------------fin init FFTW Holo-------------------------------------------------
    if(is_readable(chemin_acquis+"/Intensite_ref.pgm" )==1)
         charger_image2D_OCV(intensite_ref,chemin_acquis+"/Intensite_ref.pgm", coin, dimROI);
         else cout<<"/!\\  fichier intensité référence absent, création intensité unité"<<endl;


    for(size_t cpt=0;cpt<NbPixROI2d;cpt++)
    {
        if(intensite_ref[cpt]!=0)
        ampli_ref[cpt]=sqrt(intensite_ref[cpt]);
        else
        ampli_ref[cpt]=1;
    }
    ///------------init_ampli_ref------------------------

    ///----------------------------------------------
    ///------Init variable hologramme hors axe
    const int  NbPixUBorn=4*m1.NXMAX*m1.NXMAX;
    cout<<"m1.NXMAX="<<m1.NXMAX<<endl;
    Var2D dim2DHA= {2*m1.NXMAX,2*m1.NXMAX},coinHA= {m1.circle_cx-m1.NXMAX,m1.circle_cy-m1.NXMAX},decal2DHA= {m1.NXMAX,m1.NXMAX},posSpec= {0,0};
    cout<<"dimension decoupe dimDHA="<<dim2DHA.x<<endl;
    const int NbAngle=m1.NbAngle;//m1.Num_Angle_final-m1.premier_plan;
    vector<complex<double>> TF_UBornTot(NbPixUBorn*NbAngle);///variable stockant les N champs complexes decoupés depuis la zone 1024 (pour utilsier wisdom en 1024)
    vector<double>  TF_champMod(NbPixUBorn);//module du champ
    vector<double> centre(NbPixUBorn);///centre pour controler balayage
    size_t  NbAngleOk=0;
    cout<<"NbAngle="<<NbAngle<<endl;
    FILE* test_existence;//tester l'existence des fichiers

     auto start_decoupeHA = std::chrono::system_clock::now();///démarrage chrono
    ///Charger les acqusitions
    for(size_t cptAngle=0; cptAngle<NbAngle; cptAngle++){
            if((cptAngle-100*(cptAngle/100))==0)
                cout<<cptAngle<<endl;
            sprintf(charAngle,"%03i",cptAngle);
            string nomFichierHolo=m1.chemin_acquis+"/i"+charAngle+".pgm";
            test_existence = fopen(nomFichierHolo.c_str(), "rb");
            //cout<<nomFichierHolo<<endl;
            if(test_existence!=NULL) {
                fclose(test_existence);
                charger_image2D_OCV(holo1,nomFichierHolo, coin, dimROI);
                  for(size_t cpt=0;cpt<NbPixROI2d;cpt++){
                  holo1[cpt]=holo1[cpt]/ampli_ref[cpt];
                  }
               // holo2TF_UBorn_old(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA,NbAngleOk, masqueTukeyHolo);
               ///calculer la TF des hologrammes et la découper de dimROI à 2NXMAX
                holo2TF_UBorn(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA,NbAngleOk, masqueTukeyHolo, in_holo,out_holo,p_forward_holo);
                //holo2TF_UBorn_INPLACE(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA,NbAngleOk, masqueTukeyHolo, in_out_holo,p_forward_holo);
                NbAngleOk++;
            }
             else cout<<"fichier "<<cptAngle<<" inexistant\n";
    }
    auto end_decoupeHA = std::chrono::system_clock::now();
    auto elapsed = end_decoupeHA - start_decoupeHA;
    std::cout <<"Temps pour FFT holo+découpe Spectre= "<< elapsed.count()/(pow(10,9)) << '\n';
    ///--------------libérer allocation FFTW holo-------------------------------------
        fftw_destroy_plan(p_forward_holo);
        fftw_free(in_holo);
        fftw_free(out_holo);
        fftw_forget_wisdom();
    ///--------------Init FFTW Hors axe-------------------------------------------------
        fftw_plan p_backward_HA,p_forward_HA;
        //fftw_complex *in_out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimROI.x*dimROI.y);//
        fftw_complex *in_HA=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixUBorn);//in=out pour transformation "inplace".
        fftw_complex *out_HA=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixUBorn);//
        p_backward_HA=fftw_plan_dft_2d(dim2DHA.x, dim2DHA.y, in_HA, out_HA, FFTW_BACKWARD,FFTW_ESTIMATE);
        p_forward_HA=fftw_plan_dft_2d(dim2DHA.x, dim2DHA.y, in_HA, out_HA, FFTW_FORWARD,FFTW_ESTIMATE);
        float alpha=0.1;//coeff pour le masque de tuckey
        vector<double> masqueTukeyHA(tukey2D(dim2DHA.x,dim2DHA.y,alpha));
     ///---------fin prepare fftw HA-----------------

    ///initialiser les variables de champ, phase, amplitude etc.
    //vector<complex<double>> TF_UBorn(TF_UBornTot.begin(),TF_UBornTot.begin()+NbPixUBorn);
    vector<complex<double>> TF_UBorn(NbPixUBorn);
    vector<complex<double>> UBorn(NbPixUBorn);
    vector<double> phase_2Pi_vec(NbPixUBorn);
    vector<double> phase_deroul_vec(NbPixUBorn);
    vector<double> UnwrappedPhase(NbPixUBorn);

    double *poly_aber=new double[NbPixUBorn];
    double *UnwrappedPhase_herraez=new double[NbPixUBorn];

    vector<double> PhaseFinal(NbPixUBorn);
    Mat src=Mat(1, ampli_ref.size(), CV_64F, ampli_ref.data());
    ///Charger masque aberration
    Mat mask_aber=init_mask_aber(Chemin_mask,dim2DHA);
    mask_aber.convertTo(mask_aber, CV_8U);
    int NbPtOk=countM(mask_aber);

    vector<complex<double>> UBornFinal(NbPixUBorn);
    vector<double> UBornAmpFinal(NbPixUBorn);
    vector<double> UBornAmp(NbPixUBorn);
    vector<complex<double>> UBornFinalDecal(NbPixUBorn);
    vector<complex<double>> TF_UBorn_norm(NbPixUBorn);
    vector<double> tabPosSpec(NbAngle*2);  ///stockage des speculaires pour exportation vers reconstruction

    cout<<"\n#########################Calcul champs cplx 2D Uborn/Rytov + eventuelle Correction aberrations#############################"<<endl;

    for(size_t cpt_angle=0; cpt_angle<NbAngleOk; cpt_angle++){ //boucle sur tous les angles : correction aberrations
            ///Récupérer la TF2D dans la pile de spectre2D
            TF_UBorn.assign(TF_UBornTot.begin()+NbPixUBorn*cpt_angle, TF_UBornTot.begin()+NbPixUBorn*(cpt_angle+1));
            //SAVCplx(TF_UBorn,"Re","/home/mat/tomo_test/TF_Uborn_iterateur.raw",t_float,"a+b");
            ///Recherche de la valeur maximum du module dans ref non centré-----------------------------------------
                int cpt_max=coordSpec(TF_UBorn, TF_champMod,decal2DHA);
                double  max_part_reel = TF_UBorn[cpt_max].real(),///sauvegarde de la valeur cplx des  spéculaires
                        max_part_imag = TF_UBorn[cpt_max].imag(),
                            //max_module = TF_champMod[cpt_max];
                        max_module = sqrt(TF_UBorn[cpt_max].imag()*TF_UBorn[cpt_max].imag()+TF_UBorn[cpt_max].real()*TF_UBorn[cpt_max].real());

                int kxmi=cpt_max%(2*m1.NXMAX), kymi=cpt_max/(2*m1.NXMAX);
                posSpec= {kxmi,kymi}; ///coord informatique speculaire
                tabPosSpec[cpt_angle]=(double)posSpec.x;
                tabPosSpec[cpt_angle+NbAngle]=(double)posSpec.y;
                centre[kxmi*2*m1.NXMAX+kymi]=cpt_angle;



            if(m1.b_CorrAber==true){
                    calc_Uborn(TF_UBorn,UBorn,dim2DHA,posSpec,in_HA,out_HA,p_backward_HA);///--/!\ recale le spectre dans Uborn!
                    //SAVCplx(UBorn,"Re",chemin_result+"/ampli_UBorn_debut_extract.raw",t_float,"a+b");
                    ///----------Calcul phase + deroulement--------------------------------
                    calcPhase_mpi_pi_atan2(UBorn,phase_2Pi_vec); ///fonction atan2
                    //SAV2(phase_2Pi_vec,chemin_result+"/phasePI_atan2.raw",t_float,"a+b");
                    //phase2pi(UBorn, dim2DHA,phase2Pi);//asin

                    if(m1.b_Deroul==true){
                            if(m1.b_volkov==0)
                            {
                                phaseUnwrapping_Mat(dim2DHA, phase_2Pi_vec, UnwrappedPhase_herraez);
                                for(int cpt=0;cpt<NbPixUBorn;cpt++)
                                UnwrappedPhase[cpt]=UnwrappedPhase_herraez[cpt];///plutôt passer pointeur ?
                            }
                            else{

                                deroul_volkov(in_HA, out_HA,phase_2Pi_vec,UnwrappedPhase, p_forward_HA, p_backward_HA);
                            }
                        }
                    else UnwrappedPhase=phase_2Pi_vec;

                    //SAV2(UnwrappedPhase,chemin_result+"/phase_deroul_volkov.raw",t_float,"a+b");
                    ///-------------Correction aberration phase-------------------------------
                    src=Mat(dim2DHA.x,dim2DHA.y,CV_64F, UnwrappedPhase.data());
                    Mat Phase_corr(aberCorr(src, mask_aber,poly_aber,4,  NbPtOk));
                   // SAV2(poly_aber,NbPixUBorn,chemin_result+"/poly_aber.raw",t_float,"a+b");
                    ///---------------Correction amplitude----------------------------------------
                    for(int cpt=0; cpt<(NbPixUBorn); cpt++)
                    UBornAmp[cpt]=abs(UBorn[cpt]);
                    Mat srcAmp=Mat(dim2DHA.x, dim2DHA.y, CV_64F, UBornAmp.data());
                    //Mat UBornAmp_corr(dim2DHA.x, dim2DHA.y, CV_64F);
                    Mat UBornAmp_corr(ampliCorr(srcAmp, mask_aber,poly_aber,4,  NbPtOk));

                    ///Fin Correction amplitude----------------------------------------*/
                    // reconstruire l'onde complexe/Recalculate the complex field
                    for(size_t x=0; x<dim2DHA.x; x++){
                            for(size_t y=0; y<dim2DHA.y; y++){
                                    size_t cpt=x+y*dim2DHA.x;
                                   // size_t cpt3D=x+y*dim2DHA.x+NbPixUBorn*(cpt_angle-1);
                                    PhaseFinal[cpt]=Phase_corr.at<double>(y,x);//copie opencV->Tableau
                                    UBornAmpFinal[cpt]=UBornAmp_corr.at<double>(y,x);
                                    if(m1.b_Born==true){ //UBORN=U_tot-U_inc=u_tot_norm-1
                                            UBornFinal[cpt].real((UBornAmpFinal[cpt]-1)*cos(PhaseFinal[cpt]));///correction amplitude
                                            UBornFinal[cpt].imag((UBornAmpFinal[cpt]-1)*sin(PhaseFinal[cpt]));
                                        }
                                    else{ //RYTOV URytov = log a_t/a_i=log a_t
                                            UBornFinal[cpt].real(log(sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])));///racine(U^2) pour éliminer les éventueles amplitudes négatives
                                            UBornFinal[cpt].imag(PhaseFinal[cpt]);
                                        }
                                }
                        }
                    //SAV2(PhaseFinal,chemin_result+"/phase_finale.raw",t_float,"a+b");
                    //SAV2(UBornAmp,chemin_result+"/UbornAmp_final.bin",t_float,"a+b");
                    //SAV2(UBornAmpFinal,chemin_result+"/UBornAmpFinal.raw",t_float,"a+b");
                    ///Recalculer la TF décalée pour le programme principal.
                    Var2D recal= {kxmi,kymi};
                    decal2DCplxGen(UBornFinal,UBornFinalDecal, dim2DHA,decal2DHA);
                    TF2Dcplx_vec(in_HA,out_HA,UBornFinalDecal,TF_UBorn_norm,p_forward_HA);
                 //  SAVCplx(UBornFinal,"Re", chemin_result+"/UBornfinal_Re.raw",t_double, "a+b");

                    SAVCplx(UBornFinal,"Re", chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw", t_double, "a+b");
                  //  SAVCplx(UBornFinal,"Im", chemin_result+"/UBornfinal_Im.raw",t_double, "a+b");
                     SAVCplx(UBornFinal,"Im", chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw", t_double, "a+b");

                }
            else{ ///sauvegarde onde avec aberration

                    if(m1.b_Born==0)
                        cout<<"Attention : Rytov nécessite un déroulement"<<endl;
                    //cout<<"normalise spec"<<endl;
                    for(int cpt=0; cpt<(4*m1.NXMAX*m1.NXMAX); cpt++){//correction phase à l'ordre zéro et normalisatoin amplitude par ampli_spec=ampli_inc*ampli_ref
                           TF_UBorn_norm[cpt].real((TF_UBorn[cpt].real()*max_part_reel+TF_UBorn[cpt].imag()*max_part_imag)/max_module);
                        TF_UBorn_norm[cpt].imag((TF_UBorn[cpt].imag()*max_part_reel-TF_UBorn[cpt].real()*max_part_imag)/max_module);
                        }

                     //SAVCplx(TF_UBorn,"Re", chemin_result+"/TF_Uborn_Re.raw", t_double, "a+b");
                    calc_Uborn(TF_UBorn_norm,UBorn,dim2DHA,posSpec,in_HA,out_HA,p_backward_HA);///--/!\ recale le spectre dans support Uborn!

                    SAVCplx(UBorn,"Re", chemin_result+"/UBornfinal_Re"+m1.dimImg+".raw", t_double, "a+b");
                    SAVCplx(UBorn,"Im", chemin_result+"/UBornfinal_Im"+m1.dimImg+".raw", t_double, "a+b");
                }
        }//fin de boucle for sur tous les angles
    delete[] UnwrappedPhase_herraez;
    ///------------SAUVER LES PARAMETRES UTILES À RECONSTRUCTION ou au contrôle
    cout<<"NXMAX sauve="<<m1.NXMAX<<endl;
    vector<double> param{m1.NXMAX,NbAngle,m1.rayon,dimROI.x,m1.tailleTheoPixelHolo};
    SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");
    SAV2(tabPosSpec,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");
   // SAV2(centre,m1.chemin_result+"/centres.bin", t_int,"wb");
    SAV_Tiff2D(centre,m1.chemin_result+"/centres.tif",m1.NA/m1.NXMAX); //exportation des spéculaires en "unité NA"
//        SAV_Tiff2D2(centre, chemin_result+"/centres.tif", 2*NXMAX );
cout<<"Fin prétraitement"<<endl;
return 0;

}
