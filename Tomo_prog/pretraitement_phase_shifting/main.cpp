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
    m1.init();//initialiser les constantes de la manip
    string chemin_result=m1.chemin_result;
    string chemin_acquis=m1.chemin_acquis;
    string Chemin_mask=m1.chemin_acquis+"Image_mask.pgm";
    Var2D dimROI= {ROI_X,ROI_Y};
    size_t NbPixROI2d=dimROI.x*dimROI.y;
    vector<double> holo1(NbPixROI2d);
    char charAngle[4+1];
    char charNumImgPS[1+1];

    ///--------------Init FFTW Holo-------------------------------------------------
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        fftw_plan_with_nthreads(4);
        ///prepare fftw plan+tableaux-----------------
        fftw_plan p_forward_holo;
        fftw_complex *in_out_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//in=out pour transformation "inplace".
        fftw_complex *in_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//
        fftw_complex *out_holo=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixROI2d);//
        p_forward_holo=fftw_plan_dft_2d(dimROI.x, dimROI.y, in_holo, out_holo, FFTW_FORWARD,FFTW_ESTIMATE);
        vector<double> masqueTukeyHolo(NbPixROI2d);
        masqueTukeyHolo=tukey2D(dimROI.x,dimROI.y,0.05);
    ///--------------fin init FFTW Holo-------------------------------------------------

    ///------Init variable hologramme hors axe
    const int  NXMAX=m1.NXMAX,NbPixUBorn=4*NXMAX*NXMAX;
    Var2D dim2DUBorn= {m1.dim2DUBorn,m1.dim2DUBorn},coinROI= {0,0},decal2DUBorn= {NXMAX,NXMAX},posSpec= {0,0};

    const int NbAngle=m1.Num_Angle_final-m1.premier_plan;//m1.Num_Angle_final-m1.premier_plan;
    cout<<"Nb_angle="<<NbAngle<<endl;

    vector<complex<double>> TF_UBornTot(NbPixUBorn*NbAngle);///variable stockant les N champs complexes decoupés depuis la zone 1024 (pour utilsier wisdom en 1024)
    vector<double>  TF_champMod(NbPixUBorn);//module du champ
    vector<double> centre(NbPixUBorn);

    vector<double> holoPS1(dimROI.x*dimROI.y);
    vector<double> holoPS2(dimROI.x*dimROI.y);
    vector<double> holoPS3(dimROI.x*dimROI.y);
    vector<double> holoPS4(dimROI.x*dimROI.y);
    vector<double> modRef(NbPixROI2d);
    vector<double> ampRef(NbPixROI2d);

    vector<complex <double>> holoPSCplx(dimROI.x*dimROI.y);

    size_t  NbAngleOk=0;
    FILE* test_existence;//tester l'existence des fichiers
     auto start = std::chrono::system_clock::now();

   // int nbpix=1024*1024;
    Mat image;

    vector<double> imagePS(NbPixROI2d); //vector stockant une image de PS
    //vector<vector<double> > tab_image(4,imagePS);//tableau pour le stockage des 4 images de 1024x1024 pour calculer

    string cheminRef=m1.chemin_acquis+"modref.bmp";
    if(fopen(cheminRef.c_str(), "rb")!=NULL){

    charger_image2D_OCV(modRef,cheminRef, coinROI, dimROI);
    for(int cpt=0;cpt<NbPixROI2d;cpt++)
        if(modRef[cpt]!=0)
        ampRef[cpt]=sqrt(modRef[cpt]);
    else
        ampRef[cpt]=1;
    }
    else{
         for(int cpt=0;cpt<NbPixROI2d;cpt++)
         {
           modRef[cpt]=1;
         }

    }


    for(size_t cptAngle=0; cptAngle<NbAngle; cptAngle++){
            if((cptAngle-100*(cptAngle/100))==0)

            cout<<cptAngle<<endl;
            sprintf(charAngle,"%01i",cptAngle+1);//+1 car pas d'image i000.

            string StrAngle=charAngle;
            //cout<<"charangle="<<charAngle<<endl;

            for(int numImgPS=0;numImgPS<4;numImgPS++)
            {
                //sprintf(charNumImgPS,"%03i",numImgPS+1);
               // sprintf(charNumImgPS,"%03i",numImgPS+1);
                sprintf(charNumImgPS,"%03i",numImgPS+1);
               // string nomFichierHolo=m1.chemin_acquis+"i"+StrAngle+"-"+charNumImgPS+".pgm";
               string nomFichierHolo=m1.chemin_acquis+"session09032004-record"+StrAngle+"-"+charNumImgPS+".bmp";
                test_existence = fopen(nomFichierHolo.c_str(), "rb");
               // cout<<nomFichierHolo<<endl;


                if(test_existence!=NULL) {
                    fclose(test_existence);

                      if(numImgPS==0) charger_image2D_OCV(holoPS1,nomFichierHolo, coinROI, dimROI);
                      if(numImgPS==1)  charger_image2D_OCV(holoPS2,nomFichierHolo, coinROI, dimROI);
                      if(numImgPS==2)  charger_image2D_OCV(holoPS3,nomFichierHolo, coinROI, dimROI);
                      if(numImgPS==3)  charger_image2D_OCV(holoPS4,nomFichierHolo, coinROI, dimROI);
                }
                else cout<<"fichier "<<cptAngle<<" inexistant\n";
            }

          //  SAV_Tiff2D(holoPS1,m1.chemin_result+"/holops.tif",1);
            for(int cpt=0;cpt<NbPixROI2d;cpt++)
                    {
                        holoPSCplx[cpt].real((holoPS1[cpt]-holoPS3[cpt]));
                        holoPSCplx[cpt].imag((holoPS4[cpt]-holoPS2[cpt]));
                    }
                    SAVCplx(holoPSCplx,"Re", m1.chemin_result+"/holo_extract.bin",t_float,"a+b");
                   holo2TF_UBorn_PS(holoPSCplx,TF_UBornTot,NbAngleOk, masqueTukeyHolo, in_holo,out_holo,p_forward_holo,m1);
                    //holo2TF_UBorn(holo1,TF_UBornTot,dimROI,dim2DHA,coinHA,NbAngleOk, masqueTukeyHolo, in_holo,out_holo,p_forward_holo);
                    NbAngleOk++;
    }

//SAVCplx(TF_UBornTot,"Re",m1.chemin_result+"/TF_UBornTot.raw",t_float,"a+b");
    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    std::cout <<"Découpe hors axe= "<< elapsed.count()/(pow(10,9)) << '\n';
    ///--------------libérer allocation FFTW holo-------------------------------------
        fftw_destroy_plan(p_forward_holo);
        fftw_free(in_holo);
        fftw_free(out_holo);
        fftw_forget_wisdom();
    ///--------------Init FFTW Hors axe-------------------------------------------------
        fftw_plan p_backward_PS,p_forward_PS;
        //fftw_complex *in_out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimROI.x*dimROI.y);//
        fftw_complex *in_PS=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixUBorn);//in=out pour transformation "inplace".
        fftw_complex *out_PS=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPixUBorn);//
        p_backward_PS=fftw_plan_dft_2d(dim2DUBorn.x, dim2DUBorn.y, in_PS, out_PS, FFTW_BACKWARD,FFTW_ESTIMATE);
         p_forward_PS=fftw_plan_dft_2d(dim2DUBorn.x, dim2DUBorn.y, in_PS, out_PS, FFTW_FORWARD,FFTW_ESTIMATE);

        float alpha=0.1;//coeff pour le masque de tuckey
        vector<double> masqueTukeyHA(tukey2D(dim2DUBorn.x,dim2DUBorn.y,alpha));
     ///---------fin prepare fftw HA-----------------
     ///----------------init  vectors stockants les champs, spectres, phase etc.--------------------
    //vector<complex<double>> TF_UBorn(TF_UBornTot.begin(),TF_UBornTot.begin()+NbPixUBorn);
    vector<complex<double>> TF_UBorn(NbPixUBorn);
    vector<complex<double>> UBorn(NbPixUBorn);
    vector<double> phase_2Pi_vec(NbPixUBorn);
    vector<double> phase_deroul_vec(NbPixUBorn);
    vector<double> UnwrappedPhase(NbPixUBorn);
    vector<double> ampli_ref(NbPixUBorn);
    double *poly_aber=new double[NbPixUBorn];
    double *UnwrappedPhase_herraez=new double[NbPixUBorn];

    vector<double> PhaseFinal(NbPixUBorn);
    Mat src=Mat(1, ampli_ref.size(), CV_64F, ampli_ref.data());
    ///Charger masque aberration
    Mat mask_aber=init_mask_aber(Chemin_mask,dim2DUBorn);
    mask_aber.convertTo(mask_aber, CV_8U);
    int NbPtOk=countM(mask_aber);
    ///initialiser variable champ et spectre complexe
    vector<complex<double>> UBornFinal(NbPixUBorn);
    vector<double> UBornAmpFinal(NbPixUBorn);
    vector<double> UBornAmp(NbPixUBorn);
    vector<complex<double>> UBornFinalDecal(NbPixUBorn);
    vector<complex<double>> TF_UBorn_norm(NbPixUBorn);
    vector<double> tabPosSpec(NbAngle*2);  ///stockage des speculaires pour exportation vers reconstruction

    cout<<"\n#########################Calcul champs cplx 2D Uborn/Rytov + eventuelle Correction aberrations#############################"<<endl;
    cout<<"Nb angle OK="<<NbAngleOk<<endl;

    for(int cpt_angle=0; cpt_angle<NbAngleOk; cpt_angle++){ //boucle sur tous les angles : correction aberrations
            ///Récupérer la TF2D dans la pile de spectre2D. Place le spectre à cpt angle dans TF_UBorn
            TF_UBorn.assign(TF_UBornTot.begin()+NbPixUBorn*cpt_angle, TF_UBornTot.begin()+NbPixUBorn*(cpt_angle+1));
            //cout<<"cpt_angle="<<cpt_angle<<endl;
            SAVCplx(TF_UBorn,"Re",m1.chemin_result+"/TF_Uborn_iterateur.raw",t_float,"a+b");


            ///Recherche de la valeur maximum du module dans ref non centré-----------------------------------------
                int cpt_max=coordSpec(TF_UBorn, TF_champMod,decal2DUBorn);
                double  max_part_reel = TF_UBorn[cpt_max].real(),///sauvegarde de la valeur cplx des  spéculaires
                        max_part_imag = TF_UBorn[cpt_max].imag(),
                            //max_module = TF_champMod[cpt_max];
                        max_module = sqrt(TF_UBorn[cpt_max].imag()*TF_UBorn[cpt_max].imag()+TF_UBorn[cpt_max].real()*TF_UBorn[cpt_max].real());

                int kxmi=cpt_max%(2*NXMAX), kymi=cpt_max/(2*NXMAX);
                posSpec= {kxmi,kymi}; ///coord informatique speculaire
                tabPosSpec[cpt_angle]=(double)posSpec.x;
                tabPosSpec[cpt_angle+NbAngle]=(double)posSpec.y;
                centre[kxmi*2*NXMAX+kymi]=cpt_angle;


            if(m1.b_CorrAber==true){
                    calc_Uborn(TF_UBorn,UBorn,dim2DUBorn,posSpec,in_PS,out_PS,p_backward_PS);///--/!\ recale le spectre dans Uborn!
                    //SAVCplx(UBorn,"Re",chemin_result+"/ampli_UBorn_debut_extract.raw",t_float,"a+b");
                    ///----------Calcul phase + deroulement--------------------------------
                    calcPhase_mpi_pi_atan2(UBorn,phase_2Pi_vec); ///fonction atan2
                    SAV2(phase_2Pi_vec,chemin_result+"/phasePI_atan2.raw",t_float,"a+b");
                    //phase2pi(UBorn, dim2DHA,phase2Pi);//asin

                    if(m1.b_Deroul==true){
                            phaseUnwrapping_Mat(dim2DUBorn, phase_2Pi_vec, UnwrappedPhase_herraez);

                            for(int cpt=0;cpt<NbPixUBorn;cpt++)
                                UnwrappedPhase[cpt]=UnwrappedPhase_herraez[cpt];
                           // deroul_volkov(in_HA, out_HA,phase_2Pi_vec,UnwrappedPhase, p_forward_HA, p_backward_HA);
                        }
                    else UnwrappedPhase=phase_2Pi_vec;

                    //SAV2(UnwrappedPhase,chemin_result+"/phase_deroul_volkov.raw",t_float,"a+b");
                    ///-------------Correction aberration phase-------------------------------
                    src=Mat(dim2DUBorn.x,dim2DUBorn.y,CV_64F, UnwrappedPhase.data());
                    Mat Phase_corr(aberCorr(src, mask_aber,poly_aber,4,  NbPtOk));
                    SAV2(poly_aber,NbPixUBorn,chemin_result+"/poly_aber.raw",t_float,"a+b");
                    ///---------------Correction amplitude----------------------------------------
                    for(int cpt=0; cpt<(NbPixUBorn); cpt++)
                    UBornAmp[cpt]=abs(UBorn[cpt]);
                    Mat srcAmp=Mat(dim2DUBorn.x, dim2DUBorn.y, CV_64F, UBornAmp.data());
                    //Mat UBornAmp_corr(dim2DHA.x, dim2DHA.y, CV_64F);
                    Mat UBornAmp_corr(ampliCorr(srcAmp, mask_aber,poly_aber,4,  NbPtOk));

                    ///Fin Correction amplitude----------------------------------------*/
                    // reconstruire l'onde complexe/Recalculate the complex field
                    for(size_t x=0; x<dim2DUBorn.x; x++){
                            for(size_t y=0; y<dim2DUBorn.y; y++){
                                    size_t cpt=x+y*dim2DUBorn.x;
                                    size_t cpt3D=x+y*dim2DUBorn.x+NbPixUBorn*(cpt_angle-1);
                                    PhaseFinal[cpt]=Phase_corr.at<double>(y,x);//copie opencV->Tableau
                                    UBornAmpFinal[cpt]=UBornAmp_corr.at<double>(y,x);
                                    if(m1.b_Born==true){ //UBORN=U_tot-U_inc=u_tot_norm-1
                                            UBornFinal[cpt].real((UBornAmpFinal[cpt]-1)*cos(PhaseFinal[cpt]));///correction amplitude
                                            UBornFinal[cpt].imag((UBornAmpFinal[cpt]-1)*sin(PhaseFinal[cpt]));
                                        }
                                    else{ //RYTOV URytov = log a_t/a_i=log a_t
                                            UBornFinal[cpt].real(log(sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt])));///Rytov nécessite correction amplitude
                                            UBornFinal[cpt].imag(PhaseFinal[cpt]);
                                        }
                                }
                        }

                    SAV2(PhaseFinal,chemin_result+"/phase_finale.raw",t_float,"a+b");
                    //SAV2(UBornAmp,chemin_result+"/UbornAmp_final.bin",t_float,"a+b");
                    SAV2(UBornAmpFinal,chemin_result+"/UBornAmpFinal.raw",t_float,"a+b");
                    ///Recalculer la TF décalée pour le programme principal.
                    Var2D recal= {kxmi,kymi};
                    decal2DCplxGen(UBornFinal,UBornFinalDecal, dim2DUBorn,decal2DUBorn);
                    TF2Dcplx_vec(in_PS,out_PS,UBornFinalDecal,TF_UBorn_norm,p_forward_PS);
                    SAVCplx(UBornFinal,"Re", chemin_result+"/UBornfinal_Re.raw",t_double, "a+b");
                    SAVCplx(UBornFinal,"Im", chemin_result+"/UBornfinal_Im.raw",t_double, "a+b");

                }
            else{ ///sauvegarde onde avec aberration

                    if(m1.b_Born==0)
                        cout<<"Attention : Rytov nécessite un déroulement"<<endl;
                    //cout<<"normalise spec"<<endl;
                    for(int cpt=0; cpt<(4*NXMAX*NXMAX); cpt++){//correction phase à l'ordre zéro et normalisatoin amplitude par ampli_spec=ampli_inc*ampli_ref
                           TF_UBorn_norm[cpt].real((TF_UBorn[cpt].real()*max_part_reel+TF_UBorn[cpt].imag()*max_part_imag)/max_module);
                        TF_UBorn_norm[cpt].imag((TF_UBorn[cpt].imag()*max_part_reel-TF_UBorn[cpt].real()*max_part_imag)/max_module);
                        }

                     //SAVCplx(TF_UBorn,"Re", chemin_result+"/TF_Uborn_Re.raw", t_double, "a+b");
                    calc_Uborn(TF_UBorn_norm,UBorn,dim2DUBorn,posSpec,in_PS,out_PS,p_backward_PS);///--/!\ recale le spectre dans support Uborn!

                    SAVCplx(UBorn,"Re", chemin_result+"/UBornfinal_Re.raw", t_double, "a+b");
                    SAVCplx(UBorn,"Im", chemin_result+"/UBornfinal_Im.raw", t_double, "a+b");
                }
        }//fin de boucle for sur tous les angles

    delete[] UnwrappedPhase_herraez;
    ///------------SAUVER LES PARAMETRES UTILES À RECONSTRUCTION ou au contrôle
    vector<double> param{NXMAX,NbAngle,m1.rayon,dimROI.x,m1.tailleTheoPixelHolo};
    SAV2(param,m1.chemin_result+"/parametres.raw", t_double, "wb");
    SAV2(tabPosSpec,m1.chemin_result+"/tab_posSpec.raw",t_double,"wb");
    SAV_Tiff2D(centre,m1.chemin_result+"/centres.tif",m1.NA/(m1.NXMAX*100));//exportation centre en NA
//        SAV_Tiff2D2(centre, chemin_result+"/centres.tif", 2*NXMAX );
cout<<"Fin prétraitement"<<endl;
return 0;

}
