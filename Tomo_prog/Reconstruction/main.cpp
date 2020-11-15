#include "projet.h"
#include "fonctions.h"
//#include "deroulement.h"
//#include "Correction_aberration.h"
//#include "CImg.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <octave-2.9.9/octave/oct.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <chrono>
//#include <Magick++.h>
#include <fftw3.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <omp.h>
#include "manip.h"

#include "FFT_fonctions.h"
#define PI M_PI
//using namespace cimg_library;
using namespace std;
using namespace chrono;
//using namespace Magick;

/* -------------------------------------------------------------------------- */
// Usage
/* -------------------------------------------------------------------------- */

/*static void
usage(int argc, char **argv)
{
        if ((argc - 1) == 0) {
                printf("Programme de reconstruction en tomographie hors-axe\n");
                printf("Usage: %s <paramètres obligatoires> <paramètres optionnels>\n", argv[0]);
                printf("Paramètres obligatoires: \n");
                printf("-i <répertoire>: répertoire des images acquises à traiter \n");
                printf("-c <cx> <cy>: centre du cercle sur l'image 1024x1024 (orig°: top-left corner) \n");

                exit(EXIT_FAILURE);
        }

}*/
/* --------------------------------------------------------------------------- */
// Main
/* --------------------------------------------------------------------------- */

//int main(int argc, char *argv[])
int main()
{
    manip m1;
    int NXMAX=m1.NXMAX, NbAngle=m1.NbAngle;
    Var2D NMAX= {NXMAX,NXMAX},dimChpCplx= {2*NXMAX,2*NXMAX};
    Point2D dim2DUBorn(2*NXMAX,2*NXMAX,2*NXMAX);
    string dimImg=to_string(dimChpCplx.x)+"x"+to_string(dimChpCplx.y)+"x"+to_string(NbAngle);

    ///------Récupérer spéculaire dans le fichier binaire---
    //lire_bin(chemin_fichier,variable de stockage, type de donnée en bits (32,64...), nombre d'images 2D)
    double *TabPosSpec=new double[2*m1.NbAngle];
    vector<double> centres(dimChpCplx.x*dimChpCplx.y);
    cout<<"chemin tomo test ="<<m1.chemin_result+"/tab_posSpec.raw"<<endl;
    cout<<"lecture posSpec"<<endl;
    lire_bin(m1.chemin_result+"/tab_posSpec.raw",TabPosSpec,64,2*m1.NbAngle);

    for(size_t holo_numero=0; holo_numero<m1.NbAngle; holo_numero++){
        int cpt2D=round(TabPosSpec[holo_numero+m1.NbAngle])*dimChpCplx.x+round(TabPosSpec[holo_numero]);
        centres[cpt2D]=holo_numero;//save centres in a image file, for quick visualisation
    }
    SAV_Tiff2D(centres,m1.chemin_result+"centre_reconstruction.tif",1);

    ///--------------FFTW2D init----
    FFT_encaps  tf2D(dim2DUBorn);
    ///---Variable de dimensions (taille pixel, volume, zoom)-------------------

    const  int dim_final=m1.dim_final;//peut être différent de 4*NXMAX, mais l'image final sera (dé)zoomée;
    cout<<"Dimension forcée à "<<dim_final<<endl;
    float tailleTheoPixelTomo=m1.tailleTheoPixelTomo;
    double zoom=double(dim_final)/double(4*NXMAX);//double(dim_final)/double(4*n0*TpCam*dimROI.x/(Gt*lambda)*NA/n0);//dim_final/4NXMAX
    if(zoom!=1)
    {
        cout<<"dim_final forcée à "<<dim_final<<" au lieu de "<<4*NXMAX<<"-->facteur de zoom numérique="<<zoom<<endl;
    }
    ///---------------Calcul de quelques constantes, afin d'éviter de les calculer dans la boucle--------------
    Var3D  decal3DTF= {(int)round(dim_final/2),(int)round(dim_final/2),(int)round(dim_final/2)};
    fflush(stdout);
    cout.flush();
    ///---------------Recharger le champ complexe depuis sav prétraitement-----------------------------------------------------
    Var2D dim2DHA= {2*NXMAX,2*NXMAX};//decal2D= {NXMAX,NXMAX};//dim2DHA=dimUBorn
    const int NbPixU_Born=dim2DHA.x*dim2DHA.y;

    double *UBornFinal3D_Re=new double[NbPixU_Born*NbAngle];
    double *UBornFinal3D_Im=new double[NbPixU_Born*NbAngle];
    vector<double> mask_tukey2D(NbPixU_Born);
    vector<complex<double>> UBornFinal3D(NbPixU_Born*NbAngle);
    cout<<"lecture Uborn Re et Im"<<endl;

    lire_bin(m1.chemin_result+"/UBornfinal_Im"+dimImg+".raw",UBornFinal3D_Im,64,NbPixU_Born*NbAngle);
    lire_bin(m1.chemin_result+"/UBornfinal_Re"+dimImg+".raw",UBornFinal3D_Re,64,NbPixU_Born*NbAngle);
    #pragma omp parallel for
    for(int cpt=0; cpt<NbPixU_Born*NbAngle; cpt++){
        UBornFinal3D[cpt].real(UBornFinal3D_Re[cpt]);
        UBornFinal3D[cpt].imag(UBornFinal3D_Im[cpt]);
    }

    // SAV3D_Tiff(UBornFinal3D,dim_stack_holo,"Im",m1.chemin_result+"/UBornfinal3D_Im.tif",tailleTheoPixelUborn*1000);
    delete[] UBornFinal3D_Re;
    delete[] UBornFinal3D_Im;


    vector<complex <double>> UBornFinal2D(NbPixU_Born),  UBorn2DFinalDecal(NbPixU_Born), TF_UBorn_normC(NbPixU_Born), TF_UBorn_normI(NbPixU_Born);

    ///---------Réservation variables espace3D final --------------------------------------------------
    const size_t N_tab=dim_final*dim_final*dim_final;//Nombre de pixel dans l'esapce final (indice max du tableau 3D+1)
    vector <complex <double>> TF3D_PotObj(N_tab);
    vector<double> sup_redon(N_tab);
    ///-----Réservation chrono----------------
    clock_t
    temps_depart  = clock(), /*temps de départ de programme */
    temps_initial = clock (), /* temps initial en micro-secondes */
    temps_final=0,/* temps d'arrivée du programme */
    temps_arrivee;   /* temps final en micro-secondes */
    float temps_cpu=0,     /* temps total en secondes */
          temps_total=0;


    float alpha=0.1;//coeff pour le masque de tuckey
    mask_tukey2D=tukey2D(dim2DHA.x,dim2DHA.y,alpha); //  fenetre de Tukey,

    printf("*******************************************\n");
    printf("remplissage de l'espace réciproque\n");
    printf("*******************************************\n");
    printf("  \\,`//\n");
    printf(" _).. `_\n");
    printf("( __  -\\ \n");
    printf("    '`.\n");
    printf("   ( \\>\n");
    printf("   _||_ \n");
    const int premier_plan=0;
    Var2D posSpec= {0,0}, recal= {0,0};

    auto start_part1 = std::chrono::system_clock::now();
    cout<<"Reconstruction..."<<endl;

    for(int cpt_angle=premier_plan; cpt_angle<NbAngle; cpt_angle++) //boucle sur tous les angles    {
    {
        posSpec= {(int)TabPosSpec[cpt_angle],(int)TabPosSpec[NbAngle+cpt_angle]};  //récupérer spéculaire puis champ cplx depuis sauvegarde prétraitement
        for(int cpt=0; cpt<NbPixU_Born; cpt++)//retrieve complex fields in the stack
        {
            UBornFinal2D[cpt].real(UBornFinal3D[cpt+cpt_angle*NbPixU_Born].real()*mask_tukey2D[cpt]);
            UBornFinal2D[cpt].imag(UBornFinal3D[cpt+cpt_angle*NbPixU_Born].imag()*mask_tukey2D[cpt]);
        }

        decal2DCplxGen(UBornFinal2D,UBorn2DFinalDecal, dim2DHA,NMAX);
        TF2Dcplx(UBorn2DFinalDecal,TF_UBorn_normI, tf2D, m1.tailleTheoPixelUborn);

        recal= {posSpec.x,posSpec.y};
        decal2DCplxGen(TF_UBorn_normI,TF_UBorn_normC,dim2DHA,recal);  ///recaler le spectre à la position du spéculaire (la variable est attendue ainsi par la fonction retropropag)
        //SAVCplx(TF_UBorn_normC,"Re",m1.chemin_result+"TF_UBorn_250x250_normC.raw",t_float,"a+b");
        if((cpt_angle-100*(cpt_angle/100))==0)
            printf("cpt_angle=%i\n",cpt_angle);

        retroPropag_Born(TF3D_PotObj, TF_UBorn_normC, sup_redon, dim_final, posSpec, decal3DTF, NMAX, m1.rayon, m1);  ///--Mapping 3D=retropropagation
    }//fin de boucle for sur tous les angles

    auto end_part1 = std::chrono::system_clock::now();
    auto elapsed_part1 = end_part1 - start_part1;
    std::cout <<"Temps total fft+retroprop= "<< elapsed_part1.count()/(pow(10,9)) << '\n';
    delete[] TabPosSpec;
    printf("temps_total proj : %f \n",temps_total);

    ///---------------------mesure du temps de calcul--------------------------

    temps_final = clock ();
    temps_cpu = (temps_final - temps_initial) * 1e-6/NbAngle;
    printf("temps moyen pour 1 angle: %f\n",temps_cpu);
    temps_cpu = (temps_final - temps_initial) * 1e-6;
    printf("temps total pour %i angle(s): %f\n", NbAngle, temps_cpu);
    temps_initial=clock();//enclenchement du chronometre

    ///------------------moyennage par sup_redon------------------------------------------------
    for(size_t cpt=0; cpt<N_tab; cpt++)
    {
        if (sup_redon[cpt]==0) {  ////////////////remplace les 0 de sup_redon par des 1---> evite la division par 0
            sup_redon[cpt]=1;
        }
        TF3D_PotObj[cpt]= TF3D_PotObj[cpt]/sup_redon[cpt];//moyennage par sup_redon
    }

    vector<double>().swap(sup_redon);//forcer la libération mémoire de sup_redon
    //interp_lin3D(TF3D_PotObj);
    temps_final = clock ();
    temps_cpu = (temps_final - temps_initial) * 1e-6;
    printf("temps apres normalisation : %lf\n",temps_cpu);
    //SAV(TF3D_PotObj_Re, N_tab, "/home/aziz/Projet_tomo/Tomo_Images/TF2d_apres_masquage/TF3D_PotObj_Re_norm_avant.bin", float,"wb");

    //////////////////////////////papillon binarisé
    vector<double> papillon_masque(N_tab,0);
    for(size_t compteur=0; compteur<N_tab; compteur++)
    {
        if(TF3D_PotObj[compteur].real()!=0)
            papillon_masque[compteur]=1;
    }

    if(m1.b_Export_OTF==1)
        SAV3D_Tiff(papillon_masque,m1.chemin_result+"/OTF3D.tif",tailleTheoPixelTomo);

    printf("*******************************************\n");
    printf("circshift avant TF3D \n");
    printf("*******************************************\n");

    //////////// tableaux qui recoivent le circshift
    // nbCplx *TF3D_PotObj_shift=new nbCplx[N_tab];
    vector <complex<double >> TF3D_PotObj_shift(N_tab);

    temps_initial = clock();//enclenchement chronometre

    fftshift3D(TF3D_PotObj,TF3D_PotObj_shift);
    vector<complex<double>>().swap(TF3D_PotObj);///libérer mémoire
    temps_final = clock ();
    temps_cpu = (temps_final - temps_initial) * 1e-6;
    printf("temps apres fftshift 3D: %f\n",temps_cpu);

    ///------------------------ TF3D-------------------------------------

    vector <complex<double >> PotObj_shift(N_tab);
    Point3D dim3D(dim_final,dim_final,dim_final,dim_final);


    string chemin_wisdom=m1.chemin_config_defaut+"Wisdom/wisdom3d_512_c2c_dbl_back_EXHAUSTIVE_inplace_6Thrd_Ryzen3600x.txt";


    cout<<"dim3D="<<dim3D.x<<endl;
    high_resolution_clock::time_point t1v = high_resolution_clock::now();//temps vrai (pas CPU)

    printf("*******************************************\n");
    printf("TF3D...\n");
    printf("*******************************************\n");
    printf("2s calcul de la TF3D\n");
    printf("     .\"\".    .\"\".\n");
    printf("     |  |   /  /\n");
    printf("     |  |  /  /\n");
    printf("     |  | /  /\n");
    printf("     |  |/  ;-.__ \n");
    printf("     |  ` _/  /  /\n");
    printf("     |  /` ) /  /\n");
    printf("     | /  /_/\\_/ \n");
    printf("     |/  /      |\n");
    printf("     (  ' \\ '-  |\n");
    printf("      \\    `.  /\n");
    printf("       |      |\n");
    printf("       |      |\n");

    FFT_encaps tf3D_out(dim3D,m1.nbThreads,0);
    temps_initial = clock();//enclenchement chronometre

    //  FFT_encaps tf3D_out(dim3D,m1.nbThreads,0);
    TF3Dcplx_INV(tf3D_out.in, tf3D_out.out,TF3D_PotObj_shift, PotObj_shift, tf3D_out.p_backward_OUT,m1.Delta_fUborn);

    // FFT_encaps tf3D_IN(dim3D,m1.nbThreads,1);
    //TF3Dcplx_Inplace_INV(TF3D_PotObj_shift , PotObj_shift,tf3D_IN,m1.Delta_fUborn);

    vector<complex<double>>().swap(TF3D_PotObj_shift);
    high_resolution_clock::time_point t2v = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2v - t1v ).count();
    cout << "temps TF3D="<< duration/1000<<" ms"<<endl;

    ///--------------------------------circshift apres TF3D
    /////////////////////creation des tableaux qui recoivent le circshift

    vector <complex<double >> PotObj3D(N_tab);
    // circshift3DCplx(PotObj_shift, PotObj3D, dimFinal, decal3DTF);
    fftshift3D(PotObj_shift,PotObj3D);
    vector<complex<double>>().swap(PotObj_shift);
    ///--------------------------------chrono et écriture.
    temps_final = clock ();
    temps_cpu = (temps_final - temps_initial) * 1e-6;
    printf("temps apres circshift 3D final et calcul du module: %f\n",temps_cpu);
    temps_initial = clock();//enclenchement chronometre

    vector <complex<double >> indice_cplx(N_tab);
    double k_v=2*PI/m1.lambda0;
    complex<double>  ctePot2Ind(-1/(2*m1.n0*k_v*k_v),0);//toujours correct (et identique dans tous les papiers)

    for(size_t cpt=0; cpt<N_tab; cpt++)
    {
        indice_cplx[cpt]=ctePot2Ind*PotObj3D[cpt];
    }
    printf("ecriture de PotObj3D.Re, PotObj3D.Im \n");
    auto start_tiff = std::chrono::system_clock::now();

    SAV3D_Tiff(indice_cplx,"Re", m1.chemin_result+"/indice.tif",tailleTheoPixelTomo);
    SAV3D_Tiff(indice_cplx,"Im", m1.chemin_result+"/absorption.tif",tailleTheoPixelTomo);

    auto end_tiff = std::chrono::system_clock::now();
    auto elapsed_tiff = end_tiff - start_tiff;
    std::cout <<"Temps total ecriture tiff= "<< elapsed_tiff.count()/(pow(10,9)) << '\n';

    temps_arrivee = clock ();
    temps_cpu = (temps_arrivee-temps_depart )/CLOCKS_PER_SEC;
    printf("temps total: %f\n",temps_cpu);
    cout<<"#######-- Reconstruction terminée --#######"<<endl;
    return 0;
}

