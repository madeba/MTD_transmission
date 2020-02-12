#include "manip.h"
#include <iostream>
#include "fonctions.h"
#include "projet.h"
#include <iomanip>
using namespace std;



manip::manip(int dimROI)
{
    Var2D dim2D= {dimROI,dimROI};
    nbHolo=200;
    cout<<"\n##################### INFO MANIP ##################\n"<<endl;
    n0=1.515,//indice de l'huile
    NACond=1.4, NAObj=1.4,	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
    lambda_v=0.532*pow(10,-6), Gt=120; //longueur d'onde dans le vide, grandissement total microscope
    TpCam=5.5*pow(10,-6);//taille pixel caméra
    Tp_holo=TpCam/Gt;
    Delta_f_holo=1/(dimROI*Tp_holo);//echantillonnage dans Fourier = constante car aucun bourrage dans l'espace direct
    theta_max=asin(NAObj/n0);
    //NXMAX=dimROI*Tp_Uborn*n0/(lambda_v);
    NXMAX=NAObj/(lambda_v*Delta_f_holo);//pas de Gt car inclu dans Delta_f via Tp_holo
    NXMAX_cond=NACond/(lambda_v*Delta_f_holo);
    dim_final=4*NXMAX;
    dim_Uborn=2*NXMAX;

    cout<<"Rayon Ewald calcule via lambda="<<n0/(lambda_v*Delta_f_holo)<<endl;

    //double tailleTheoPixelUborn=Tp_holo*dimROI/(2*NXMAX);
  //  double tailleTheoPixelTomo=tailleTheoPixelUborn*(2*NXMAX)/dim_final;

    R_EwaldPix=round(NXMAX/sin(theta_max));
    R_EwaldMet=round(NXMAX/sin(theta_max)*Delta_f_holo);
    Tp_Uborn=Tp_holo*dimROI/(2*NXMAX);
    Delta_f_Uborn=1/(2*NXMAX*Tp_Uborn);

    cout<<"Rayon Ewald classique="<<dim_Uborn*Tp_Uborn*n0/(lambda_v)<<endl;
    Tp_Tomo=Tp_Uborn*(2*NXMAX)/dim_final;
    Delta_f_tomo=1/(Tp_Tomo*dim_final);
    string home=getenv("HOME");
    string sav_param=home+"/tomo_test/SAV_param_manip.txt";
    cout<<sav_param<<endl;
    ofstream fichier_sav_parametre(sav_param);

    cout<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    cout<<"+----------------+----------------+"<<endl;
    cout<<"|    Grandeur    |    Valeur      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|      Gt        |     x"<<Gt<<"       |"<<endl;
    cout<<"|      λ_v       |     "<<lambda_v*pow(10,9)<<" nm     |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|Tp Caméra       |"<<TpCam*pow(10,6)<<"µm           |"<<endl;
    cout<<"|Tp_holo         |     "<<Tp_holo*pow(10,9)<<"nm  |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|champ holo pixel|     "<<dimROI          <<" pix    |"<<endl;
    cout<<"|Champ holo metre|     "<<dimROI*Tp_holo*pow(10,6)<<" µm |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp UBorn    |     "<<Tp_Uborn*pow(10,9)<<" nm  |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Delta_f_holo|     "<<round(Delta_f_holo)<<" µm-1 |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Delta_f_tomo|     "<<round(Delta_f_tomo)<<" µm-1 |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     R_Ewald    |     "<<R_EwaldPix<<" pixels  |"<<endl;
    cout<<"|     R_Ewald m  |    "<< setprecision (2) <<round(R_EwaldMet)<<" µm-1|"<<setprecision (4)<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     NXMAX      |     "<<round(NXMAX)<<" pix     |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"| chp cplx mini  |     "<<2*round(R_EwaldPix*NAObj/n0)<<" pix    |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     Dim_finale |     "<<dim_final<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp Tomo theo|     "<<Tp_Tomo*pow(10,9)<<" nm   |"<<endl;
    cout<<"+---------------------------------+"<<endl;
    cout<<"\n#######################################\n"<<endl;

    //fichier_sav_parametre<<"Fichiers d'acquisition : "<<chemin_acquis<<endl;
  /*  fichier_sav_parametre<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    fichier_sav_parametre<<"+----------------+----------------+"<<endl;
    fichier_sav_parametre<<"|    Grandeur    |    Valeur      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"| champ en pixel |     "<<dimROI          <<" pix     |"<<endl;
    fichier_sav_parametre<<"|    Tp Holo     |     "<<Tp_Uborn<<" nm            |"<<endl;
    fichier_sav_parametre<<"|    Delta_f     |     "<<Delta_f<<" µm-1      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     Champ      |     "<<dimROI*Tp_Uborn<<" µm      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     R_Ewald    |     "<<R_EwaldPix<<" pixels |"<<endl;
    fichier_sav_parametre<<"|     R_Ewald m  |     "<<R_EwaldMet<<" µm-1    |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     NXMAX_theo |     "<<round(R_EwaldPix*NAObj/n0)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"| Taille chp cplx|     "<<2*round(R_EwaldPix*NAObj/n0)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Chp cplx |     "<<2*round(R_EwaldPix*NAObj/n0)*Tp_Uborn<<" nm     |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Tomo theo|     "<<(2*NXMAX)/dim_final*Tp_Tomo<<" nm      |"<<endl;
    fichier_sav_parametre<<"+---------------------------------+"<<endl;
    fichier_sav_parametre<<"\n##################### FIN INFO Reconstruction ##################\n"<<endl;

    fichier_sav_parametre.close();*/


}
