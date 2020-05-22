#include "manip.h"
#include <iostream>
#include "fonctions.h"
#include "projet.h"
#include <iomanip>
using namespace std;



manip::manip(unsigned short int dimROI)
{
    Var2D dim2D= {dimROI_Cam,dimROI};
    dimROI_Cam=dimROI;

    string home=getenv("HOME");
    //string sav_param=home+"/tomo_test/SAV_param_manip.txt";
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    string fic_cfg_recon=repertoire_config+"/recon.txt";
    cout<<"fichier recon="<<fic_cfg_recon<<endl;
    string fic_cfg_manip=repertoire_config+"/config_manip.txt";
    cout<<"chemin config="<<fic_cfg_manip<<endl;

    nbHolo=extract_val("NB_HOLO",fic_cfg_manip);
    cout<<"\n##################### INFO MANIP ##################\n"<<endl;
    n0=extract_val("N0",fic_cfg_manip);//indice de l'huile
    NACond=extract_val("NA_COND",fic_cfg_manip);
    NAObj=extract_val("NA_OBJ",fic_cfg_manip);	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
    lambda_v=extract_val("LAMBDA",fic_cfg_manip);
    f_tube=extract_val("F_TUBE",fic_cfg_manip);
    f_obj=extract_val("F_OBJ",fic_cfg_manip);
    Rf=extract_val("RF",fic_cfg_manip);
    Gt=f_tube/(f_obj*Rf); //longueur d'onde dans le vide, grandissement total microscope
    TpCam=extract_val("TPCAM",fic_cfg_manip);//taille pixel caméra
    Tp_holo=TpCam/Gt;
    Delta_f_holo=1/(dimROI_Cam*Tp_holo);//echantillonnage dans Fourier = constante car aucun bourrage dans l'espace direct
    theta_max=asin(NAObj/n0);
    //NXMAX=dimROI*Tp_Uborn*n0/(lambda_v);
    NXMAX=NAObj/(lambda_v*Delta_f_holo);//pas de Gt car inclu dans Delta_f via Tp_holo
    NXMAX_cond=NACond/(lambda_v*Delta_f_holo);
    dim_final=4*NXMAX;
    dim_Uborn=2*NXMAX;

    cout<<"Rayon Ewald calcule via lambda="<<n0/(lambda_v*Delta_f_holo)<<endl;

    //double tailleTheoPixelUborn=Tp_holo*dimROI/(2*NXMAX);
  //  double tailleTheoPixelTomo=tailleTheoPixelUborn*(2*NXMAX)/dim_final;

    R_EwaldPix=round(NXMAX/sin(theta_max));//R_Ewald pixels
    R_EwaldMet=round(NXMAX/sin(theta_max)*Delta_f_holo);//R_Ewald metric
    Tp_Uborn=Tp_holo*dimROI_Cam/(2*NXMAX);//pixel size of complex field
    Delta_f_Uborn=1/(2*NXMAX*Tp_Uborn);//sampling size of complex field spectrum

    cout<<"Rayon Ewald classique="<<dim_Uborn*Tp_Uborn*n0/(lambda_v)<<endl;
    Tp_Tomo=Tp_Uborn*(2*NXMAX)/dim_final;
    Delta_f_tomo=1/(Tp_Tomo*dim_final);
    dimImg=to_string(dim_Uborn)+"x"+to_string(dim_Uborn)+"x"+to_string(nbHolo);

    string sav_param=chemin_result+"/SAV_param_manip.txt";
    cout<<sav_param<<endl;
    ofstream fichier_sav_parametre(sav_param);

    cout<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    cout<<"+----------------+------------------+"<<endl;
    cout<<"|    Grandeur    |    Valeur        |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|      Gt        |     x"<<Gt<<"         |"<<endl;
    cout<<"|      λ_v       |     "<<lambda_v*pow(10,9)<<" nm       |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|Tp Caméra       |"<<TpCam*pow(10,6)<<"µm             |"<<endl;
    cout<<"|Tp_holo         |     "<<Tp_holo*pow(10,9)<<"nm    |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|champ holo pixel|     "<<dimROI_Cam          <<" pix      |"<<endl;
    cout<<"|Champ holo metre|     "<<dimROI_Cam*Tp_holo*pow(10,6)<<" µm   |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|    Tp UBorn    |     "<<Tp_Uborn*pow(10,9)<<" nm    |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|    Delta_f_holo|     "<<round(Delta_f_holo)<<" µm-1   |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|    Delta_f_tomo|     "<<round(Delta_f_tomo)<<" µm-1   |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|     R_Ewald    |     "<<R_EwaldPix<<" pixels    |"<<endl;
    cout<<"|     R_Ewald m  |    "<< setprecision (2) <<round(R_EwaldMet)<<" µm-1  |"<<setprecision (4)<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|     NXMAX      |     "<<round(NXMAX)<<" pix       |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"| chp cplx mini  |     "<<2*round(R_EwaldPix*NAObj/n0)<<" pix      |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|     Dim_finale |     "<<dim_final<<" pixels   |"<<endl;
    cout<<"|-----------------------------------|"<<endl;
    cout<<"|    Tp Tomo theo|     "<<Tp_Tomo*pow(10,9)<<" nm     |"<<endl;
    cout<<"+-----------------------------------+"<<endl;
    cout<<"\n#######################################\n"<<endl;

    fichier_sav_parametre<<"Fichiers d'acquisition : "<<chemin_acquis<<endl;
    fichier_sav_parametre<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    fichier_sav_parametre<<"+----------------+-------------------+"<<endl;
    fichier_sav_parametre<<"|    Grandeur    |    Valeur         |"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"| champ en pixel |     "<<dimROI_Cam         <<" pix       |"<<endl;
    fichier_sav_parametre<<"|    Tp Holo     |     "<<Tp_Uborn<<" nm |"<<endl;
    fichier_sav_parametre<<"|    Delta_f     |     "<<Delta_f_tomo<<" µm-1  |"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"|     Champ      |     "<<dimROI_Cam*Tp_Uborn<<" µm|"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"|     R_Ewald    |     "<<R_EwaldPix<<" pixels     |"<<endl;
    fichier_sav_parametre<<"|     R_Ewald m  |     "<<R_EwaldMet<<" µm-1|"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"|     NXMAX_theo |     "<<round(R_EwaldPix*NAObj/n0)<<" pixels     |"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"| Taille chp cplx|     "<<2*round(R_EwaldPix*NAObj/n0)<<" pixels    |"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Chp cplx |     "<<2*round(R_EwaldPix*NAObj/n0)*Tp_Uborn<<" nm|"<<endl;
    fichier_sav_parametre<<"|------------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Tomo theo|     "<<Tp_Tomo*pow(10,9)<<" nm    |"<<endl;
    fichier_sav_parametre<<"+------------------------------------+"<<endl;
    fichier_sav_parametre<<"\n##################### FIN INFO Reconstruction ##################\n"<<endl;

    fichier_sav_parametre.close();


}
