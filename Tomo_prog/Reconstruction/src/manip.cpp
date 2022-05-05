#include "manip.h"
#include <iostream>
#include "fonctions.h"
#include "projet.h"

using namespace std;


manip::manip()
{
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    cout<<"chemin resultat="<<chemin_result<<endl;
    string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
    chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);
    string fic_cfg_recon=repertoire_config+"/recon.txt";
    cout<<"fichier recon="<<fic_cfg_recon<<endl;
    string fic_cfg_manip=repertoire_config+"/config_manip.txt";
    cout<<"chemin config="<<fic_cfg_manip<<endl;
    chemin_config=fic_cfg_manip;
    chemin_config_defaut=extract_string("CHEMIN_CONFIG_PC_ACQUIS",home+fin_chemin_gui_tomo)+"/";


        cout<<"##################"<<endl;
        cout<<"# RECONSTRUCTION #"<<endl;
         cout<<"##################"<<endl;

    ///-----Recharger les paramètres--------///
    const int nbParam=5;
    double parametre[nbParam];
    lire_bin(chemin_result+"/parametres.raw",parametre,64,nbParam);
    NXMAX=parametre[0],NYMAX=parametre[0],NbAngle=parametre[1],rayon=parametre[2],dimROI=parametre[3],tailleTheoPixelHolo=parametre[4];

    cout<<"NbAngle="<<NbAngle<<endl;
    cout<<"NXMAX="<<NXMAX<<endl;
    cout<<"rayon param=="<<rayon<<endl;
    cout<<"dimROIparam="<<dimROI<<endl;
    cout<<"\n##################### INFO MANIP ##################\n"<<endl;
    n0=extract_val("N0",fic_cfg_manip);	//indice de l'huile
    NA=extract_val("NA",fic_cfg_manip);	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
    lambda0=extract_val("LAMBDA",fic_cfg_manip);
    f_tube=extract_val("F_TUBE",fic_cfg_manip), ///focale lentille tube
    f_obj=extract_val("F_OBJ",fic_cfg_manip),///focale objectif
    G=f_tube/f_obj,	//grossissement telan+objectif
    TpCam=extract_val("TPCAM",fic_cfg_manip),//cam. photon focus
    //NXMAX_OBJ=extract_val("NXMAX_OBJ",fic_cfg_manip);
    Rf=extract_val("RF",fic_cfg_manip),//1.2;=1/facteur grossissement
    Gt=G/Rf;
    cout<<"Gt="<<Gt<<endl;
    tailleTheoPixelHolo=TpCam/Gt;//Pour info, taille des pixels sur un hologramme=Tpcam/GT
    cout<<"taille theorique pixel holo="<<tailleTheoPixelHolo*pow(10,9)<<"nm"<<endl;
    theta=asin(NA/n0);
 //   circle_cx=extract_val("CIRCLE_CX",fic_cfg_manip);
   // circle_cy=extract_val("CIRCLE_CY",fic_cfg_manip);
  //  NXMAX=extract_val("NXMAX",fic_cfg_manip);

    if(NXMAX==0)
        cout<<"----- ERREUR : NXMAX absent du fichier de configuration" <<fic_cfg_manip<<endl;
    dim_final=extract_val("DIM_FINAL",fic_cfg_recon);
    if(dim_final<4*NXMAX) cout<<"Attention, dim_final<4*NXMAX"<<endl;
    tailleTheoPixelUborn=tailleTheoPixelHolo*dimROI/(2*NXMAX);
    Delta_fUborn=1/(tailleTheoPixelUborn*2*NXMAX);//echantillonnage fréquentiel Uborn, utile pour repasser le rayon en métrique dans le volume 3D final
    tailleTheoPixelTomo=tailleTheoPixelUborn*(2*NXMAX)/dim_final;

    premier_plan=extract_val("PREMIER_ANGLE",fic_cfg_recon),
    Num_Angle_final=extract_val("NB_HOLO",fic_cfg_manip),//
    nbThreads=extract_val("NB_THREADS",fic_cfg_recon);
    cout<<"\n##################### Options de RECONSTRUCTION ##################\n"<<endl;
    b_CorrAber=extract_val("C_ABER",fic_cfg_recon);///corriger les aberrations?
    b_Deroul=extract_val("DEROUL",fic_cfg_recon);///Dérouler la phase?
    b_Born=extract_val("BORN",fic_cfg_recon);///Born vrai ? Sinon Rytov
    if(b_Born==true)
        cout<<"RYTOV=0"<<endl;
    else
        cout<<"RYTOV=1"<<endl;

    b_Export_OTF=extract_val("EXPORT_OTF",fic_cfg_recon);///Exporter OTF ?
    cout<<"\n##########################################################"<<endl;

      rayon=round(NXMAX*n0/NA);//calcul du rayon à partir de la fréquence NXMAX defini pare l'utlisateur
      double R_Ewald=IMAGE_DIMX*tailleTheoPixelHolo*n0/(lambda0); //vraie valeur de R_Ewald.
      double NXMAX_theo=R_Ewald*NA/n0;
    string sav_param=chemin_result+"/SAV_param_manip.txt";
    cout<<sav_param<<endl;
    ofstream fichier_sav_parametre(sav_param);

    cout<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    cout<<"+----------------+----------------+"<<endl;
    cout<<"|    Grandeur    |    Valeur      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp Holo     |     "<<tailleTheoPixelHolo*pow(10,9)<<" nm      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Lambda      |     "<<lambda0*pow(10,9)<<" nm     |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Gt          |     "<<"X "<<Gt<<"      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     Champ holo |     "<<round(IMAGE_DIMX*tailleTheoPixelHolo*pow(10,6))<<" µm      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     R_Ewald    |     "<<round(R_Ewald)<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     NXMAX_theo |     "<<round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"| chp UBorn pix  |     "<<2*round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp UBorn    |     "<<(tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI)*pow(10,9)<<" nm     |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"| chp UBorn µm   |     "<<2*round(R_Ewald*NA/n0)*(tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI)*pow(10,6) <<" µm    |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"| Chp tomo pixel |     "<<dim_final<<" pix^3  |" <<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|  Dela_fUborn   |   "<<Delta_fUborn<<" m-1  |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp Tomo theo|     "<<(2*NXMAX_theo)/dim_final*tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI*pow(10,9)<<" nm |"<<endl;
    cout<<"+---------------------------------+"<<endl;

    cout<<"\n##################### FIN INFO Reconstruction ##################\n"<<endl;

    fichier_sav_parametre<<"Fichiers d'acquisition : "<<chemin_acquis<<endl;
    fichier_sav_parametre<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    fichier_sav_parametre<<"+----------------+----------------+"<<endl;
    fichier_sav_parametre<<"|    Grandeur    |    Valeur      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Holo     |     "<<tailleTheoPixelHolo*pow(10,9)<<" nm      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Gt          |     "<<"X "<<Gt<<"      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     Champ holo |     "<<round(IMAGE_DIMX*tailleTheoPixelHolo*pow(10,6))<<" µm      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     R_Ewald    |     "<<round(R_Ewald)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     NXMAX_theo |     "<<round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"| chp UBorn pix  |     "<<2*round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp UBorn    |     "<<(tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI)*pow(10,9)<<" nm     |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"| chp UBorn µm   |     "<<2*round(R_Ewald*NA/n0)*(tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI)*pow(10,6) <<" µm    |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"| Chp tomo pixel |     "<<dim_final<<" pix^3  |" <<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|  Delta_fUborn   |   "<<Delta_fUborn<<" m-1  |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Tomo theo|     "<<(2*NXMAX_theo)/dim_final*tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI*pow(10,9)<<" nm |"<<endl;
    fichier_sav_parametre<<"+---------------------------------+"<<endl;

    fichier_sav_parametre<<"\n##################### FIN INFO Reconstruction ##################\n"<<endl;

    fichier_sav_parametre.close();
}
