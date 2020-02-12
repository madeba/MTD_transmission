#include "manip.h"
#include <iostream>
#include "fonctions.h"
#include "projet.h"

using namespace std;

manip::manip()
{

}

void manip::init()
{
    Var2D dimROI= {WINDOW_X, WINDOW_Y};
    cout<<"dans la classe manip"<<endl;
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);

    string fic_cfg_recon=repertoire_config+"recon.txt";
    cout<<"fichier recon="<<fic_cfg_recon<<endl;
    string fic_cfg_manip=repertoire_config+"config_manip.txt";
    cout<<"chemin config="<<fic_cfg_manip<<endl;

    cout<<"\n##################### INFO MANIP ##################\n"<<endl;
    n0=extract_val("N0",fic_cfg_manip);	//indice de l'huile
    NA=extract_val("NA",fic_cfg_manip);	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
    lambda0=extract_val("LAMBDA",fic_cfg_manip);
    f_tube=extract_val("F_TUBE",fic_cfg_manip), ///focale lentille tube
    f_obj=extract_val("F_OBJ",fic_cfg_manip),///focale objectif
    G=f_tube/f_obj,	//grossissement telan+objectif
      TpCam=extract_val("TPCAM",fic_cfg_manip),//cam. photon focus
    Rf=extract_val("RF",fic_cfg_manip),//1.2;=1/facteur grossissement
    Gt=G/Rf;
    cout<<"Gt="<<Gt<<endl;
    size_t dim_final=extract_val("DIM_FINAL",fic_cfg_recon);
    //K=lambda*G/(2*NA*TpCam*Rf),// Facteur d'echelle
    //Tps=TpCam,//*K; //Taille pixel apres mise à l'echelle important si RF<>1
    tailleTheoPixelHolo=TpCam/Gt*pow(10,9);//Pour info, taille des pixels sur un hologramme=Tpcam/GT
    cout<<"taille theorique pixel holo="<<tailleTheoPixelHolo<<endl;
    theta=asin(NA/n0);
    circle_cx=extract_val("CIRCLE_CX",fic_cfg_manip);
    circle_cy=extract_val("CIRCLE_CY",fic_cfg_manip);
    NXMAX=extract_val("NXMAX",fic_cfg_manip);

    if(NXMAX==0)
        cout<<"----- ERREUR : NXMAX absent du fichier de configuration" <<fic_cfg_manip<<endl;
    dim_final=extract_val("DIM_FINAL",fic_cfg_recon);
    tailleTheoPixelUborn=tailleTheoPixelHolo*dimROI.x/(2*NXMAX);
    tailleTheoPixelTomo=tailleTheoPixelUborn*(2*NXMAX)/dim_final;

    premier_plan=extract_val("PREMIER_ANGLE",fic_cfg_recon),
    Num_Angle_final=extract_val("NB_HOLO_RECON",fic_cfg_recon),//

    cout<<"\n##################### Options de RECONSTRUCTION ##################\n"<<endl;
    b_CorrAber=extract_val("C_ABER",fic_cfg_recon);///corriger les aberrations?
    b_Deroul=extract_val("DEROUL",fic_cfg_recon);///Dérouler la phase?
    b_Born=extract_val("BORN",fic_cfg_recon);///Born vrai ? Sinon Rytov
    if(b_Born==true)
        cout<<"RYTOV=0"<<endl;
    else
        cout<<"RYTOV=1"<<endl;

    cout<<"\n##########################################################"<<endl;

        //Effacer précédent résultats.


        string tampon="UBornfinal_Im.raw";
        string result=chemin_result+tampon;
        if( remove(result.c_str())== 0 )
        perror( "Fichier Ubornfinal_Im impossible à effacer" );
        tampon="UBornfinal_Re.raw";
        result=chemin_result+tampon;
        if( remove(result.c_str()) == 0 )
        perror( "Fichier UBornfinal_Re impossible à effacer" );

      rayon=round(NXMAX*n0/NA);//calcul du rayon à partir de la fréquence NXMAX defini pare l'utlisateur
      double R_Ewald=IMAGE_DIMX*tailleTheoPixelHolo*n0/(lambda0*pow(10,9)); //vraie valeur de R_Ewald.
      double NXMAX_theo=R_Ewald*NA/n0;
    string sav_param="/home/mat/tomo_test/SAV_param_manip.txt";
    cout<<sav_param<<endl;
    ofstream fichier_sav_parametre(sav_param);

    cout<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    cout<<"+----------------+----------------+"<<endl;
    cout<<"|    Grandeur    |    Valeur      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp Holo     |     "<<tailleTheoPixelHolo<<" nm      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     Champ      |     "<<round(IMAGE_DIMX*tailleTheoPixelHolo/1000)<<" µm      |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     R_Ewald    |     "<<round(R_Ewald)<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|     NXMAX_theo |     "<<round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"| Taille chp cplx|     "<<2*round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp Chp cplx |     "<<round(tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI.x)<<" nm     |"<<endl;
    cout<<"|---------------------------------|"<<endl;
    cout<<"|    Tp Tomo theo|     "<<(2*NXMAX_theo)/dim_final*tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI.x<<" nm      |"<<endl;
    cout<<"+---------------------------------+"<<endl;
    cout<<"\n##################### FIN INFO Reconstruction ##################\n"<<endl;

    fichier_sav_parametre<<"Fichiers d'acquisition : "<<chemin_acquis<<endl;
    fichier_sav_parametre<<"\n##################### INFO Reconstruction ##################\n"<<endl;
    fichier_sav_parametre<<"+----------------+----------------+"<<endl;
    fichier_sav_parametre<<"|    Grandeur    |    Valeur      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Holo     |     "<<tailleTheoPixelHolo<<" nm      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     Champ      |     "<<round(IMAGE_DIMX*tailleTheoPixelHolo/1000)<<" µm      |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     R_Ewald    |     "<<round(R_Ewald)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|     NXMAX_theo |     "<<round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"| Taille chp cplx|     "<<2*round(R_Ewald*NA/n0)<<" pixels |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Chp cplx |     "<<round(tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI.x)<<" nm     |"<<endl;
    fichier_sav_parametre<<"|---------------------------------|"<<endl;
    fichier_sav_parametre<<"|    Tp Tomo theo|     "<<(2*NXMAX_theo)/dim_final*tailleTheoPixelHolo/(2*NXMAX_theo)*dimROI.x<<" nm      |"<<endl;
    fichier_sav_parametre<<"+---------------------------------+"<<endl;
    fichier_sav_parametre<<"\n##################### FIN INFO Reconstruction ##################\n"<<endl;

    fichier_sav_parametre.close();


}
