#include "manip.h"
#include <iostream>
#include "fonctions.h"
#include "projet.h"

using namespace std;

manip::manip()
{



    cout<<"dans la classe manip"<<endl;
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    string repertoire_config=extract_string("CHEMIN_CONFIG_PC_ACQUIS",home+fin_chemin_gui_tomo);
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);
//    dimROI=extract_val("DIM_ROI", chemin_config_manip);
    string fic_cfg_recon=repertoire_config+"recon.txt";
    cout<<"fichier recon="<<fic_cfg_recon<<endl;
    string fic_cfg_manip=repertoire_config+"config_manip.txt";
    cout<<"chemin config manip="<<fic_cfg_manip<<endl;

    // Var2D dimROI= {extract_val("DIM_ROI", fic_cfg_manip),extract_val("DIM_ROI", fic_cfg_manip)};
    dimROI=extract_val("DIM_ROI", fic_cfg_manip);
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
    tailleTheoPixelUborn=tailleTheoPixelHolo*dimROI/(2*NXMAX);
    tailleTheoPixelTomo=tailleTheoPixelUborn*(2*NXMAX)/dim_final;


    rayon=round(NXMAX*n0/NA);//calcul du rayon à partir de la fréquence NXMAX defini pare l'utlisateur
    int R_Ewald=dimROI*tailleTheoPixelHolo*n0/(lambda0*pow(10,9)); //vraie valeur de R_Ewald.
    cout<<"R_Ewald="<<R_Ewald<<endl;
    R_th_Ewald=R_Ewald;
    NXMAX_theo=round(R_Ewald*NA/n0);
    cout<<"NXMAX_theo="<<NXMAX_theo<<endl;
    coord_porteuse=3*NXMAX_theo/sqrt(2);//distance entre porteuse et axes abscisses/ordonnées


}
