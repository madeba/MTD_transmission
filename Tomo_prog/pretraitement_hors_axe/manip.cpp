#include "manip.h"
#include <iostream>
#include "fonctions.h"
#include "projet.h"

using namespace std;
///Initialize paramaters from the set-up + path for config files
manip::manip()
{

    Var2D dimROI= {WINDOW_X, WINDOW_Y};
    cout<<"dans la classe manip"<<endl;
    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);

    string fic_cfg_recon=repertoire_config+"/recon.txt";
    cout<<"fichier recon="<<fic_cfg_recon<<endl;
    string fic_cfg_manip=repertoire_config+"/config_manip.txt";
    cout<<"chemin config="<<fic_cfg_manip<<endl;

    cout<<"\n##################### INFO MANIP ##################\n"<<endl;
    n0=extract_val("N0",fic_cfg_manip);	//indice de l'huile
    NA=extract_val("NA",fic_cfg_manip);	/// NA=Numerical aperture of the objective///ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
    lambda0=extract_val("LAMBDA",fic_cfg_manip);///lambda0 = laser wavelength
    f_tube=extract_val("F_TUBE",fic_cfg_manip), ///f_tube = tube lens focal length
    f_obj=extract_val("F_OBJ",fic_cfg_manip),///f_obj=focale objectif
    G=f_tube/f_obj,	///magnitude of (telan+objectif)
    TpCam=extract_val("TPCAM",fic_cfg_manip),///Tpcam = pixel size of the camera
    Rf=extract_val("RF",fic_cfg_manip),//1.2;=1/facteur grossissement
    Gt=G/Rf;///Gt : total magnitude of the optical system. grandissement total
    cout<<"Gt="<<Gt<<endl;
    size_t dim_final=extract_val("DIM_FINAL",fic_cfg_recon);
    tailleTheoPixelHolo=TpCam/Gt*pow(10,9);//Pour info, taille des pixels sur un hologramme=Tpcam/GT
    cout<<"taille theorique pixel holo="<<tailleTheoPixelHolo<<endl;
    theta=asin(NA/n0);
    circle_cx=extract_val("CIRCLE_CX",fic_cfg_manip);///coordinates of the off-axis career
    circle_cy=extract_val("CIRCLE_CY",fic_cfg_manip);
    NXMAX=extract_val("NXMAX",fic_cfg_manip);
    cout<<"m1.NXMAX="<<NXMAX<<endl;
    if(NXMAX==0)
        cout<<"----- ERREUR : NXMAX absent du fichier de configuration" <<fic_cfg_manip<<endl;
    dim_final=extract_val("DIM_FINAL",fic_cfg_recon);
    tailleTheoPixelUborn=tailleTheoPixelHolo*dimROI.x/(2*NXMAX);
    tailleTheoPixelTomo=tailleTheoPixelUborn*(2*NXMAX)/dim_final;

    premier_plan=extract_val("PREMIER_ANGLE",fic_cfg_recon),
    Num_Angle_final=extract_val("NB_HOLO",fic_cfg_manip),//
    NbAngle=Num_Angle_final-premier_plan;
    nbThreads=extract_val("NB_THREADS",fic_cfg_recon),
    CamDimROI=extract_val("CCD_ROIX",fic_cfg_recon),

    cout<<"\n##################### Options de RECONSTRUCTION ##################\n"<<endl;
    ///boolean used to choose the type of reconstruction/phase unwrapping
    b_CorrAber=extract_val("C_ABER",fic_cfg_recon);///corriger les aberrations?
    b_Deroul=extract_val("DEROUL",fic_cfg_recon);///Dérouler la phase?
    b_Born=extract_val("BORN",fic_cfg_recon);///Born vrai ? Sinon Rytov
    b_volkov=extract_val("VOLKOV",fic_cfg_recon);///VOLKOV vrai ? Sinon Herraez
    if(b_Born==true)
        cout<<"RYTOV=0"<<endl;
    else
        cout<<"RYTOV=1"<<endl;

    cout<<"\n##########################################################"<<endl;

        ///Effacer précédent résultatscar fonction d'enregistrement cumulative ("a+b")

        Var2D dim2DHA={2*NXMAX,2*NXMAX};
        dimImg=to_string(dim2DHA.x)+"x"+to_string(dim2DHA.y)+"x"+to_string(NbAngle);

        //string tampon="/UBornfinal_Im.raw";
        string tampon="/UBornfinal_Im"+dimImg+".raw";
        string result=chemin_result+tampon;
        cout<<"Effacement Uborn : "<<result<<endl;
        if( remove(result.c_str())== 0 )
        perror( "Fichier Ubornfinal_Im impossible à effacer" );
       // tampon="/UBornfinal_Re.raw";
        tampon="/UBornfinal_Re"+dimImg+".raw";
        result=chemin_result+tampon;
        if( remove(result.c_str()) == 0 )
        perror( "Fichier UBornfinal_Re impossible à effacer" );

      rayon=round(NXMAX*n0/NA);//calcul du rayon à partir de la fréquence NXMAX defini pare l'utlisateur
      double R_Ewald=IMAGE_DIMX*tailleTheoPixelHolo*n0/(lambda0*pow(10,9)); //vraie valeur de R_Ewald.
      double NXMAX_theo=R_Ewald*NA/n0;
    ///-------------enregistrer les paramètres dans un fichier log.--------------------------------------
    string sav_param=chemin_result+"/SAV_param_manip.txt";
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
