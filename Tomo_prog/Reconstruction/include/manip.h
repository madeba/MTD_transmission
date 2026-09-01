#ifndef DEF_MANIP// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_MANIP

#include <cstdlib>
#include <iostream>
class manip {
public :
        double NA_obj;//objective numerical aperture
        double coef_NA_obj_limit;
        double Tp;//pixel size
        double n0,nM;//oil index and mounting medium (background)
        double rayon;
        double lambda0;
        int dimROI; //region of interest dimension cropped on camera
        size_t NXMAX, NYMAX; //maxium frequency for measured complex field
        //size_t NXMAX_OBJ;
        size_t premier_plan;
        size_t Num_Angle_final;
        size_t NbAngle;//peut différer de premier_angle-Num_angle_final  si des angles ont été enlevés dans le prétraitement
        size_t nbThreads;
        double theta;
        double f_tube; //Tube lens focal length
        double f_obj;//objective focal length
        double G_obj,Gt,Rf;
        bool b_CorrAber=false;
        bool b_Deroul=false;
        bool b_Born=true;

        bool b_Export_OTF=true;
        bool b_polar=false;
        bool b_reflex=false;
        size_t dim_final;
        double TpCam;//pixel size (camera)
        double R_th;
        double tailleTheoPixelHolo;
        double tailleTheoPixelUborn;
        double Delta_fUborn;
        double tailleTheoPixelTomo;
        size_t circle_cx,circle_cy;
        std::string chemin_result;
        std::string chemin_acquis;
        std::string chemin_racine; //polar only
        std::string chemin_config;
        std::string chemin_config_defaut;
        manip(std::string str_config_manip,std::string etat_polar, bool b_polar);//You can guess that the param 2 & 3 are polar only

       //void init();
} ;

#endif // DEF_MANIP
