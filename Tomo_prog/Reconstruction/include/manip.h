#ifndef DEF_MANIP// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_MANIP

#include <cstdlib>
#include <iostream>
struct manip {
        double NA_obj;
        double coef_NA_obj_limit;
        double Tp;
        double n0;
        double rayon;
        double lambda0;
        int dimROI;
        size_t NXMAX, NYMAX;
        //size_t NXMAX_OBJ;
        size_t premier_plan;
        size_t Num_Angle_final;
        size_t NbAngle;//peut différer de premier_angle-Num_angle_final  si des angles ont été enlevés dans le prétraitement
        size_t nbThreads;
        double theta;
        double f_tube;
        double f_obj;
        double G,Gt,Rf;
        bool b_CorrAber=false;
        bool b_Deroul=false;
        bool b_Born=true;
        bool b_Export_OTF=true;
        size_t dim_final;
        double TpCam;
        double R_th;
        double tailleTheoPixelHolo;
        double tailleTheoPixelUborn;
        double Delta_fUborn;
        double tailleTheoPixelTomo;
        size_t circle_cx,circle_cy;
        std::string chemin_result;
        std::string chemin_acquis;
        std::string chemin_config;
        std::string chemin_config_defaut;
        manip();
       //void init();
} ;

#endif // DEF_MANIP
