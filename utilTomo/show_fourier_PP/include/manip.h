#ifndef DEF_MANIP// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_MANIP

#include <cstdlib>
#include <iostream>
struct manip {
        int dimROI;
        double NA;
        double Tp;
        double n0;
        double rayon;
        double lambda0;
        size_t NXMAX;
        unsigned int NXMAX_theo;
        size_t premier_plan;
        size_t Num_Angle_final;
        double theta;
        double f_tube;
        double f_obj;
        double G,Gt,Rf;
        bool b_volkov=false;
        bool b_CorrAber=false;
        bool b_Deroul=false;
        bool b_Born=true;
        size_t dim_final;
        double TpCam;
        double R_th_Ewald;
        double tailleTheoPixelHolo;
        double tailleTheoPixelUborn;
        double tailleTheoPixelTomo;
        int coord_porteuse;
        double delta_f_Holo;
        size_t circle_cx,circle_cy;
        std::string chemin_result;
        std::string chemin_acquis;
        manip();
       // void init();
} ;

#endif // DEF_MANIP
