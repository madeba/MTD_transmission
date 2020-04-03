#ifndef DEF_MANIP// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_MANIP

#include <cstdlib>
#include <iostream>
#include "Point2D.h"
//#include "struct.h"
struct manip {

        double NA;
        double Tp;
        double n0;
        double rayon;
        double lambda0;
        Point2D dimROI;
        size_t NXMAX;
        size_t dim2DUBorn;
        size_t premier_plan;
        size_t Num_Angle_final;
        double theta;
        double f_tube;
        double f_obj;
        double G,Gt,Rf;
        bool b_CorrAber=false;
        bool b_Deroul=false;
        bool b_Born=true;
        bool b_Volkov=true;
        size_t dim_final;
        double TpCam;
        double R_th;
        double tailleTheoPixelHolo;
        double tailleTheoPixelUborn;
        double tailleTheoPixelTomo;
        size_t circle_cx,circle_cy;

        std::string chemin_result;
        std::string chemin_acquis;
        manip();
        void init();
} ;

#endif // DEF_MANIP
