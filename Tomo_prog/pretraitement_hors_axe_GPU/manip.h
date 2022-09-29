#ifndef DEF_MANIP// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_MANIP

#include <cstdlib>
#include <iostream>
#include "struct.h"
class manip {
public :
        double NA_obj, Tp;
        double n0;
        double rayon;
        double lambda0;
        size_t NXMAX;

        size_t premier_plan;
        size_t Num_Angle_final;
        size_t NbAngle;
        size_t nbThreads;
        size_t CamDimROI;
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
        double R_th;
        double tailleTheoPixelHolo;
        double tailleTheoPixelUborn;
        double tailleTheoPixelTomo;
        size_t circle_cx,circle_cy;
        Var2D fPort, fPortShift;
        std::string dimImg;//string containing dimension of the field, and put in the name of the saved file
        std::string chemin_result;
        std::string chemin_acquis;

        manip();

} ;

#endif // DEF_MANIP
