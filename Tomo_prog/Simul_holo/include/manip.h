#ifndef DEF_MANIP// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_MANIP

#include <cstdlib>
#include <iostream>
struct manip {

        double NAObj,NACond; //ouverture numérique objectif et condenseur
        double n0;         //indice du milieu d'immersion
        double R_EwaldPix; //Rayon de la sphere d'Ewald en pixel
        double R_EwaldMet; //Rayon de la sphere d'Ewald métrique
        double lambda_v;   //longueur d'onde dans levide
        int NXMAX;      //fréquence spatiale maximum enregistrée par le miscrocope en pixel
        int NXMAX_cond;
        size_t nbHolo;
        double theta_max;  //angle de collection max de l'objectif
        double Gt;         //grandissemen total
        size_t dim_final,dim_Uborn; //dimension epsace objet tomo et dimension champ complexe
        double TpCam;               //---échantillonnage
        double Tp_Tomo;
        double Tp_holo;
        double Tp_Uborn;
        double Delta_f_holo;
        double Delta_f_tomo;//a priori inutile car egal à Delta_f_Uborn
        double Delta_f_Uborn;
        std::string dimImg;//string containing dimension of the field, and put in the name of the saved file
        std::string chemin_result;
        std::string chemin_acquis;
        manip(int dimROI);
} ;

#endif // DEF_MANIP
