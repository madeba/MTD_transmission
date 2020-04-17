#ifndef VECTEUR_INCLUDED
#define VECTEUR_INCLUDED
#include <iostream>
class vecteur {
    private :
        double x;
        double y;
        double z;

    public :
        vecteur();
        vecteur(double vx,double vy, double vz);
        void set_coord(double vx,double vy, double vz);

        double getx();
        double gety();
        double getz();
        void setx(double vx);
        void sety(double vy);
        void setz(double vy);
        double norm();
        double calc_angle(vecteur const &vec2);

        double operator*(vecteur const &vec2);
        vecteur operator*(double scalaire);
        friend vecteur operator*(double a, vecteur &v);//fonctoion amie pour la commutativité de *
       // friend vecteur operator^(vecteur v1, vecteur v2);//fonctoion amie pour la commutativité de *
        vecteur operator+(vecteur const &vec2);
        vecteur operator-(vecteur const &vec2);
        vecteur operator^(vecteur  vec2); //produit vectoriel
        //void record2D(vecteur vec, std::string axe, std::string chemin, std::string options, int nbPix)
} ;


#endif // VECTEUR_INCLUDED

