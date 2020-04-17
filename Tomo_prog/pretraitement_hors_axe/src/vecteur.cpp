#include "vecteur.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <string>

vecteur :: vecteur (){
    this->x=0;
    this->y=0;
    this->z=0;
}
vecteur :: vecteur(double vx,double vy, double vz){
    this->x=vx;
    this->y=vy;
    this->z=vz;
}

void vecteur ::set_coord(double vx,double vy, double vz){
    this->x=vx;
    this->y=vy;
    this->z=vz;
}
double vecteur :: getx(){
    return this->x;
}
double vecteur :: gety(){
    return this->y;
}
double vecteur :: getz(){
    return this->z;
}

void vecteur :: setx(double vx){
    this->x=vx;
}
void vecteur :: sety(double vy){
    this->y=vy;
}
void vecteur :: setz(double vz){
    this->z=vz;
}

double vecteur :: norm(){
    return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
}

double vecteur::calc_angle(vecteur const &vec2){
    //double prod_scal=this->x*vec2.x+this->y*vec2.y+this->z*vec2.z;
//   double prod_scal=(*this)*vec2;//this->x*vec2.x+this->y*vec2.y+this->z*vec2.z;
double angle = atan2(  (*this^vec2).norm(), (*this)*vec2  );
return angle;
   // angle = atan(norm(cross(a,b)), dot(a,b))
   // std::cout<<"norme vec 1="<<this->norm()<<std::endl;
    //std::cout<<"norm vec2="<<vec2.norm()<<std::endl;
    //double angle=acos(prod_scal/(this->norm()*vec2.norm()));

}

///--------------Opérateurs---
///surcharge : produit scalaire
double vecteur::operator*(vecteur const &vec2)
{//std::cout<<"this->x="<<this->x<<std::endl;
    return this->x*vec2.x+this->y*vec2.y+this->z*vec2.z;
}
///multiplication par un nombre à droite
vecteur vecteur::operator*(double scalaire)
{//std::cout<<"this->x="<<this->x<<std::endl;
    vecteur result(scalaire*this->x,scalaire*this->y,scalaire*this->z);
    return result;
}

vecteur operator*(double scalaire, vecteur &v){ //fonction amie pour la commutativité de *
    return v*scalaire;
}

vecteur vecteur::operator+(vecteur const &vec2){
    vecteur result;
    result.x=this->x+vec2.x;
    result.y=this->y+vec2.y;
    result.z=this->z+vec2.z;
    return result;
}

vecteur vecteur::operator-(vecteur const &vec2){
vecteur result;
result.x=this->x-vec2.x;
result.y=this->y-vec2.y;
result.z=this->z-vec2.z;
return result;
}


/*vecteur operator^(vecteur v,vecteur w) //produit vectoriel
{vecteur z(
 v.gety()*w.getz()-w.gety()*v.getz() ,
 v.getz()*w.getx()-w.getz()*v.getx() ,
 v.getx()*w.gety()-w.getx()*v.gety()
 );
 return(z);
}*/


vecteur vecteur::operator^(vecteur w) //produit vectoriel
{vecteur z(
 this->gety()*w.getz()-w.gety()*this->getz() ,
 this->getz()*w.getx()-w.getz()*this->getx() ,
 this->getx()*w.gety()-w.getx()*this->gety()
 );
 return(z);
}





