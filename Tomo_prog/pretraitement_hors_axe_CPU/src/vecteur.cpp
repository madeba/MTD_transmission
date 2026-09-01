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
void vecteur :: setNorm(int newNormValue)
{
this->norm=newNormValue;
}
double vecteur :: calc_norm(){
     if(this->z==0){
        std::cout<<"warning : kz=0"<<std::endl;
     }
    return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
}

void vecteur::calc_angle(){

    if(this->norm==0.0)
        {
        throw std::runtime_error("Vecteur nul : impossible de calculer les angles.");
    }
//kz is known inside the class, but is caculated again to avoid rounded error due to component being integer
//it may create differences betwwen float and integer calculations.
double kz=sqrt(this->norm*this->norm-this->x*this->x-this->y*this->y);
theta=acos(kz/this->norm);
this->phi=atan2(this->y,this->x);
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





