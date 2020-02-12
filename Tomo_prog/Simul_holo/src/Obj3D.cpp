#include "Obj3D.h"
#include <vector>
#include <complex>
#include <fftw3.h>
#include "Point3D.h"

//Obj3D::Obj3D(int dim):m
Obj3D::Obj3D(int dim)
{
    nbPix=dim*dim*dim;
//    dim3D(dim,dim,dim);
    for(size_t cpt=0;cpt<nbPix;cpt++)
    {
      Valeur[cpt].real(0);
      Valeur[cpt].imag(0);
    }
    //ctor
}



Obj3D::~Obj3D()
{
    //dtor
}


double Obj3D::energie()
{
    double E=0;
    for(int cpt=0;cpt<dim3D.x*dim3D.y*dim3D.z;cpt++)
    {
        E=E+pow(Valeur[cpt].real(),2)+pow(Valeur[cpt].imag(),2);
    }
    return E;
}
