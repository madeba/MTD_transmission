#ifndef SPECIMEN_H_INCLUDED
#define SPECIMEN_H_INCLUDED

#include "Obj3D.h"
#include <complex>

class Specimen:public Obj3D
{

    Specimen(int dim,std::complex<double> n);
};

#endif // SPECIMEN_H_INCLUDED
