#ifndef __OBJ3D__
#define __OBJ3D__

#include <fftw3.h>
#include <vector>
#include <complex>
#include "Point3D.h"
#include "manip.h"
#include "struct.h"

class Obj3D{
    private:



    public:
        //Obj3D();
        Obj3D(int dim);
        virtual ~Obj3D();
        double energie();


        std::vector<std::complex<double>> Valeur;
        Point3D dim3D;
        Point2D dim2D;
        size_t nbPix;

      //  manip m;

    protected:

};

#endif // OBJ3D_H
