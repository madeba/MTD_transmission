#ifndef POINT3D_H
#define POINT3D_H


#include <fftw3.h>
#include <vector>
#include <complex>
#include <iostream>
#include "Point2D.h"

class Point3D{
    private:



    public:
        double x,y,z;
        int dim3D;
        //dim3d int sinon problème multipliccation par unsigned

        Point3D();
        Point3D(double x, double y, double z);

        Point3D(double x, double y, double z,int m_dim3D);
        Point3D(Point2D Pt2D, double z, int dim3D);
        virtual ~Point3D();
        size_t cpt3D();
        void decale(int decalx, int decaley, int decalez);
        Point3D coordI();
        Point3D coordH();
        int get_dim3D();
        void set_coord3D(double x, double y, double z);
        Point3D sym_xoz();
        Point3D operator-(Point3D const &Pt2);///const & permet d'éviter une copie mais rend non modifiable

    protected:

};

#endif // CHAMPOPT_H
