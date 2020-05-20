#ifndef Point2D_H
#define Point2D_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <complex>


class Point2D{
    private:



    public:
        unsigned short int dim2D;
        double x,y;



        Point2D();
        Point2D(double x, double y,int dim);
        Point2D(int x, int y,int m_dim2D);

        virtual ~Point2D();
        size_t cpt2D();
        void decale(int decalx, int decaley);
        Point2D coordI();
        Point2D coordI2coordH();
        void set_xy(double x, double y);
        void setDim(unsigned short int dim);

    protected:
};

#endif //
