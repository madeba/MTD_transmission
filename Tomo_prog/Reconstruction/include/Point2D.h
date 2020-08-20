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

        double x,y;
        int dim2D;
        Point2D();
        Point2D(double a_x, double a_y,int a_dim):x{a_x},y{a_y},dim2D{a_dim}{}; //liste d'init pour éviter l'appel au constructeur par défaut : plus performant.
        virtual ~Point2D();
        size_t cpt2D();
        void decale(int decalx, int decaley);
        Point2D coordI();
        Point2D coordI2coordH();
        void set_xy(double x, double y);

    protected:
};

#endif // CHAMPOPT_H
