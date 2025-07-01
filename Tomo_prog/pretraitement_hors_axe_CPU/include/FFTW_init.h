#ifndef __FFTW_INIT__
#define __FFTW_INIT__
#include "Point3D.h"
#include "Point2D.h"
#include <fftw3.h>
#include "struct.h"
class FFTW_init{

private:

public:
    size_t m_Nthread;
    fftw_plan p_forward_IN, p_backward_IN, p_forward_OUT, p_backward_OUT;
    int fftwThreadInit;

    fftw_complex *in=nullptr,*out=nullptr;
    double *in_double=nullptr;// fftw r2c


    FFTW_init(Point3D dim,size_t nbThread);
    FFTW_init(std::vector<double> const &entree, Point2D dim, size_t nbThread);//init r2c (real input)
    //init r2c or c2c
    FFTW_init(std::vector<double>const &entree, std::string str_geometry,size_t nbThread);

    FFTW_init(std::vector<std::complex<double>> const &entree,size_t nbThread);///init c2r for complex input
    FFTW_init(std::vector<double> const &entree, size_t nbThread);///init c2r for real input
    FFTW_init(Var2D dim, size_t nbThread);
    FFTW_init(Var2D dim2D,bool b_inplace, size_t nbThread);///INIT  INPLACE or OUTPLACE with boolean
    ~FFTW_init();
 void prepare_wisdom2D(Var2D dim, std::string chemin);

};

#endif
