#ifndef __FFTW_INIT__
#define __FFTW_INIT__
#include "Point3D.h"
#include "Point2D.h"
#include <fftw3.h>
class FFTW_init{

private:

public:
    unsigned int m_Nthread;
    fftw_plan p_forward_IN, p_backward_IN, p_forward_OUT, p_backward_OUT;
    int fftwThreadInit;
    fftw_complex *in,*out;

    FFTW_init(Point3D dim);
    FFTW_init(Point2D dim);
    FFTW_init(Point3D dim,size_t nbThread,bool b_inplace, unsigned int plan_type);
    ~FFTW_init();


};

#endif
