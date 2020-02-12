#ifndef __FFT_ENCAPS__
#define __FFT_ENCAPS__
#include "Point3D.h"
#include "Point2D.h"
#include "manip.h"
#include <fftw3.h>
class FFT_encaps{

private:

public:
    unsigned int m_Nthread;
    fftw_plan p_forward_IN, p_backward_IN, p_forward_OUT, p_backward_OUT;
    int fftwThreadInit;
    fftw_complex *in,*out;

    FFT_encaps(Point3D dim,manip m1);
    FFT_encaps(Point2D dim,manip m1);
    ~FFT_encaps();


};

#endif
