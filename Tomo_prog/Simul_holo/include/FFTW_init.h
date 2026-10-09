#ifndef __FFTW_INIT__
#define __FFTW_INIT__
#include "Point3D.h"
#include "Point2D.h"
#include <fftw3.h>
class FFTW_init{

private:


public:
    unsigned int nbPix;//number of pixel used to reserve memory (i.e size of the image)

    fftw_plan p_forward_IN, p_backward_IN, p_forward_OUT, p_backward_OUT;

    fftw_complex *in,*out;
   /* fftw_complex *in = nullptr;      // ← Initialisation par défaut
    fftw_complex *out = nullptr;
    fftw_plan p_forward_IN = nullptr;
    fftw_plan p_backward_IN = nullptr;
    fftw_plan p_forward_OUT = nullptr;
    fftw_plan p_backward_OUT = nullptr;*/
    unsigned int m_Nthread = 4;
    int fftwThreadInit = 0;
    unsigned int getFFTSize();
    FFTW_init(Point3D dim);
    FFTW_init(Point2D dim);
    FFTW_init(std::vector<double>const &entree,size_t nbThread);
    FFTW_init(Point3D dim,size_t nbThread,bool b_inplace, unsigned int plan_type);
    ~FFTW_init();


};

#endif
