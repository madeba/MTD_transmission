#include "FFT_encaps.h"
using namespace std;
///init fftw pour FFT3D
FFT_encaps::FFT_encaps(Point3D dim)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;
    m_Nthread=4;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
    p_backward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_BACKWARD, FFTW_ESTIMATE );
}

///init fftw pour FFT2D
FFT_encaps::FFT_encaps(Point2D dim)
{
    unsigned int nbPix=dim.x*dim.y;
    m_Nthread=4;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_2d(dim.x, dim.y,  in, out,FFTW_FORWARD, FFTW_MEASURE );
    p_backward_OUT=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_BACKWARD, FFTW_MEASURE);
}

///destructeur
FFT_encaps::~FFT_encaps(){
fftw_free(in);
fftw_free(out);
fftw_destroy_plan(p_forward_OUT);
fftw_destroy_plan(p_backward_OUT);
}

