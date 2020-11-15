#include "FFTW_init.h"
using namespace std;
///init fftw pour FFT3D
FFTW_init::FFTW_init(Point3D dim)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;
    m_Nthread=4;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_FORWARD, FFTW_ESTIMATE );
    p_backward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_BACKWARD, FFTW_ESTIMATE );
}

///init fftw pour FFT2D
FFTW_init::FFTW_init(Point2D dim)
{
    unsigned int nbPix=dim.x*dim.y;
    m_Nthread=4;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_2d(dim.x, dim.y,  in, out,FFTW_FORWARD, FFTW_ESTIMATE );
    p_backward_OUT=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE );
}

///init3D fftw /outinplace, nb_threads, fftw_measure/fftw_wisdom
FFTW_init::FFTW_init(Point3D dim,size_t nbThread,bool b_inplace, unsigned int plan_type)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;

    m_Nthread=nbThread;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    if(b_inplace==1){
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);

    p_forward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in, FFTW_FORWARD, plan_type );
    p_backward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in, FFTW_BACKWARD, plan_type);
    }
    else{
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out, FFTW_FORWARD, plan_type );
    p_backward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out, FFTW_BACKWARD, plan_type);
    }
}
///destructeur
FFTW_init::~FFTW_init(){
fftw_free(in);
fftw_free(out);
fftw_destroy_plan(p_forward_OUT);
fftw_destroy_plan(p_backward_OUT);
}

