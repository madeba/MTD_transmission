#include "FFT_encaps.h"
using namespace std;
///init fftw pour FFT3D
FFT_encaps::FFT_encaps(Point3D dim,size_t nbThreads)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;
    m_Nthread=nbThreads;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
    p_backward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_BACKWARD, FFTW_ESTIMATE );
}
FFT_encaps::FFT_encaps(Point3D dim,size_t nbThreads, bool b_inPlace)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;
    m_Nthread=nbThreads;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    if(b_inPlace==0){
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
    p_backward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_BACKWARD, FFTW_ESTIMATE );
    }
    else{
    int isWisdomOK=import_wisdom("/home/mat/Dropbox/projet_c/2020/Projet_tomo/Tomo_config/Wisdom/wisdom3D_512_c2c_double_backward_inplace_i5-3550.txt");
    cout<<"isWisdomOk="<<isWisdomOK<<endl;
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    //out=nullptr;
    p_forward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_FORWARD, FFTW_WISDOM_ONLY);
    p_backward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_BACKWARD, FFTW_WISDOM_ONLY);
    }
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
    p_forward_OUT=fftw_plan_dft_2d(dim.x, dim.y,  in, out,FFTW_FORWARD, FFTW_MEASURE);
    p_backward_OUT=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_BACKWARD, FFTW_MEASURE);
}

bool
FFT_encaps::import_wisdom(const char *filename)
{
  FILE *wisdom_file;
  wisdom_file = fopen(filename, "r");

  if ( wisdom_file )
    {
      bool value = fftw_import_wisdom_from_file(wisdom_file);
      fclose(wisdom_file);
      return value;
    }

  return false;
}
///destructeur
FFT_encaps::~FFT_encaps(){
//fftw_free(in);
//fftw_free(out);
//fftw_destroy_plan(p_forward_OUT);
//fftw_destroy_plan(p_backward_OUT);
}

