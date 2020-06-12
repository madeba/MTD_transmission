#include "FFTW_init.h"
#include "struct.h"
#include <fftw3.h>
using namespace std;
///Different constructor to initialize useful variables for fftw


/*void FFTW_init::prepare_wisdom2D(Var2D dim, string chemin)
{
        fftw_plan_with_nthreads(4);
        int N=dim.x*dim.y;

        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_BACKWARD, FFTW_EXHAUSTIVE);
        fftw_export_wisdom_to_filename(chemin.c_str());
        fftw_destroy_plan(p);

}*/


FFTW_init::FFTW_init(Point3D dim,size_t nbThreads, bool b_inPlace)
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
    int isWisdomOK=import_wisdom("/home/mat/Dropbox/projet_c/2020/Projet_tomo/Tomo_config/Wisdom/wisdom3D_512_c2c_double_backward_EXHAUSTIVE_inplace_4Thrd_i5-3550.txt");
   // int isWisdomOK=import_wisdom("wisdom3D_512_c2c_double_backward_EXHAUSTIVE_inplace_4Thrd_i5-3550.txt");
    cout<<"isWisdomOk="<<isWisdomOK<<endl;
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=nullptr;
    p_forward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_FORWARD, FFTW_MEASURE);///problème, wisdom calculée en backward !!->MEASURE! possiblité d'avoir 2 wisdom ?
    p_backward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_BACKWARD, FFTW_MEASURE);
    //p_backward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_BACKWARD, FFTW_WISDOM_ONLY);
   // p_forward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_FORWARD, FFTW_ESTIMATE);///problème, wisdom calculée en backward !!
   // p_backward_IN=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, in,FFTW_BACKWARD, FFTW_ESTIMATE);
    }
}
///init3D fftw inplace
FFTW_init::FFTW_init(Point3D dim,size_t nbThread,bool b_inplace, unsigned int plan_type)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;

    m_Nthread=nbThread;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    if(b_inplace=1){
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

///init3D fftw inplace
FFTW_init::FFTW_init(vector<double> const &entree,size_t nbThread,bool b_inplace, unsigned int plan_type)
{
    unsigned int nbPix=entree.size();
    size_t dim=pow(entree.size(),-3);
    m_Nthread=nbThread;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    if(b_inplace=1){
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);

    p_forward_IN=fftw_plan_dft_3d(dim, dim, dim, in, in, FFTW_FORWARD, plan_type );
    p_backward_IN=fftw_plan_dft_3d(dim, dim, dim, in, in, FFTW_BACKWARD, plan_type);
    }
    else{
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_IN=fftw_plan_dft_3d(dim, dim, dim, in, out, FFTW_FORWARD, plan_type );
    p_backward_IN=fftw_plan_dft_3d(dim, dim, dim, in, out, FFTW_BACKWARD, plan_type);
    }
}



///init3D fftw outplace
FFTW_init::FFTW_init(Point3D dim,size_t nbThread)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;

    m_Nthread=nbThread;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
    p_backward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
}




bool FFTW_init::import_wisdom(const char *filename)
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
FFTW_init::~FFTW_init(){
//fftw_free(in);
//fftw_free(in_double);
//fftw_free(out);
//fftw_destroy_plan(p_forward_OUT);
//fftw_destroy_plan(p_backward_OUT);
}

