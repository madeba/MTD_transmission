#include "FFTW_init.h"
#include "struct.h"
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

///init fftw for FFT3D
FFTW_init::FFTW_init(Point3D dim,size_t nbThread)
{
    unsigned int nbPix=dim.x*dim.y*dim.z;

    m_Nthread=nbThread;

    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_FORWARD, FFTW_ESTIMATE );
    p_backward_OUT=fftw_plan_dft_3d(dim.x, dim.y, dim.z, in, out,FFTW_BACKWARD, FFTW_ESTIMATE );
}

///init fftw for fft2D r2c
FFTW_init::FFTW_init(vector<double>const &entree, string str_geometry,size_t nbThread)
{
    size_t dim=sqrt(entree.size());
    Var2D dim2D={dim,dim};
    m_Nthread=nbThread;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
  //  in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    if(str_geometry=="r2c"){
    unsigned int nbPixCplx=dim2D.x*(dim2D.y/2+1);
    unsigned int nbPix = dim2D.x*dim2D.y;


    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPixCplx);//Attention out est défini avec 2 fois moins de pixel !
    in_double=(double*) fftw_malloc(sizeof(double) * nbPix);
    p_forward_OUT=fftw_plan_dft_r2c_2d(dim2D.x, dim2D.y, in_double, out,FFTW_MEASURE);
    }
}

//choisir inplace ou outplace (str  IN || OUT),

///init fftw c2c : input is a square complex matrix
FFTW_init::FFTW_init(vector<complex<double>>const &entree,size_t nbThread)
{
    size_t nbPix=entree.size();
    size_t dim=sqrt(nbPix);
    m_Nthread=nbThread;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_2d(dim, dim, in, out,FFTW_FORWARD , FFTW_MEASURE);///FFTW_MEASURE= 0.15s, FFTW_PATIENT=8s
    p_backward_OUT=fftw_plan_dft_2d(dim, dim, in, out,FFTW_BACKWARD, FFTW_MEASURE);
}

///init fftw c2c : input is a square complex matrix
FFTW_init::FFTW_init(Var2D dim, size_t nbThread)
{

    size_t nbPix=dim.x*dim.y;
    m_Nthread=nbThread;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_FORWARD , FFTW_MEASURE);
    p_backward_OUT=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_BACKWARD, FFTW_MEASURE);
}

///init fftw c2c : input is a square matrix
/*FFTW_init::FFTW_init(Var2D dim, size_t nbThread, string str_inplace)
{

    size_t nbPix=dim.x*dim.y;
    m_Nthread=nbThread;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_IN=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_FORWARD , FFTW_MEASURE);
   // p_backward_OUT=fftw_plan_dft_2d(dim.x, dim.y, in, out,FFTW_BACKWARD, FFTW_MEASURE);
}*/

///init fftw c2c : input is a square REAL matrix
FFTW_init::FFTW_init(vector<double>const &entree,size_t nbThread)
{

    size_t nbPix=entree.size();
     size_t dim=sqrt(nbPix);
    m_Nthread=nbThread;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
    p_forward_OUT=fftw_plan_dft_2d(dim, dim, in, out,FFTW_FORWARD , FFTW_MEASURE);
    p_backward_OUT=fftw_plan_dft_2d(dim, dim, in, out,FFTW_BACKWARD, FFTW_MEASURE);
   //  fftw_export_wisdom_to_filename("planner_patient.txt");
}

///-----------Init INPLACE (ou OUTPLACE) C2C-------------------
FFTW_init::FFTW_init(Var2D dim2D,bool b_inplace,size_t nbThread)
{
   // unsigned int nbPixCplx=dim2D.x*(dim2D.y/2+1);
    unsigned int nbPix = dim2D.x*dim2D.y;
    m_Nthread=nbThread;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(m_Nthread);

   if(b_inplace==false){
    out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);//Attention out est défini avec 2 fois moins de pixel !
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
//fftw_plan forwardA = fftw_plan_dft_r2c_2d(nrows,ncols,  Aa,         Af,FFTW_FORWARD | FFTW_ESTIMATE);
    p_forward_OUT=fftw_plan_dft_2d(dim2D.x, dim2D.y, in, out,FFTW_FORWARD,FFTW_MEASURE);
    p_backward_OUT=fftw_plan_dft_2d(dim2D.x, dim2D.y, in, out,FFTW_BACKWARD,FFTW_MEASURE);
   }
   else
   {
    in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbPix);
//fftw_plan forwardA = fftw_plan_dft_r2c_2d(nrows,ncols,  Aa,         Af,FFTW_FORWARD | FFTW_ESTIMATE);
    p_forward_IN=fftw_plan_dft_2d(dim2D.x, dim2D.y, in, in,FFTW_FORWARD,FFTW_MEASURE);
    p_backward_IN=fftw_plan_dft_2d(dim2D.x, dim2D.y, in, in,FFTW_BACKWARD,FFTW_MEASURE);
   }
}
///destructeur
FFTW_init::~FFTW_init(){
//fftw_free(in);
//fftw_free(in_double);
//fftw_free(out);
//fftw_destroy_plan(p_forward_OUT);
//fftw_destroy_plan(p_backward_OUT);
}

