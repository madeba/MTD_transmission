#ifndef __FFT_ENCAPS__
#define __FFT_ENCAPS__
#include "Point3D.h"
#include "Point2D.h"
#include <fftw3.h>
class FFT_encaps{

private:

public:
    unsigned int m_Nthread;
    fftw_plan p_forward_IN, p_backward_IN, p_forward_OUT, p_backward_OUT;
    int fftwThreadInit;
    fftw_complex *in,*out;
    FFT_encaps();
    FFT_encaps(Point3D dim,size_t nbThread);
    FFT_encaps(Point3D dim,size_t nbThread,bool inplace);
    FFT_encaps(Point3D dim,size_t nbThreads, bool b_inPlace,std::string chemin_wisdom);
    explicit FFT_encaps(Point2D dim);//explicit evite les conversion à l'arrache d'un paramètre d'entree different de Point2D

    bool import_wisdom(const char *filename);
    ~FFT_encaps();

};

#endif
