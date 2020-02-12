#ifndef CHAMPOPT_H
#define CHAMPOPT_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <complex>

class champOpt : public std::complex<double>{
    private:

    unsigned int dimx,dimy;


    public:
        champOpt();
        virtual ~champOpt();
        void fft2D(fftw_complex *in_out,std::vector<complex<double> > sortie,fftw_plan p);
        void fftshift();

    protected:

};

#endif // CHAMPOPT_H
