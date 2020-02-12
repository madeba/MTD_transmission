#ifndef SPECTRE_H
#define SPECTRE_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include "Point3D.h"

class Spectre{
    private:



    public:
        Spectre();
        virtual ~Spectre();

        int dim;
        int length;
       // std::vector<std::complex<double>> valeur;
        fftw_plan p_forward,p_backward;

        Spectre(int m_dim3D);
        void prepare_fftw();

    protected:

};

#endif //
