#ifndef OTF_H
#define OTF_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include "Point3D.h"
#include  "Obj3D.h"

class OTF{
    private:


    public:
        OTF(manip m1);
        //OTF(size_t rayon_ewald);


        //const std::vector<std::complex<double>>& GetV() const{return valeur;}

        virtual ~OTF();

        manip manipOTF;


        std::vector<std::complex<double>> Valeur;
       // std::complex<double> n0;

        size_t nbPixEff;
        size_t nbPixRedon;

      //  std::vector<std::complex<double>> valeur;
        void retropropag(Point2D spec);
       // void bFleur();

        std::vector<Var2D> bFleur();
        void bSpiral();
        void bSpiralNU();
        void bDblSpiral();
        void bCercle(int pourcentage_NA);
        void bGrille();
        void bMultiCercleUNI(int nb_cercle);

        void symetrize_xoy();

    protected:

};

#endif //
