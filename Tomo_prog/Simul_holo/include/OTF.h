#ifndef OTF_H
#define OTF_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include "Point2D.h"
#include  "manip.h"

class OTF{
    private:


    public:
        OTF(manip m1);
        OTF(const manip &m1, bool b_reflex);//boolean to override manip->b_reflex value (to adapt simulation)
        //OTF(size_t rayon_ewald);
        //const std::vector<std::complex<double>>& GetV() const{return valeur;}

        virtual ~OTF();
        manip manipOTF;
        std::vector<std::complex<double>> Valeur;
       // std::complex<double> n0;

        size_t nbPixEff;
        size_t nbPixRedon;
        bool b_Reflex;
      //  std::vector<std::complex<double>> valeur;
        void retropropag(Point2D spec);
       // void bFleur();

       // std::vector<Var2D> bFleur();

        std::vector<Point2D> bFleur(short unsigned int nbAxes);
        void bFleur(std::vector<Point2D> &CoordSpec);
        void bFleur(std::vector<Point2D> &CoordSpec, size_t const nbAxes);
        void bSpiral();
        void bSpiralNU();
        void bDblSpiral();
        void bCercle(int pourcentage_NA);
        void bGrille();
        void bMultiCercleUNI(int nb_cercle);
        void bFermat(int nbHolo);
        void symetrize_xoy();
        void symetrize_central();
        void scan_uniform3D(std::vector<Point2D> &CoordSpec, float coef_limit);//specular coordinate and coeffcient to limit the scanning (to avoid pixel on the  border NXMAX
    protected:

};

#endif //
