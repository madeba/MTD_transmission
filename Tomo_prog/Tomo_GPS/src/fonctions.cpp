#include <fstream>//ifstream
#include <vector>
#include <chrono>
#include <complex>
#include <fftw3.h>
#include "struct.h"
//#include <cv.h>
//#include <highgui.h>//imread
//#include "projet.h"
#include "FFT_fonctions.h"
#include "fonctions.h"
#include "IO_fonctions.h"
using namespace std;
using namespace cv;
///free functions used outside classes

///translate normal coordinate to shifted coordinate
///allows to work with shifted version of the spectrum (to avoir fftshift and speedup spectrum crop)


Var2D coord_to_coordShift(Var2D coord2D, Var2D dimROI)
{
Var2D coord2D_shift;
    if(coord2D.x-dimROI.x/2>0)
        coord2D_shift.x=coord2D.x-dimROI.x/2;
    else coord2D_shift.x=dimROI.x/2+coord2D.x;
     if(coord2D.y-dimROI.y/2>0)
        coord2D_shift.y=coord2D.y-dimROI.y/2;
    else coord2D_shift.y=dimROI.y/2+coord2D.y;
    return coord2D_shift;
}
///calculate the wrapped phase from -pi to pi
void calcPhase_mpi_pi_atan2(vector<complex<double>> const &cplxField, vector<double> &phaseMod2pi)///calcul phase -PI-PI
{
for(int cpt=0;cpt<cplxField.size();cpt++)
phaseMod2pi[cpt]=atan2(cplxField[cpt].imag(),cplxField[cpt].real());
}
///--------Initialize reference
vector<double> initRef(string chemin_ref, Var2D coin, Var2D dimROI){
size_t nbPixROI2d=dimROI.x*dimROI.y;
vector<double> ampli_ref(nbPixROI2d);
vector<double> intensite_ref(nbPixROI2d);

   if(is_readable(chemin_ref)==1)
         charger_image2D_OCV(intensite_ref,chemin_ref, coin, dimROI);
         else cout<<"/!\\  fichier intensité référence absent, création intensité unité"<<endl;
    for(size_t cpt=0;cpt<nbPixROI2d;cpt++){
        if(intensite_ref[cpt]!=0)
        ampli_ref[cpt]=sqrt(intensite_ref[cpt]);
        else
        ampli_ref[cpt]=1;
    }
return ampli_ref;
}
///---arbitrary shift (not fftshift)
void   decal2DCplxGen2(vector<complex<double>> const &entree,vector<complex<double>>  &result, Var2D   &decalGen){
    size_t nbPix2D=entree.size();
    unsigned short int dim=sqrt(nbPix2D);
    Var2D const dim2D={dim,dim};
    decalGen.y=decalGen.y%dim2D.y;
    decalGen.x=decalGen.x%dim2D.x;
    size_t pixel=0,pixel_shift=0;
    if(decalGen.x<0)
      decalGen.x=dim2D.x+decalGen.x;
    if(decalGen.y<0)
      decalGen.y=dim2D.y+decalGen.y;
    Var2D const decal={decalGen.x,decalGen.y};
    unsigned int yi;
        for(yi=0; yi<dim2D.y-decal.y; yi++){
            copy(entree.begin()+yi*dim2D.x,   entree.begin()+yi*dim2D.x+dim2D.x-decal.x,   result.begin()+(yi+decal.y)*dim2D.x+decal.x);
            copy(entree.begin()+yi*dim2D.x+dim2D.x-decal.x,   entree.begin()+yi*dim2D.x+dim2D.x,   result.begin()+(yi+decal.y)*dim2D.x);
        }
        for(yi=dim2D.y-decal.y; yi<dim2D.y; yi++){
            copy(entree.begin()+yi*dim2D.x,   entree.begin()+yi*dim2D.x+dim2D.x-decal.x,   result.begin()+(-dim2D.y+yi+decal.y)*dim2D.x+decal.x);
            copy(entree.begin()+yi*dim2D.x+dim2D.x-decal.x,   entree.begin()+yi*dim2D.x+dim2D.x,   result.begin()+(-dim2D.y+yi+decal.y)*dim2D.x);
        }
}

///@parameters PosSpec : position of the specular beam
int coordSpec(vector<complex<double>> const &TF_UBorn, vector<double> &TF_champMod,Var2D NMAX)
 {
    int cpt_max=0;
    TF_champMod[0]=pow(TF_UBorn[0].real(),2)+pow(TF_UBorn[0].imag(),2);

    for(int cpt=1; cpt<(4*NMAX.x*NMAX.y); cpt++) {
       // TF_champMod[cpt]=sqrt(pow(TF_UBorn[cpt].real(),2)+pow(TF_UBorn[cpt].imag(),2));
        TF_champMod[cpt]=abs(TF_UBorn[cpt]);
        if(TF_champMod[cpt]>TF_champMod[cpt_max]) {
        cpt_max=cpt;
        }
    }
    return cpt_max;
 }
///crop src2D[0:dim_src,0:dim_src] to dest3D(coin.x:coin.x+dim_dest,coin.y+dim_dest), human
void coupeCplx(vector<complex<double>> const &src, vector<complex<double>> &dest, Var2D dim_src, Var2D dim_dest, Var2D coin, size_t NumAngle)
{
size_t nbPixSrc=src.size();
//Var2D dim_src={sqrt(nbPixSrc),sqrt(nbPixSrc)};
size_t X_dest,Y_dest, cpt_dest1D,
X_src, Y_src, cpt_src1D, cpt_Z_dest;

cpt_Z_dest=(dim_dest.x*dim_dest.y)*NumAngle;
        for(Y_dest=0; Y_dest<dim_dest.y; Y_dest++){
          size_t num_lgn= Y_dest*dim_dest.x;
          for(X_dest=0; X_dest<dim_dest.x; X_dest++){
            cpt_dest1D=cpt_Z_dest+X_dest+num_lgn;///coord 1D destination

            //coordonnées de découpe dans la source
            X_src=coin.x+X_dest;///coord X src
            Y_src=coin.y+Y_dest;///coord Y src
            cpt_src1D=X_src+Y_src*dim_src.x;///coord 1D source
            //copie src->dest
            dest[cpt_dest1D]=src[cpt_src1D];
            //dest[cpt_dest1D]->imag=src[cpt_src1D].imag;
        }
    }

}


