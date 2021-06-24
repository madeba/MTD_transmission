
#include <vector>
#include <complex>
#include "fonctions.h"
#include "FFT_fonctions.h"
#include "vecteur.h"
#include "deroulement_volkov3.h"
#include <chrono>
#define M_2PI 2*M_PI
using namespace std;

///main function
void deroul_volkov3(vector<double> const &phase_enroul,vector<double> &phase_deroul,vector<vecteur>  &kvect_shift,  FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=phase_enroul.size();
  vector<complex<double>> gradx_enroul_fft(nbPix), grady_enroul_fft(nbPix);
  vector<complex<double>> Z(nbPix);
  vector<complex<double>> Gradx_Z_fft(nbPix), Grady_Z_fft(nbPix);
  vector<double> gradx_IntM(nbPix), grady_IntM(nbPix);

  gradient_fft3(phase_enroul, gradx_enroul_fft,grady_enroul_fft, kvect_shift,param_c2c);

  for(size_t cpt=0;cpt<nbPix;cpt++){
    Z[cpt].real(cos(phase_enroul[cpt]));
    Z[cpt].imag(sin(phase_enroul[cpt]));
  }
  gradient_fft3(Z, Gradx_Z_fft,Grady_Z_fft, kvect_shift, param_c2c);
  //SAVCplx(TF_gradx,"Re",nbPix,"/home/mat/tomo_test/TF_gradx.bin",t_float,"w+b");
  /// Calcul du champ des entiers de déroulement
  complex<double> ax,ay;
  for(size_t cpt=0;cpt<nbPix;cpt++){
    ax=-I*(Gradx_Z_fft[cpt]/Z[cpt]);
    ay=-I*(Grady_Z_fft[cpt]/Z[cpt]);
    gradx_IntM[cpt]=(ax.real()-gradx_enroul_fft[cpt].real())/(M_2PI);//*(-0.159);//-1/2pi
    grady_IntM[cpt]=(ay.real()-grady_enroul_fft[cpt].real())/(M_2PI);
  }
  vector<complex<double>> IntM(nbPix);
  integ_grad3(gradx_IntM,grady_IntM,IntM,kvect_shift,param_c2c);
  ///calcul phase dépliée
  for(size_t cpt=0;cpt<nbPix;cpt++){
    phase_deroul[cpt]=phase_enroul[cpt]+M_2PI*IntM[cpt].real();
  }
}

///gradient par fft, entrée réelle
void gradient_fft3(vector<double> const &entree, vector<complex<double>> &gradx, vector<complex<double>> &grady, vector<vecteur> &kvect_shift, FFTW_init &param_c2c)
{
     complex<double> I(0,1);
     unsigned int nbPix=entree.size();
     unsigned int dim=sqrt(nbPix);
     vector<complex<double>> tamponx(nbPix), tampony(nbPix);
     vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);
     //SAV2(entree,nbPix,"/home/mat/tomo_test/phase.bin",t_float,"w+b");
    // TF2D_vec(in,out, entree, spectre, p_forward);
    TF2Dcplx(entree,spectre,param_c2c);

     for(size_t cpt=0;cpt<nbPix;cpt++){
        tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
        tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
    }
    TF2Dcplx_INV(tamponx, gradx, param_c2c);
    TF2Dcplx_INV(tampony, grady, param_c2c);
}
///gradient par fft, entrée complexe
void gradient_fft3(vector<complex<double>> const&entree, vector<complex<double>> &gradx, vector<complex<double>> &grady,vector<vecteur>  &kvect_shift, FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=entree.size();
  unsigned int dim=sqrt(nbPix);
  vector<complex<double>> tamponx(nbPix);
  vector<complex<double>> tampony(nbPix);
  vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);

  TF2Dcplx(entree, spectre, param_c2c);

  for(size_t cpt=0;cpt<nbPix;cpt++){
    tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
    tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
  }

  TF2Dcplx_INV(tamponx, gradx, param_c2c);
  TF2Dcplx_INV(tampony, grady, param_c2c);
}
///intégration du gradient par fft
void integ_grad3(vector<double> const &gradx, vector<double> const& grady, vector<complex<double>> &sortie,std::vector<vecteur> &kvect_shift, FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=gradx.size(),dim=sqrt(nbPix);
  vector<complex<double>> TF_gradx(nbPix), TF_grady(nbPix), tampon(nbPix);
  vector<double> kvect_mod_sq(nbPix);
//  vector<vecteur> kvect(nbPix),kvect_shift(nbPix);

  TF2Dcplx(gradx, TF_gradx, param_c2c);
  TF2Dcplx(grady, TF_grady, param_c2c);

  //SAV_vec3D(kvect_shift,"x","/home/mat/tomo_test/kvect_shift_x.bin","wb",nbPix);
  //SAV_vec3D(kvect_shift,"y","/home/mat/tomo_test/kvect_shift_y.bin","wb",nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
    kvect_mod_sq[cpt]=pow(kvect_shift[cpt].getx(),2)+pow(kvect_shift[cpt].gety(),2);
    //if(isfinite(kvect_mod_sq[cpt])==0)
        //cout<<"division par zéro="<<kvect_mod_sq[cpt]<<endl;
    tampon[cpt]=-I*(TF_gradx[cpt]*kvect_shift[cpt].getx()+TF_grady[cpt]*kvect_shift[cpt].gety())/(M_2PI*kvect_mod_sq[cpt]);
  }
  //kvect_mod_sq[0]=(kvect_mod_sq[nbPix]+kvect_mod_sq[1])/2;
  //cout<<"kvect_mod[0]"<<kvect_mod_sq[0]<<endl;
  //SAV2(kvect_mod_sq,nbPix,"/home/mat/tomo_test/kvect_mod.bin",t_float,"wb");
  //SAVCplx(tampon,"Re",nbPix,"/home/mat/tomo_test/tampon_Re.bin",t_float,"w+b");
  tampon[0]=0;//(tampon[nbPix]+tampon[1])*0.5;//moyennage de la fréquence zéro pour éviter NAN
  TF2Dcplx_INV(tampon,sortie,param_c2c);
}
///calculate the kvector field (array whose value a simply  kx and ky)
std::vector<vecteur> init_kvect_shift(Var2D dim2DHA)
{
  size_t nbPix=dim2DHA.x*dim2DHA.y;
  vector<vecteur> kvect(nbPix),kvect_shift(nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
    kvect[cpt].setx((cpt%dim2DHA.x-round(dim2DHA.x/2))/(dim2DHA.x));
    kvect[cpt].sety((cpt/dim2DHA.y-round(dim2DHA.y/2))/(dim2DHA.y));
  }
  kvect_shift=fftshift2D(kvect);
  return kvect_shift;
}
///calculate the kvector square modulus field (array whose value a simply  kx^2+ky^2)
std::vector<double> init_kvect_mod2Shift(vector<vecteur>  &kvect_shift)
{
  size_t nbPix=kvect_shift.size();
  vector<double> kvect_mod_sq_shift(nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
        kvect_mod_sq_shift[cpt]=pow(kvect_shift[cpt].getx(),2)+pow(kvect_shift[cpt].gety(),2);
  }

  return kvect_mod_sq_shift;
}
