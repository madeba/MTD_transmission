
#include <vector>
#include <complex>
#include "fonctions.h"
#include "FFT_fonctions.h"
#include "src/vecteur.h"
#include "deroulement_volkov2.h"
#define M_2PI 2*M_PI
using namespace std;


void deroul_volkov2(vector<double> const &phase_enroul,vector<double> &phase_deroul , FFTW_init param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=phase_enroul.size();
  vector<complex<double>> gradx_enroul_fft(nbPix), grady_enroul_fft(nbPix);
  vector<complex<double>> Z(nbPix);
  vector<complex<double>> Gradx_Z_fft(nbPix), Grady_Z_fft(nbPix);
  vector<double> gradx_IntM(nbPix), grady_IntM(nbPix);
  //string repertoire_sav="/home/mat/tomo_test/";
  //SAV2_vec(phase_enroul,nbPix,repertoire_sav+"phase_enroul_dans_deroul_volkov.bin",t_float,"wb");

  gradient_fft2(phase_enroul, gradx_enroul_fft,grady_enroul_fft, param_c2c);
  //SAVCplx(gradx_enroul_fft,"Re",nbPix,"/home/mat/tomo_test/gradx_enroul_re.bin",t_float,"w+b");
  //SAVCplx(gradx_enroul_fft,"Im",nbPix,"/home/mat/tomo_test/gradx_enroul_Im.bin",t_float,"w+b");
  for(size_t cpt=0;cpt<nbPix;cpt++){
    Z[cpt].real(cos(phase_enroul[cpt]));
    Z[cpt].imag(sin(phase_enroul[cpt]));
  }
  gradient_fft2(Z, Gradx_Z_fft,Grady_Z_fft, param_c2c);
  //SAVCplx(TF_gradx,"Re",nbPix,"/home/mat/tomo_test/TF_gradx.bin",t_float,"w+b");
  /// Calcul du champ des entiers de déroulement
  complex<double> ax,ay;
  for(size_t cpt=0;cpt<nbPix;cpt++){
    ax=-I*(Gradx_Z_fft[cpt]/Z[cpt]);
    ay=-I*(Grady_Z_fft[cpt]/Z[cpt]);
    gradx_IntM[cpt]=(ax.real()-gradx_enroul_fft[cpt].real())/(2*M_PI);//*(-0.159);//-1/2pi
    grady_IntM[cpt]=(ay.real()-grady_enroul_fft[cpt].real())/(2*M_PI);
  }
  vector<complex<double>> IntM(nbPix);
  integ_grad2(gradx_IntM,grady_IntM,IntM,param_c2c);

  for(size_t cpt=0;cpt<nbPix;cpt++){
    phase_deroul[cpt]=phase_enroul[cpt]+2*3.1415*IntM[cpt].real();
  }

}


//surcharge avec entree <double>

void gradient_fft2(vector<double> const &entree, vector<complex<double>> &gradx, vector<complex<double>> &grady, FFTW_init &param_c2c)
{
     // string repertoire_sav="/home/mat/tomo_test/";
     complex<double> I(0,1);
     unsigned int nbPix=entree.size();
     unsigned int dim=sqrt(nbPix);
     vector<complex<double>> tamponx(nbPix);
     vector<complex<double>> tampony(nbPix);
     vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);
     //SAV2(entree,nbPix,"/home/mat/tomo_test/phase.bin",t_float,"w+b");
    // TF2D_vec(in,out, entree, spectre, p_forward);
    TF2Dcplx(entree,spectre,param_c2c);

     vector<vecteur> kvect(nbPix),kvect_shift(nbPix);

     //spectre_shift=fftshift2D(spectre);
    // SAVCplx(spectre_shift,"Im",nbPix,repertoire_sav+"spectre_pyramide_shift.bin",t_float,"w+b");
     for(size_t cpt=0;cpt<nbPix;cpt++)
     {
       kvect[cpt].setx((cpt%dim-round(dim/2))/(dim));
       kvect[cpt].sety((cpt/dim-round(dim/2))/(dim));
     }

     kvect_shift=fftshift2D(kvect);
     //SAV_vec3D(kvect, "x",repertoire_sav+"kx.bin","w+b",nbPix);
     //SAV_vec3D(kvect_shift, "x",repertoire_sav+"kx_shift.bin","w+b",nbPix);

     for(size_t cpt=0;cpt<nbPix;cpt++){
        tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
        tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
    }
    TF2Dcplx_INV(tamponx, gradx, param_c2c);
    TF2Dcplx_INV(tampony, grady, param_c2c);
}

void gradient_fft2(vector<complex<double>> const&entree, vector<complex<double>> &gradx, vector<complex<double>> &grady,FFTW_init &param_c2c)
{
// string repertoire_sav="/home/mat/";
  complex<double> I(0,1);
  unsigned int nbPix=entree.size();
  unsigned int dim=sqrt(nbPix);
  vector<complex<double>> tamponx(nbPix);
  vector<complex<double>> tampony(nbPix);
  vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);
  //SAV2(entree,nbPix,"/home/mat/tomo_test/phase.bin",t_float,"w+b");
  TF2Dcplx(entree, spectre, param_c2c);
  vector<vecteur> kvect(nbPix),kvect_shift(nbPix);
  //SAVCplx(spectre,"Re",nbPix,repertoire_sav+"spectre_pyramide.bin",t_float,"w+b");
  //spectre_shift=fftshift2D(spectre);
  //SAVCplx(spectre_shift,"Im",nbPix,repertoire_sav+"spectre_pyramide_shift.bin",t_float,"w+b");
  for(size_t cpt=0;cpt<nbPix;cpt++){
      kvect[cpt].setx((cpt%dim-round(dim/2))/(dim));
      kvect[cpt].sety((cpt/dim-round(dim/2))/(dim));
  }
  kvect_shift=fftshift2D(kvect);
  //SAV_vec3D(kvect, "x",repertoire_sav+"kx.bin","w+b",nbPix);
  //SAV_vec3D(kvect_shift, "x",repertoire_sav+"kx_shift.bin","w+b",nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
    tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
    tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
  }

  TF2Dcplx_INV(tamponx, gradx, param_c2c);
  TF2Dcplx_INV(tampony, grady, param_c2c);
}

void integ_grad2(vector<double> const &gradx, vector<double> const& grady, vector<complex<double>> &sortie,FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=gradx.size(),dim=sqrt(nbPix);
  vector<complex<double>> TF_gradx(nbPix), TF_grady(nbPix), tampon(nbPix);
  vector<double> kvect_mod_sq(nbPix);
  vector<vecteur> kvect(nbPix),kvect_shift(nbPix);

  TF2Dcplx(gradx, TF_gradx, param_c2c);
  TF2Dcplx(grady, TF_grady, param_c2c);
  //SAVCplx(TF_gradx,"Re","/home/mat/tomo_test/TF_gradx_re.bin",t_float,"a+b");
  //SAVCplx(TF_gradx,"Im","/home/mat/tomo_test/TF_gradx_im.bin",t_float,"a+b");

  //ne pas recalculer à chaque fois?
  for(size_t cpt=0;cpt<nbPix;cpt++)
  {
   kvect[cpt].setx((cpt%dim-round(dim/2))/(dim));
   kvect[cpt].sety((cpt/dim-round(dim/2))/(dim));
  }
  kvect_shift=fftshift2D(kvect);
  //SAV_vec3D(kvect_shift,"x","/home/mat/tomo_test/kvect_shift_x.bin","wb",nbPix);
  //SAV_vec3D(kvect_shift,"y","/home/mat/tomo_test/kvect_shift_y.bin","wb",nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
    kvect_mod_sq[cpt]=pow(kvect_shift[cpt].getx(),2)+pow(kvect_shift[cpt].gety(),2);
    //if(isfinite(kvect_mod_sq[cpt])==0)
        //cout<<"division par zéro="<<kvect_mod_sq[cpt]<<endl;
    tampon[cpt]=-I*(TF_gradx[cpt]*kvect_shift[cpt].getx()+TF_grady[cpt]*kvect_shift[cpt].gety())/(2*3.14159*kvect_mod_sq[cpt]);
  }
  //kvect_mod_sq[0]=(kvect_mod_sq[nbPix]+kvect_mod_sq[1])/2;
  //cout<<"kvect_mod[0]"<<kvect_mod_sq[0]<<endl;
  //SAV2(kvect_mod_sq,nbPix,"/home/mat/tomo_test/kvect_mod.bin",t_float,"wb");
  //SAVCplx(tampon,"Re",nbPix,"/home/mat/tomo_test/tampon_Re.bin",t_float,"w+b");
  tampon[0]=0;//(tampon[nbPix]+tampon[1])*0.5;//moyennage de la fréquence zéro pour éviter NAN
  TF2Dcplx_INV(tampon,sortie,param_c2c);
}
