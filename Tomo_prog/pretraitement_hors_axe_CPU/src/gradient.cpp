#include "gradient.h"
#include "FFTW_init.h"
#include "FFT_fonctions.h"
using namespace std;
///gradient par fft, entrée réelle
///calculate gradient by fft, real input
void gradient_fft4(vector<double> &entree, vector<complex<double>> &gradx, vector<complex<double>> &grady, vector<vecteur> &kvect_shift, FFTW_init &param_c2c)
{
     complex<double> I(0,1);
     unsigned int nbPix=entree.size();
//     unsigned int dim=sqrt(nbPix);
     vector<complex<double>> tamponx(nbPix), tampony(nbPix);
     vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);
     //SAV2(entree,nbPix,"/home/mat/tomo_test/phase.bin",t_float,"w+b");
    // TF2D_vec(in,out, entree, spectre, p_forward);

    TF2Dcplx((entree),spectre,param_c2c);
     for(size_t cpt=0;cpt<nbPix;cpt++){
        tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
        tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
    }
    TF2Dcplx_INV(tamponx, gradx, param_c2c);
    TF2Dcplx_INV(tampony, grady, param_c2c);
}


///gradient par fft, entrée complexe
///calculate gradient by fft, complex input (overloaded function)
void gradient_fft4(vector<complex<double>>   &entree, vector<complex<double>> &gradx, vector<complex<double>> &grady,vector<vecteur>  &kvect_shift, FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=entree.size();
  //unsigned int dim=sqrt(nbPix);
  vector<complex<double>> tamponx(nbPix);
  vector<complex<double>> tampony(nbPix);
  vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);

  TF2Dcplx((entree), spectre, param_c2c);

  for(size_t cpt=0;cpt<nbPix;cpt++){
    tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
    tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
  }
  TF2Dcplx_INV(tamponx, gradx, param_c2c);
  TF2Dcplx_INV(tampony, grady, param_c2c);

}



///-----///////////////////////calcul du gradient par simple difference dans l'espace image; pour comparaison---------------------
void gradient_central(const std::vector<double> src, std::vector<double> &grad,string direction)
{
const unsigned int dim=sqrt(src.size());
//const unsigned int nbPix=dim*dim;
//vector<double> decalx(nbPix);
//vector<double> decaly(nbPix);
if(direction!="x" && direction !="y")
        cout<<"direction de gradient inconnue. Choix possibles : x ou y"<<endl;
size_t cpt,cpt_decal_avt,cpt_decal_arr;
    if(direction=="x")
    {
        for(size_t x=1;x<dim-1;x++)
            for(size_t y=0;y<dim;y++){
                cpt=x+y*dim;
                cpt_decal_avt=x+1+y*dim;
                cpt_decal_arr=x-1+y*dim;
                grad[cpt]=(src[cpt_decal_avt]-src[cpt_decal_arr])/2;
    }
    }

    if(direction=="y"){
      for(size_t x=0;x<dim;x++)
        for(size_t y=1;y<dim-1;y++){
             cpt=x+y*dim;
             cpt_decal_avt=x+y*dim+dim;
             cpt_decal_arr=x+y*dim-dim;
            grad[cpt]=(src[cpt_decal_avt]-src[cpt_decal_arr])/2;
            }
    }
}

void gradient_central(const std::vector<complex<double>> src, std::vector<complex<double>> &grad,string direction)
{
const unsigned int dim=sqrt(src.size());
//const unsigned int nbPix=dim*dim;
//vector<double> decalx(nbPix);
//vector<double> decaly(nbPix);
if(direction!="x" && direction !="y")
        cout<<"direction de gradient inconnue. Choix possibles : x ou y"<<endl;
size_t cpt,cpt_decal_avt,cpt_decal_arr;
    if(direction=="x")
    {
        for(size_t x=1;x<dim-1;x++)
            for(size_t y=0;y<dim;y++){
                cpt=x+y*dim;
                cpt_decal_avt=x+1+y*dim;
                cpt_decal_arr=x-1+y*dim;
                grad[cpt].real((src[cpt_decal_avt].real()-src[cpt_decal_arr].real())/2);
                grad[cpt].imag((src[cpt_decal_avt].imag()-src[cpt_decal_arr].imag())/2);
    }
    }

    if(direction=="y"){
      for(size_t x=0;x<dim;x++)
        for(size_t y=1;y<dim-1;y++){
             cpt=x+y*dim;
             cpt_decal_avt=x+y*dim+dim;
             cpt_decal_arr=x+y*dim-dim;
            grad[cpt].real((src[cpt_decal_avt].real()-src[cpt_decal_arr].real())/2);
            grad[cpt].imag((src[cpt_decal_avt].imag()-src[cpt_decal_arr].imag())/2);
            }
    }
}

void gradient_back(const std::vector<complex<double>> src, std::vector<complex<double>> &grad,string direction)
{
const unsigned int dim=sqrt(src.size());
//const unsigned int nbPix=dim*dim;
//vector<double> decalx(nbPix);
//vector<double> decaly(nbPix);
if(direction!="x" && direction !="y")
        cout<<"direction de gradient inconnue. Choix possibles : x ou y"<<endl;
size_t cpt,cpt_decal_avt,cpt_decal_arr;
    if(direction=="x")
    {
        for(size_t x=1;x<dim-1;x++)
            for(size_t y=0;y<dim;y++){
                cpt=x+y*dim;
                cpt_decal_avt=x+y*dim;
                cpt_decal_arr=x-1+y*dim;
                grad[cpt].real(src[cpt_decal_avt].real()-src[cpt_decal_arr].real());
                grad[cpt].imag(src[cpt_decal_avt].imag()-src[cpt_decal_arr].imag());
    }
    }

    if(direction=="y"){
      for(size_t x=0;x<dim;x++)
        for(size_t y=1;y<dim-1;y++){
             cpt=x+y*dim;
             cpt_decal_avt=x+y*dim;
             cpt_decal_arr=x+y*dim-dim;
            grad[cpt].real(src[cpt_decal_avt].real()-src[cpt_decal_arr].real());
            grad[cpt].imag(src[cpt_decal_avt].imag()-src[cpt_decal_arr].imag());
            }
    }
}


void gradient_back(const std::vector<double> src, std::vector<double> &grad,string direction)
{
const unsigned int dim=sqrt(src.size());
//const unsigned int nbPix=dim*dim;
//vector<double> decalx(nbPix);
//vector<double> decaly(nbPix);
if(direction!="x" && direction !="y")
        cout<<"direction de gradient inconnue. Choix possibles : x ou y"<<endl;
size_t cpt,cpt_decal_avt,cpt_decal_arr;
    if(direction=="x")
    {
        for(size_t x=1;x<dim-1;x++)
            for(size_t y=0;y<dim;y++){
                cpt=x+y*dim;
                cpt_decal_avt=x+y*dim;
                cpt_decal_arr=x-1+y*dim;
                grad[cpt]=(src[cpt_decal_avt]-src[cpt_decal_arr]);

    }
    }

    if(direction=="y"){
      for(size_t x=0;x<dim;x++)
        for(size_t y=1;y<dim-1;y++){
             cpt=x+y*dim;
             cpt_decal_avt=x+y*dim;
             cpt_decal_arr=x+y*dim-dim;
            grad[cpt]=(src[cpt_decal_avt]-src[cpt_decal_arr]);

            }
    }
}
