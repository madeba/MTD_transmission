
//#include <vector>
//#include <complex>
#include "fonctions.h"
//#include "FFT_fonctions.h"
//#include "vecteur.h"
#include "deroulement_volkov3_GPU.h"
#include "GPU_fonctions.h"
#include <chrono>
//#define M_2PI 2*M_PI

using namespace std;

///main function
void deroul_volkov3_GPU2(af::array const &phase_enroul, af::array  &phase_deroul,af::array  const &kvect_shiftX,af::array const &kvect_shiftY)
{   unsigned int dim=phase_enroul.dims(0);
    af::cdouble i_cdouble = { 0, 1 };//cdouble, creates an object of class std::complex, its use does not create an actual ArrayFire object.
    af::array i = constant(i_cdouble,dim,dim, c64);/* imaginary unit */

///array complex double
    af::array gradx_enroul_fft(dim,dim,c64), grady_enroul_fft(dim,dim,c64),  Gradx_Z_fft(dim,dim,c64), Grady_Z_fft(dim,dim,c64);
    ///array double
    gradient_fft3_GPU2(phase_enroul, gradx_enroul_fft,grady_enroul_fft, kvect_shiftX, kvect_shiftY);

   // SAV_Tiff2D(real(grady_enroul_fft),"/home/mat/tomo_test/grady_enroul_fft_GPU.tif",1);
  //   SAV_Tiff2D(real(gradx_enroul_fft),"/home/mat/tomo_test/gradx_enroul_fft_GPU.tif",1);

    // display(imag(gradx_enroul_fft));
    af::array Z=af::complex(cos(phase_enroul),sin(phase_enroul));

    gradient_fft3_GPU2(Z, Gradx_Z_fft,Grady_Z_fft, kvect_shiftX,kvect_shiftY);
    // SAV_Tiff2D(real(Gradx_Z_fft),"/home/mat/tomo_test/gradx_Z_fft_GPU.tif",1);
      //    SAV_Tiff2D(real(Grady_Z_fft),"/home/mat/tomo_test/grady_Z_fft_GPU.tif",1);

    //SAVCplx(TF_gradx,"Re",nbPix,"/home/mat/tomo_test/TF_gradx.bin",t_float,"w+b");
    /// Calcul du champ des entiers de déroulement
    af::array ax=-i*(Gradx_Z_fft/Z);///c64
    af::array ay=-i*(Grady_Z_fft/Z);//c64
    // af::array gradx_IntM(dim,dim,f64), grady_IntM(dim,dim,f64);
    af::array gradx_IntM=(af::real(ax)-real(gradx_enroul_fft))/(2*af::Pi);  //f64
    af::array grady_IntM=(af::real(ay)-real(grady_enroul_fft))/(2*af::Pi); //f64
    //SAV_Tiff2D(real(gradx_IntM),"/home/mat/tomo_test/gradx_IntM_GPU.tif",1);
       // SAV_Tiff2D(real(grady_IntM),"/home/mat/tomo_test/grady_IntM_GPU.tif",1);
    af::array IntM=af::array(dim,dim, c64);
    integ_grad3_GPU2(gradx_IntM,grady_IntM,IntM,kvect_shiftX,kvect_shiftY);
    ///calcul phase dépliée
    phase_deroul=phase_enroul+2*af::Pi*(real(IntM));
//SAV_Tiff2D(real(IntM),"/home/mat/tomo_test/Int_M_GPU.tiff",1);
  }


///gradient par fft, pas de typage pour l'entrée car array
void gradient_fft3_GPU2(af::array const &entree, af::array &gradx, af::array &grady, af::array const &kvect_shiftX,af::array const &kvect_shiftY) {
af::cdouble i_cdouble = { 0, 1 };//cdouble, creates an object of class std::complex, its use does not create an actual ArrayFire object.
af::array i = constant(i_cdouble,1, c64);/* imaginary unit */

af::array spectre=fft2(entree);
af::array tamponX=2*af::Pi*kvect_shiftX*spectre*i_cdouble;
af::array tamponY=2*af::Pi*kvect_shiftY*spectre*i_cdouble;
//cout<<"coucou"<<endl;
gradx=af::ifft2(tamponX);
grady=af::ifft2(tamponY);
}


///intégration du gradient par fft
void integ_grad3_GPU2(af::array const &gradx, af::array const& grady, af::array &sortie,af::array const &kvect_shiftX, af::array const &kvect_shiftY)
{
    af::cdouble i_cdouble = { 0, 1 };//cdouble, creates an object of class std::complex, its use does not create an actual ArrayFire object.
    af::array i = constant(i_cdouble,1, c64);/* imaginary unit */
    int nbPix=gradx.elements();

    af::array  TF_gradx=af::fft2(gradx), TF_grady=af::fft2(grady), kvect_mod_sq=pow(kvect_shiftX,2)+pow(kvect_shiftY,2);

  //  SAV_Tiff2D(kvect_mod_sq,"/home/mat/tomo_test/kvect_mod_GPU.tif",1);

   af::array tampon=-1*i_cdouble*(TF_gradx*kvect_shiftX+TF_grady*kvect_shiftY)/(2*af::Pi*kvect_mod_sq);
 //  SAV_Tiff2D(kvect_mod_sq,"/home/mat/tomo_test/kvect_mod_sq_GPU.tif",1);
   //   SAV_Tiff2D(imag(tampon),"/home/mat/tomo_test/tampon_imag_GPU.tif",1);
    //SAV_vec3D(kvect_shift,"x","/home/mat/tomo_test/kvect_shift_x.bin","wb",nbPix);
    //SAV_vec3D(kvect_shift,"y","/home/mat/tomo_test/kvect_shift_y.bin","wb",nbPix);
    //kvect_mod_sq[0]=(kvect_mod_sq[nbPix]+kvect_mod_sq[1])/2;
    //cout<<"kvect_mod[0]"<<kvect_mod_sq[0]<<endl;
    //SAV2(kvect_mod_sq,nbPix,"/home/mat/tomo_test/kvect_mod.bin",t_float,"wb");
    //SAVCplx(tampon,"Re",nbPix,"/home/mat/tomo_test/tampon_Re.bin",t_float,"w+b");
    tampon(0,0)=0;
   // tampon(tampon.dims(0)/2,tampon.dims(0)/2)=0;//(tampon[nbPix]+tampon[1])*0.5;//moyennage de la fréquence zéro pour éviter NAN
    sortie=ifft2(tampon);
    // SAV_Tiff2D(real(tampon),"/home/mat/tomo_test/tampon_real_GPU.tif",1);
}
/*
af::array init_kvect_shiftX_GPU(Var2D dim2DHA)
{
 af::array  kvect_shiftX=af::array(dim2DHA.x,dim2DHA.y);
 kvect_shiftX(af::seq(0,dim2DHA.x))=af::seq(0,dim2DHA.x);
 display(kvect_shiftX);
 return kvect_shiftX;
}*/


