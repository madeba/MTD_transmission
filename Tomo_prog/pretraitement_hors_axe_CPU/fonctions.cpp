#include <fstream>//ifstream
#include <vector>
#include <chrono>
#include <complex>
#include <fftw3.h>
#include "struct.h"
//#include <cv.h>
//#include <highgui.h>//imread
#include "projet.h"
#include "FFT_fonctions.h"
#include "fonctions.h"
#include "IO_fonctions.h"
using namespace std;
using namespace cv;
///free functions used outside classes

void sav_param2D(string texte,string file_path){
    ///open an ofstream to save preprocessing informations
    ofstream fichier_sav_parametre;
    fichier_sav_parametre.open(file_path, std::ios_base::app);//append to file
    fichier_sav_parametre<<texte<<endl;
    fichier_sav_parametre.close();
}
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

///coupe dans le repère informatique puis copie vers pile 3D.
///crop on a fftshifted spectru, then copy to 3D stack
//fonction très particulière à n'utiliser qu'après un fftshift, typiquement sur un spectre centre en repère informatique
/*
void coupe2D_RefI_to_I3D(vector<complex<double>> const &src, vector<complex<double>> &dest, Var2D dim_dest, unsigned short int numAngle)
{
size_t nbPixSrc=src.size();

Var2D dim_src={sqrt(nbPixSrc),sqrt(nbPixSrc)};
//cout<<"dim_src.x="<<dim_src.x<<endl;
//Var2D dim_dest={sqrt(nbPixDest),sqrt(nbPixDest)};
Var2D Nmax={dim_dest.x/2,dim_dest.y/2};
//cout<<"Nmax="<<Nmax.x<<endl;
unsigned short int X_dest,Y_dest, X_src, Y_src;
size_t cpt_src,cpt_dest,Z_dest;

Z_dest=(dim_dest.x*dim_dest.y)*numAngle;
//partie haute (A&B)
for(Y_src=0;Y_src<Nmax.y;Y_src++){
  Y_dest=Y_src;///coord Y src
  size_t num_ligne_dest=Y_dest*dim_dest.x,
                     num_ligne_src=Y_src*dim_src.x;
  for(X_src=0;X_src<Nmax.y;X_src++){
     X_dest=X_src;///coord X src

     cpt_src=X_src+num_ligne_src;///coord 1D source
     cpt_dest=X_dest+num_ligne_dest+Z_dest;
     dest[cpt_dest]=src[cpt_src];
  }
  for(X_src=dim_src.x-Nmax.x;X_src<dim_src.x;X_src++){
     X_dest=X_src-(dim_src.x-2*Nmax.x);///coord X src
     cpt_src=X_src+num_ligne_src;///coord 1D source
     cpt_dest=X_dest+num_ligne_dest+Z_dest;
     dest[cpt_dest]=src[cpt_src];
  }
}
  //partie basse (C&D)
 for(Y_src=dim_src.y-Nmax.y;Y_src<dim_src.y;Y_src++){
    Y_dest=Y_src-(dim_src.y-2*Nmax.y);///coord X src
      size_t num_ligne_dest=Y_dest*dim_dest.x,
                     num_ligne_src=Y_src*dim_src.x;
   for(X_src=0;X_src<Nmax.y;X_src++){
     X_dest=X_src;///coord X src
     cpt_src=X_src+num_ligne_src;///coord 1D source
     cpt_dest=X_dest+num_ligne_dest+Z_dest;
     dest[cpt_dest]=src[cpt_src];
  }
  for(X_src=dim_src.x-Nmax.x;X_src<dim_src.x;X_src++){
     X_dest=X_src-(dim_src.x-2*Nmax.x);///coord X src
     cpt_src=X_src+num_ligne_src;///coord 1D source
     cpt_dest=X_dest+num_ligne_dest+Z_dest;
     dest[cpt_dest]=src[cpt_src];
   }
  }
}*/


///crop dans le repère informatique vers repère humain . Crop src into dest, human-centered (zero=middle of the image)
//fonction très particulière à n'utiliser qu'après un fftshift, typiquement sur un spectre centre en repère informatique
void coupe2D_I_to_H3D(vector<complex<double>> const &src2D, vector<complex<double>> &dest3D,Var2D dim_dest2D, unsigned short int numAngle)
{
  size_t nbPixSrc=src2D.size();
  Var2D dim_src={sqrt(nbPixSrc),sqrt(nbPixSrc)};
  Var2D Nmax={dim_dest2D.x/2,dim_dest2D.y/2};
  unsigned short int X_dest,Y_dest, X_src, Y_src;
  size_t cpt_src,cpt_dest,Z_dest;

  Z_dest=(dim_dest2D.x*dim_dest2D.y)*numAngle;
  //source partie haute (A&B) ver partie basse dest
  for(Y_src=0;Y_src<Nmax.y;Y_src++){
    Y_dest=Y_src+Nmax.x;
    size_t num_ligne_dest=Y_dest*dim_dest2D.x,  num_ligne_src=Y_src*dim_src.x;

    for(X_src=0;X_src<Nmax.x;X_src++){//A ver sA'
       X_dest=X_src+Nmax.x;///coord X src
       cpt_src=X_src+num_ligne_src;///coord 1D source
       cpt_dest=X_dest+num_ligne_dest+Z_dest;
       dest3D[cpt_dest]=src2D[cpt_src];
    }
    for(X_src=dim_src.x-Nmax.x;X_src<dim_src.x;X_src++){//B vers B'
       X_dest=X_src-(dim_src.x-Nmax.x);///coord X src
       cpt_src=X_src+num_ligne_src;///coord 1D source
       cpt_dest=X_dest+num_ligne_dest+Z_dest;
       dest3D[cpt_dest]=src2D[cpt_src];
    }
  }
  for(Y_src=dim_src.y-Nmax.y;Y_src<dim_src.y;Y_src++){
    Y_dest=Y_src-(dim_src.y-Nmax.y);
    size_t num_ligne_dest=Y_dest*dim_dest2D.x,
                     num_ligne_src=Y_src*dim_src.x;
    for(X_src=0;X_src<Nmax.x;X_src++){//D vers D'
        X_dest=X_src+Nmax.x;
        cpt_src=X_src+num_ligne_src;///coord 1D source
        cpt_dest=X_dest+num_ligne_dest+Z_dest;
        dest3D[cpt_dest]=src2D[cpt_src];
    }
    for(X_src=dim_src.x-Nmax.x;X_src<dim_src.x;X_src++){//C vers C'
      X_dest=X_src-(dim_src.x-Nmax.x);
      cpt_src=X_src+num_ligne_src;///coord 1D source
      cpt_dest=X_dest+num_ligne_dest+Z_dest;
      dest3D[cpt_dest]=src2D[cpt_src];
    }
  }
}

///r2c symetric to 2D, the hologram is fftshifted, but the spectrum is not inverse-fftshifted. The  shifted spectrum is  (cropped @ coin_shifted and send to stack) by the function coupeCplx.
void holo2TF_UBorn2_shift(vector<double>  &holo1,vector<complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA_shift, size_t NbAngleOk, vector<double> const &tukeyHolo,FFTW_init  &param_fftw2DHolo)
{
    size_t NbPixROI2d=holo1.size();
    vector<complex<double>> TF_Holo(NbPixROI2d);
    for(size_t pixel=0; pixel<NbPixROI2d; pixel++)
      holo1[pixel]=(double)holo1[pixel]*tukeyHolo[pixel];

    TF2D_r2c_symetric(fftshift2D2(holo1),TF_Holo,param_fftw2DHolo);
//SAVCplx(TF_Holo,"Im","/home/mat/tmp/Tfholo_1024x1024x599x32.bin",t_float,"a+b");

    coupeCplx(TF_Holo, TF_UBornTot, dimROI, dim2DHA, coinHA_shift, NbAngleOk);///Découpe à [-Nxmax,+NXmax]dans repère humain-lisible +envoi dans pile3D
 //   SAVCplx(TF,"Im","/home/mat/tmp/Tfholo_220x220x60.bin",t_float,"a+b");
}

///r2c non symmetrized to 3D stack, fastest method
void holo2TF_UBorn2_shift_r2c(vector<double>  &holo1,vector<complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA_shift, size_t NbAngleOk, vector<double> const &tukeyHolo,FFTW_init  &param_fftw2D_r2c_Holo)
{
    size_t NbPixROI2d=holo1.size();
    vector<complex<double>> TF_Holo(NbPixROI2d);
    for(size_t pixel=0; pixel<NbPixROI2d; pixel++)
      holo1[pixel]=(double)holo1[pixel]*tukeyHolo[pixel];

   // TF2D_r2c_symetric(fftshift2D2(holo1),TF_Holo,param_fftw2DHolo);
     TF2D_r2c_coupeHA_to_stack(fftshift2D2(holo1), TF_UBornTot, dim2DHA, coinHA_shift,  NbAngleOk, param_fftw2D_r2c_Holo);///warning, fftshift for the 1st argument
   // SAVCplx(TF_Holo,"Im","/home/mat/tmp/Tfholo_1024x1024.bin",t_float,"a+b");
   // coupeCplx(TF_Holo, TF_UBornTot, dimROI, dim2DHA, coinHA_shift, NbAngleOk);///Découpe à [-Nxmax,+NXmax]dans repère humain-lisible +envoi dans pile3D
 //   SAVCplx(TF,"Im","/home/mat/tmp/Tfholo_220x220x60.bin",t_float,"a+b");
}

///r2c but symmetrized.
///Cons : The input hologram must be fftshifted,  then the spectrum must be back-fftshifted. Finally the image is  cropped and send to the stack of complex fields  by the function "CoupeCplx"
///pros : slower but easier to understand, because the fft is complete.
void holo2TF_UBorn2(vector<double>  &holo1,vector<complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, vector<double> const &tukeyHolo,FFTW_init  &param_fftw2DHolo)
{
    size_t NbPixROI2d=holo1.size();
    vector<complex<double>> TF_Holo(NbPixROI2d);
    for(size_t pixel=0; pixel<NbPixROI2d; pixel++)
      holo1[pixel]=(double)holo1[pixel]*tukeyHolo[pixel];

    TF2D_r2c_symetric(fftshift2D2(holo1),TF_Holo,param_fftw2DHolo);

    coupeCplx(fftshift2D2(TF_Holo), TF_UBornTot, dimROI, dim2DHA, coinHA, NbAngleOk);///Découpe à [-Nxmax,+NXmax]dans repère humain-lisible +envoi dans pile3D
   // coupeCplx(TF_Holo, TF_UBornTot, dimROI, dim2DHA, coinHA, NbAngleOk);///Découpe à [-Nxmax,+NXmax]dans repère humain-lisible +envoi dans pile3D
   //   SAVCplx(TF,"Im","/home/mat/tmp/Tfholo_220x220x60.bin",t_float,"a+b");
}

///ancienne fonction, lente, avec plan calculé à l'extérieur
void holo2TF_UBorn(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukey_holo, fftw_complex *in,fftw_complex *out,fftw_plan p_forward_holo)
{
    ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
    size_t dimx=sqrt(NbPix2dROI);

    size_t NbPixROI2d=holo1.size();
    vector<double> holo_shift(NbPixROI2d);
    vector<complex<double>> TF_Holo(NbPixROI2d);
    vector<complex<double>> TFHoloCentre(NbPixROI2d);
    nbCplx *TF_UBorn_A=new nbCplx[NbPixROI2d];

    for(size_t pixel=0; pixel<NbPixROI2d; pixel++){ holo1[pixel]=(double)holo1[pixel]*tukey_holo[pixel]; }///multiply by Tukey windows

    ///--------Circshift et TF2D HOLOGRAMME------
    holo_shift=fftshift2D(holo1);
    //SAV2(holo1, "/home/mat/tomo_test/holo_shift_extract_holo.bin",t_float,"a+b");

    TF2Dcplx_vec(in,out,holo_shift, TF_Holo,p_forward_holo);
    TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
    //  SAVCplx(TFHoloCentre,"Re","/home/mat/tomo_test/TFHoloCentre.raw",t_float,"a+b");
    coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]
    //SAVCplx(TFHoloCentre,"Re","/home/mat/TFHoloCentre.raw",t_float,"a+b");
    ///--------Découpe hors axée------------------
    // coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]
}

///ancienne fonction, avec plans calculé dans la fonction, donc autonome.
void holo2TF_UBorn_old(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukey_holo)
{
     ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
    size_t dimx=sqrt(NbPix2dROI);
    int fftwThreadInit;
    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(4);
    ///prepare fftw plan+tableaux-----------------
    fftw_plan p_forward_holo, p_backward_holo;
    //fftw_complex *in_out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * dimROI.x*dimROI.y);//in=out pour transformation "inplace".
    fftw_complex *in=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPix2dROI);//in=out pour transformation "inplace".
    fftw_complex *out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPix2dROI);//in=out pour transformation "inplace".
    p_forward_holo=fftw_plan_dft_2d(dimROI.x, dimROI.y, in, out, FFTW_FORWARD,FFTW_ESTIMATE);
    p_backward_holo=fftw_plan_dft_2d(dimROI.x, dimROI.y, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        size_t NbPixROI2d=holo1.size();
        vector<double> holo_shift(NbPixROI2d);
        vector<complex<double>> TF_Holo(NbPixROI2d);
        vector<complex<double>> TFHoloCentre(NbPixROI2d);
        nbCplx *TF_UBorn_A=new nbCplx[NbPixROI2d];

        for(size_t pixel=0; pixel<NbPixROI2d; pixel++) {
                holo1[pixel]=(double)holo1[pixel]*tukey_holo[pixel];
        }

        ///--------Circshift et TF2D HOLOGRAMME------
        holo_shift=fftshift2D(holo1);
        //SAV2(holo1, "/home/mat/tomo_test/holo_shift_extract_holo.bin",t_float,"a+b");

        TF2Dcplx_vec(in,out,holo_shift, TF_Holo,p_forward_holo);
        TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
        //  SAVCplx(TFHoloCentre,"Re","/home/mat/tomo_test/TFHoloCentre.raw",t_float,"a+b");

        coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]
        //SAVCplx(TFHoloCentre,"Re","/home/mat/tomo_test/TFHoloCentre.raw",t_float,"a+b");
        ///--------Découpe hors axée------------------
       // coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]
}

///@parameters PosSpec : position of the specular beam //surcharge FFTW_init
void calc_Uborn2(vector<complex<double>> const &TF_UBorn,vector<complex<double>> &UBorn,Var2D dim2DHA,Var2D PosSpec,FFTW_init &param_c2c)
{
    Var2D recalUBorn={-PosSpec.x,-PosSpec.y},DecalU_Born={dim2DHA.x/2,dim2DHA.y/2};
    size_t NbPixUBorn=dim2DHA.x*dim2DHA.y;

    vector<complex<double>> TF_UBorn_I(NbPixUBorn);

    //save shifted image (with fringes)
    /*   {
    vector<complex<double>> UBorn_I2(NbPixUBorn);
    TF2Dcplx_INV(fftshift2D((TF_UBorn)), UBorn_I2, param_c2c);

   // decal2DCplxGen2(UBorn_I,UBorn_I2,DecalU_Born);
    SAVCplx(fftshift2D(UBorn_I2),"Re","/home/mat/tomo_test/Uborn_decal.raw",t_float,"a+b");
    }*/


    decal2DCplxGen2(TF_UBorn,TF_UBorn_I,recalUBorn);

    vector<complex<double>> UBorn_I(NbPixUBorn);
    TF2Dcplx_INV(TF_UBorn_I, UBorn_I, param_c2c);

//SAVCplx(fftshift2D(UBorn_I),"Re","/home/mat/tomo_test/Uborn_I.raw",t_float,"a+b");

    decal2DCplxGen2(UBorn_I,UBorn,DecalU_Born);

}


