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
#include "symetrisation.h"
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

///calculate the wrapped phase from -pi to pi. useless now, please used atan2 in C++ standard library
void calcPhase_mpi_pi_atan2(vector<complex<double>> const &cplxField, vector<double> &phaseMod2pi)///calcul phase -PI-PI
{
for(int cpt=0;cpt<cplxField.size();cpt++)
phaseMod2pi[cpt]=atan2(cplxField[cpt].imag(),cplxField[cpt].real());
}

///--------Initialize reference
vector<double> initRef(string chemin_ref, Var2D coin, Var2D dimROI, Var2D dim2DHA){
size_t nbPixROI2d=dimROI.x*dimROI.y;
size_t nbPixHA=dim2DHA.x*dim2DHA.y;
vector<double> ampli_ref(nbPixROI2d);
vector<double> intensite_ref(nbPixROI2d);

   if(is_readable(chemin_ref)==1){
         charger_image2D_OCV_UNI(intensite_ref,chemin_ref, coin, dimROI,dim2DHA);//reduction de la taille aux dimension hors axe
             //  SAV_Tiff2D(intensite_ref,"/home/mat/tomo_test/intensite_ref.pgm",1);
         }
         else cout<<"/!\\  fichier intensité référence absent, création intensité unité"<<endl;

    for(size_t cpt=0;cpt<nbPixHA;cpt++){
        if(intensite_ref[cpt]!=0){
        ampli_ref[cpt]=sqrt(intensite_ref[cpt]);
        }
        else
        ampli_ref[cpt]=1;
    }

   // SAV_Tiff2D(ampli_ref,"/home/mat/tomo_test/ampli_ref.pgm",1);
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
///crop src2D[0:dim_src,0:dim_src] to dest3D(coin.x:coin.x+dim_dest,coin.y+dim_dest), human
//cette fonction est compliquée ! il faut 3 indices : un dans l'image 1024x204, 1 dans le stack 3D 220*220*nbAngle et 1 pour l'image 2D en cours
void coupeCplxTukey(vector<complex<double>> const &src, vector<complex<double>> &dest, Var2D dim_src, Var2D dim_dest, Var2D coin, size_t NumAngle, vector<double>  &tukeyHA)
{
    vector<complex<double>> imgCrop(dim_dest.x*dim_dest.y);
    size_t X_dest,Y_dest, cpt_dest1D,
           X_src, Y_src, cpt_src1D, cpt_Z_dest;
    size_t NA_pix_carre=(dim_dest.x/2)*(dim_dest.y/2);
    cpt_Z_dest=(dim_dest.x*dim_dest.y)*NumAngle;
    for(Y_dest=0; Y_dest<dim_dest.y; Y_dest++)
    {
        size_t num_lgn= Y_dest*dim_dest.x;
        for(X_dest=0; X_dest<dim_dest.x; X_dest++)
        {
            cpt_dest1D=cpt_Z_dest+X_dest+num_lgn;///coord 1D destination

            //coordonnées de découpe dans la source
            X_src=coin.x+X_dest;///coord X src
            Y_src=coin.y+Y_dest;///coord Y src
            cpt_src1D=X_src+Y_src*dim_src.x;///coord 1D source
            size_t idx = X_dest + Y_dest * dim_dest.x;//coord for tukey windows
            //copie src->dest
            //dest[cpt_dest1D]=src[cpt_src1D];
            if((X_dest-dim_dest.x/2)*(X_dest-dim_dest.x/2)+(Y_dest-dim_dest.y/2)*(Y_dest-dim_dest.y/2)<=NA_pix_carre)
            {

                //imgCrop[idx]=src[cpt_src1D]*tukeyHA[idx];
                imgCrop[idx]=src[cpt_src1D];
            }
            else
            {
                imgCrop[idx]=0;
            }
              dest[cpt_dest1D].real(imgCrop[idx].real());
            dest[cpt_dest1D].imag(imgCrop[idx].imag());
           // dest[cpt_dest1D].real(tukeyHA[idx]*imgCrop[idx].real());
           // dest[cpt_dest1D].imag(tukeyHA[idx]*imgCrop[idx].imag());
        }
    }

  //  SAVCplx(imgCrop,"Im","/home/mat/tomo_test/imgCrop_Im_90x90x534x32.raw",t_float,"a+b");
    // SAVCplx(imgCrop,"Re","/home/mat/tomo_test/imgCrop_Re_208x208x500x32.raw",t_float,"a+b");
}

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
void holo2TF_UBornTukeyHA_r2c(vector<double>  &holo1,vector<complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, vector<double> const &tukeyHA,FFTW_init  &param_fftw2DHolo)
{
    size_t NbPixROI2d=holo1.size();
    vector<complex<double>> TF_Holo(NbPixROI2d);
    for(size_t pixel=0; pixel<NbPixROI2d; pixel++)
      holo1[pixel]=(double)holo1[pixel];//*tukeyHolo[pixel];

    TF2D_r2c_symetric(fftshift2D2(holo1),TF_Holo,param_fftw2DHolo);

    coupeCplx(fftshift2D2(TF_Holo), TF_UBornTot, dimROI, dim2DHA, coinHA, NbAngleOk);///Découpe à [-Nxmax,+NXmax]dans repère humain-lisible +envoi dans pile3D
    //coupeCplxTukey(fftshift2D2(TF_Holo), TF_UBornTot, dimROI, dim2DHA, coinHA, NbAngleOk, tukeyHA);
    //coupeCplx(TF_Holo, TF_UBornTot, dimROI, dim2DHA, coinHA, NbAngleOk);///Découpe à [-Nxmax,+NXmax]dans repère humain-lisible +envoi dans pile3D
   // SAVCplx(TF_Holo,"Im","/home/mat/tomo_test/Tfholo_208x208x60.bin",t_float,"a+b");
}

///ancienne fonction, lente, mais avec plan calculé à l'extérieur
void holo2TF_UBorn(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukey_holo, fftw_complex *in,fftw_complex *out,fftw_plan p_forward_holo)
{
    ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
   // size_t dimx=sqrt(NbPix2dROI);

    size_t NbPixROI2d=holo1.size();
    vector<double> holo_shift(NbPixROI2d);
    vector<complex<double>> TF_Holo(NbPixROI2d);
    vector<complex<double>> TFHoloCentre(NbPixROI2d);

   // for(size_t pixel=0; pixel<NbPixROI2d; pixel++){ holo1[pixel]=(double)holo1[pixel]*tukey_holo[pixel]; }///multiply by Tukey windows

    ///--------Circshift et TF2D HOLOGRAMME------
    holo_shift=fftshift2D(holo1);
    //SAV2(holo1, "/home/mat/tomo_test/holo_shift_extract_holo.bin",t_float,"a+b");

    TF2Dcplx_vec(in,out,holo_shift, TF_Holo,p_forward_holo);
    TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
    //  SAVCplx(TFHoloCentre,"Re","/home/mat/tomo_test/TFHoloCentre.raw",t_float,"a+b");
    coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]

    SAVCplx(TFHoloCentre,"Re","/home/mat/TFHoloCentre.raw",t_float,"a+b");
    ///--------Découpe hors axée------------------
    // coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]
}

void holo2TF_UBornTukeyHA(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukeyHA, fftw_complex *in,fftw_complex *out,fftw_plan p_forward_holo)
{
    ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
   // size_t dimx=sqrt(NbPix2dROI);
    size_t NbPixROI2d=holo1.size();
    vector<double> holo_shift(NbPixROI2d);
    vector<complex<double>> TF_Holo(NbPixROI2d);
    vector<complex<double>> TFHoloCentre(NbPixROI2d);
   // for(size_t pixel=0; pixel<NbPixROI2d; pixel++){ holo1[pixel]=(double)holo1[pixel]*tukey_holo[pixel]; }///multiply by Tukey windows
    ///tukey on hologramms
  /*  vector<double>  masqueTukeyHolo=tukey2D(dimROI.x,dimROI.y,0.01);
    for(int cpt=0;cpt<holo.size();cpt++)    holo[cpt]=holo[cpt]*tukeyHA[cpt];*/

    ///--------Circshift et TF2D HOLOGRAMME------
    holo_shift=fftshift2D(holo1);
    //SAV2(holo1, "/home/mat/tomo_test/holo_shift_extract_holo.bin",t_float,"a+b");

    TF2Dcplx_vec(in,out,holo_shift, TF_Holo,p_forward_holo);

   //  SAVCplx(TF_Holo,"Re","/home/mat/tomo_test/TFHolo_1024x1024.raw",t_float,"a+b");
    TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)

     // SAVCplx(TFHoloCentre,"Re","/home/madeba/tomo_test/log_TFHoloCentre.raw",t_float,"a+b");
    //coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]
    coupeCplxTukey(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle,tukeyHA);
    ///--------Découpe hors axée------------------
    // coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]

}

void holo2TF_UBornSym(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukeyHA, fftw_complex *in,fftw_complex *out,fftw_plan p_forward_holoSym)
{
    ///--------------Init FFTW-------------------------------------------------
    size_t NbPixROI2d=holo1.size();
   // size_t dimx=sqrt(NbPix2dROI);
    vector<double> holoSym(4*NbPixROI2d);
    vector<double> holoSymShift(4*NbPixROI2d);
    vector<complex<double>> TF_HoloSym(4*NbPixROI2d);
    vector<complex<double>> TFHoloCentreSym(4*NbPixROI2d);
   // for(size_t pixel=0; pixel<NbPixROI2d; pixel++){ holo1[pixel]=(double)holo1[pixel]*tukey_holo[pixel]; }///multiply by Tukey windows
    ///tukey on hologramms
  /*  vector<double>  masqueTukeyHolo=tukey2D(dimROI.x,dimROI.y,0.01);
    for(int cpt=0;cpt<holo.size();cpt++)    holo[cpt]=holo[cpt]*tukeyHA[cpt];*/
   Symetrise_mirror2(holo1, holoSym);
   //SAV2(holoSym, "/home/madeba/tomo_test/holo_sym_1024x1024.bin",t_float,"a+b");
    ///--------Circshift et TF2D HOLOGRAMME------
    holoSymShift=fftshift2D(holo1);
    //SAV2(holoSymShift, "/home/madeba/tomo_test/holo_shift_sym_1024x1024.bin",t_float,"a+b");

    TF2Dcplx_vec(in,out,holoSymShift, TF_HoloSym,p_forward_holoSym);

    SAVCplx(TF_HoloSym,"Re","/home/madeba/tomo_test/TFHolo_1024x1024.raw",t_float,"a+b");
    TFHoloCentreSym=fftshift2D(TF_HoloSym);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)

     // SAVCplx(TFHoloCentreSym,"Re","/home/madeba/tomo_test/TFHoloCentre.bin",t_float,"a+b");
    //coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]
//    coupeCplxTukey(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle,tukeyHA);
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

///@parameters PosSpec : position of the specular beam //surcharge FFTW_init*
//export TF with 0 at (0,0) to calculate later grad u /u (grad u need spectrum)
vector<complex<double>> calc_Uborn2exportTF(vector<complex<double>> const &TF_UBorn,vector<complex<double>> &UBorn, Var2D dim2DHA,Var2D PosSpec,FFTW_init &param_c2c)
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
    return TF_UBorn_I;
}

///@parameters PosSpec : position of the specular beam //surcharge FFTW_init*
//export TF with 0 at (0,0) to calculate later grad u /u (grad u need spectrum)
vector<complex<double>> Symcalc_Uborn2exportTF(vector<complex<double>> const &TF_UBorn,vector<complex<double>> &UBorn, Var2D dim2DHA,Var2D PosSpec,FFTW_init &doubled_param_c2c)
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
    vector<complex<double>> Sym_TF_UBorn_I(4*NbPixUBorn);
    vector<complex<double>> Sym_UBorn_I(4*NbPixUBorn);
    Symetrise_mirror2(TF_UBorn_I,Sym_TF_UBorn_I);
    SAVCplx(TF_UBorn_I,"Re","/home/mat/tomo_test/TF_Uborn_I_Re_208x208x500.raw",t_float,"a+b");
    SAVCplx(Sym_TF_UBorn_I,"Re","/home/mat/tomo_test/Sym_TF_Uborn_I_Re_416x416x500.raw",t_float,"a+b");
    TF2Dcplx_INV(Sym_TF_UBorn_I, Sym_UBorn_I, doubled_param_c2c);

SAVCplx(Sym_UBorn_I,"Re","/home/mat/tomo_test/Sym_Uborn_I_416x416.raw",t_float,"a+b");
SAVCplx(Sym_UBorn_I,"Im","/home/mat/tomo_test/Sym_Uborn_I_416x416.raw",t_float,"a+b");
    UBorn_I=cut_quad4_cplx(Sym_UBorn_I);
    SAVCplx(UBorn_I,"Im","/home/mat/tomo_test/Uborn_I_Im_208x208x500.raw",t_float,"a+b");
    SAVCplx(UBorn_I,"Re","/home/mat/tomo_test/Uborn_I_Re_208x208x500.raw",t_float,"a+b");
    decal2DCplxGen2(UBorn_I,UBorn,DecalU_Born);
    return TF_UBorn_I;
}
