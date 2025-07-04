
#include <vector>
#include <complex>
#include "fonctions.h"
#include "FFT_fonctions.h"
#include "vecteur.h"
#include "deroulement_volkov4.h"
#include <chrono>
#define M_2PI 2*M_PI
using namespace std;


//uniquement symétrie miroir, déroulement parfait à un très faible piston près.
void deroul_volkov4_total_sym_paire(vector<double>  &phase_enroul,vector<double> &phase_deroul,vector <vecteur> double_kvect_shift,FFTW_init &param_c2c_double)
{
  complex<double> I(0,1);
  unsigned int nbPix=phase_enroul.size();
  Var2D dimROI={sqrt(nbPix),sqrt(nbPix)};

 // vector<complex<double>> gradx_enroul_fft(nbPix), grady_enroul_fft(nbPix);//attntion il s'git du gradient de la phase enroulée calculé par fft (et non le spectre)
  //vector<complex<double>> Z(nbPix);
 // vector<complex<double>> Gradx_Z_fft(nbPix), Grady_Z_fft(nbPix);

///-------variable 4 fois plus grande pour symétrie
///+ init opérateur et variables pour le calcul du gradient symétrisé
//vector<vecteur>  double_kvect_shift(4*nbPix);
vector<double> phase_enroul_sym(4*nbPix);
vector<complex<double>> gradx_enroul_fft_sym(4*nbPix), grady_enroul_fft_sym(4*nbPix);
//FFTW_init param_c2c_double(phase_enroul_sym,5);//initialiser la fftw avec une taille espace 4 fois plus grande
//double_kvect_shift=init_kvect_shift({2*dimROI.x,2*dimROI.y});
//------------------calcul gradient de la phase enroulé: on symétrize d'abord (symétrie paire)
Symetrise_mirror(phase_enroul,phase_enroul_sym);
//SAV_Tiff2D(phase_enroul,"/ramdisk/phase_enroul.tiff",1);

 // SAV_Tiff2D(phase_enroul_sym,"/ramdisk/phase_enroul_sym.tiff",1);

gradient_fft4(phase_enroul_sym, gradx_enroul_fft_sym,grady_enroul_fft_sym, double_kvect_shift,param_c2c_double);
///----------------------------------------------------
//SAV_Tiff2D(gradx_enroul_fft_sym,"Re","/ramdisk/gradx_enroul_fft_sym.tiff",1);
//SAV_Tiff2D(grady_enroul_fft_sym,"Re","/ramdisk/grady_enroul_fft_sym.tiff",1);

  vector<complex<double>> Z_sym(4*nbPix);
  vector<complex<double>> Gradx_Z_fft_sym(4*nbPix), Grady_Z_fft_sym(4*nbPix);


for(size_t cpt=0;cpt<4*nbPix;cpt++){
    Z_sym[cpt].real(cos(phase_enroul_sym[cpt]));
    Z_sym[cpt].imag(sin(phase_enroul_sym[cpt]));
  }
//SAV_Tiff2D(gradx_enroul_fft,"Re","/home/mat/tomo_test/gradx_enroul.tiff",1);
 //SAV_Tiff2D(Z_sym,"Re","/ramdisk/Z_sym.tiff",1);


 gradient_fft4(Z_sym, Gradx_Z_fft_sym,Grady_Z_fft_sym, double_kvect_shift, param_c2c_double);

 //SAV_Tiff2D(Gradx_Z_fft_sym,"Re","/ramdisk/gradx_Z_fft_sym.tiff",1);
  //SAVCplx(Z_sym,"Re",nbPix,"/home/mat/tomo_test/TF_gradx.bin",t_float,"w+b");
  /// Calcul du champ des entiers de déroulement

  complex<double> ax,ay;
  vector<double> gradx_IntM_sym(4*nbPix),grady_IntM_sym(4*nbPix);

  for(size_t cpt=0;cpt<4*nbPix;cpt++){
    ax=-I*(Gradx_Z_fft_sym[cpt]/Z_sym[cpt]);
    ay=-I*(Grady_Z_fft_sym[cpt]/Z_sym[cpt]);
    gradx_IntM_sym[cpt]=(ax.real()-gradx_enroul_fft_sym[cpt].real())/(M_2PI);//*(-0.159);//-1/2pi
    grady_IntM_sym[cpt]=(ay.real()-grady_enroul_fft_sym[cpt].real())/(M_2PI);
  }
//  SAV_Tiff2D(gradx_IntM_sym,"/ramdisk/gradx_IntM_sym.tiff",1);
   // SAV_Tiff2D(grady_IntM_sym,"/ramdisk/grady_IntM_sym.tiff",1);

 vector<complex<double>> double_IntM(4*nbPix);

integ_grad4(gradx_IntM_sym,grady_IntM_sym,double_IntM,double_kvect_shift,param_c2c_double);
 // SAV_Tiff2D(double_IntM,"Re","/ramdisk/double_Int_M.tiff",1);

  vector<double> IntM(nbPix);
  IntM=cut_quad4(double_IntM);

 for(size_t cpt=0;cpt<nbPix;cpt++){
        //if(IntM[cpt]>0)     IntM[cpt]=ceil(100*IntM[cpt])/100;
       //if(IntM[cpt]<0.01)     IntM[cpt]=0;
       //IntM[cpt]=round(IntM[cpt]*100)/100;
     // double residu=(round(IntM[cpt])-IntM[cpt])*2*M_PI;
    phase_deroul[cpt]=phase_enroul[cpt]+M_2PI*(IntM[cpt]);//+residu;
  }
//  SAV_Tiff2D(IntM,"/home/mat/tomo_test/Int_M.tiff",1);
  ///calcul phase dépliée en ajoutant 2\pi*carte des entiers
//SAV_Tiff2D(IntM,"/ramdisk/Int_M.tiff",1);
//SAV_Tiff2D(IntM,"Re","/home/mat/tomo_test/Int_M.tiff",1);

}



///gradient par fft, entrée réelle
///calculate gradient by fft, real input
void gradient_fft4(vector<double>   &entree, vector<complex<double>> &gradx, vector<complex<double>> &grady, vector<vecteur> &kvect_shift, FFTW_init &param_c2c)
{
     complex<double> I(0,1);
     unsigned int nbPix=entree.size();
     unsigned int dim=sqrt(nbPix);
     vector<complex<double>> tamponx(nbPix), tampony(nbPix);
     vector<complex<double>> spectre(nbPix);//spectre_shift(nbPix);
     //SAV2(entree,nbPix,"/home/mat/tomo_test/phase.bin",t_float,"w+b");
    // TF2D_vec(in,out, entree, spectre, p_forward);

    TF2Dcplx((entree),spectre,param_c2c);
//-----------------------------------------------------------------------------------

  /* vector<complex<double>> spectre_shift=fftshift2D(spectre);

   for(int cpt=0;cpt<entree.size();cpt++)
    {
        int cptXi=cpt%dim, cptYi=cpt/dim;
        int cptXc=cptXi-dim/2, cptYc=cptYi-dim/2;
       if(sqrt(cptXc*cptXc+cptYc*cptYc)>round(0.1*dim/2))//fixer ouverture numérique à 0.9*dim
       {
            spectre_shift[cpt].real(0);
            spectre_shift[cpt].imag(0);
        }

    }*/
//     SAV_Tiff2D(spectre_shift, "Re","/home/mat/tomo_test/spectre_shift_volkov4.tiff", 1);
   // spectre=fftshift2D(spectre_shift);
   // SAV_Tiff2D(spectre, "Re","/home/mat/tomo_test/spectre_volkov4.tiff", 1);
//-----------------------------------------------------------------------------------


     for(size_t cpt=0;cpt<nbPix;cpt++){
        tamponx[cpt]=2*M_PI*kvect_shift[cpt].getx()*spectre[cpt]*I;
        tampony[cpt]=2*M_PI*kvect_shift[cpt].gety()*spectre[cpt]*I;
    }
    TF2Dcplx_INV(tamponx, gradx, param_c2c);
    TF2Dcplx_INV(tampony, grady, param_c2c);
}


///gradient par fft, entrée complexe
///calculate gradient by fft, complex input (overloaded function)
void gradient_fft4(vector<complex<double>>  &entree, vector<complex<double>> &gradx, vector<complex<double>> &grady,vector<vecteur>  &kvect_shift, FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=entree.size();
  unsigned int dim=sqrt(nbPix);
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

///intégration du gradient par fft, entrée réelle
void integ_grad4(vector<double> const &gradx, vector<double> const& grady, vector<complex<double>> &sortie,std::vector<vecteur> &kvect_shift, FFTW_init &param_c2c)
{
  complex<double> I(0,1);
  unsigned int nbPix=gradx.size(),dim=sqrt(nbPix);
  vector<complex<double>> TF_gradx(nbPix), TF_grady(nbPix), tampon(nbPix);
  vector<double> kvect_mod_sq(nbPix);
//  vector<vecteur> kvect(nbPix),kvect_shift(nbPix);

  TF2Dcplx((gradx), TF_gradx, param_c2c);
  TF2Dcplx((grady), TF_grady, param_c2c);
   TF_gradx=(TF_gradx);
  TF_grady=(TF_grady);
 // SAV_Tiff2D(TF_gradx,"Re","/home/mat/tomo_test/TF_gradx.tiff",1);

//  SAV_vec3D(kvect_shift,"x","/home/mat/tomo_test/kvect_shift_x.bin","wb",nbPix);
  //SAV_vec3D(kvect_shift,"y","/home/mat/tomo_test/kvect_shift_y.bin","wb",nbPix);
  for(size_t cpt=0;cpt<nbPix;cpt++){
    kvect_mod_sq[cpt]=pow(kvect_shift[cpt].getx(),2)+pow(kvect_shift[cpt].gety(),2);
   // cout<<kvect_shift[cpt].getx()<<endl;
    if(isfinite(kvect_mod_sq[cpt])==0)
        cout<<"division par zéro="<<kvect_mod_sq[cpt]<<endl;
    tampon[cpt]=-I*(TF_gradx[cpt]*kvect_shift[cpt].getx()+TF_grady[cpt]*kvect_shift[cpt].gety())/(M_2PI*kvect_mod_sq[cpt]);
  }
  //kvect_mod_sq[0]=(kvect_mod_sq[nbPix]+kvect_mod_sq[1])/2;
  //cout<<"kvect_mod[0]"<<kvect_mod_sq[0]<<endl;
  //SAV2(kvect_mod_sq,nbPix,"/home/mat/tomo_test/kvect_mod.bin",t_float,"wb");
  //SAV_Tiff2D(tampon,"Re","/home/mat/tomo_test/tampon_Re.tiff",1);
  tampon[0]=(tampon[nbPix]+tampon[1])*0.5;//moyennage de la fréquence zéro pour éviter NAN
  TF2Dcplx_INV(tampon,sortie,param_c2c);
}



///calculate the kvector field (array whose value a simply  kx and ky)
/*std::vector<vecteur> init_kvect_shift(Var2D dim2DHA)
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
*/


void Symetrise_mirror(vector<double> const &monImg, vector<double> &monImgSymetric){
    size_t cpt_orig,cpt0_final,cpt1_final,cpt2_final,cpt3_final;
    size_t dim=sqrt(monImg.size());
    for(size_t x=0;x<dim;x++)
    for(size_t y=0;y<dim;y++){
        cpt_orig=x+y*dim;
        cpt0_final=dim-1-x+(dim-1-y)*(2*dim);
        monImgSymetric[cpt0_final]=monImg[cpt_orig];//-dxW(-x,-y)


        cpt1_final=dim+x+(dim-y-1)*(2*dim);
        monImgSymetric[cpt1_final]=monImg[cpt_orig];//dxW(-x,y)


        cpt2_final=dim-x-1+(dim+y)*(2*dim);
        monImgSymetric[cpt2_final]=monImg[cpt_orig];//dxW(x,-y)


        cpt3_final=dim+x+(dim+y)*(2*dim);
        monImgSymetric[cpt3_final]=monImg[cpt_orig];//dxW(x,y)//image originale

    }
//     return monImgSymetric;
}

//découpe le 4eme cadrant pour récupérer l'image dans une image symétrisée
vector<double> cut_quad4(vector<double> const &monImg4Quad)
{
    size_t cpt_orig,cpt4_final;
    size_t dim=sqrt(monImg4Quad.size()/4);//on balaye avec les coef de la petite image finale
    vector<double> monImgCut(dim*dim);
    for(size_t x=0; x<dim; x++)
        for(size_t y=0; y<dim; y++)
            {
                cpt_orig=x+y*dim;
                cpt4_final=dim+x+(dim+y)*(2*dim);
                monImgCut[cpt_orig]=monImg4Quad[cpt4_final];//dxW(x,y) //image originale

            }
    return monImgCut;
}
//découpe le 4eme cadrant pour récupérer l'image dans une image symétrisée (surcharge complexe)
vector<double> cut_quad4(vector<complex<double>> const &monImg4Quad)
{
    size_t cpt_orig,cpt4_final;
    size_t dim=sqrt(monImg4Quad.size()/4);//on balaye avec les coef de la petite image finale
    vector<double> monImgCut(dim*dim);
    for(size_t x=0; x<dim; x++)
        for(size_t y=0; y<dim; y++)
            {
                cpt_orig=x+y*dim;
                cpt4_final=dim+x+(dim+y)*(2*dim);
                monImgCut[cpt_orig]=monImg4Quad[cpt4_final].real();//dxW(x,y) //image originale
            }
    return monImgCut;
}

void gradient_central(const std::vector<double> src, std::vector<double> &grad,string direction)
{
const unsigned int dim=sqrt(src.size());
const unsigned int nbPix=dim*dim;
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
const unsigned int nbPix=dim*dim;
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
const unsigned int nbPix=dim*dim;
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
const unsigned int nbPix=dim*dim;
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
