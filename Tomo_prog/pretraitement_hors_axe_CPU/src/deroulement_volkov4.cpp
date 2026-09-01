///Functions used to unwrap phase, with symetrisation before fft
#include <vector>
#include <complex>
#include "fonctions.h"
#include "FFT_fonctions.h"
#include "vecteur.h"
#include "deroulement_volkov4.h"
#include <chrono>
#include "integration.h"
#include "symetrisation.h"
#include "regularized_division.h"
#define M_2PI 2*M_PI
using namespace std;
//Bisou

// symétrie paire sur la mesure + déroulement global i.e. calcul de grad phi =Im((grad uBorn)/uBorn)
//this function is called "exact solution" in Volkov paper, but need a division by uborn, which must be properly  done (see damped_division function)
vector<double> deroul_volkov5_total_sym_paire_gradu(vector<complex<double>> & UBorn, vector<vecteur> double_kvect_shift,FFTW_init &param_c2c_double, double alpha_damp)
{
    complex<double> I(0,1);
    size_t nbPix=UBorn.size();
    vector<complex<double>> Sym_UBorn(4*nbPix);

    Symetrise_mirror2(UBorn,Sym_UBorn);
   // SAVCplx(Sym_UBorn,"Re","/home/mat/tomo_test/Sym_UBorn_Re_416x416x500.raw",t_float,"a+b");

    vector<complex<double>> Sym_gradxUBorn(4*nbPix),Sym_gradyUBorn(4*nbPix);
    vector<double> Sym_gradxPhi(4*nbPix),Sym_gradyPhi(4*nbPix);

    gradient_fft4(Sym_UBorn,Sym_gradxUBorn,Sym_gradyUBorn,double_kvect_shift,param_c2c_double);

    Sym_gradxPhi=gradU_U_dampedDiv(Sym_UBorn,Sym_gradxUBorn, alpha_damp);
    Sym_gradyPhi=gradU_U_dampedDiv(Sym_UBorn,Sym_gradyUBorn, alpha_damp);

    ///Symetrise_mirror_gradient(Sym_gradxPhi,Sym_gradyPhi,Sym_gradxPhi,Sym_gradyPhi);
    // SAV2(Sym_gradxPhi,"/home/mat/tomo_test/Sym_gradxPHI_416x416x500.raw",t_float,"a+b");

    vector<complex<double>> Sym_phase_deroul(4*nbPix);
    //le déroulement est implicitement résolu par l’intégration (pas de champ d'entiers)
    //unwrapping is  implicitly solved through integration (no integer field)
    integ_grad4(Sym_gradxPhi,Sym_gradyPhi,Sym_phase_deroul,double_kvect_shift,param_c2c_double);
    // SAVCplx(Sym_phase_deroul,"Re","/home/mat/tomo_test/Sym_phase_deroul_Re_416x416x500.raw",t_float,"a+b");
    return cut_quad4(Sym_phase_deroul);//get back the right quadrant in the symmetrized image
}

///antisymétrie avant intégration + calcul de grad phi= (grad uBorn)/uBorn
vector<double> deroul_volkov5_sym_paire_gradu(vector<complex<double>> & UBorn, vector<vecteur> kvect_shift, vector<vecteur> double_kvect_shift,FFTW_init &param_c2c,FFTW_init &param_c2c_double)
{
    complex<double> I(0,1);
    size_t nbPix=UBorn.size();
    vector<complex<double>> gradxUBorn(nbPix),gradyUBorn(nbPix);
    vector<double> gradxPhi(nbPix),gradyPhi(nbPix);

    gradient_fft4(UBorn,gradxUBorn,gradyUBorn,kvect_shift,param_c2c);
        /*  SAVCplx(gradxUBorn,"Im","/home/mat/tomo_test/gradxUBorn_Im_208x208x500.raw",t_float,"a+b");*/

    for(size_t cpt=0;cpt<nbPix;cpt++){
    gradxPhi[cpt]=std::imag(gradxUBorn[cpt]/UBorn[cpt]);
    gradyPhi[cpt]=std::imag(gradyUBorn[cpt]/UBorn[cpt]);
    }
    //  SAV2(gradxPhi,"/home/mat/tomo_test/gradxPHI_Im_208x208x500.raw",t_float,"a+b");

    vector<double> Sym_gradxPhi(4*nbPix), Sym_gradyPhi(4*nbPix);

    Symetrise_mirror_gradient(gradxPhi,gradyPhi,Sym_gradxPhi,Sym_gradyPhi);
    //SAV2(Sym_gradxPhi,"/home/mat/tomo_test/Sym_gradxPHI_416x416x500.raw",t_float,"a+b");
    // SAV2(Sym_gradyPhi,"/home/mat/tomo_test/Sym_gradyPHI_416x416x500.raw",t_float,"a+b");

    vector<complex<double>> Sym_phase_deroul(4*nbPix);
    //le déroulement est implicitement résolu par l’intégration (pas de champ d'entiers)
    integ_grad4(Sym_gradxPhi,Sym_gradyPhi,Sym_phase_deroul,double_kvect_shift,param_c2c_double);
   //SAVCplx(Sym_phase_deroul,"Re","/home/mat/tomo_test/Sym_phase_deroul_Re_416x416x500.raw",t_float,"a+b");
   //  SAVCplx(Sym_phase_deroul,"Im","/home/mat/tomo_test/Sym_phase_deroul_Im_416x416x500.raw",t_float,"a+b");
       return cut_quad4((Sym_phase_deroul));
}


///----------------déroulement avec symétrie miroir sur la phase enroulée d'entrée+calcul champ d'entiers, déroulement quaasi parfait à un très faible piston près.
//Volkov mentionne 0.5 à 2% de bruit lié aux discontinuités du champ d'entier M
void deroul_volkov4_total_sym_paire(vector<double>  &phase_enroul,vector<double> &phase_deroul,vector <vecteur> double_kvect_shift,FFTW_init &param_c2c_double)
{
    complex<double> I(0,1);
    unsigned int nbPix=phase_enroul.size();

    ///-------variable 4 fois plus grande pour symétrie
    ///+ init opérateur et variables pour le calcul du gradient symétrisé
    vector<double> phase_enroul_sym(4*nbPix);
    vector<complex<double>> gradx_enroul_fft_sym(4*nbPix), grady_enroul_fft_sym(4*nbPix);
    //------------------calcul gradient de la phase enroulé: on symétrize d'abord (symétrie paire)
    Symetrise_mirror2(phase_enroul,phase_enroul_sym);
    gradient_fft4(phase_enroul_sym, gradx_enroul_fft_sym,grady_enroul_fft_sym, double_kvect_shift,param_c2c_double);
    ///----------------------------------------------------
    //SAV_Tiff2D(gradx_enroul_fft_sym,"Re","/ramdisk/gradx_enroul_fft_sym.tiff",1);

    vector<complex<double>> Z_sym(4*nbPix);
    vector<complex<double>> Gradx_Z_fft_sym(4*nbPix), Grady_Z_fft_sym(4*nbPix);

    for(size_t cpt=0; cpt<4*nbPix; cpt++)
    {
        Z_sym[cpt].real(cos(phase_enroul_sym[cpt]));
        Z_sym[cpt].imag(sin(phase_enroul_sym[cpt]));
    }
    gradient_fft4(Z_sym, Gradx_Z_fft_sym,Grady_Z_fft_sym, double_kvect_shift, param_c2c_double);
/// Calcul du champ des entiers de déroulement

    complex<double> ax,ay;
    vector<double> gradx_IntM_sym(4*nbPix),grady_IntM_sym(4*nbPix);

    for(size_t cpt=0; cpt<4*nbPix; cpt++)
    {
        ax=-I*(Gradx_Z_fft_sym[cpt]/Z_sym[cpt]);
        ay=-I*(Grady_Z_fft_sym[cpt]/Z_sym[cpt]);
        gradx_IntM_sym[cpt]=(ax.real()-gradx_enroul_fft_sym[cpt].real())/(M_2PI);//*(-0.159);//-1/2pi
        grady_IntM_sym[cpt]=(ay.real()-grady_enroul_fft_sym[cpt].real())/(M_2PI);
    }

    vector<complex<double>> double_IntM(4*nbPix);

    integ_grad4(gradx_IntM_sym,grady_IntM_sym,double_IntM,double_kvect_shift,param_c2c_double);

    vector<double> IntM(nbPix);
    IntM=cut_quad4(double_IntM);

    for(size_t cpt=0; cpt<nbPix; cpt++)
    {
        phase_deroul[cpt]=phase_enroul[cpt]+M_2PI*(IntM[cpt]);//+residu;
    }
}


