#include "FFT_fonctions.h"
#include "fonctions.h"
#include "omp.h"

using namespace std;


void fftshift3D(vector<complex<double>> const &entree, vector<complex<double>> &result)
{

    size_t nbPix=entree.size();
    size_t dimFinale=round(std::pow(entree.size(), 1.0/3.0));
    size_t nbPix2D=dimFinale*dimFinale;//taille d'un plan 2D
    Var3D decal={round(dimFinale/2),round(dimFinale/2),round(dimFinale/2)}, dim={dimFinale,dimFinale,dimFinale};
//    vector<complex<double>> result(nbPix);
    size_t yi=0,xi=0,zi=0;
    int nbPix_z=0, nbPix_zdecal=0,nbPixy=0;
  // #pragma omp parallel for
    for(zi=0;zi<dim.z/2;zi++){
        nbPix_z=zi*nbPix2D;
        nbPix_zdecal=(zi+decal.z)*nbPix2D;

        for(yi=0; yi<dim.y/2; yi++) {
                decal.x=dim.x/2;
                decal.y=dim.y/2;
                //#pragma omp parallel for private(xi)
                for(xi=0; xi<dim.x/2; xi++){///cadran 1<->8
                    int pixel=nbPix_z+yi*dim.x+xi;
                    int pixel_shift=nbPix_zdecal+(yi+decal.y)*dim.x+xi+decal.x;
                    // 1e quadrant vers 8e

                    result[pixel_shift].real(entree[pixel].real());
                    result[pixel_shift].imag(entree[pixel].imag());
                  //  result[pixel_shift]=entree[pixel];

                    // 8eme quadrant vers 1er
                    result[pixel].real(entree[pixel_shift].real());
                    result[pixel].imag(entree[pixel_shift].imag());
                   // result[pixel]=entree[pixel_shift];
                }
                decal.x=-dim.x/2;
                for(xi=dim.x/2; xi<dimFinale; xi++){
                    int pixel=nbPix_z+yi*dim.x+xi;
                    int pixel_shift=nbPix_zdecal+(yi+decal.y)*dim.x+xi+decal.x;
                    //2e quadrant vers 7 eme
                    result[pixel_shift].real(entree[pixel].real());
                    result[pixel_shift].imag(entree[pixel].imag());
                    // 7eme quadrant vers 2e
                    result[pixel].real(entree[pixel_shift].real());
                    result[pixel].imag(entree[pixel_shift].imag());
                }
        }

        for(yi=dim.y/2; yi<dim.y; yi++){
             decal.y=-dim.y/2;
             decal.x=dim.x/2;
            for(xi=0; xi<dim.x/2; xi++){///cadran 3<->6
                int pixel=zi*nbPix2D+yi*dim.x+xi;
                int pixel_shift=(zi+decal.z)*nbPix2D+(yi+decal.y)*dim.x+xi+decal.x;
                //3e quadrant vers 6e
                result[pixel_shift].real(entree[pixel].real());
                result[pixel_shift].imag(entree[pixel].imag());
                // 6eme vers 3e
                result[pixel].real(entree[pixel_shift].real());
                result[pixel].imag(entree[pixel_shift].imag());
            }
                decal.x=-dim.x/2;
            for(xi=dim.x/2; xi<dimFinale; xi++){///cadran 4<->5
                int pixel=zi*nbPix2D+yi*dim.x+xi;
                int pixel_shift=(zi+decal.z)*nbPix2D+(yi+decal.y)*dim.x+xi+decal.x;
                //2e quadrant vers 7 eme
                result[pixel_shift].real(entree[pixel].real());
                result[pixel_shift].imag(entree[pixel].imag());
                // 7eme quadrant vers 2e
                result[pixel].real(entree[pixel_shift].real());
                result[pixel].imag(entree[pixel_shift].imag());
            }
        }
    }
}


vector<complex<double>> fftshift3D(vector<complex<double>> const &entree)
{
    static vector<complex<double>> result(entree.size());
    size_t nbPix=entree.size();
    size_t dimFinale=round(std::pow(entree.size(), 1.0/3.0));
    size_t nbPix2D=dimFinale*dimFinale;//taille d'un plan 2D
    Var3D decal={round(dimFinale/2),round(dimFinale/2),round(dimFinale/2)}, dim={dimFinale,dimFinale,dimFinale};
//    vector<complex<double>> result(nbPix);
    size_t yi=0,xi=0,zi=0;
    int nbPix_z=0, nbPix_zdecal=0,nbPixy=0;
  // #pragma omp parallel for
    for(zi=0;zi<dim.z/2;zi++){
        nbPix_z=zi*nbPix2D;
        nbPix_zdecal=(zi+decal.z)*nbPix2D;

        for(yi=0; yi<dim.y/2; yi++) {
                decal.x=dim.x/2;
                decal.y=dim.y/2;
                //#pragma omp parallel for private(xi)
                for(xi=0; xi<dim.x/2; xi++){///cadran 1<->8
                    int pixel=nbPix_z+yi*dim.x+xi;
                    int pixel_shift=nbPix_zdecal+(yi+decal.y)*dim.x+xi+decal.x;
                    // 1e quadrant vers 8e

                    result[pixel_shift].real(entree[pixel].real());
                    result[pixel_shift].imag(entree[pixel].imag());
                  //  result[pixel_shift]=entree[pixel];

                    // 8eme quadrant vers 1er
                    result[pixel].real(entree[pixel_shift].real());
                    result[pixel].imag(entree[pixel_shift].imag());
                   // result[pixel]=entree[pixel_shift];
                }
                decal.x=-dim.x/2;
                for(xi=dim.x/2; xi<dimFinale; xi++){
                    int pixel=nbPix_z+yi*dim.x+xi;
                    int pixel_shift=nbPix_zdecal+(yi+decal.y)*dim.x+xi+decal.x;
                    //2e quadrant vers 7 eme
                    result[pixel_shift].real(entree[pixel].real());
                    result[pixel_shift].imag(entree[pixel].imag());
                    // 7eme quadrant vers 2e
                    result[pixel].real(entree[pixel_shift].real());
                    result[pixel].imag(entree[pixel_shift].imag());
                }
        }

        for(yi=dim.y/2; yi<dim.y; yi++){
             decal.y=-dim.y/2;
             decal.x=dim.x/2;
            for(xi=0; xi<dim.x/2; xi++){///cadran 3<->6
                int pixel=zi*nbPix2D+yi*dim.x+xi;
                int pixel_shift=(zi+decal.z)*nbPix2D+(yi+decal.y)*dim.x+xi+decal.x;
                //3e quadrant vers 6e
                result[pixel_shift].real(entree[pixel].real());
                result[pixel_shift].imag(entree[pixel].imag());
                // 6eme vers 3e
                result[pixel].real(entree[pixel_shift].real());
                result[pixel].imag(entree[pixel_shift].imag());
            }
                decal.x=-dim.x/2;
            for(xi=dim.x/2; xi<dimFinale; xi++){///cadran 4<->5
                int pixel=zi*nbPix2D+yi*dim.x+xi;
                int pixel_shift=(zi+decal.z)*nbPix2D+(yi+decal.y)*dim.x+xi+decal.x;
                //2e quadrant vers 7 eme
                result[pixel_shift].real(entree[pixel].real());
                result[pixel_shift].imag(entree[pixel].imag());
                // 7eme quadrant vers 2e
                result[pixel].real(entree[pixel_shift].real());
                result[pixel].imag(entree[pixel_shift].imag());
            }
        }
    }
    return result;
}

vector<double> tukey2D(int dimx,int dimy, float alpha)
{
        int N=dimx;
        vector<double> tuk2D(dimx*dimy);
        vector<double> tuk1Dx(dimx*dimy);
        vector<double> tuk1Dy(dimx*dimy);


        int borne1=round(alpha*(N-1)/2);
        int borne2=round((N-1)*(1-alpha/2));

        for(int cpt=0; cpt<borne1+1; cpt++)
                tuk1Dx[cpt]=0.5*(1+cos(3.1415*(2*cpt/(alpha*(N-1))-1)));
        for(int cpt=borne1+1; cpt<borne2+1; cpt++)
                tuk1Dx[cpt]=1;
        for(int cpt=borne2+1; cpt<N; cpt++)
                tuk1Dx[cpt]=0.5*(1+cos(3.1415*(2*cpt/(alpha*(N-1))-2/alpha+1)));

        for(int cpt=0; cpt<N*N; cpt++) {
                int cptx=cpt%(N);
                int cpty=cpt/(N);
                tuk2D[cpt]=tuk1Dx[cptx]*tuk1Dx[cpty];
        }
        return tuk2D;
}






/////////----------------------------TF FFTW_INIT-------------------------------------------------------------------------
// surcharge fftw_init+entree <double>+calcul fft
void TF2Dcplx(vector<double> const &entree, vector<complex<double>> &sortie, FFTW_init &tf2D_c2r)
{
    size_t nbPix=entree.size();
    for(size_t cpt=0; cpt<nbPix; cpt++) {
        tf2D_c2r.in[cpt][0]=entree[cpt];
        tf2D_c2r.in[cpt][1]=0;//
    }
//in = reinterpret_cast<fftw_complex*>(&entree);
    fftw_execute(tf2D_c2r.p_forward_OUT);

    for(size_t cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(tf2D_c2r.out[cpt][0]/nbPix); //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
        sortie[cpt].imag(tf2D_c2r.out[cpt][1]/nbPix);
    }
}

// surcharge fftw_init+entree <complex double>+calcul fft
///FFT2D complex c2c
void TF2Dcplx(vector<complex<double>> const &entree, vector<complex<double>> &sortie,FFTW_init &param_c2c)
{
    size_t nbPix=entree.size();
    for(size_t cpt=0; cpt<nbPix; cpt++) {
        param_c2c.in[cpt][0]=entree[cpt].real();
        param_c2c.in[cpt][1]=entree[cpt].imag();
    }
//in = reinterpret_cast<fftw_complex*>(&entree);
//fftw_complex * in = reinterpret_cast<fftw_complex*>(entree.data().begin());
//param_c2c.in = reinterpret_cast<fftw_complex*>(entree[0]);
    fftw_execute(param_c2c.p_forward_OUT);

    for(size_t cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(param_c2c.out[cpt][0]/nbPix); //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
        sortie[cpt].imag(param_c2c.out[cpt][1]/nbPix);
    }
}
///calcul TF vraie
void TF2Dcplx(vector<complex<double>> const &entree, vector<complex<double> > &sortie, FFTW_init &tf2D, double delta_x){
    size_t nbPix=entree.size();
//    size_t dim=sqrt(nbPix);
    //double Coef_norm=dim*m1.Tp_Uborn;
    double Coef_norm=pow(delta_x,2);
    for(size_t cpt=0; cpt<nbPix; cpt++){
            tf2D.in[cpt][0]=entree[cpt].real();
            tf2D.in[cpt][1]=entree[cpt].imag();
        }
    // in = reinterpret_cast<fftw_complex*>(&entree);
    fftw_execute(tf2D.p_forward_OUT);

    for(size_t cpt=0; cpt<(nbPix); cpt++){
            sortie[cpt].real(tf2D.out[cpt][0]*Coef_norm);
            sortie[cpt].imag(tf2D.out[cpt][1]*Coef_norm);
        }
}

void FFT2Dcplx_INV(vector<complex<double>> const &entree, vector<complex<double> > &sortie, FFTW_init &tf2D)
{
    size_t nbPix=entree.size();
//    int dim=sqrt(nbPix);

   // cout<<"delta_f="<<delta_f<<endl;
    size_t cpt=0;
    for(size_t cpt=0; cpt<nbPix; cpt++) {
        tf2D.in[cpt][0]=entree[cpt].real();
        tf2D.in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(tf2D.p_backward_OUT);

    for(cpt=0; cpt<(nbPix); cpt++) {///FFT inverse, no normalization
        sortie[cpt].real(tf2D.out[cpt][0]); //division par N (dim*dim) pour normaliser l'énergie
        sortie[cpt].imag(tf2D.out[cpt][1]);
    }
}
//Surcharge avec echantillonnage pour calcul exact de la transformée fourier
void TF2Dcplx_INV(vector<complex<double>> const &entree, vector<complex<double> > &sortie, FFTW_init &tf2D, double delta_f)
{
    size_t nbPix=entree.size();
//    int dim=sqrt(nbPix);
    double Coef_norm=pow(delta_f,2);
   // cout<<"delta_f="<<delta_f<<endl;
    size_t cpt=0;
    for(size_t cpt=0; cpt<nbPix; cpt++) {
        tf2D.in[cpt][0]=entree[cpt].real();
        tf2D.in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(tf2D.p_backward_OUT);

    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(tf2D.out[cpt][0]*Coef_norm); //division par N (dim*dim) pour normaliser l'énergie
        sortie[cpt].imag(tf2D.out[cpt][1]*Coef_norm);
    }
}

void FFT3Dcplx_Inplace_INV(vector<complex<double>> const & entree, vector<complex<double> > &sortie, FFTW_init &param_Inplace_c2c)
{
    int nbPix=entree.size();
   // double Coef_norm=sqrt((double)nbPix*m1.Tp_Uborn);
    //double Coef_norm=pow(delta_f,3);
    size_t cpt=0;
    for(int cpt=0; cpt<nbPix; cpt++) {
        param_Inplace_c2c.in[cpt][0]=entree[cpt].real();
        param_Inplace_c2c.in[cpt][1]=entree[cpt].imag();
    }
    fftw_execute(param_Inplace_c2c.p_backward_IN);
    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(param_Inplace_c2c.in[cpt][0]);
        sortie[cpt].imag(param_Inplace_c2c.in[cpt][1]);
    }
}

void FFT3Dcplx_Inplace(vector<complex<double>> const & entree, vector<complex<double> > &sortie, FFTW_init &param_Inplace_c2c)
{
    int nbPix=entree.size();
    double Coef_norm=nbPix;
    //double Coef_norm=pow(delta_f,3);
    size_t cpt=0;
    for(int cpt=0; cpt<nbPix; cpt++) {
        param_Inplace_c2c.in[cpt][0]=entree[cpt].real();
        param_Inplace_c2c.in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(param_Inplace_c2c.p_forward_IN);

    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(param_Inplace_c2c.in[cpt][0]/nbPix);
        sortie[cpt].imag(param_Inplace_c2c.in[cpt][1]/nbPix);
    }
}
