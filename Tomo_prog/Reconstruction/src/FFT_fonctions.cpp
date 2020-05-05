#include "FFT_fonctions.h"
#include "fonctions.h"
#include "manip.h"
#include "omp.h"
using namespace std;
///Set of functions to calculate FFT with fftw : c2c, r2c, r2c symetric, fftshift
 vector<complex<double>> fftshift3D(vector<complex<double>> &entree)
{       //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo

    int nbPix=entree.size();
    int dimFinale=round(std::pow(entree.size(), 1.0/3.0));
    int nbPix2D=dimFinale*dimFinale;//taille d'un plan 2D
    Var3D  decal={round(dimFinale/2),round(dimFinale/2),round(dimFinale/2)}, dim={dimFinale,dimFinale,dimFinale};
    vector<complex<double>> result(nbPix);
    size_t yi=0,xi=0,zi=0;
    int nbPix_z=0, nbPix_zdecal=0,nbPixy=0;
   #pragma omp parallel for private(zi)
    for(zi=0;zi<dim.z/2;zi++){
        nbPix_z=zi*nbPix2D;
        nbPix_zdecal=(zi+decal.z)*nbPix2D;
        //#pragma omp parallel for private(yi)
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
                    // 8eme quadrant vers 1er
                    result[pixel].real(entree[pixel_shift].real());
                    result[pixel].imag(entree[pixel_shift].imag());
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
void fftshift3D(vector<complex<double>> const &entree, vector<complex<double>> &result)
{       //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo

    int nbPix=entree.size();
    int dimFinale=round(std::pow(entree.size(), 1.0/3.0));
    int nbPix2D=dimFinale*dimFinale;//taille d'un plan 2D
    Var3D decal={round(dimFinale/2),round(dimFinale/2),round(dimFinale/2)}, dim={dimFinale,dimFinale,dimFinale};
//    vector<complex<double>> result(nbPix);
    size_t yi=0,xi=0,zi=0;
    int nbPix_z=0, nbPix_zdecal=0,nbPixy=0;
   #pragma omp parallel for
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

 vector<complex<double>> fftshift2D(vector<complex<double>> &entree)///
{
    int nbPix=entree.size();
    int dimFinale=round(std::pow(entree.size(), 0.5));
    Var2D decal={round(dimFinale/2),round(dimFinale/2)}, dim={dimFinale,dimFinale};
    vector<complex<double>> result(nbPix);
    int yi=0,xi=0;
        //#pragma omp parallel for private(yi)
       for(yi=0; yi<decal.y; yi++){
            size_t nbPixY=yi*dim.x;
            size_t nbpixY_Decal=(yi+decal.y)*dim.x;
                decal.x=dim.x/2;
                for(xi=0; xi<decal.x; xi++)
                {
                      int pixel=nbPixY+xi;
                      int pixel_shift=nbpixY_Decal+xi+decal.x;
                      //1er quadrant vers 4 eme
                      result[pixel_shift].real(entree[pixel].real());
                      result[pixel_shift].imag(entree[pixel].imag());
                      //4 eme quadrant vers 1er
                      result[pixel].real(entree[pixel_shift].real());
                      result[pixel].imag(entree[pixel_shift].imag());
                }
                decal.x=-dim.x/2;
                for(xi=dim.x/2; xi<dim.x; xi++){
                    int pixel=nbPixY+xi;
                    int pixel_shift=nbpixY_Decal+xi+decal.x;
                    //1er quadrant vers 4 eme
                    result[pixel_shift].real(entree[pixel].real());
                    result[pixel_shift].imag(entree[pixel].imag());
                    //4 eme quadrant vers 1er
                    result[pixel].real(entree[pixel_shift].real());
                    result[pixel].imag(entree[pixel_shift].imag());
                }
        }
return result;
}/*
vector<vecteur>  fftshift2D(vector<vecteur> &entree)
{
unsigned int dim=sqrt(entree.size());
vector<vecteur > result(dim*dim);
int decal=dim/2;
size_t yi=0;
       size_t xi=0;
        //#pragma omp parallel for private(yi)
       for(yi=0; yi<decal; yi++) {
            size_t num_ligne=yi*dim;
                for(xi=0; xi<decal; xi++)
                {
                      int pixel=num_ligne+xi;
                      //  cout<<"pixel="<<pixel<<endl;
                      // cout<<"result[pixel]="<<result[pixel];
                      int pixel_shift=(yi+decal)*dim+xi+decal;
                      //  cout<<"pixel_shift="<<pixel_shift<<endl;
                      //1er quadrant vers 4 eme
                      result[pixel_shift].setx(entree[pixel].getx());
                      result[pixel_shift].sety(entree[pixel].gety());

                      //4 eme quadrant vers 1er
                      result[pixel].setx(entree[pixel_shift].getx());
                      result[pixel].sety(entree[pixel_shift].gety());
                      //2eme vers 3eme
                      result[(yi+decal)*dim+xi].setx(entree[pixel+decal].getx());
                      result[(yi+decal)*dim+xi].sety(entree[pixel+decal].gety());
                      //3eme vers 2eme
                      result[pixel+decal].setx(entree[(yi+decal)*dim+xi].getx());
                      result[pixel+decal].sety(entree[(yi+decal)*dim+xi].gety());
                }
        }
        return result;
}*/
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

/*vector<complex<double>>  fftshift2D(vector<complex<double> > &entree)
{
unsigned int dim=sqrt(entree.size());
vector<complex<double> > result(dim*dim);
int decal=dim/2;
size_t yi=0;
       size_t xi=0;
        //#pragma omp parallel for private(yi)
       for(yi=0; yi<decal; yi++) {
            size_t num_ligne=yi*dim;
                for(xi=0; xi<decal; xi++)
                {
                      int pixel=num_ligne+xi;
                      //  cout<<"pixel="<<pixel<<endl;
                      // cout<<"result[pixel]="<<result[pixel];
                      int pixel_shift=(yi+decal)*dim+xi+decal;
                      //  cout<<"pixel_shift="<<pixel_shift<<endl;
                      //1er quadrant vers 4 eme
                      result[pixel_shift]=entree[pixel];

                      //4 eme quadrant vers 1er
                      result[pixel]=entree[pixel_shift];

                      //2eme vers 3eme
                      result[(yi+decal)*dim+xi]=entree[pixel+decal];

                      //3eme vers 2eme
                      result[pixel+decal]=entree[(yi+decal)*dim+xi];
                }
        }
        return result;
}*/

vector<double>  fftshift2D(vector<double> &entree)
{
unsigned int dim=sqrt(entree.size());
vector<double> result(dim*dim);
int decal=dim/2;
size_t yi=0;
       size_t xi=0;
        //#pragma omp parallel for private(yi)
       for(yi=0; yi<decal; yi++) {
            size_t num_ligne=yi*dim;
                for(xi=0; xi<decal; xi++)
                {
                      int pixel=num_ligne+xi;
                      //  cout<<"pixel="<<pixel<<endl;
                      // cout<<"result[pixel]="<<result[pixel];
                      int pixel_shift=(yi+decal)*dim+xi+decal;
                      //  cout<<"pixel_shift="<<pixel_shift<<endl;
                      //1er quadrant vers 4 eme
                      result[pixel_shift]=entree[pixel];
                      //4 eme quadrant vers 1er
                      result[pixel]=entree[pixel_shift];
                      //2eme vers 3eme
                      result[(yi+decal)*dim+xi]=entree[pixel+decal];

                      //3eme vers 2eme
                      result[pixel+decal]=entree[(yi+decal)*dim+xi];
                }
        }
        return result;
}

void TF2D(vector<double> entree, vector<complex<double> > &sortie, FFT_encaps &tf2D, double delta_x){

    int nbPix=entree.size();
    //int dim=sqrt(nbPix*m1.Tp_Uborn);
   // double Coef_norm=dim*m1.Tp_Uborn;
    double Coef_norm=pow(delta_x,2);
    for(int cpt=0; cpt<nbPix; cpt++){
            tf2D.in[cpt][0]=entree[cpt];
            tf2D.in[cpt][1]=0;
        }//
    fftw_execute(tf2D.p_forward_OUT);

    for(int cpt=0; cpt<(nbPix); cpt++){
            sortie[cpt].real(tf2D.out[cpt][0]*Coef_norm);
            sortie[cpt].imag(tf2D.out[cpt][1]*Coef_norm);
        }
}


void TF2Dcplx(vector<complex<double>> entree, vector<complex<double> > &sortie, FFT_encaps &tf2D, double delta_x){
    int nbPix=entree.size();
    int dim=sqrt(nbPix);
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


void TF2Dcplx_INV(vector<complex<double> > entree, vector<complex<double> > &sortie, FFT_encaps &tf2D, double delta_f)
{
    int nbPix=entree.size();
    int dim=sqrt(nbPix);
    double Coef_norm=pow(delta_f,2);
   // cout<<"delta_f="<<delta_f<<endl;
    size_t cpt=0;
    for(int cpt=0; cpt<nbPix; cpt++) {
        tf2D.in[cpt][0]=entree[cpt].real();
        tf2D.in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(tf2D.p_backward_OUT);

    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(tf2D.out[cpt][0]*Coef_norm); //division par N (dim*dim) pour normaliser l'énergie
        sortie[cpt].imag(tf2D.out[cpt][1]*Coef_norm);
    }
}


void TF3Dcplx(fftw_complex *in, fftw_complex *out, vector<complex<double> > entree, vector<complex<double> > &sortie, fftw_plan p3d_forward, double delta_x)
{
    int nbPix=entree.size();

   //double Coef_norm=sqrt((double)nbPix);
    double Coef_norm=pow(delta_x,3);
    size_t cpt=0;
    for(int cpt=0; cpt<nbPix; cpt++) {
        in[cpt][0]=entree[cpt].real();
        in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(p3d_forward);

    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(out[cpt][0]*Coef_norm);
        sortie[cpt].imag(out[cpt][1]*Coef_norm);
    }
}


void TF3Dcplx_INV(fftw_complex *in, fftw_complex *out, vector<complex<double> > entree, vector<complex<double> > &sortie, fftw_plan p3d_back, double delta_f)
{

    int nbPix=entree.size();
   // double Coef_norm=sqrt((double)nbPix*m1.Tp_Uborn);
    double Coef_norm=pow(delta_f,3);
    size_t cpt=0;
    for(int cpt=0; cpt<nbPix; cpt++) {
        in[cpt][0]=entree[cpt].real();
        in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(p3d_back);

    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(out[cpt][0]*Coef_norm); //division par N (dim^3) pour normaliser l'énergie (Parseval)
        sortie[cpt].imag(out[cpt][1]*Coef_norm);//division par N (dim^3) pour normaliser l'énergie (Parseval)
    }
}
void TF3Dcplx_Inplace_INV(vector<complex<double>> const & entree, vector<complex<double> > &sortie, FFT_encaps &param_Inplace_c2c, double delta_f)
{

    int nbPix=entree.size();
   // double Coef_norm=sqrt((double)nbPix*m1.Tp_Uborn);
    double Coef_norm=pow(delta_f,3);
    size_t cpt=0;
    for(int cpt=0; cpt<nbPix; cpt++) {
        param_Inplace_c2c.in[cpt][0]=entree[cpt].real();
        param_Inplace_c2c.in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(param_Inplace_c2c.p_backward_IN);
    for(cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(param_Inplace_c2c.in[cpt][0]*Coef_norm); //division par N (dim^3) pour normaliser l'énergie (Parseval)
        sortie[cpt].imag(param_Inplace_c2c.in[cpt][1]*Coef_norm);//division par N (dim^3) pour normaliser l'énergie (Parseval)
    }
}






