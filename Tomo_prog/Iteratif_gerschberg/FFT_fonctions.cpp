#include "FFT_fonctions.h"
#include "fonctions.h"


using namespace std;

/*vector<vecteur>  fftshift2D(vector<vecteur> &entree)
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

vector<complex<double> >  fftshift2D(vector<complex<double> > &entree)
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
}

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
void TF2D_vec(fftw_complex *in,fftw_complex *out, vector<double> entree, vector<complex<double> > &sortie, fftw_plan p){

    int nbPix=entree.size();
    int dim=sqrt(nbPix);
    for(int cpt=0; cpt<nbPix; cpt++){
            in[cpt][0]=entree[cpt];
            in[cpt][1]=0;
        }//
    fftw_execute(p);

    for(int cpt=0; cpt<(nbPix); cpt++){
            sortie[cpt].real(out[cpt][0]/nbPix); //division par N (dim*dim) si FORWARD pour normaliser la fftw qui n'est pas normalisée
            sortie[cpt].imag(out[cpt][1]/nbPix);
        }
}
///FFT2D entree=vector complex
void TF2Dcplx_vec(fftw_complex *in, fftw_complex *out, vector<complex<double> > entree, vector<complex<double> > &sortie,fftw_plan p)
{
    unsigned int nbPix=entree.size();
    for(int cpt=0; cpt<nbPix; cpt++) {
        in[cpt][0]=entree[cpt].real();
        in[cpt][1]=entree[cpt].imag();
    }
//in = reinterpret_cast<fftw_complex*>(&entree);
    fftw_execute(p);

    for(int cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(out[cpt][0]/nbPix); //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
        sortie[cpt].imag(out[cpt][1]/nbPix);
    }
}

void TF2Dcplx_vec_INV(fftw_complex *in,fftw_complex *out, vector<complex<double> > entree, vector<complex<double> > &sortie, fftw_plan p)
{
    unsigned int nbPix=entree.size();

    for(int cpt=0; cpt<nbPix; cpt++){
        in[cpt][0]=entree[cpt].real();
        in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(p);

    for(int cpt=0; cpt<(nbPix); cpt++){
        sortie[cpt].real(out[cpt][0]);
        sortie[cpt].imag(out[cpt][1]);
    }

}
///FFT2D entree=TF2D vector double
void TF2Dcplx_vec(fftw_complex *in,fftw_complex *out, vector<double> entree, vector<complex<double> > &sortie, fftw_plan p){
    int nbPix=entree.size();
    int dim=sqrt(nbPix);
    for(size_t cpt=0; cpt<nbPix; cpt++){
            in[cpt][0]=entree[cpt];
            in[cpt][1]=0;
        }
    // in = reinterpret_cast<fftw_complex*>(&entree);
    fftw_execute(p);

    for(size_t cpt=0; cpt<(nbPix); cpt++){
            sortie[cpt].real(out[cpt][0]/nbPix); //division par N (dim*dim) si FORWARD pour normaliser la fftw qui n'est pas normalisée
            sortie[cpt].imag(out[cpt][1]/nbPix);
        }
}

///FFT2D entree=TF2D vector double
void TF2Dcplx_vec_INPLACE(fftw_complex *in_out,vector<double> entree, vector<complex<double> > &sortie, fftw_plan p){

    int nbPix=entree.size();
    int dim=sqrt(nbPix);

    for(int cpt=0; cpt<nbPix; cpt++)
        {
            in_out[cpt][0]=entree[cpt];
            in_out[cpt][1]=0;
        }
    // in = reinterpret_cast<fftw_complex*>(&entree);
// SAV2_vec(entree,nbPix,"/home/mat/tomo_test/entree_dans_tf2D.bin",t_float,"a+b");

    fftw_execute(p);

    for(int cpt=0; cpt<(nbPix); cpt++)
        {
            sortie[cpt].real(in_out[cpt][0]/nbPix); //division par N (dim*dim) si FORWARD pour normaliser la fftw qui n'est pas normalisée
            sortie[cpt].imag(in_out[cpt][1]/nbPix);
        }
    //  SAVCplx(sortie,"Im",nbPix,"/home/mat/tomo_test/spectre_dans_tf2D.bin",t_float,"a+b");
}

void TF2Dcplx_vec_INV(fftw_complex *in, fftw_complex *out, vector<double> entree, vector<complex<double> > &sortie, fftw_plan p){
    unsigned int nbPix=entree.size();
    for(int cpt=0; cpt<nbPix; cpt++){
        in[cpt][0]=entree[cpt];
        in[cpt][1]=0;
    }

    fftw_execute(p);

        for(int cpt=0; cpt<(nbPix); cpt++){
            sortie[cpt].real(out[cpt][0]);
            sortie[cpt].imag(out[cpt][1]);
        }
}


