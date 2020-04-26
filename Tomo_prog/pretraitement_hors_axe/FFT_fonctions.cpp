#include "FFT_fonctions.h"
#include "fonctions.h"
#include "omp.h"

using namespace std;

vector<vecteur>  fftshift2D(vector<vecteur> &entree)
{
unsigned int dim=sqrt(entree.size());
vector<vecteur > result(dim*dim);
size_t decal=dim/2;
size_t yi=0;
size_t xi=0;
       // #pragma omp parallel for
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

vector<complex<double> >  fftshift2D(vector<complex<double>> const &entree)
{
unsigned int dim=sqrt(entree.size());
vector<complex<double> > result(dim*dim);
size_t decal=dim/2;
size_t yi=0;
       size_t xi=0;
    //#pragma omp parallel for
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


void  fftshift2D(vector<complex<double>> const &entree,vector<complex<double>>  &result)
{
unsigned int dim=sqrt(entree.size());
//vector<complex<double> > result(dim*dim);
size_t decal=dim/2;
size_t yi=0;
       size_t xi=0;
       // #pragma omp parallel for
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
                      result[pixel_shift]=entree[pixel];                      //4 eme quadrant vers 1er
                      result[pixel]=entree[pixel_shift];
                      //2eme vers 3eme
                      result[(yi+decal)*dim+xi]=entree[pixel+decal];
                      //3eme vers 2eme
                      result[pixel+decal]=entree[(yi+decal)*dim+xi];
                }
        }
     //   return result;
}
//surcharge avec passage par paramètre
void  fftshift2D(vector<double> const &entree,vector<double>  &result)
{
unsigned int dim=sqrt(entree.size());
//vector<complex<double> > result(dim*dim);
size_t decal=dim/2;
size_t yi=0;
       size_t xi=0;
       // #pragma omp parallel for
       for(yi=0; yi<decal; yi++) {
            size_t num_ligne=yi*dim, num_lgn_decal=(yi+decal)*dim;
                for(xi=0; xi<decal; xi++)
                {
                      int pixel=num_ligne+xi;
                      //  cout<<"pixel="<<pixel<<endl;
                      // cout<<"result[pixel]="<<result[pixel];
                      int pixel_shift=num_lgn_decal+xi+decal;
                      //  cout<<"pixel_shift="<<pixel_shift<<endl;
                      //1er quadrant vers 4 eme
                      result[pixel_shift]=entree[pixel];                      //4 eme quadrant vers 1er
                      result[pixel]=entree[pixel_shift];
                      //2eme vers 3eme
                      result[num_lgn_decal+xi]=entree[pixel+decal];
                      //3eme vers 2eme
                      result[pixel+decal]=entree[num_lgn_decal+xi];
                }
        }
     //   return result;
}
//surcharge avec passage par paramètre

vector<double>  fftshift2D(vector<double> const &entree)
{
unsigned int dim=sqrt(entree.size());
vector<double> result(dim*dim);
size_t decal=dim/2;
size_t yi=0;
       size_t xi=0;
    //  #pragma omp parallel for
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
//    int dim=sqrt(nbPix);
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
    size_t nbPix=entree.size();
    for(size_t cpt=0; cpt<nbPix; cpt++) {
        in[cpt][0]=entree[cpt].real();
        in[cpt][1]=entree[cpt].imag();
    }
//in = reinterpret_cast<fftw_complex*>(&entree);
    fftw_execute(p);

    for(size_t cpt=0; cpt<(nbPix); cpt++) {
        sortie[cpt].real(out[cpt][0]/nbPix); //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
        sortie[cpt].imag(out[cpt][1]/nbPix);
    }
}

///FFT2D entree=vector complex



void TF2Dcplx_vec_INV(fftw_complex *in,fftw_complex *out, vector<complex<double> > entree, vector<complex<double> > &sortie, fftw_plan p_backward)
{
    size_t nbPix=entree.size();

    for(size_t cpt=0; cpt<nbPix; cpt++){
        in[cpt][0]=entree[cpt].real();
        in[cpt][1]=entree[cpt].imag();
    }

    fftw_execute(p_backward);

    for(size_t cpt=0; cpt<(nbPix); cpt++){//BACKWARD, no normalization
        sortie[cpt].real(out[cpt][0]);
        sortie[cpt].imag(out[cpt][1]);
    }

}
///FFT2D entree=TF2D vector double
void TF2Dcplx_vec(fftw_complex *in,fftw_complex *out, vector<double> entree, vector<complex<double> > &sortie, fftw_plan p_forward){
    size_t nbPix=entree.size();
    //size_t dim=sqrt(nbPix);
    for(size_t cpt=0; cpt<nbPix; cpt++){
            in[cpt][0]=entree[cpt];
            in[cpt][1]=0;
        }
    // in = reinterpret_cast<fftw_complex*>(&entree);
    fftw_execute(p_forward);

    for(size_t cpt=0; cpt<(nbPix); cpt++){
            sortie[cpt].real(out[cpt][0]/nbPix); //division par N (dim*dim) si FORWARD pour normaliser la fftw qui n'est pas normalisée
            sortie[cpt].imag(out[cpt][1]/nbPix);
        }
}

///FFT2D entree=TF2D vector double+inplace
void TF2Dcplx_vec_INPLACE(fftw_complex *in_out,vector<double> entree, vector<complex<double> > &sortie, fftw_plan p){

    size_t nbPix=entree.size();
//    int dim=sqrt(nbPix);

    for(size_t cpt=0; cpt<nbPix; cpt++)
        {
            in_out[cpt][0]=entree[cpt];
            in_out[cpt][1]=0;
        }
    // in = reinterpret_cast<fftw_complex*>(&entree);
// SAV2_vec(entree,nbPix,"/home/mat/tomo_test/entree_dans_tf2D.bin",t_float,"a+b");

    fftw_execute(p);

    for(size_t cpt=0; cpt<(nbPix); cpt++)
        {
            sortie[cpt].real(in_out[cpt][0]/nbPix); //division par N (dim*dim) si FORWARD pour normaliser la fftw qui n'est pas normalisée
            sortie[cpt].imag(in_out[cpt][1]/nbPix);
        }
    //  SAVCplx(sortie,"Im",nbPix,"/home/mat/tomo_test/spectre_dans_tf2D.bin",t_float,"a+b");
}

void TF2Dcplx_vec_INV(fftw_complex *in, fftw_complex *out, vector<double> entree, vector<complex<double> > &sortie, fftw_plan p_backward){
    size_t nbPix=entree.size();
    for(size_t cpt=0; cpt<nbPix; cpt++){
        in[cpt][0]=entree[cpt];
        in[cpt][1]=0;
    }

    fftw_execute(p_backward);

        for(size_t cpt=0; cpt<(nbPix); cpt++){//BACKWARD : no normalization
            sortie[cpt].real(out[cpt][0]);
            sortie[cpt].imag(out[cpt][1]);
        }
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

void TF2Dcplx_INV(vector<complex<double>> const &entree, vector<complex<double> > &sortie, FFTW_init &tf2D)
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


///-------TF2D purement reelle
///TF2D avec entrée purement réelle, utilise r2c et exporte le demi-espace.
void TF2D_r2c(vector<double> const &entree, vector<complex<double> > &sortie, FFTW_init const &tf2D_Re, double delta_x){

  size_t nbPix=entree.size(),   dim=sqrt(nbPix);
  Var2D dimROI{dim,dim};
  size_t nbPixCplx=dim*(dim/2+1);
  // double Coef_norm=dim*m1.Tp_Uborn;
  double Coef_norm=nbPixCplx;//ow(delta_x,2);
  // #pragma omp parallel
  for(size_t  cpt=0; cpt<nbPix; cpt++){
    tf2D_Re.in_double[cpt]=entree[cpt];
  }
  fftw_execute(tf2D_Re.p_forward_OUT);
  size_t cpt=0;
  for( cpt=0; cpt<(nbPixCplx); cpt++){
    sortie[cpt].real(tf2D_Re.out[cpt][0]/Coef_norm);
    sortie[cpt].imag(tf2D_Re.out[cpt][1]/Coef_norm);
  }
}
///TF2D avec entrée purement réelle, utlise r2c, mais resymétrise la sortie
void TF2D_r2c_symetric(vector<double> const &entree, vector<complex<double> > &sortie, FFTW_init  &tf2D_Re){
  size_t dim=sqrt(entree.size());
  Var2D dimROI{dim,dim};
  double nbPix=entree.size();
  // double Coef_norm=dim*m1.Tp_Uborn;
 // double Coef_norm=pow(delta_x,2);
  double Coef_norm=1/(nbPix);
 // cout<<"Coef_norm"<<Coef_norm<<endl;
// #pragma omp parallel
  for(int cpt=0; cpt<nbPix; cpt++){
    tf2D_Re.in_double[cpt]=entree[cpt];
  }

int    fftwThreadInit=fftw_init_threads();
    fftw_plan_with_nthreads(tf2D_Re.m_Nthread);
  fftw_execute(tf2D_Re.p_forward_OUT);

  //#pragma omp parallel for num_threads(2)
  for(int yA=0;yA<dimROI.y;yA++){///copie de A vers  la partie gauche (B), avec l'axe de symétrie (donc jusque dim.x/2+1)
    size_t num_lgnA=yA*(dimROI.x/2+1);
    size_t num_lgnB=yA*(dimROI.x);
    for(int xA=0;xA<dimROI.x/2+1;xA++){
      int cptA=xA+num_lgnA;
      int cptB=xA+num_lgnB;
      sortie[cptB].real(tf2D_Re.out[cptA][0]*Coef_norm);
      sortie[cptB].imag(tf2D_Re.out[cptA][1]*Coef_norm);
    }
  }
  for(int yA=1;yA<dimROI.y;yA++){///copie de A* vers  la partie droite (C), sans l'axe de symétrie (donc partie A jusque dim.x/2-1)
    size_t num_lgnA=yA*(dimROI.x/2+1);
    size_t num_lgnB=(dimROI.y-yA)*(dimROI.x);
    for(size_t xA=1;xA<dimROI.x/2;xA++){
      size_t cptA=xA+num_lgnA;
      size_t cptC=(dimROI.x-xA)+num_lgnB;
      sortie[cptC].real(tf2D_Re.out[cptA][0]*Coef_norm);
      sortie[cptC].imag(-tf2D_Re.out[cptA][1]*Coef_norm);
    }
  }
  //copie (conjuguée) de la demi ligne de y=0, x=[0,dim/2-1]
  for(size_t xA=1;xA<dimROI.x/2;xA++){
    size_t cptA=xA;
    size_t cptD=(dimROI.x-xA);
    sortie[cptD].real(tf2D_Re.out[cptA][0]*Coef_norm);
    sortie[cptD].imag(-tf2D_Re.out[cptA][1]*Coef_norm);
  }
}
/*
void TF2D_r2c_symetric(vector<double> const &entree, vector<complex<double> > &sortie, FFTW_init  &tf2D_Re){
  size_t dim=sqrt(entree.size()),  nbPix=entree.size();
  Var2D dimROI{dim,dim};
 // double Coef_norm=pow(delta_x,2);
  double Coef_norm=1;///nbPix; //car durant le hors, il ya un aller retour avec le même nombre de dimension (mails il ya découpe ?)
  // #pragma omp parallel
 // SAV2(entree,"/home/mat/tmp/Holo_shift.raw",t_float,"a+b");
  for(size_t cpt=0; cpt<nbPix; cpt++){
    tf2D_Re.in_double[cpt]=entree[cpt];
  }
  fftw_execute(tf2D_Re.p_forward_OUT);

  size_t num_lgnA=0,num_lgnB=0, cptA=0,cptB=0;
  //#pragma omp parallel for num_threads(2)
  for(size_t yA=0;yA<dimROI.y;yA++){///copie de A vers  la partie gauche (B), avec l'axe de symétrie (donc jusque dim.x/2+1)
    num_lgnA=yA*(dimROI.x/2+1);
    num_lgnB=yA*(dimROI.x);
    for(size_t xA=0;xA<dimROI.x/2+1;xA++){
      cptA=xA+num_lgnA;
      cptB=xA+num_lgnB;
      sortie[cptB].real(tf2D_Re.out[cptA][0]*Coef_norm);
      sortie[cptB].imag(tf2D_Re.out[cptA][1]*Coef_norm);
      //cout<<tf2D_Re.out[cptA][0]<<endl;
    }
  }
  for(int yA=1;yA<dimROI.y;yA++){///copie de A* vers  la partie droite (C), sans l'axe de symétrie (donc partie A jusque dim.x/2-1)
    size_t num_lgnA=yA*(dimROI.x/2+1);
    size_t num_lgnB=(dimROI.y-yA)*(dimROI.x);
    for(size_t xA=1;xA<dimROI.x/2;xA++){
      size_t cptA=xA+num_lgnA;
      size_t cptC=(dimROI.x-xA)+num_lgnB;
      sortie[cptC].real(tf2D_Re.out[cptA][0]*Coef_norm);
      sortie[cptC].imag(-tf2D_Re.out[cptA][1]*Coef_norm);
    }
  }
  //copie (conjuguée) de la demi ligne de y=0, x=[0,dim/2-1]
  for(size_t xA=1;xA<dimROI.x/2;xA++){
    size_t cptA=xA;
    size_t cptD=(dimROI.x-xA);
    sortie[cptD].real(tf2D_Re.out[cptA][0]*Coef_norm);
    sortie[cptD].imag(-tf2D_Re.out[cptA][1]*Coef_norm);
  }

}*/


