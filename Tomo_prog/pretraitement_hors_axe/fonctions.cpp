#include <fstream>//ifstream
#include <vector>
#include <chrono>
#include <complex>
#include <fftw3.h>
#include "struct.h"
#include <cv.h>
#include <highgui.h>//imread
#include "projet.h"
#include "FFT_fonctions.h"
#include "fonctions.h"

using namespace std;
using namespace cv;

void calcPhase_mpi_pi_atan2(vector<complex<double>> obj, vector<double> &phaseMod2pi)///calcul phase -PI-PI
{
for(int cpt=0;cpt<obj.size();cpt++)
phaseMod2pi[cpt]=atan2(obj[cpt].imag(),obj[cpt].real());
}

void circshift2DCplx(vector<complex<double>> entree, vector<complex<double>> &result, Var2D dim,Var2D decal)///___/!\ ne fonctionne qu'avec des demi espace??
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo

        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;
       // cout<<"decal.x="<<decal.x<<"; decal.y="<<decal.y<<endl;
       size_t yi=0;
       size_t xi=0;
        //#pragma omp parallel for private(yi)
       for(yi=0; yi<decal.y; yi++) {
            size_t num_ligne=yi*dim.x;
                for(xi=0; xi<decal.x; xi++)
                {
                        int pixel=num_ligne+xi;
                      //  cout<<"pixel="<<pixel<<endl;
                       // cout<<"result[pixel]="<<result[pixel];
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                      //  cout<<"pixel_shift="<<pixel_shift<<endl;
                        //1er quadrant vers 4 eme
                        result[pixel_shift]=entree[pixel];

                        //4 eme quadrant vers 1er
                        result[pixel]=entree[pixel_shift];

                        //2eme vers 3eme
                        result[(yi+decal.y)*dim.x+xi]=entree[pixel+decal.x];

                        //3eme vers 2eme
                        result[pixel+decal.x]=entree[(yi+decal.y)*dim.x+xi];
                }
        }
       //  #pragma omp barrier
}

void decal2DCplxGen(vector<complex<double>> entree, vector<complex<double>> &result, Var2D dim,Var2D decal)
{
            //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
            //cout<<"decal.x,y="<<decal.x<<","<<decal.y<<endl;
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;
        if(decal.x<0)
        decal.x=dim.x+decal.x;
        if(decal.y<0)
        decal.y=dim.y+decal.y;
       size_t yi,xi;
       //#pragma omp parallel for private(yi)
       for(yi=0; yi<dim.y-decal.y; yi++) {

                for(xi=0; xi<dim.x-decal.x; xi++)
                { //cout<<"xi,yi="<<xi<<","<<yi<<endl;
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                       // result[pixel_shift].real()=entree[pixel].real();

                        result[pixel_shift].real(entree[pixel].real());
                        result[pixel_shift].imag(entree[pixel].imag());
                }


                for(xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift].real(entree[pixel].real());
                        result[pixel_shift].imag(entree[pixel].imag());
                       // result[pixel_shift]=entree[pixel];
                }

        }
      // #pragma omp barrier
         //   #pragma omp parallel for private(yi)
              for(int yi=dim.y-decal.y; yi<dim.y; yi++) {

                for(int xi=0; xi<dim.x-decal.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift].real(entree[pixel].real());
                        result[pixel_shift].imag(entree[pixel].imag());
                        //result[pixel_shift]=entree[pixel];
                }

                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift].real(entree[pixel].real());
                        result[pixel_shift].imag(entree[pixel].imag());
                        //result[pixel_shift]=entree[pixel];
                }

        }
      //  #pragma omp barrier
}



void calc_Uborn(vector<complex<double>> const &TF_UBorn,vector<complex<double>> &UBorn,Var2D dim2DHA,Var2D PosSpec,fftw_complex *in,fftw_complex *out,fftw_plan p)
{

    Var2D recalUBorn={-PosSpec.x,-PosSpec.y},DecalU_Born={dim2DHA.x/2,dim2DHA.y/2};
   // cout<<"recalUbonrxy"<<recalUBorn.x<<","<<recalUBorn.y<<endl;
    //Var2D recal_filtrage={recalUBorn.x+DecalU_Born.x, recalUBorn.y+DecalU_Born.y};
    size_t NbPixUBorn=dim2DHA.x*dim2DHA.y;

  //SAVCplx(TF_UBorn, "Re", "/home/mat/tomo_test/TFUborn_extract_holo_re.raw", t_float, "a+b");

    vector<complex<double>> TF_UBorn_I(NbPixUBorn);

    decal2DCplxGen(TF_UBorn,TF_UBorn_I,dim2DHA,recalUBorn);/// ///
  //  SAVCplx(TF_UBorn_I, "Re", "/home/mat/tomo_test/TFUbornI_extract_holo_re.raw", t_float, "a+b");

    vector<complex<double>> UBorn_I(NbPixUBorn);
    TF2Dcplx_vec_INV(in,out,TF_UBorn_I, UBorn_I, p);
 //SAVCplx(UBorn_I, "Re", "/home/mat/tomo_test/UbornI_extract_holo_re.raw", t_float, "a+b");


    circshift2DCplx(UBorn_I,UBorn,dim2DHA,DecalU_Born);
   // SAVCplx(UBorn, "Re", "/home/mat/tmp/Uborn_re_extract_holo.raw", t_float, "a+b");

}
int coordSpec(vector<complex<double>> TF_UBorn, vector<double> &TF_champMod,Var2D NMAX)
 {
    int cpt_max=0;
    TF_champMod[0]=pow(TF_UBorn[0].real(),2)+pow(TF_UBorn[0].imag(),2);

    for(int cpt=1; cpt<(4*NMAX.x*NMAX.y); cpt++) {
        TF_champMod[cpt]=sqrt(pow(TF_UBorn[cpt].real(),2)+pow(TF_UBorn[cpt].imag(),2));
        if(TF_champMod[cpt]>TF_champMod[cpt_max]) {
        cpt_max=cpt;
        }
    }
   /* double  max_part_reel = TF_UBorn[cpt_max].Re,///sauvegarde de la valeur cplx des  spéculaires
    max_part_imag = TF_UBorn[cpt_max].Im,
    max_module = TF_champMod[cpt_max];*/

    //int kxmi=cpt_max%(2*NMAX.x), kymi=cpt_max/(2*NMAX.y);
    //Var2D posSpec={kxmi,kymi};///coord informatique speculaire
    //recalUBorn={-kxmi,-kymi};
    return cpt_max;
 }

void coupeCplx(vector<complex<double>> src, vector<complex<double>> &dest, Var2D dim_src, Var2D dim_dest, Var2D coin, size_t NumAngle)
{
size_t nbPixSrc=src.size();
//Var2D dim_src={sqrt(nbPixSrc),sqrt(nbPixSrc)};
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

///--------------------------------Sauver-Charger---------------------------------------------
void holo2TF_UBorn(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukey_holo, fftw_complex *in,fftw_complex *out,fftw_plan p_forward_holo)
{

            ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
    size_t dimx=sqrt(NbPix2dROI);

        size_t NbPixROI2d=holo1.size();
        vector<double> holo_shift(NbPixROI2d);
        vector<complex<double>> TF_Holo(NbPixROI2d);
        vector<complex<double>> TFHoloCentre(NbPixROI2d);
       // nbCplx *TF_UBorn_A=new nbCplx[NbPixROI2d];

        for(size_t pixel=0; pixel<NbPixROI2d; pixel++) {
                holo1[pixel]=(double)holo1[pixel]*tukey_holo[pixel];
        }

        ///--------Circshift et TF2D HOLOGRAMME------
        holo_shift=fftshift2D(holo1);
        //SAV2(holo1, "/home/mat/tomo_test/holo_shift_extract_holo.bin",t_float,"a+b");

        TF2Dcplx_vec(in,out,holo_shift, TF_Holo,p_forward_holo);
        TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
        SAVCplx(TF_Holo,"Re","/home/mat/tmp/TFHolo_Re_classique_1024x1024x59.raw",t_float,"a+b");

        coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]

  //  SAVCplx(TFHoloCentre,"Re","/home/mat/tmp/TFHoloCentre_Re_1024x1024x59.raw",t_float,"a+b");

        ///--------Découpe hors axée------------------
       // coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]
}
void holo2TF_UBorn2(vector<double>  &holo1,vector<complex<double>> &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, vector<double> const &tukeyHolo,FFTW_init  &param_fftw2DHolo)
{
    size_t NbPix2dROI=holo1.size(), dimx=sqrt(NbPix2dROI);

        size_t NbPixROI2d=holo1.size();
        vector<double> holo_shift(NbPixROI2d);
        vector<complex<double>> TF_Holo(NbPixROI2d), TFHoloCentre(NbPixROI2d);
    //    nbCplx *TF_UBorn_A=new nbCplx[NbPixROI2d];
    // auto start_tukey = std::chrono::system_clock::now();

        for(size_t pixel=0; pixel<NbPixROI2d; pixel++) {
                holo1[pixel]=(double)holo1[pixel]*tukeyHolo[pixel];
        }
                /*  auto end_tukey = std::chrono::system_clock::now();
           auto elapsed_tukey = end_tukey - start_tukey;
    std::cout <<"Temps tukey= "<< elapsed_tukey.count()/(pow(10,9)) << '\n';*/
        ///--------Circshift et TF2D HOLOGRAMME------
             auto start_fftshift = std::chrono::system_clock::now();
       // holo_shift=fftshift2D(holo1);
        fftshift2D(holo1, holo_shift);

     /*     auto end_fftshift = std::chrono::system_clock::now();
           auto elapsed_fftshift = end_fftshift - start_fftshift;
    std::cout <<"Temps fftshift= "<< elapsed_fftshift.count()/(pow(10,9)) << '\n';*/

   //    SAV2(holo_shift, "/home/mat/tmp/holo_shift_extract_holo.bin",t_float,"a+b");

        TF2D_r2c_symetric(holo_shift,TF_Holo,param_fftw2DHolo);

//        TF2Dcplx_vec(in,out,holo_shift, TF_Holo,p_forward_holo);
  // SAVCplx(TF_Holo,"Re","/home/mat/tmp/TFHolo_Re_r2c_1024x1024x59.raw",t_float,"a+b");
        TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
      //  SAVCplx(TFHoloCentre,"Re","/home/mat/tmp/TFHoloCentre_Re_r2c_1024x1024x59.raw",t_float,"a+b");
 //auto start_coupe= std::chrono::system_clock::now();
        coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NbAngleOk);///Découpe à [-Nxmax,+NXmax]
/*auto end_coupe= std::chrono::system_clock::now();
  auto elapsed_coupe = end_coupe - start_coupe;
    std::cout <<"Temps pour coupe FFT r2c= "<< elapsed_coupe.count()/(pow(10,9)) << '\n';*/
}
/*void holo2TF_UBorn_INPLACE(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukey_holo, fftw_complex *in_out,fftw_plan p_forward_holo)
{

            ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
    size_t dimx=sqrt(NbPix2dROI);

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

        TF2Dcplx_vec_INPLACE(in_out,holo_shift, TF_Holo,p_forward_holo);
        TFHoloCentre=fftshift2D(TF_Holo);//Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
      //  SAVCplx(TFHoloCentre,"Re","/home/mat/tomo_test/TFHoloCentre.raw",t_float,"a+b");

        coupeCplx(TFHoloCentre, TF_UBornTot, dimROI, dim2DHA, coinHA, NumAngle);///Découpe à [-Nxmax,+NXmax]

        //SAVCplx(TFHoloCentre,"Re","/home/mat/tomo_test/TFHoloCentre.raw",t_float,"a+b");

        ///--------Découpe hors axée------------------
       // coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]
}*/
///--------------------------------Sauver-Charger---------------------------------------------
/*void holo2TF_UBorn_old(vector<double> holo1, vector<complex<double>> &TF_UBornTot,Var2D dimROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, vector<double> tukey_holo)
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
}*/

void charger_image2D_OCV(std::vector<double> &imgTab, string imgFile, Var2D coin, Var2D dimROI)
{
        //Mat img(taille.x, taille.y,CV_8UC1);
        Mat img=imread(imgFile, 0);//0=grayscale
        if(! img.data )  // Check for invalid input
       {
              cout <<  "##################### /!\\ ##############################"<<endl;
              cout <<  "Impossible d'ouvrir" <<imgFile<< endl ;
              cout <<  "#########################################################"<<endl;
       }
        Var2D taille={img.cols,img.rows};
        Rect myROI(coin.x, coin.y, dimROI.x, dimROI.y);
        Mat imgCrop = img(myROI);
       /*  for(size_t y=0;y<dimROI.y;y++){
        for(size_t x=0;x<dimROI.x;x++){
                    imgTab[taille.x*y+x]=(double)imgCrop.at<uchar>(y,x);//openCv->tableau
                }
        }*/
        imgTab.assign(imgCrop.begin<uchar>(), imgCrop.end<uchar>());

}



string extract_string(std::string token,  std::string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    string valeur;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
            valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
            cout<<motcle<<"="<<valeurMot<<endl;
            valeur=valeurMot.c_str();
            }
        }
    }
    if(valeur.empty())
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
    fichier.close();
    return valeur;
}

float extract_val(string token,  string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    float valeur=0;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
            valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
            cout<<motcle<<"="<<valeurMot<<endl;
            valeur=atof(valeurMot.c_str());
            }
        }
    }
    fichier.close();
    return valeur;
}

void SAVCplx2(std::vector<complex<double> > var_sav, string partie, std::string chemin, enum PRECISION2 precision, char options[])
{        //double* var_sav = &v[0];

        size_t NbPix2D=var_sav.size();
        unsigned int cpt;
        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision){
        case t_double:{ //64 bit
                double tampon=0;
                if(partie=="Re"||partie=="re"){
                    for(cpt=0; cpt<NbPix2D; cpt++) {
                        tampon=var_sav[cpt].real();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                if(partie=="Im"|| partie=="im"){
                     for(cpt=0; cpt<NbPix2D; cpt++){
                        tampon=var_sav[cpt].imag();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                break;
                }
        case t_float:{//32 bits float
                float tampon=0;

                if(partie=="Re"|| partie=="re"){
                    for(cpt=0; cpt<NbPix2D; cpt++) {
                        tampon=var_sav[cpt].real();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                if(partie=="Im"|| partie=="im"){
                     for(cpt=0; cpt<NbPix2D; cpt++){
                        tampon=var_sav[cpt].imag();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                break;
            }
        }
        fclose(fichier_ID);
}


void SAVCplx(std::vector<complex<double> > var_sav, string partie, std::string chemin, enum PRECISION2 precision, char options[])
{        //double* var_sav = &v[0];

        size_t NbPix2D=var_sav.size();
        unsigned int cpt;
        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision){
        case t_double:{ //64 bit
                double tampon=0;
                if(partie=="Re"||partie=="re"){
                    for(cpt=0; cpt<NbPix2D; cpt++) {
                        tampon=var_sav[cpt].real();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                if(partie=="Im"|| partie=="im"){
                     for(cpt=0; cpt<NbPix2D; cpt++){
                        tampon=var_sav[cpt].imag();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                break;
                }
        case t_float:{//32 bits float
                float tampon=0;

                if(partie=="Re"|| partie=="re"){
                    for(cpt=0; cpt<NbPix2D; cpt++) {
                        tampon=var_sav[cpt].real();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                if(partie=="Im"|| partie=="im"){
                     for(cpt=0; cpt<NbPix2D; cpt++){
                        tampon=var_sav[cpt].imag();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                break;
            }
        }
        fclose(fichier_ID);
}


void SAV2(vector<double> v, std::string chemin, enum PRECISION2 precision, char options[])
{
        size_t NbPix2D=v.size();
        double* var_sav = &v[0];
        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case t_double: //64 bit

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        double tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case t_float://32 bits float

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        float tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case t_int: //32 bit signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case t_uint://32 bit non signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        unsigned int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);

                }
                break;
        case t_char: //8 bits

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        char tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}

void SAV2(double *var_sav,int NbPix2D, std::string chemin, enum PRECISION2 precision, char options[])
{

        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case t_double: //64 bit

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        double tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case t_float://32 bits float

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        float tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case t_int: //32 bit signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case t_uint://32 bit non signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        unsigned int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);

                }
                break;
        case t_char: //8 bits

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        char tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}



void SAV_Tiff2D(double *var_sav, string chemin, int dim)
{
    TIFF *tif= TIFFOpen(chemin.c_str(), "a");
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    tsize_t strip_size = TIFFStripSize (tif);
    tstrip_t strips_num = TIFFNumberOfStrips (tif);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt];
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}

void SAV_Tiff2D(std::vector<double> var_sav, string chemin, double taille_pixel)
{   int dim=pow(var_sav.size(),0.5);
    float xres, yres;
    uint16  res_unit;
    TIFF *tif_out= TIFFOpen(chemin.c_str(), "w");//"w"->ecraser, "a"->ajouter
    TIFFSetField (tif_out, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif_out, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif_out, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif_out, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif_out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif_out, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif_out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    //fixer la résolution
    xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre)=0.01m/taille_pixel en metre
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(tif_out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(tif_out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(tif_out, TIFFTAG_RESOLUTIONUNIT, res_unit);

    tsize_t strip_size = TIFFStripSize (tif_out);
    tstrip_t strips_num = TIFFNumberOfStrips (tif_out);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt];


        }
        TIFFWriteEncodedStrip (tif_out, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif_out);
    TIFFClose(tif_out);
}


void SAV_Tiff2D(std::vector<complex<double>> var_sav, string partie, string chemin, double taille_pixel)
{   int dim=pow(var_sav.size(),0.5);
    float xres, yres;
    uint16  res_unit;
    TIFF *tif_out= TIFFOpen(chemin.c_str(), "w");//"w"->ecraser, "a"->ajouter
    TIFFSetField (tif_out, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif_out, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif_out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif_out, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif_out, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif_out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif_out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif_out, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif_out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif_out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    //fixer la résolution
    xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre)=0.01m/taille_pixel en metre
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(tif_out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(tif_out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(tif_out, TIFFTAG_RESOLUTIONUNIT, res_unit);

    tsize_t strip_size = TIFFStripSize (tif_out);
    tstrip_t strips_num = TIFFNumberOfStrips (tif_out);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            if(partie=="Re" || partie=="re")
            strip_buf[col]=(float)var_sav[cpt].real();
            else
            {
                if(partie=="Im"|| partie=="im")
                strip_buf[col]=(float)var_sav[cpt].imag();
                else
                cout<<"La chaine fixant la partie doit être Re,re, Im ou im"<<endl;
            }
        }
        TIFFWriteEncodedStrip (tif_out, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif_out);
    TIFFClose(tif_out);
}

void SAV3D_Tiff(vector<complex <double>> var_sav, string partie, string chemin, double taille_pixel)
{
    int dim=round(std::pow(var_sav.size(), 1.0/3.0));
    uint32 image_width, image_height, dimz;
    float xres, yres;
    uint16 spp, bpp, photo, res_unit,zpage;
    TIFF *out;
    int x, y;
    float buffer2D[dim * dim];
    out = TIFFOpen(chemin.c_str(), "w");
    if (!out)
            fprintf (stderr, "Can't open  for writing\n");

    image_width = dim;
    image_height = dim;
    dimz=dim;
    spp = 1; /* Samples per pixel */
    bpp = 32; /* Bits per sample */
   // photo = PHOTOMETRIC_MINISBLACK;

    for (int num_page = 0; num_page < dim; num_page++)//z=page
    {
        for (y = 0; y < dim; y++)
            for(x = 0; x < dim; x++){
                    if(partie=="Re" || partie=="re"){
                      buffer2D[y * dim + x] = (float)var_sav[y * dim + x+num_page*dim*dim].real();
                    }
                    else{
                        if(partie=="Im" || partie=="im"){
                        buffer2D[y * dim + x] = (float)var_sav[y * dim + x+num_page*dim*dim].imag();
                        }

                        else
                            cout<<"Partie non identifiée : Re, re, Im ou im"<<endl;
                    }

            }


        TIFFSetField(out, TIFFTAG_IMAGEWIDTH, image_width / spp);
        TIFFSetField(out, TIFFTAG_IMAGELENGTH, image_height);
        TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bpp);
        TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, spp);
        TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
       // TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photo);
        TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
        TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); //image en Floating point
        /* It is good to set resolutions too (but it is not nesessary) */
        xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre, on multiplie par 0.01 pour tout passer en mètre)
        res_unit = RESUNIT_CENTIMETER;
        TIFFSetField(out, TIFFTAG_XRESOLUTION, xres);
        TIFFSetField(out, TIFFTAG_YRESOLUTION, yres);
        TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);
        /* We are writing single page of the multipage file */
        TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        /* Set the page number */
        TIFFSetField(out, TIFFTAG_PAGENUMBER, num_page, dimz);

        for (y = 0; y < image_height; y++)
            TIFFWriteScanline(out, &buffer2D[y * image_width], y, 0);

        TIFFWriteDirectory(out);
    }

    TIFFClose(out);
}
bool is_readable( const std::string & file )
{
    std::ifstream fichier( file.c_str() );
    return !fichier.fail();
}
