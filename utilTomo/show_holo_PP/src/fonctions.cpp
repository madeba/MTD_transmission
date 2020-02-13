#include <fstream>//ifstream
#include <vector>
#include <complex>
#include <fftw3.h>
#include "struct.h"
//#include <cv.h>
//#include <highgui.h>//imread
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

void decal2DCplxGen(vector<complex<double>> entree, vector<complex<double>> &result, Var2D decal)
{
            //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
            //cout<<"decal.x,y="<<decal.x<<","<<decal.y<<endl;
        Var2D dim={sqrt(entree.size()),sqrt(entree.size())};
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
void crop2DCplx(vector<complex<double>> src, vector<complex<double>> &dest, Var2D coin)
{
    size_t nbPixSrc=src.size();
    Var2D dim_src={sqrt(nbPixSrc),sqrt(nbPixSrc)};
    Var2D dim_dest={sqrt(dest.size()),sqrt(dest.size())};

    //Var2D dim_src={sqrt(nbPixSrc),sqrt(nbPixSrc)};
    size_t X_dest,Y_dest, cpt_dest1D,
            X_src, Y_src, cpt_src1D;
    //balyage destination en 2NXMAX*2NXMAX

    for(X_dest=0; X_dest<dim_dest.x; X_dest++){
        for(Y_dest=0; Y_dest<dim_dest.y; Y_dest++){
            cpt_dest1D=X_dest+Y_dest*dim_dest.x;;///coord 1D destination
            X_src=coin.x+X_dest;///coord X src
            Y_src=coin.y+Y_dest;///coord Y src
           //   cout<<"X,yrsc=("<<X_src<<","<<Y_src<<")"<<endl;
            cpt_src1D=X_src+Y_src*dim_src.x;///coord 1D source
            dest[cpt_dest1D]=src[cpt_src1D];///copie src->dest
        }
    }
}
///--------------------------------Sauver-Charger---------------------------------------------
void holo2TF_UBorn(vector<double>   holo1, vector<complex<double>> &TF_UBorn,Var2D coinHA, FFT_encaps &tf2dROI)
{

    ///--------------Init FFTW-------------------------------------------------
    size_t NbPix2dROI=holo1.size();
    size_t dimx=sqrt(NbPix2dROI);
    int dim_Uborn=sqrt(TF_UBorn.size());
    Var2D dimROI= {dimx,dimx};
    Var2D dim2DHA= {dim_Uborn,dim_Uborn};

    size_t NbPixROI2d=holo1.size();
  //  vector<double> holo_shift(NbPixROI2d);
    vector<complex<double>> TF_Holo(NbPixROI2d);
    vector<complex<double>> TFHoloCentre(NbPixROI2d);

    ///--------Circshift et TF2D HOLOGRAMME------
    TF2D(fftshift2D(holo1), TF_Holo,tf2dROI,1);


    crop2DCplx(fftshift2D(TF_Holo), TF_UBorn,  coinHA);///Découpe à [-Nxmax,+NXmax]
    ///--------Découpe hors axée------------------
}

///--------------------------------Sauver-Charger---------------------------------------------


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
        for(size_t x=0;x<dimROI.x;x++){
                  for(size_t y=0;y<dimROI.y;y++){
                    imgTab[taille.x*y+x]=(double)imgCrop.at<uchar>(y,x);//openCv->tableau
                }
        }

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

void SAVCplx2(std::vector<complex<double> > var_sav, string partie, std::string chemin, enum PRECISION precision, char options[])
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


void SAVCplx(std::vector<complex<double> > var_sav, string partie, std::string chemin, enum PRECISION precision, char options[])
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


void SAV2(vector<double> v, std::string chemin, enum PRECISION precision, char options[])
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

void SAV2(double *var_sav,int NbPix2D, std::string chemin, enum PRECISION precision, char options[])
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

bool is_readable( const std::string & file )
{
    std::ifstream fichier( file.c_str() );
    return !fichier.fail();
}
