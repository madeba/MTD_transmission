#include "IO_fonctions.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>

using namespace std;
using namespace cv;
void charger_image2D_OCV(std::vector<double> &imgTab, string imgFile, Var2D coin, Var2D dimROI)
{
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
       /*  for(size_t y=0;y<dimROI.y;y++)
        for(size_t x=0;x<dimROI.x;x++){
                    imgTab[taille.x*y+x]=(double)imgCrop.at<uchar>(y,x);//openCv->tableau
                }    */
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
    if(valeurMot.empty())
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
    fichier.close();

    return valeur;
}

void SAVCplx(std::vector<complex<double>> const &var_sav, string partie, std::string chemin, enum PRECISION2 precision, char options[])
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


void SAV2(vector<double> &v, std::string chemin, enum PRECISION2 precision, char options[])
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

void SAV_Tiff2D(std::vector<double> const &var_sav, string chemin, double taille_pixel)
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


void SAV_Tiff2D(std::vector<complex<double>> const & var_sav, string partie, string chemin, double taille_pixel)
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

void SAV3D_Tiff(vector<complex <double>> const &var_sav, string partie, string chemin, double taille_pixel)
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
