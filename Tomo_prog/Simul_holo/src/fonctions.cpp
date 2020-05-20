
#include "fonctions.h"
#include <tiffio.h>
#include <opencv/cxcore.h>
#include <chrono>
//#include <highgui.h>
#include <opencv/cv.h>
#define PI 3.14159
#include "projet.h"
#include "OTF.h"

using namespace cv;
using namespace std;
//calcule le spectre d'hologramme à partir d'une sphere d'ewald de centre (Point2D spec) et du spectre de l'objet.
void calcPhase_mpi_pi_atan2(vector<complex<double>> obj, vector<double> &phaseMod2pi)///calcul phase -PI-PI
{
for(int cpt=0;cpt<obj.size();cpt++)
phaseMod2pi[cpt]=atan2(obj[cpt].imag(),obj[cpt].real());
}
///prend un point 2D en entrée et retourne le suivant le long de la rosace//
///pas optimale :  refait le calcul de tous les paramètres (dont dtheta) à chaque appel à la fonction...
Point2D maj_fleur(Point2D Vin, float rho, int nbHolo, double *theta, manip m1)//tension x,y/rayon/nbHolo,/angle
{   int dim2D=m1.dim_Uborn;
    Point2D Vout(0,0,dim2D);
    float t=*theta;//transfer theta->t
    float LongTot=17.15*rho, delta_curv=LongTot/nbHolo;//longueur totale, delta_abs_curv "régulier"

    float delta_theta=(1/rho*delta_curv)/(sqrt(15*sin(4*t)*sin(4*t)+1));

    t=delta_theta+t; ///abscisse curv  point suivant

    Vout.y=rho*cos(4*t)*sin(t);
    Vout.x=rho*cos(4*t)*cos(t);
    *theta=t;///maj theta prog principal
    return Vout;
}
///extraction d'une calotte centrée sur (spec.x,spec.y) d'Ewald dans le spectre 3D de l'objet
void calcHolo(Point2D spec,std::vector<std::complex<double>> const &TF_vol3D,std::vector<std::complex<double>> &TF_hologramme, manip m1)
{
int Nmax=m1.NXMAX,rayon=m1.R_EwaldPix, dim2D=m1.dim_Uborn, dim3D=m1.dim_final;
double fmcarre=Nmax*Nmax;///frequence max
double rcarre=rayon*rayon;///rayon Ewald au carré
double kv=2*PI/m1.lambda_v, k0=kv*m1.n0;
double sdz=0;
complex<double> cteInd2Pot(-2*kv*kv*m1.n0,0);
complex<double> ctePot2UBorn(0,-PI/k0);
//complex <double> coef_global(0,-2*kv*PI);

complex<double> cteNormalisation(-1/(2*PI),0);
//std::vector<std::complex<double>> support_holo(dim*dim*dim);
    Point3D fi(spec,sqrt(rcarre-spec.x*spec.x-spec.y*spec.y),dim3D);
   // cout<<"fix="<<fi.x<<", fiy="<<fi.y<<", fiz="<<fi.z<<endl;
    Point3D fobj(0,0,0,dim3D);
    Point3D fd(0,0,0,dim3D);//
    Point2D fLat(0,0,dim2D);//kx, ky

    for(fd.x=-Nmax; fd.x<=Nmax; fd.x++){
        for(fd.y=-Nmax; fd.y<=Nmax; fd.y++){
                sdz=(double)sqrt(rcarre-(fd.x)*(fd.x)-(fd.y)*(fd.y))/(double)rayon;
                fd.z=round(sqrt(rcarre-(fd.x)*(fd.x)-(fd.y)*(fd.y)));

                if(fd.x*fd.x+fd.y*fd.y<fmcarre){
                    fobj=fd-fi;
                    fLat.set_xy(fd.x,fd.y);
                    TF_hologramme[fLat.coordI().cpt2D()]=cteNormalisation*ctePot2UBorn*cteInd2Pot/(sdz)*TF_vol3D[fobj.coordI().cpt3D()];
            }
        }
    }
}

void SAV2D_Tiff(std::vector<complex<double>> var_sav, string partie, string chemin, double taille_pixel)
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
                        if(partie=="Mod"|| partie=="mod")
                            strip_buf[col]=(float)abs(var_sav[cpt]);
                            else
                                cout<<"La chaine fixant la partie doit être Re/re, Im/im, Mod/mod"<<endl;
            }
        }
        TIFFWriteEncodedStrip (tif_out, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif_out);
    TIFFClose(tif_out);
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


void SAV3D_Tiff(vector<complex <double>> var_sav, string partie, string chemin, double taille_pixel)
{

    int dim=round(std::pow(var_sav.size(), 1.0/3.0));
    //cout<<"dim savtiff3d ======"<<dim<<endl;
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
        xres = yres = 0.01/taille_pixel; //nbpixel par resunit (par centimetre)
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

void circshift3DCplx(vector<complex<double>> volume3D, vector<complex<double>> volume3D_shift, Var3D dimFinal3D, Var3D decal3D)//(entrée, sortie_decalee, dim, decalage)
{
        decal3D.x=decal3D.x%dimFinal3D.x;//éliminer les "modulos"
        decal3D.y=decal3D.y%dimFinal3D.y;
        decal3D.z=decal3D.z%dimFinal3D.y;

        unsigned short int xi,yi,zi=0;
        short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
        const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;

        for(zi=0; zi<dimFinal3D.z; zi++) {
                if(zi+decal3D.z>dimFinal3D.z-1) { //dépassement à droite
                        z2=zi+decal3D.z-dimFinal3D.z;
                } else {
                        if(zi+decal3D.z<0) { //dépassement à gauche
                                z2=dimFinal3D.z+(decal3D.z+zi);
                        } else {
                                z2=zi+decal3D.z;
                        }
                }
                int nb_pixelz_decal=z2*taille_plan;
                unsigned int nb_pixelz=zi*taille_plan;
                for(yi=0; yi<dimFinal3D.y; yi++) {
                        if(yi+decal3D.y>dimFinal3D.y-1) { //dépassement à droite
                                y2=yi+decal3D.y-dimFinal3D.y;
                        } else {
                                if(yi+decal3D.y<0) { //dépassement à gauche
                                        y2=dimFinal3D.y+(decal3D.y+yi);
                                } else {
                                        y2=yi+decal3D.y;
                                }
                        }
                        int nb_lignes=yi*dimFinal3D.x;
                        int nb_lignes_decal=y2*dimFinal3D.x;
                        //#pragma omp parallel for private(xi)
                        for(xi=0; xi<dimFinal3D.x; xi++) {
                                if(xi+decal3D.x>dimFinal3D.x-1) { //dépassement à droite
                                        x2=xi+decal3D.x-dimFinal3D.x;
                                } else {
                                        if(xi+decal3D.x<0) { //dépassement à gauche
                                                x2=dimFinal3D.x+(decal3D.x+xi);
                                        } else {
                                                x2=xi+decal3D.x;
                                        }
                                }
                            volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].real(volume3D[nb_pixelz+nb_lignes+xi].real() );
                            volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].imag(volume3D[nb_pixelz+nb_lignes+xi].imag() );
                        }
                       // #pragma omp barrier
                }
        }
}

//#############################################################"
void genere_bille(vector<complex<double>> &vol_bille,  Point3D centre, size_t rayon_bille, complex<double> indiceObj, complex<double> indiceImmersion, size_t dim_espace)
{
double zmax=centre.z+rayon_bille;
double zmin=centre.z-rayon_bille;
double xmax=centre.x+rayon_bille;
double xmin=centre.x-rayon_bille;
double ymax=centre.y+rayon_bille;
double ymin=centre.y-rayon_bille;
int nbPixBille=0;
for(size_t cpt=0;cpt<vol_bille.size();cpt++){
     vol_bille[cpt].real(indiceImmersion.real());
     vol_bille[cpt].imag(indiceImmersion.imag());
}

int rayon_carre=pow(rayon_bille,2);

for(double z=zmin;z<zmax;z++)
    for(double y=ymin;y<ymax;y++){
         for(double x=xmin;x<xmax;x++){
               // Point3D mon_point(x,y,z);
            int cpt=z*dim_espace*dim_espace+y*dim_espace+x;

            if((x-centre.x)*(x-centre.x)+(y-centre.y)*(y-centre.y)+(z-centre.z)*(z-centre.z)<rayon_carre){
                vol_bille[cpt].real(indiceObj.real());
                vol_bille[cpt].imag(indiceObj.imag());
                nbPixBille++;
            }
         }
        }
       // cout<<"nbPixBille="<<nbPixBille<<endl;

}





void circshift3DCplx(vector<complex<double>> volume3D, vector<complex<double> > &volume3D_shift, Var3D decal3D)//(entrée, sortie_decalee, dim, decalage)
{
        unsigned int dim=pow(volume3D.size(),1.0/3.0);//récuperer les dimensions de l'espace 3D
        Var3D dimFinal3D={dim,dim,dim};
        decal3D.x=decal3D.x%dimFinal3D.x;//éliminer les "modulos"
        decal3D.y=decal3D.y%dimFinal3D.y;
        decal3D.z=decal3D.z%dimFinal3D.y;

        unsigned short int xi,yi,zi=0;
        short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
        const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;

        for(zi=0; zi<dimFinal3D.z; zi++) {
                if(zi+decal3D.z>dimFinal3D.z-1) { //dépassement à droite
                        z2=zi+decal3D.z-dimFinal3D.z;
                } else {
                        if(zi+decal3D.z<0) { //dépassement à gauche
                                z2=dimFinal3D.z+(decal3D.z+zi);
                        } else {
                                z2=zi+decal3D.z;
                        }
                }
                int nb_pixelz_decal=z2*taille_plan;
                unsigned int nb_pixelz=zi*taille_plan;
                for(yi=0; yi<dimFinal3D.y; yi++) {
                        if(yi+decal3D.y>dimFinal3D.y-1) { //dépassement à droite
                                y2=yi+decal3D.y-dimFinal3D.y;
                        } else {
                                if(yi+decal3D.y<0) { //dépassement à gauche
                                        y2=dimFinal3D.y+(decal3D.y+yi);
                                } else {
                                        y2=yi+decal3D.y;
                                }
                        }
                        int nb_lignes=yi*dimFinal3D.x;
                        int nb_lignes_decal=y2*dimFinal3D.x;
                        //#pragma omp parallel for private(xi)
                        for(xi=0; xi<dimFinal3D.x; xi++) {
                                if(xi+decal3D.x>dimFinal3D.x-1) { //dépassement à droite
                                        x2=xi+decal3D.x-dimFinal3D.x;
                                } else {
                                        if(xi+decal3D.x<0) { //dépassement à gauche
                                                x2=dimFinal3D.x+(decal3D.x+xi);
                                        } else {
                                                x2=xi+decal3D.x;
                                        }
                                }
                            volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].real(volume3D[nb_pixelz+nb_lignes+xi].real());
                            volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].imag(volume3D[nb_pixelz+nb_lignes+xi].imag());
                        }
                       // #praSAVCplxgma omp barrier
                }
        }
}




Mat vector2mat(std::vector<double> &src, unsigned short int dim)
{
Mat dst(dim,dim,CV_32FC1);
size_t cpt=0;
 for(size_t x=0;x<dst.rows;x++)
        for(size_t y=0;y<dst.cols;y++){
            cpt=x+y*dst.rows;
            dst.at<float>(y, x)=src[cpt];
        }
 return dst;
}

void mat2vector(Mat src, std::vector<double> &dst)
{
 for(size_t x=0;x<src.rows;x++)
        for(size_t y=0;y<src.cols;y++){
            int cpt=x+y*src.rows;
            dst[cpt]=src.at<float>(y, x);
        }
}
void decal2DCplxGen(vector<complex<double>> &entree, vector<complex<double>> &result, Var2D decal)
{
    Var2D dim={sqrt(entree.size()),sqrt(entree.size())};
    //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;
        if(decal.x<0)
        decal.x=dim.x+decal.x;
        if(decal.y<0)
        decal.y=dim.y+decal.y;
       for(int yi=0; yi<dim.y-decal.y; yi++) {
                for(int xi=0; xi<dim.x-decal.x; xi++)
                { //cout<<"xi,yi="<<xi<<","<<yi<<endl;
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }
              for(int yi=dim.y-decal.y; yi<dim.y; yi++) {
                for(int xi=0; xi<dim.x-decal.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }
}




void SAV2(vector<double> v, std::string chemin, enum PRECISION precision, char options[])
{

        unsigned int NbPix2D=v.size();

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

void SAVCplx(std::vector<complex<double> > var_sav, string partie, std::string chemin, enum PRECISION precision, char options[])
{        //double* var_sav = &v[0];
        unsigned int cpt;
        unsigned int NbPix=var_sav.size();
//        cout<<"taille volume="<<NbPix<<endl;
        FILE *fichier_ID;
        fichier_ID= fopen(chemin.c_str(), options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision){
        case t_double:{ //64 bit
                double tampon=0;
                if(partie=="Re"|| partie=="re"){
                    for(cpt=0; cpt<NbPix; cpt++) {
                        tampon=var_sav[cpt].real();

                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                if(partie=="Im"|| partie=="im"){
                     for(cpt=0; cpt<NbPix; cpt++){
                        tampon=var_sav[cpt].imag();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                break;
                }
        case t_float:{//32 bits float
                float tampon=0;

                if(partie=="Re"|| partie=="re"){
                    for(cpt=0; cpt<NbPix; cpt++) {
                        tampon=var_sav[cpt].real();
                          //cout<<"tampon="<<tampon<<endl;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                if(partie=="Im"|| partie=="im"){
                     for(cpt=0; cpt<NbPix; cpt++){
                        tampon=var_sav[cpt].imag();
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                        }
                }
                break;
            }
        }
        fclose(fichier_ID);
}


string type2str(int type) {
  string r;

  uchar depth = type & CV_MAT_DEPTH_MASK;
  uchar chans = 1 + (type >> CV_CN_SHIFT);

  switch ( depth ) {
    case CV_8U:  r = "8U"; break;
    case CV_8S:  r = "8S"; break;
    case CV_16U: r = "16U"; break;
    case CV_16S: r = "16S"; break;
    case CV_32S: r = "32S"; break;
    case CV_32F: r = "32F"; break;
    case CV_64F: r = "64F"; break;
    default:     r = "User"; break;
  }

  r += "C";
  r += (chans+'0');

  return r;
}

///#########lecture propre d'un fichier binaire  image carré ou cubique, connaissant sa taille et sa largeur.
//type de donnée possible : float ou double.

void lire_bin_vector(string chemin, vector<double> &dst, unsigned short int dim, unsigned  int nbPix)
{
       unsigned long lTaille;
       size_t nb_elmnt_lu;
      //unsigned short int dimData=precision/8;//taille en octet d'un element.

       FILE* pFichier = NULL;

       pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire
        cout<<pFichier<<endl;


       if(pFichier==NULL){
           fputs("Impossible d'ouvrir le fichier\n",stderr);
           cout<<chemin<<endl;
           exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
       }
       else{
           fseek(pFichier,0,SEEK_END);//trouver la fin de fichier
           lTaille = ftell (pFichier);//retourne la position courante (en octet) du curseur de fichier : ici, position de la fin du fichier
          // printf("taille trouvée en octet par ftell %li, taille estimée : %i\n",lTaille,NbPix*dimData);//
           rewind(pFichier);

           int type_donnee=lTaille*8/(nbPix);
            cout<<"type_donnee"<<type_donnee<<endl;
            if(type_donnee==32){
                float resultat;
                size_t taille_case=type_donnee/8;

                for(int cpt=0;cpt<nbPix;cpt++){
                nb_elmnt_lu = fread (&resultat,1,taille_case,pFichier);//lecture
                dst[cpt]=(double)resultat;
                }
            }
            if(type_donnee==64){
                double resultat;//pas important?
                size_t taille_case=type_donnee/8;

                for(int cpt=0;cpt<nbPix;cpt++){
                nb_elmnt_lu = fread(&resultat,1,taille_case,pFichier);//lecture
                dst[cpt]=resultat;
                }
            }

        }

            fclose(pFichier);
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


//convolution par une sphere d'Ewald de centre spec
void Conv_Ewald(Point2D spec,std::vector<std::complex<double>> TF_vol3D,std::vector<std::complex<double>> &TF_conv3D, manip m1)
{
int Nmax=m1.NXMAX,rayon=m1.R_EwaldPix, dim=m1.dim_Uborn;
double n0=m1.n0;
double fmcarre=Nmax*Nmax;///frequence max
double rcarre=rayon*rayon;///rayon Ewald au carré


complex<double> cteInd2pot(-2*(2*PI)*(2*PI)*n0/(m1.lambda_v*m1.lambda_v),0);//constante pour passer de l'indice 3D vers le potentiel objet 3D
complex<double> ctePot2UBorn(0,-m1.lambda_v/(2*m1.n0));

std::vector<std::complex<double>> support_holo(dim*dim*dim);
    Point3D fi(spec,sqrt(rcarre-spec.x*spec.x-spec.y*spec.y),dim);
    Point3D fobj(0,0,0,dim);
    Point3D fd(0,0,0,dim);//dim espace erronée mais sinon problème soustraction
    Point2D fLat(0,0,dim);//kx, ky
    double sdz=0;
    //double kiz_d=sqrt(rayon*rayon-spec.x*spec.x-spec.y*spec.y);
   // cout<<"kxix2+kiy2="<<sqrt(spec.x*spec.x+spec.y*spec.y)<<endl;
   // cout<<"Nmax="<<-Nmax<<endl;
    for(fd.x=-Nmax; fd.x<=Nmax; fd.x++){
        for(fd.y=-Nmax; fd.y<=Nmax; fd.y++){
                fd.x=round(fd.x);
                fd.y=round(fd.y);
                fd.z=round(sqrt(rcarre-(fd.x)*(fd.x)-(fd.y)*(fd.y)));
                sdz=(double)fd.z/rayon;
              // cout<<round(kd.x*kd.x)+round(kd.y*kd.y)<<endl;
               if(round(fd.x*fd.x)+round(fd.y*fd.y)<fmcarre){
                fobj=fd-fi;
                fLat.set_xy(fobj.x,fobj.y);
                if(fd.z==0)
                    cout<<"fd.z=0!!"<<endl;
                TF_conv3D[fobj.coordI().cpt3D()]=cteInd2pot*ctePot2UBorn/sdz*TF_vol3D[fobj.coordI().cpt3D()];
            }
        }
    }
}


