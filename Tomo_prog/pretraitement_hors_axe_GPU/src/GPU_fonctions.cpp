#include "GPU_fonctions.h"
#include "fonctions.h"
#include "struct.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>//imread

using namespace af;
using namespace std;
using namespace cv;

///normalement dans fonctions.h
af::array init_kvect_shiftX_GPU(int dim2DHA)
{
  af::array  kvect_shiftX=af::array(dim2DHA,dim2DHA,f64);
//for(int y=0;y<dim2DHA-1;y++)
gfor(seq y,dim2DHA)
 kvect_shiftX(af::seq(0,dim2DHA-1),y)=(af::seq(0,dim2DHA-1)-round(dim2DHA/2))/dim2DHA;
 kvect_shiftX=shift(kvect_shiftX,dim2DHA/2,dim2DHA/2);
 return kvect_shiftX;
}
///normalement dans fonctions.h
af::array init_kvect_shiftY_GPU(int dim2DHA)
{
  af::array  kvect_shiftY=af::array(dim2DHA,dim2DHA,f64);
for(int x=0;x<dim2DHA-1;x++)
//gfor(seq x,dim2DHA)
 kvect_shiftY(x,af::seq(0,dim2DHA-1))=(af::seq(0,dim2DHA-1)-round(dim2DHA/2))/dim2DHA;
 kvect_shiftY=shift(kvect_shiftY,dim2DHA/2,dim2DHA/2);
 return kvect_shiftY;
}
///dans IO_fonctions. Tres lent, voir chargeviaOCV (10 fois + rapide)
void chargeImageGPU(af::array &img, string nomImg, Var2D coin){
    Var2D dimCrop{img.dims(0),img.dims(1)};//crop @ dimension of the hologram

  img=crop(loadImage(nomImg.c_str(), false),dimCrop,coin);
    //img=af::rotate(img, 90, 1,AF_INTERP_NEAREST);
   // return img;
}

///dans IO_fonctions
void chargeImageGPU_via_OCV(af::array &img_array, string nomImg, Var2D coin){

    Var2D dimCrop{img_array.dims(0),img_array.dims(1)};//crop @ dimension of the hologram
    std::vector<double> imgTab(dimCrop.x*dimCrop.y);
    Mat img=imread(nomImg, 0);//0=grayscale
    if(! img.data )   // Check for invalid input
    {
        cout <<  "##################### /!\\ ##############################"<<endl;
        cout <<  "Impossible d'ouvrir" <<nomImg<< endl ;
        cout <<  "#########################################################"<<endl;
    }
    Var2D taille={img.cols,img.rows};
    Rect myROI(coin.x, coin.y, dimCrop.x, dimCrop.y);
    Mat imgCrop = img(myROI);
    /*  for(size_t y=0;y<dimROI.y;y++)
     for(size_t x=0;x<dimROI.x;x++){
                 imgTab[taille.x*y+x]=(double)imgCrop.at<uchar>(y,x);//openCv->tableau
             }    */
    imgTab.assign(imgCrop.begin<uchar>(), imgCrop.end<uchar>());
    img_array=af::array(img.cols,img.rows,imgTab.data());
}

void display(af::array img){
    const static int width = img.dims(0), height = img.dims(1);
Window window(width, height, "Title");
do{
 window.image(img.as(f32)/max<float>(img));
//drawing functions here
} while( !window.close() );
}
/*void SAV_afArray(af::array myArray){
   vector<std::complex<double>> TF_holo_vctr(myArray.dims(0)*myArray.dims(1));
    myArray.host(TF_holo_vctr.data());
    SAVCplx(TF_holo_vctr,"/home/mat/tomo_test/Tf_holo_GPU.raw",t_float,"a+b");
}*/

///------------------------GPU functions----------------------------------

///dans FFT_fonctions
//calculate 2D Tukey window whose width is controled by alpha (to be applied before fft)
af::array tukey2D_GPU(size_t dimx,size_t dimy, float alpha){
    size_t N=dimx;
    size_t nbPix2D=dimx*dimy;
    af::array tuk2D(dimx,dimy);
    af::array tuk1Dx=constant(1,dimx);

    ///partie gauche de la courbe
    size_t borne1=round(alpha*(N-1)/2);
    af::array     cpt1=seq(0,borne1);//ne fonctionne pas avec range"
    tuk1Dx(seq(0,borne1))=0.5*(1+cos(3.1415*(2*seq(0,borne1)/(alpha*(N-1))-1)));
    ///partie droite de la courbe
    size_t borne2=round((N-1)*(1-alpha/2));
    af::array cpt2=seq(borne2,dimx);
    tuk1Dx(cpt2)=0.5*(1+cos(3.1415*(2*cpt2/(alpha*(N-1))-2/alpha+1)));

    gfor(seq y,dimx)      //  for(size_t x=0; x<dimx; x++)
            tuk2D(span,y)=tuk1Dx(span)*tuk1Dx(y);

    return tuk2D;//segfault si pas de return !
}
af::array crop(af::array img_orig,Var2D dim_crop, Var2D coord_coin){
    int cx=coord_coin.x,cy=coord_coin.y;
    af::array img_crop(dim_crop.x,dim_crop.y,f64);
    img_crop=img_orig(seq(cx,cx+dim_crop.x-1),seq(cy,cy+dim_crop.y-1));
    return img_crop;
}
void holo2TF_UBorn_GPU(af::array  &holo1,af::array &TF_UBornTot,Var2D dimROI,Var2D dim2DHA,Var2D coinHA, size_t NbAngleOk, af::array  const &tukeyHolo, af::array const &tukeyBorn)
{
  //   timer start1=timer::start();
    const int width=holo1.dims(0), height=holo1.dims(1), FFT_NORM=width*height;
    af::array TF_Holo=af::array(width, height,c64);
    holo1=holo1*tukeyHolo;
    TF_Holo=fft2(shift(holo1,width/2,height/2))/FFT_NORM;
    TF_Holo.eval();
    TF_Holo=shift(TF_Holo,width/2,height/2);
    TF_UBornTot(span,span,NbAngleOk)=crop(TF_Holo,dim2DHA,coinHA)*tukeyBorn;//<- fenêtre de tukey peu utile car appliquée surtout dans le bruit, à R>NA, partie éliminée lors de la projection sur sphere d'Ewald
   //   printf("elapsed seconds: %g\n", timer::stop(start1));
}
