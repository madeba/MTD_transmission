//Manuscript Title: Tomographic diffracitive microscopy and multiview profilometry
//                  with flexible aberration correction
//
//Authors: H. Liu, J. Bailleul, B. Simon, M. Debailleul, B. Colicchio, and O. Haeberlé
//Program title: CorrectionAberration
//Distribution format: Zip
//Download link: to come
//Licensing provisions: CeCILL-B
//Journal reference: Applied optics, Vol. 53, issue 4, p748-755 (2014)
//Programming language: C++
//Computer: General computers
//Operating Synstem: All
//Keywords: Aberration compensation; Holography; Tomographic imaging; Fourier optics and signal processing;
//External routines: Opencv
//
//Nature of problem: Automatic aberration compensation in holographic imaging
//Solution method: Polynomial approximation of the aberration function,
//                 by using a random point selection in the background
//The object is excluded thanks to a mask, which can be specified as "Image_mask.pgm" in the data folder
//Running time: 0,3s with an Intel quadcore i7 processor and 16 Gbytes of memory
///warning ; the correction is not automatic in this version, you must specify a mask "Image_mask.pgm" in the data folder
///without this mask, the object will be included in the "background" calculation
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include <chrono>
//using namespace cv;
using namespace std;

/// Function header
#include "Correction_aberration.h"
#include "FFT_fonctions.h"

///function, used to correct only for reference amplitude, useless with a blank acquisition
void correctAmpliRef(std::vector<std::complex<double>> &UBorn,std::vector<double> const &ampli_ref)
{   size_t nbPix=UBorn.size();
    for(size_t cpt=0;cpt<nbPix;cpt++)
    {
        UBorn[cpt]=UBorn[cpt]/ampli_ref[cpt];
    }

}

///Pre-calculate power of x and y in the polynom , which is  a constant table, to avoid calculation in a loop
void initCorrAber(std::string Chemin_mask, Mat const & mask_aber, size_t degre_poly, Var2D dim2DHA, Mat &polynome_to_fit, Mat &polynomeUs_to_fit)
{
CalcPolyUs_xy(degre_poly, mask_aber, dim2DHA, polynomeUs_to_fit);
CalcPoly_xy(degre_poly, dim2DHA, polynome_to_fit);
}

///load object mask used to correct aberrations
Mat init_mask_aber(string Chemin_mask, string Chemin_acquis, Var2D dim2DHA)
{
    Mat mask = imread(Chemin_mask, 0);
    ///Charger masque aberration

    if(! mask.data ){
        cout <<  "Masque : "<<Chemin_mask <<" non ---trouvé création masque unité" << std::endl ;
        mask=255*Mat::ones(dim2DHA.x,dim2DHA.y, CV_8UC1);
      /*  vector<double> masqueTukeyHA(tukey2D(dim2DHA.x,dim2DHA.y,0.05));//on utilise alpha/2 sinon bord du masque un peu épais
        for(int x=0; x<dim2DHA.x; x++)
            for(int y=0; y<dim2DHA.y; y++){
                int cpt=x+y*dim2DHA.x;
                if(mask.at<uchar>(y,x)*masqueTukeyHA[cpt]<255)//appliquer le masque de tukey puis binariser
                    mask.at<uchar>(y,x)=0;
        }
     imwrite(Chemin_acquis+"/Image_mask.pgm",mask);*/
    }
    else{
        cout<<"chargement du masque"<<Chemin_mask<<endl;
    }
    if(mask.rows!=dim2DHA.x){
        cout<<"Problème aberrations : masque "<<Chemin_mask<< " de largeur "<<mask.rows<<", mais image de largeur "<<dim2DHA.x<<endl;
    }

    mask.convertTo(mask, CV_8U);  ///Chargement du masque pour correction aberration/ampli

    return mask;
}

/// Return the size of polynomial. Used to calculate the number of columns for the vector "polynom_to_fit" (ex [1,x,x^2,xy,y,y^]->size=6], degre=3, size=10 etc.
///note that the general formula is just (deg+1)(deg+2)/2
int sizePoly2D(int deg){
    int j = 0, size = 0;
    for (int i = 0; i <= deg; i ++){
        while ((i+j) <= deg){
            size ++;
            j ++;
        }
        j = 0;
    }
    return size;
}
/*
int sizePoly2D(int deg){
    int NbCoef=(deg+1)*(deg+2)/2
    return NbCoef;
}*/

/// Count pixels in mask old version, (include margin) and exclude 1 out of 5 pixels to speed up calculations
///the return integer is used to initliae the size of the undersampled polynome.
int countM(Mat mask){
    const int step=5;//échantillonnage divisé par step=5
    int count = 0;
    for (int i = 0; i < mask.rows; i ++){
        for (int j = 0; j < mask.cols; j ++){
            if (mask.at<uchar>(j,i) > 45){//le masque doit être créé avec des valeurs  >45 (seuil arbritaire) hors objet, 0 dans l'objet
               count ++;
            }
        }
    }
    size_t NbPtRand=0;
    for (int i = 0; i < mask.rows; i =i+step){///point hors objet
        for (int j = 0; j < mask.cols; j=j+step){
            if (mask.at<uchar>(j,i) > 45){//Si les points sont hors objet
            mask.at<uchar>(j,i)=220;// alors mettre 220 dans un point sur 5.
            NbPtRand++;
            }
        }
    }
    return NbPtRand;
}
///new version, excluding border by "margin", because artefact can arise on the borders and distort the corrective polynom
int countM(Mat mask, int margin){
    const int step=5;//échantillonnage divisé par step=5
    int count = 0;
    for (int i = 0; i < mask.rows; i ++){
        for (int j = 0; j < mask.cols; j ++){
            if (mask.at<uchar>(j,i) > 45){//le masque doit être créé avec des valeurs  >45 (seuil arbritaire) hors objet, 0 dans l'objet
               count ++;
            }
        }
    }
    size_t NbPtRand=0;
    for (int i = margin; i < mask.rows-margin; i =i+step){///point hors objet et hors bords (même sans que les brds soit exclu par le masque binaire)
        for (int j = margin; j < mask.cols-margin; j=j+step){
            if (mask.at<uchar>(j,i) > 45){//Si les points sont hors objet
            mask.at<uchar>(j,i)=220;// alors mettre 220 dans un point sur 5 (chiffre pris au hasard, ne cherchez pas d'explications).
            NbPtRand++;
            }
        }
    }
    return NbPtRand;
}
///Generate the numerical value of the vector [1,x,x^2,xy,y,y^2, etc.] for each (x,y) outside the  mask area, undersampled to accelerate calculations  ("US"=undersampled)
void CalcPolyUs_xy(int degre_poly, Mat const & mask, Var2D dimImg, Mat &polynomeUs_to_fit)
{
   size_t nbPtUs=polynomeUs_to_fit.rows;
   if(nbPtUs>9){//you need at least 9 points for a deg 3 polynome
     int num_coef, Coord1D = 0;// coordinate (1D)
     for (int y = 0; y < dimImg.y; y ++){//largeur de 262
       for (int x = 0; x < dimImg.x; x ++){//hauteur de 262
         num_coef = 0;
         if (mask.at<uchar>(y,x)==220){//220=>1 point sur 5
           for (int powX = 0; powX <= degre_poly; powX ++){
             int powY = 0;
               while((powX+powY) <= degre_poly){

                 polynomeUs_to_fit.at<double>(Coord1D,num_coef) = pow((double)x,powX) * pow((double)y,powY); //genere polynome pour l'ajustement: identique pour toutes les images
                 num_coef ++;
                 powY ++;
               }
             }
             Coord1D ++;//x+y*dimImg.x
           }
        }
    }
  }
}
///Generate the numerical value of the vector [1,x,x^2,xy,y,y^2, etc] for each (x,y) outside the mask area, but fully sampled
///cette fonction évite de recalculer les monome : xy,  xx et yy etc. pour chaque image : les valeurs des monômes sont pré-calculées une seule fois
void CalcPoly_xy(int degre_poly,Var2D dimImg, Mat &polynome_to_fit)
{
     int num_coef, Coord1D = 0;// coordinate (1D)
     for (int y = 0; y < dimImg.y; y ++){//largeur de 262
       for (int x = 0; x < dimImg.x; x ++){//hauteur de 262
         num_coef = 0;
           for (int powX = 0; powX <= degre_poly; powX ++){///boucle inutilement compliquée : on pourrait itérer sur num_coef
             int powY = 0;
               while((powX+powY) <= degre_poly){
                 polynome_to_fit.at<double>(Coord1D,num_coef) = pow((double)x,powX) * pow((double)y,powY); //genere polynome pour l'ajustement: identique pour toutes les images
                 num_coef ++;
                 powY ++;
               }
             }
             Coord1D ++;
        }
    }
}
///calculate the numerical value of the poly for *all*  (x,y)
void compuBackgr2(Mat const &coefficients, Mat const & polynome_to_fit, Mat &PolyBackgr)
{   size_t Coord1D=0;
    int poly_size=polynome_to_fit.cols;//sizePoly2D(deg);
    for (int y = 0; y < PolyBackgr.rows; y ++)///scan all the (x,y) coord.
        for (int x = 0; x < PolyBackgr.cols; x ++){
            double sum=0;
            for(int num_coef=0;num_coef<poly_size;num_coef++){///calculate the numerical value of the poly for *one* coordinate (x_num_coef,y_num_coef)
            //for(int num_coef=0;num_coef<=poly_size;num_coef++){///calculate the numerical value of the poly for *one* coordinate (x_num_coef,y_num_coef)
               sum +=(coefficients.at<double>(num_coef)) * polynome_to_fit.at<double>(Coord1D,num_coef);
               //if(abs(sum)==0) cout<<"attention"<<endl;
            }
            PolyBackgr.at<double>(y,x) = sum;//poly2DEval2(coefficients, polynome_to_fit, UsCoord1D);///calculate the numerical value of the poly for one  (x,y)
            Coord1D++;
        }
}
///--------------------//Top functions used in "main.c" containing the algorithm for phase (abercorr2) and amplitude (amplicorr2)------------------------------------------------------------------------------
Mat  aberCorr2(Mat const &image, Mat const &mask,  Mat const &polynomeUs_to_fit, Mat const &polynome_to_fit)
{
    Mat coefsolve;
    compuCoefPoly2(image, mask, coefsolve, polynomeUs_to_fit, true); /// Compute the coef of polynomial (Least Squares method)
    Mat resultatpolyBG(image.rows, image.cols, CV_64F), result_final(image.rows, image.cols, CV_64F);
    compuBackgr2(coefsolve, polynome_to_fit, resultatpolyBG);/// Compute the background image with the coef of polynomial
//    SAV2(resultatpolyBG, "/home/mat/tomo_test/poly_aber.raw", t_float, "a+b");
    result_final = image-resultatpolyBG;
    return result_final;
}
///amplitude correction, with a naive division  regularisation
Mat  ampliCorr2(Mat const & image,  Mat const &polynomeUs_to_fit, Mat const &polynome_to_fit, Mat mask)
{
    Mat coefsolve;
    compuCoefPoly2(image, mask, coefsolve, polynomeUs_to_fit, true); /// Compute the coef of polynomial (Least Squares method)
    Mat resultatpoly(image.rows, image.cols, CV_64F), result(image.rows, image.cols, CV_64F);
    compuBackgr2(coefsolve, polynome_to_fit,  resultatpoly);/// Compute the background image with the coef of polynomial

   result = image/(resultatpoly+0.00000001);

    return result;
}
///amplitude correction, with a more robust division regularisation. This one should be prefered.
Mat  ampliCorr3(Mat const & image,  Mat const &polynomeUs_to_fit, Mat const &polynome_to_fit, Mat mask, double gamma_ampli)
{
    Mat coefsolve;
    compuCoefPoly2(image, mask, coefsolve, polynomeUs_to_fit, true); /// Compute the coef of polynomial (Least Squares method)
    Mat resultatpoly(image.rows, image.cols, CV_64F), result(image.rows, image.cols, CV_64F);
    compuBackgr2(coefsolve, polynome_to_fit,  resultatpoly);/// Compute the background image with the coef of polynomial
    double minVal, maxVal;
    cv::Point minLoc, maxLoc;
    cv::minMaxLoc(resultatpoly, &minVal, &maxVal, &minLoc, &maxLoc);//find min, max and their position
    double epsilon = gamma_ampli * maxVal;  // adjust (0.001 - 0.05)
    Mat safe_poly = max(resultatpoly, epsilon);
    result = image / safe_poly;
    return result;
}

/// Compute the coef of polynomial (Least Squares method) by SVD,  i.e. solve COEF*POLYNOME_TO_FIT=BACKGROUND (with undersampled variables to speed up the process)
void compuCoefPoly2(Mat const &imagebrut, Mat const & mask, Mat& coef_polynomial, Mat const &polynomeUs_to_fit, bool method)
{
  int nbCoef=polynomeUs_to_fit.cols, nbPtPoly=polynomeUs_to_fit.rows;///
  Mat Bt(Size(nbPtPoly,nbCoef), CV_64F);//variable to stock transposed polynome_to_fit, for SVD inversion
  Mat undersampled_background(nbPtPoly, 1, CV_64F);///matrix containing  1 point out of  "step" (default_step=5) (gives NbPtOK).
  int margin = 5; // ou 5–20 pixels

  if(nbPtPoly>9){
    int UsCoord = 0;/// 1D undersampled coordinate (corresponds to the 2D undersampled coordinates (x_us,y_us))
    for (int y = 0; y < imagebrut.rows; y ++){
      for (int x = 0; x < imagebrut.cols; x ++){
        if (
        mask.at<uchar>(y,x)==220)
        {//220=>1 point sur 5
           undersampled_background.at<double>(UsCoord) = (double)imagebrut.at<double>(y,x); /// copy 1 point out of 5 from image (outside the mask area) to speed up least square
           UsCoord ++;//val max=210/5*210/5=42*42=1764
         }
     }
  }
  Mat coef(nbCoef, 1, CV_64F), D(Size(nbCoef, nbPtPoly), CV_64F), invD(Size(nbCoef, nbPtPoly), CV_64F);
  if (method){  /// Use OpenCV solve() function to solve the linear system
    cv::solve(polynomeUs_to_fit, undersampled_background, coef, DECOMP_NORMAL);//DECOMP_NORMAL->speed ++, DECOMP_SVD = robust
 // SVD svd(polynomeUs_to_fit);
//cout << "Singular values: " << svd.w << endl;
  }
  else{///alternatively, use a simple matrix inversion (less robust)
    cv::transpose(polynomeUs_to_fit, Bt);
    D = Bt * polynomeUs_to_fit;
    cv::invert(D, invD); /// Use OpenCV matrix inversion operators to solve the linear system
    coef = (invD * Bt) *undersampled_background;
  }
  coef.copyTo(coef_polynomial);
  }
   else{
       Mat coef=Mat::zeros(nbCoef, 1, CV_64F);//init with zeros
  }
}


///*-----------------------------------------------old functions, with automatic segmentation-----------------------------------
/// function threshCallback
/*void threshCallback(int thresh, void* param)
{
    RNG rng(12345);
    Mat &src=*(Mat*)param;//typecast  du ptr void vers ptr Mat
    Mat canny_output;
    Mat src_filtering;
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    imshow("Source",src);
    /// Detect edges using canny
    Canny(src, canny_output, thresh, thresh*2, 3);
    /// Detect edges using Threshold
    threshold(src, canny_output, thresh, 255, THRESH_BINARY);
    /// Dilate helps to remove potential holes between edge segments
    dilate(canny_output, canny_output, Mat(), Point(-1,-1));
    //imshow("Canny output",canny_output);
    /// Find contours
    findContours(canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
    /// Find the convex hull object for each contour
    vector<vector<Point> >hull(contours.size());
    for (int i = 0; i < contours.size(); i ++)
    {
        convexHull(Mat(contours[i]), hull[i], false);
    }
    /// Draw contours
    Mat drawing = Mat::zeros(canny_output.size(), CV_8UC3);//Uc3?
    //Mat drawing_gray;
    for (int i = 0; i < contours.size(); i ++)
    {
        Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255));
        /// hull results
        drawContours(drawing, hull, i, Scalar(255,255,255), CV_FILLED, 8, vector<Vec4i>(), 0, Point());
    }
    /// Show in a window
    imshow("Image_mask", drawing);
    imwrite("Image_mask.tif",drawing);
}*/

