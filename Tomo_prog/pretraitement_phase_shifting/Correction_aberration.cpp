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
//Running time: 0,3s with an Intel quadcore i7 processor and 16 Gbytes of memory

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include <chrono>
using namespace cv;
using namespace std;


/// Function header
#include "Correction_aberration.h"

Mat init_mask_aber(string Chemin_mask, Var2D dim2DHA)
{
    ///Charger masque aberration
    Mat mask = imread(Chemin_mask, 0);
    if(! mask.data )
    {
        cout <<  "Masque non trouvé, création masque unité" << std::endl ;
        mask=255*Mat::ones(dim2DHA.x,dim2DHA.y, CV_8UC1);
    }
    else{
        cout<<"chargement du masque"<<Chemin_mask<<endl;
    }
    if(mask.rows!=2*dim2DHA.x)
    {
        cout<<"Problème aberrations : masque "<<Chemin_mask<< " de largeur "<<mask.rows<<", mais image de largeur "<<2*dim2DHA.x<<endl;
    }

    ///Chargement du masque pour correction aberration/ampli

    mask.convertTo(mask, CV_8U);

    return mask;
}

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



/// Return the size of polynomial. Used to calculate the number of columns for the vector "polynom_to_fit"
int sizePoly2D(int deg){
    int j = 0, size = 0;
    for (int i = 0; i <= deg; i ++){
        while ((i+j) <= deg){
            size ++;
            j ++;
        }
        j = 0;
    }
    //cout<<"size="<<size<<endl;
    return size;
}

/// Count pixels in mask (for control, not used in calculations)
int countM(Mat mask){
    const int step=5;//échantillonnage divisé par step=5
    //vector<int> tableau(5);
    int count = 0;
    for (int i = 0; i < mask.rows; i ++){
        for (int j = 0; j < mask.cols; j ++){
            if (mask.at<uchar>(j,i) > 45){//le masque doit être créé avec des valeurs  >45 (seuil arbritaire) hors objet, 0 dans l'objet
               count ++;
            }
        }
    }

    const int nbPt=mask.rows*mask.cols;
    double test_pt[nbPt];
    for(size_t cpt=0;cpt<nbPt;cpt++)
    test_pt[cpt]=0;

    const int nbPointOk=count;
    // cout<<"nb pt hors objet"<<count<<endl;
    int *tabPtOk=new int[nbPointOk*2];
    size_t cpt=0;

  /*  imshow("masque",mask);
    waitKey(0);*/

    size_t NbPtRand=0;
    for (int i = 0; i < mask.rows; i =i+step)///point hors objet
    {
        for (int j = 0; j < mask.cols; j=j+step)
        {
            if (mask.at<uchar>(j,i) > 45)//Si les points sont hors objet
            {

            //cpt++;
            mask.at<uchar>(j,i)=220;// alors mettre 220 dans un point sur 5.
           // cout<<(double)mask.at<uchar>(j,i)<<endl;
            NbPtRand++;
           // cpt=262*j+i;
           // test_pt[cpt]=mask.at<uchar>(j,i);
            }
        }
    }
    //cout<<"NbPtRand="<<NbPtRand<<endl;
   // SAV(test_pt,nbPt,m1.chemin_result+"/test_pt.bin",FLOAT,"wb");

   /* cout<<"nbptrand="<<NbPtRand<<endl;;
    imshow("masque2",mask);
    waitKey(0);*/
    //return count;
    return NbPtRand;
}
///Generate the numerical value of the vector [1,x,x^2,xy,y,y^2] for each (x,y) in the image
void CalcPoly_xy(int degre_poly,int NbPtOk, Mat const & mask, Var2D dimImg,Mat &polynome_to_fit)
{
    int nbRows=NbPtOk;
    //int nbRows=polynome_to_fit.row();
   // cout<<"nbrows="<<nbRows<<endl;
   if(nbRows>9){
     int col, k = 0;
     for (int y = 0; y < dimImg.y; y ++){//largeur de 262
       for (int x = 0; x < dimImg.x; x ++){//hauteur de 262
         col = 0;
         if (mask.at<uchar>(y,x)==220){//220=>1 point sur 5
           for (int i = 0; i <= degre_poly; i ++){
             int j = 0;
               while((i+j) <= degre_poly){
                 polynome_to_fit.at<double>(k,col) = pow((double)x,i) * pow((double)y,j); //genere polynome pour l'ajustement: identique pour toutes les images
                 col ++;
                 j ++;
               }
             }
             k ++;
           }
        }
    }
  }
}


/// Compute the coef of polynomial (Least Squares method) by SVD : i.e. solve coef*polynom_to_fit=background
//void compuPoly(Mat const &imagebrut, Mat mask, Mat& coef_polynomial, int deg, bool method, int NbPtOk)
void compuPoly(Mat const &imagebrut, Mat const &mask, Mat& coef_polynomial, int deg, bool method, int NbPtOk)
{ ///auto start = std::chrono::system_clock::now();
    int nbCols = sizePoly2D(deg);//Nb coef poly
    //cout<<"nbcols="<<nbCols<<endl;
    //int nbRows = countNonZero(mask); /// Opencv function
    //int nbRows = countM(mask); //nb points hors objet
    int nbRows=NbPtOk;
   // cout<<"nbrows="<<nbRows<<endl;
   if(nbRows>9)
   {
        int col, k = 0;
        // double FOND[imagebrut.rows*imagebrut.rows];
        // for(int cpt=0;cpt<imagebrut.rows*imagebrut.rows;cpt++)
        // FOND[cpt]=0;
        Mat polynome_to_fit(Size(nbCols,nbRows), CV_64F);///Polynome to fit= function to fit. We use a polynome. we had to genrate polynome_to_fit=[1,x,x^2,xy,y^2] for each coordinate (x,y)
        Mat Bt(Size(nbRows,nbCols), CV_64F);//variable to stock transposed polynome_to_fit, for SVD inversion
        Mat undersampled_background(nbRows, 1, CV_64F);///matrice contenant 1 point sur 5 (NbPtOK). Il s'agit de l'image avec moins de points, f est donc toujours connue

        for (int y = 0; y < imagebrut.rows; y ++){//largeur de 262 Remarquer que imagebrut ne sert à rien, sauf pour les coordonnées de balayage
            for (int x = 0; x < imagebrut.cols; x ++)//hauteur de 262
            {
               col = 0;
               //if (mask.at<float>(y,x) > 210 )
              // int test=mask.at<uchar>(y,x);
               if (mask.at<uchar>(y,x)==220)//220=>1 point sur 5
                {
                    undersampled_background.at<double>(k) = (double)imagebrut.at<double>(y,x); /// copy 1 point out of 5 from image (outside the mask area) to speed up least square
                    for (int i = 0; i <= deg; i ++){
                        int j = 0;
                        while((i+j) <= deg){
                            polynome_to_fit.at<double>(k,col) = pow((double)x,i) * pow((double)y,j); //genere polynome pour l'ajustement: identique pour toutes les images
                            col ++;
                            j ++;
                        }
                    }
                    k ++;
                   //cout<<"k="<<k<<endl;
                }
            }
        }

        double t_total=0;

        Mat coef(nbCols, 1, CV_64F);
        Mat D(Size(nbCols, nbRows), CV_64F);
        Mat invD(Size(nbCols, nbRows), CV_64F);

        /// Use OpenCV solve() function to solve the linear system
        if (method){
            cv::solve(polynome_to_fit, undersampled_background, coef, DECOMP_NORMAL);//DECOMP_NORMAL->speed ++
          //  cout<<"t_total="<<(double)t_total/CLOCKS_PER_SEC<<endl;
        }///alternatively, use a simple matrix inversion (less robust)
        else{
            /// Use OpenCV matrix inversion operators to solve the linear system
            cv::transpose(polynome_to_fit, Bt);
            D = Bt * polynome_to_fit;
            cv::invert(D, invD);
            coef = (invD * Bt) *undersampled_background;
        }
        coef.copyTo(coef_polynomial);
   }
   else{
       Mat coef=Mat::zeros(nbCols, 1, CV_64F);//init with zeros
   }
  /*  auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    std::cout <<"compu_poly = "<< elapsed.count()/(pow(10,9)) <<"s"<< '\n';*/
}
/// numerically compute the value of each point of the polynom for **one** coordinate (x,y), for deg 3: a+b*x+c*x²+d*x³+e*y+f*xy+g*x²y+h*y²+i*xy²+j*y³
double poly2DEval(Mat const &coefficients, int deg, int x, int y)
{
   int k = 0;
   double sum = 0;
   for (int i = 0; i <= deg; i ++){
        int j = 0;
        while ((i+j) <= deg){
            sum +=(coefficients.at<double>(k)) * pow(x,i) * pow(y,j);/// * pow(x,i) * pow(y,j) n'a pas besoin d'être recalculé à chaque fois :c'est polynome_to_fit, déjà calculé CalcPoly_xy
            //cout<<"k="<<k<<" i= "<<i <<" j= "<<j<<endl;
            k ++;
            j ++;
        }
    }

  /*  int xc=x*x,yc=y*y;
    sum=coefficients.at<double>(0)+
    coefficients.at<double>(1)*y+
    coefficients.at<double>(2)*y*y+
    coefficients.at<double>(3)*y*y*y+
    coefficients.at<double>(4)*x+
    coefficients.at<double>(5)*x*y+
    coefficients.at<double>(6)*x*y*y+
    coefficients.at<double>(7)*x*x+
    coefficients.at<double>(8)*x*x*y+
    coefficients.at<double>(9)*x*x*x;*/
    return sum;
}
/// Compute the background polynome for *all* coordinates (x,y)
void compuBackgr(Mat const &coefficients, int deg, Mat PolyBackgr)
{
    for (int y = 0; y < PolyBackgr.rows; y ++)
        for (int x = 0; x < PolyBackgr.cols; x ++){
            PolyBackgr.at<double>(y,x) = poly2DEval(coefficients, deg, x, y);
        }
}

/// Apply the aberration correction (old version)
/*Mat  aberCorr(Mat image, Mat mask, int degpoly,  int NbPtOk){
    Mat coefsolve;
    compuPoly(image, mask, coefsolve, degpoly, true, NbPtOk); /// Compute the coef of polynomial (Least Squares method)
    Mat resultatpoly(image.rows, image.cols, CV_64F);
    Mat result(image.rows, image.cols, CV_64F);
    compuBackgr(coefsolve, degpoly, resultatpoly);/// Compute the background image with the coef of polynomial
    //SAV2((double*)resultatpoly.data,image.rows*image.cols,"/home/mat/tmp/poly_aber_phase.raw",t_float,"a+b");
    result = image-resultatpoly;
    return result;
}*/


Mat  ampliCorr(Mat const & image, Mat mask, int degpoly, int NbPtOk)
{
    Mat coefsolve;
    compuPoly(image, mask, coefsolve, degpoly, true, NbPtOk); /// Compute the coef of polynomial (Least Squares method)
    Mat resultatpoly(image.rows, image.cols, CV_64F);
    Mat result(image.rows, image.cols, CV_64F);
    compuBackgr(coefsolve, degpoly, resultatpoly);/// Compute the background image with the coef of polynomial
    //SAV2((double*)resultatpoly.data,image.rows*image.cols,"/home/mat/tmp/poly_aber_phase.raw",t_float,"a+b");
    result = image/resultatpoly;
    return result;
}

///#--------------------nouvelle fonction pour calculer poly_to_fit en dehors de la boucle------------------------------------------------------------------------------

///surcharge pour calcul poly exterieur
Mat  aberCorr2(Mat image, Mat const &mask,  Mat const &polynome_to_fit, int degpoly,  int NbPtOk)
{
    Mat coefsolve;
    /// Compute the coef of polynomial (Least Squares method)
    compuPoly2(image, mask, coefsolve, polynome_to_fit, degpoly, true, NbPtOk);
    Mat resultatpoly(image.rows, image.cols, CV_64F);
    Mat result_final(image.rows, image.cols, CV_64F);
    compuBackgr(coefsolve, degpoly, resultatpoly);/// Compute the background image with the coef of polynomial
    //SAV2((double*)resultatpoly.data,image.rows*image.cols,"/home/mat/tmp/poly_aber_phase.raw",t_float,"a+b");
    result_final = image-resultatpoly;
    return result_final;
}


/// Compute the coef of polynomial (Least Squares method) by SVD : i.e. solve coef*polynom_to_fit=background
//void compuPoly(Mat const &imagebrut, Mat mask, Mat& coef_polynomial, int deg, bool method, int NbPtOk)
//surcharge avec calciul du plynome à l extérieur
void compuPoly2(Mat const &imagebrut, Mat const & mask, Mat& coef_polynomial, Mat const &polynome_to_fit,int deg, bool method, int NbPtOk)
{
  double t_total=0;
  int nbCols = sizePoly2D(deg);
  int nbRows=NbPtOk;
  Mat Bt(Size(nbRows,nbCols), CV_64F);//variable to stock transposed polynome_to_fit, for SVD inversion
  Mat undersampled_background(nbRows, 1, CV_64F);///matrice contenant 1 point sur 5 (NbPtOK). Il s'agit de l'image avec moins de points

  if(nbRows>9){
    int col, k = 0;
    for (int y = 0; y < imagebrut.rows; y ++){//largeur de 262
      for (int x = 0; x < imagebrut.cols; x ++){//hauteur de 262
         col = 0;
         //if (mask.at<float>(y,x) > 210 )
         // int test=mask.at<uchar>(y,x);
         if (mask.at<uchar>(y,x)==220){//220=>1 point sur 5
           undersampled_background.at<double>(k) = (double)imagebrut.at<double>(y,x); /// copy 1 point out of 5 from image (outside the mask area) to speed up least square
           k ++;
         }
     }
  }
  Mat coef(nbCols, 1, CV_64F), D(Size(nbCols, nbRows), CV_64F), invD(Size(nbCols, nbRows), CV_64F);
  /// Use OpenCV solve() function to solve the linear system
  if (method){
    cv::solve(polynome_to_fit, undersampled_background, coef, DECOMP_NORMAL);//DECOMP_NORMAL->speed ++
    //  cout<<"t_total="<<(double)t_total/CLOCKS_PER_SEC<<endl;
  }///alternatively, use a simple matrix inversion (less robust)
  else{
    /// Use OpenCV matrix inversion operators to solve the linear system
    cv::transpose(polynome_to_fit, Bt);
    D = Bt * polynome_to_fit;
    cv::invert(D, invD);
    coef = (invD * Bt) *undersampled_background;
  }
  coef.copyTo(coef_polynomial);
  }
   else{
       Mat coef=Mat::zeros(nbCols, 1, CV_64F);//init with zeros
  }
}
