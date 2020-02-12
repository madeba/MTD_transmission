//Manuscript Title: Tomographic diffracitive micrscopy and multiview profilometry
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
#include <cxcore.h>
#include <highgui.h>
#include <cv.h>
#include "fonctions.h"
using namespace cv;
using namespace std;


/// Function header
#include "Correction_aberration.h"



/// function threshCallback
void threshCallback(int thresh, void* param)
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
    imshow("Canny output",canny_output);
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
}

/// Compute the value sum of the polynomial for each (x,y), for deg 3: 1+x+x²+x³+y+xy+x²y+y²+xy²+y³
double poly2DEval(Mat coefficients, int deg, int x, int y)
{
    int k = 0;
    double sum = 0;
    for (int i = 0; i <= deg; i ++)
    {
        int j = 0;
        while ((i+j) <= deg)
        {
            sum += (coefficients.at<double>(k)) * pow(x,i) * pow(y,j);
            k ++;
            j ++;
        }
    }
    return sum;
}

/// Return the size of polynomial
int sizePoly2D(int deg)
{
    int j = 0, size = 0;
    for (int i = 0; i <= deg; i ++)
    {
        while ((i+j) <= deg)
        {
            size ++;
            j ++;
        }
        j = 0;
    }
    return size;
}

/// Count pixels in mask
int countM(Mat mask)
{
    int count = 0;
    for (int i = 0; i < mask.rows; i ++)
    {
        for (int j = 0; j < mask.cols; j ++)
        {
            if (mask.at<uchar>(j,i) > 10)
                count ++;
        }
    }
   // cout<<"pixels dans le masque"<<count<<endl;
    return count;
}

/// Compute the coef of polynomial (Least Squares method)
void compuPoly(Mat imagebrut, Mat mask, Mat& polynomial, int deg, bool method)
{

    int nbCols = sizePoly2D(deg);
    //int nbRows = countNonZero(mask); /// Opencv function
    int nbRows = countM(mask);
    int col, k = 0;
   // double FOND[imagebrut.rows*imagebrut.rows];
   // for(int cpt=0;cpt<imagebrut.rows*imagebrut.rows;cpt++)
       // FOND[cpt]=0;
    Mat B(Size(nbCols,nbRows), CV_64F);
    Mat Bt(Size(nbRows,nbCols), CV_64F);
    Mat f(nbRows, 1, CV_64F);

    for (int y = 0; y < imagebrut.rows; y ++)
    {
        for (int x = 0; x < imagebrut.cols; x ++)
        {
            col = 0;
            if (mask.at<uchar>(y,x) > 10)
            {
                //FOND[x+y*imagebrut.rows]=(double)imagebrut.at<double>(y,x);

                f.at<double>(k) = (double)imagebrut.at<double>(y,x); /// uchar--> double for .bin files
                for (int i = 0; i <= deg; i ++)
                {
                    int j = 0;
                    while((i+j) <= deg)
                    {
                        B.at<double>(k,col) = pow((double)x,i) * pow((double)y,j);
                        col ++;
                        j ++;
                    }
                }
                k ++;
            }
        }
    }
//SAV(FOND,imagebrut.rows*imagebrut.rows,"/home/mat/tomo_test/fond.raw",FLOAT,"a+b");
    Mat coef(nbCols, 1, CV_64F);
    Mat D(Size(nbCols, nbRows), CV_64F);
    Mat invD(Size(nbCols, nbRows), CV_64F);
    /// Use OpenCV solve() function to solve the linear system
    if (method)  cv::solve(B, f, coef, DECOMP_NORMAL);
    else
    {
        /// Use OpenCV matrix inversion operators to solve the linear system
        cv::transpose(B, Bt);
        D = Bt * B;
        cv::invert(D, invD);
        coef = (invD * Bt) *f;
    }
    coef.copyTo(polynomial);
}

/// Compute the background image with the coef of polynomial
void compuBackgr(Mat coefficients, int deg, Mat imageBackgr)
{
    for (int y = 0; y < imageBackgr.rows; y ++)
    {
        for (int x = 0; x < imageBackgr.cols; x ++)
        {
            imageBackgr.at<double>(y,x) = poly2DEval(coefficients, deg, x, y);
        }
    }
  //  imwrite("home/mat/tomo_test/Image_background.tif", imageBackgr);
}

/// Apply the aberration correction
Mat  aberCorr(Mat image, Mat mask, double *polyAber, int degpoly)
{
    degpoly = 4;
    Mat coefsolve;
    /// Compute the coef of polynomial (Least Squares method)
    compuPoly(image, mask, coefsolve, degpoly, true);
    Mat resultatpoly(image.rows, image.cols, CV_64F);
    Mat result(image.rows, image.cols, CV_64F);
    compuBackgr(coefsolve, degpoly, resultatpoly);/// Compute the background image with the coef of polynomial
   // double minmat, maxmat;
   // polyAber = (double*)resultatpoly.data;
   for(size_t x=0; x<image.rows;x++){
    for(size_t y=0; y<image.cols;y++)
   {
        size_t cpt=x+y*image.cols;
        polyAber[cpt]=resultatpoly.at<double>(y,x);
   }
   }

  //  SAV(polyAber,image.rows*image.cols,"/home/mat/tomo_test/poly_aber_phase.raw",FLOAT,"a+b");
    result = image-resultatpoly;
    return result;
    }


Mat  ampliCorr(Mat image, Mat mask, double *polyAber,int degpoly)
{
    degpoly = 4;
    Mat coefsolve;
    /// Compute the coef of polynomial (Least Squares method)
    compuPoly(image, mask, coefsolve, degpoly, true);
    Mat resultatpoly(image.rows, image.cols, CV_64F);
    Mat result(image.rows, image.cols, CV_64F);
    compuBackgr(coefsolve, degpoly, resultatpoly);/// Compute the background image with the coef of polynomial

    //double minmat, maxmat;

     polyAber = (double*)resultatpoly.data;

     //SAV(polyAber,image.rows*image.cols,"/home/mat/tomo_test/poly_aber_ampli.raw",FLOAT,"a+b");
    result = image/resultatpoly;
    return result;
}

