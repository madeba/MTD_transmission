#ifndef __CORR_ABER__
#define __CORR_ABER__

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
using namespace cv;

/// Function header
Mat init_mask_aber(std::string Chemin_mask, Var2D dim2DHA);
int sizePoly2D(int deg);
int countM(Mat mask);
void compuBackgr2(Mat const &coefficients, Mat const & polynome_to_fit, Mat &PolyBackgr);
Mat  aberCorr2(Mat const &image, Mat const &mask, Mat const &polynomeUs_to_fit,Mat const &polynome_to_fit);
Mat  ampliCorr2(Mat const & image,  Mat const &polynomeUs_to_fit, Mat const &polynome_to_fit, Mat mask);
void compuCoefPoly2(Mat const &imagebrut, Mat const & mask, Mat& coef_polynomial, Mat const &polynome_to_fit, bool method);
void CalcPolyUs_xy(int degre_poly, Mat const & mask, Var2D dimChpCplx,Mat &polynome_to_fit);
void CalcPoly_xy(int degre_poly, Var2D dimImg, Mat &polynome_to_fit);
void initCorrAber(std::string Chemin_mask, Mat const &mask, size_t degre_poly, Var2D dim2DHA,Mat &polynome_to_fit, Mat &polynomeUs_to_fit);

//void threshCallback(int, void*);
//double poly2DEval2(Mat const &coefficients, Mat const &polynome_to_fit, size_t UsCoord1D);
//void poly2DEval(Mat const &coefficients, int deg,int x, int y);
//void compuBackgr(Mat const &coefficients, int deg, Mat imageBackgr);
//void compuPoly(Mat const &imagebrut, Mat const &mask, Mat& polynomial, int deg, bool method, int NbPtOk);
//Mat aberCorr(Mat image, Mat mask,  int degpoly, int NbPtOk);
//Mat  ampliCorr(Mat const &image, Mat mask, int degpoly, int NbPtOk);
/*class ImgSeuil
{
private:
    Mat src_gray; // non-static
    int thresh;
public:
    ///methods
    ImgSeuil();
    void threshCallback(int, void* ); // non-static
};*/
#endif
