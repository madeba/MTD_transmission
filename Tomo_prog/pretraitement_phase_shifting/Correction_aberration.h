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
//using namespace std;


/// Function header
Mat init_mask_aber(std::string Chemin_mask, Var2D dim2DHA);
//void threshCallback(int, void*);
double poly2DEval(Mat coefficients, int deg,int x, int y);
int sizePoly2D(int deg);
int countM(Mat mask);
void compuBackgr(Mat coefficients, int deg, Mat imageBackgr);
void compuPoly(Mat imagebrut, Mat mask, Mat& polynomial, int deg, bool method, int NbPtOk);
Mat aberCorr(Mat image, Mat mask, double *polyAber, int degpoly, int NbPtOk);
Mat  ampliCorr(Mat image, Mat mask, double *polyAber, int degpoly, int NbPtOk);
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

