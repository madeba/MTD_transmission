#ifndef DEF_STRUCT// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#include "struct.h"
#define DEF_STRUCT
#endif // DEF_STRUCT


//pixel information
struct PIXEL
{
	//int x;					//x coordinate of the pixel
    //int y;					//y coordinate
    int increment;			//No. of 2*pi to add to the pixel to unwrap it
    int number_of_pixels_in_group;	//No. of pixels in the pixel group
    float value;			//value of the pixel
	float reliability;
    int group;				//group No.
    int new_group;
    struct PIXEL *head;		//pointer to the first pixel in the group in the linked list
    struct PIXEL *last;		//pointer to the last pixel in the group
    struct PIXEL *next;		//pointer to the next pixel in the group
};

//the EDGE is the line that connects two pixels.
//if we have S PIXELs, then we have S horizental edges and S vertical edges
struct EDGE
{
	float reliab;			//reliabilty of the edge and it depends on the two pixels
	PIXEL *pointer_1;		//pointer to the first pixel
    PIXEL *pointer_2;		//pointer to the second pixel
    int increment;			//No. of 2*pi to add to one of the pixels to unwrap it with respect to the second
};





void phaseUnwrapping(nbCplx* obj, Var2D taille, double* WrappedImage, double* UnwrappedImage, PIXEL *pixel, EDGE *edge);
void phaseUnwrapping_Mat(Var2D taille, double* WrappedImage, double* UnwrappedImage);
void  initialisePIXELs(double *WrappedImage, PIXEL *pixel, int image_width, int image_height);
void  quick_sort(EDGE *Pointer, int size);
void    calculate_reliability(double *wrappedImage, PIXEL *pixel, int image_width, int image_height);
void  horizentalEDGEs(PIXEL *pixel, EDGE *edge, int image_width, int image_height);
void  verticalEDGEs(PIXEL *pixel, EDGE *edge, int image_width, int image_height);
void  gatherPIXELs(EDGE *edge, int image_width, int image_height);
void  unwrapImage(PIXEL *pixel, int image_width, int image_height);
void  returnImage(PIXEL *pixel, double *unwrappedImage, int image_width, int image_height);
void write_data(char *outputfile,double *Data,int length);



void phase2pi(nbCplx* obj, Var2D taille,double* WrappedImage);
