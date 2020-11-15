#ifndef __SCAN_FUNCTIONS__
#define __SCAN8FUNCTIONS__

#include <vector>


typedef struct {
float x,y;
}float2D;

typedef struct {
      int x,y;
}Var2D;

typedef struct {
  double Re,Im;
}nbCplx;

float2D maj_fleur(float rho, int nbHolo, float *theta);
void scan_fleur(float rho, int nbHolo, std::vector<float2D> &Vout_table);//tension x,y/rayon/nbHolo,/angle
//concentric circles, export pattern as a table
void scan_annular(float Vmax, float nb_circle, int nbHolo, std::vector<float2D> &Vout_table);
void scan_fermat(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dimHolo);
///generate a random cartesian coordinate between [-vmax, +vmax], return the values in Vout_table
void scan_random_cartesian(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dim);
void scan_random_polar3D(double Vmax, int nbHolo, std::vector<float2D> &Vout_table,Var2D dim);
float2D maj_spiral(float rho, int nbHolo, float *theta);//tension x,y/rayon/nbHolo,/angle
void scan_spiralOS(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, unsigned short int nbTurn);//tension x,y/rayon/nbHolo,/angle
void scan_spiral(double Vmax, int nbHolo, std::vector<float2D> &Vout_table, Var2D dim, unsigned short int nbTurn);//tension x,y/rayon/nbHolo,/angle
void scan_uniform3D(double Vmax, int nbHolo, std::vector<float2D> &Vout_table);
void scan_etoile(double Vmax, int nbHolo, int nbAxe,std::vector<float2D> &Vout_table);
#endif
