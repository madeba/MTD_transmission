#ifndef __cvPlotter__
#define __cvPlotter__


/* ============================================================================= 

   Cette classe gère le plotting dans un repère 2D par le biais d'images OpenCV.
   Sont supportés: le plotting de points

   ============================================================================= */


#include "cv.h"
#include "highgui.h"

#include "macros.h"
#include <iostream>

using namespace std;


class cvPlotter
{
  // ==================================================
 private:

  // ----------------------------------------
  float x_min, x_max, x_range;
  float y_min, y_max, y_range;

  float grid_xt, grid_yt;

  size_t dim_x, dim_y;
  CvScalar bg_color, plot_color, axis_color, grid_color;
  size_t plot_thick, axis_thick;

  IplImage* target_image;

  
  char* resident_window;
  bool resident_window_set;


  // ----------------------------------------
  // helper fonctions
  bool in_bounds(float x, float y);

  // retourne les coordonnées image d'une abcisse/ordonnée
  float project_x(float x);
  float project_y(float y);

  // dessine une ligne dans les coordonnées système
  void draw_line(float x1, float y1, float x2, float y2, CvScalar color, size_t thick); 
  void draw_dashed_line_H(float x1, float y1, float x2, float y2, CvScalar color, size_t thick, size_t gap); 
  void draw_dashed_line_V(float x1, float y1, float x2, float y2, CvScalar color, size_t thick, size_t gap); 
  

  float module(float val, float modulo);

  // ==================================================
 public:
  

  cvPlotter();
  ~cvPlotter();

  // ----------------------------------------
  // data accessors

  READWRITER(size_t, plot_thick);
  READWRITER(size_t, axis_thick);
  // grid for every xt, yt
  READWRITER(float, grid_xt);
  READWRITER(float, grid_yt);
  void set_plot_color(size_t r, size_t g, size_t b);
  void set_axis_color(size_t r, size_t g, size_t b);
  void set_grid_color(size_t r, size_t g, size_t b);


  // ----------------------------------------
  // interface functions

  void plot( float (*fun_y)(float), bool connected = true );
  void plot( float x, float y );

  void set_range(float x_min, float x_max, float y_in, float y_max);
  void resize(size_t dim_x, size_t dim_y);

  void reset();

  // offre le pointeur vers l'image interne pour permettre un affichage
  IplImage* yield_image();

  // affiche l'image interne dans une nouvelle fenête nommée à usage unique
  void display(const char* window_name);

  void resident_display(const char* window_name);
  void display_update();

};


#endif // __cvPlotter__
