#ifndef _VECTRA_PLOTTER_
#define _VECTRA_PLOTTER_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.vectraproject.com) */


// =============================================================================

// classe qui gère l'échantillonnage d'une fonction paramétrique
// n'affiche rien (cf vPlotter pour cela).

// =============================================================================




#include <iostream>
#include <string>
#include <math.h>

#include "macros.h"


using namespace std;



// =============================================================================


template <class T>
class vPlotterT
{

 private:



  // --------------------------------------------------
  // Variables d'état

  // valeur courante de t
  T t;
  // indice de la valeur courante (de 0 à nb_steps)
  size_t i;
  // valeurs de (x,y) pour t
  T x, y;
  // x and y outputs are multiplied by these factors ((T)1 by default)
  T factor_x, factor_y;


  // --------------------------------------------------
  // Variables de configuration

  // t ira de t_min à t_max en effectant nb_steps pas réguliers ou adaptatifs
  T t_min, t_max;

  bool regular_p;
  bool adaptative_p;
  bool start_zero_p; // pour ajouter un point zéro au début

  bool functions_set, range_set, ready;

  // dans le cas d'une paramétrisation régulière, nombre d'étapes
  size_t nb_steps;

  T (*fun_x)(T);
  T (*fun_y)(T);
  // dans le cas d'une paramétrisation adaptative, code de la fonction qui calcule T+1 en fonction de T
  T (*fun_updateT)(T);


  // mode spécial: tableau de données.
  // on alloue localement un tableau de taille nb_steps pour les valeurs de x
  // et un autre pour les valeurs de y.
  // compute-next -> parcours du tableau.
  // inutile alors de donner des fonctions, range, sampling, etc.

  T* data_x;
  T* data_y;
  bool data_allocated;

  bool repeat_mode_p;
  size_t repeat_steps, repeat_counter;

 public:


  // ==================================================
  // constructors

  vPlotterT(T (*fun_x)(T), T (*fun_y)(T));
  vPlotterT(); // will need to set functions

  ~vPlotterT() { };


  // ==================================================
  // configuration for normal mode
  //
  // all are mandatory, order at discretion

  // range of the parametric function
  void set_range(T t_min, T t_max);

  // only one choice
  void set_sampling_regular(size_t steps);
  void set_sampling_adaptative(T (*fun_next_t)(T));

  // only if not specified in constructor
  // both functions are fed by the same parameter
  void setFunctions(T (*fun_x)(T), T (*fun_y)(T));

  // ajoute (0,0) en début
  void add_zero()
  { start_zero_p = true; }

  void set_repeat(size_t repeats) { repeat_mode_p = true; repeat_steps = repeats; }

  // ==================================================
  // configuration for data mode

  void data_allocate(size_t nb_steps)
  {
    this -> nb_steps = nb_steps;
    ARRAY_ALLOC(data_x, nb_steps, T);
    ARRAY_ALLOC(data_y, nb_steps, T);

    data_allocated = ready = true;
  }


  // ==================================================
  // compute

  void compute_next();


  // ==================================================
  // accessors

  // retrieve computation result
  READER(T, x);
  READER(T, y);

  // factors that separately affect output x and y values
  READWRITER(T, factor_x);
  READWRITER(T, factor_y);

  // access raw data if relevant
  READER(T*, data_x);
  READER(T*, data_y);



};




// =============================================================================

#include "vPlotterT.cc"


#endif // _VECTRA_PLOTTER_


  // map a 1-arg function to matrix elements
/*   void map(T (*map_fun)(T)); */
