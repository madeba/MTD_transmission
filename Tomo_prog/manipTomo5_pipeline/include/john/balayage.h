#ifndef __BALAYAGE__
#define __BALAYAGE__


/* -------------------------------------------------------------------------- */
// Fonctions de parcours du tip-tilt
/* -------------------------------------------------------------------------- */



// --------------------------------------------------
// balayages en cercle


float
circle_x(float t);

float
circle_y(float t);


// --------------------------------------------------
// balayage en fleur


float
flower_x(float t);


float
flower_y(float t);



// --------------------------------------------------
// balayage angulaire


void
fill_angular_scan(float* data_x, float* data_y, size_t nb_theta, size_t nb_phi, float N1, float NA);


#endif //__BALAYAGE__
