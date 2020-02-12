#ifndef __FONCTIONS__
#define __FONCTIONS__

/* --------------------------------------------------------------------------- */
// incl
/* --------------------------------------------------------------------------- */
#include "struct.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <strstream>
//#include <Magick++.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <vector>

//using namespace Magick;
using namespace std;

#include "projet.h"
typedef struct{
clock_t init, fin;
float   total;
}temps;

/* --------------------------------------------------------------------------- */
// Prototypes
/* --------------------------------------------------------------------------- */

float extract_val(string token,  string chemin_fic);
float ecrire_val(string token, float valeur_token, string chemin_fic);
string extract_string(string token,  string chemin_fic);

#endif

