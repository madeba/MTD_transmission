#ifndef __IO_FONCTIONS__
#define __IO_FONCTIONS__


#include <vector>
#include <complex>
#include <tiffio.h>
#include "projet.h"


void init_tab_val(string chemin_fic, vector<string> &tab_val);
void modif_tab_val(string token,string string_valeur_token,vector<string> &tab_val);
void sav_val(string chemin_fic,vector<string> &tab_val);
#endif
