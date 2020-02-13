#ifndef __MANIP_TOMO__
#define __MANIP_TOMO__


#include "Camera.h"
#include "Ljack.h"
#include "Ljack_DAC.h"
#include "vPlotterT.h"


/* -------------------------------------------------------------------------- */
// types internes
/* -------------------------------------------------------------------------- */


// liste des types d'acquisition possibles
enum ExpType {
  _decal = 0,                     // décalage de phase classique
  _decal_hdr = 1,                 // décalage + HDR
  _hors_axe = 2,                  // hors-axe sans décalage
  _decal_shutter = 3              // décalage de phase avec shutter
};


// liste des types de balayage possibles
enum ScanType {
  _flower = 0,
  _circle = 1,
  _angular = 2,
  _segment = 3
};


typedef struct _Instruments{
  Camera* GEV_Camera;
  Ljack_DAC* DAC_tiptilt;
  Ljack_DAC* DAC_phase;
} Instruments;


typedef struct _Parametres{
  size_t nb_iterations;
  char* images_dir;
  float exposure_ms;
  float modphase_increment;
  float tiptilt_factor;
  // 
  float offset_vx;
  float offset_vy;
  // temporisations
  size_t delay_tilt_ms;
  size_t delay_phase_ms;
  size_t delay_shutter_ms;
  size_t delay_camera_ms;
  // flags
  bool show_acquis_p;
  bool show_tiptilt_p;
} Parametres;


/* -------------------------------------------------------------------------- */
// Default values
/* -------------------------------------------------------------------------- */


//#define _EXPO 10.0f
#define _PHASE_INC 1.16f



#endif // __MANIP_TOMO__
