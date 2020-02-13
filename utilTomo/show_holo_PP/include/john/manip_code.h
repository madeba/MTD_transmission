#ifndef __MANIPS__
#define __MANIPS__



#include "manipTomo.h"



// manip standard à décalage de phase, 4 clichés par point de fleur

void 
manip_std(Instruments *Instru, vPlotterT<float> *FlowerPlotter, Parametres *Param);


// hors-axe

void
manip_horsaxe(Instruments *Instru, vPlotterT<float> *FlowerPlotter, Parametres *Param);
void
manip_horsaxe_fast(Instruments *Instru, vPlotterT<float> *FlowerPlotter, Parametres *Param);


// décalage de phase + HDR ( 8 clichés par point de fleur, deux par valeur de phase: 1 à exposition spécifiée, suivi immédiatement d'un second à double exposition );

void
manip_hdr(Instruments *Instru, vPlotterT<float> *FlowerPlotter, Parametres *Param);


// manip standard à décalage de phase ET SHUTTER, 8 clichés par point de fleur
// 4 shutter ouvert, puis 4 shutter fermé

void 
manip_std_shutter(Instruments *Instru, vPlotterT<float> *FlowerPlotter, Parametres *Param);


// test du shutter

void 
manip_shutter(Instruments *Instru, vPlotterT<float> *FlowerPlotter, Parametres *Param);


#endif // __MANIPS__
