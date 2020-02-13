#ifndef __FFTW_IMAGE2__
#define __FFTW_IMAGE2__


// pas encore utilisé, mais si un programme inclut FFTW_Image.hpp avant vImg.hpp, compile pas
#include "vImg.hpp"  
// ceci résoud définitivement le problème!

#include "vChronos.h" // avant iostream
#include "macros.h"

//#include <AIR_Volume.h> 
//#include <COMPLEX_Volume.h> 

#include <fftw3.h>
#include <complex.h>
#include <math.h>
#include <iostream>



#include <omp.h>
#include <xmmintrin.h>
//madvise
#include <sys/mman.h>


using namespace std;


#ifndef OMP_THREADS
#define OMP_THREADS 4
#endif


/*
  cette classe permet de représenter une image 2D traitable par fftw.  
*/


// =============================================================================

// U est le type fftw_complex (T = double) ou fftwf_complex (T = float)

template <typename T, typename FFTW_T>
class FFTW_Image{

  // ==================================================
protected:

  size_t transform_calls;
  vChronos *timer1, *timer2, *timer_fourier;


  // ----------------------------------------
  // core 
  
  size_t x_dim;
  size_t y_dim;
  size_t array_size; 

  FFTW_T* data;
  T* data_r2c_r; // pour r2c, alloué en même temps pour simplifier
  bool allocated_p;  

  T *shift_array_r, *shift_array_i; 
  bool shift_allocated_p;
  
  void init();


  // ----------------------------------------
  // fftw core related

  // fonction de calcul générique pour calculer TF avant et arrière
  inline void fourier_perform(bool is_forward); 

  size_t nb_threads; // nombre de threads utilisés pour la TF
  // déterminé à la construction, une fois pour toutes
  inline void cleanup_threads(); 

  fftwf_plan plan_forward_f;
  fftwf_plan plan_backward_f;
  fftw_plan plan_forward_d;
  fftw_plan plan_backward_d;
  fftw_plan plan_r2c_d;
  fftwf_plan plan_r2c_f;

  // avec ou sans wisdom, un plan doit être calculé la première fois
  bool plan_fwd_set_p, plan_bwd_set_p; 
  bool plan_r2c_fwd_set_p; 
  const char* wisdom_dir; // où on doit chercher les wisdoms (. par défaut)


  // ==================================================
public:

  
  // ----------------------------------------
  // constructeurs 
   
  inline FFTW_Image(); // interdit
  inline FFTW_Image(size_t x_dim, size_t y_dim, size_t nb_threads = 3);
  ~FFTW_Image();


  // ----------------------------------------
  // memory 
  inline void allocate(); 
  inline void unallocate(); 

  bool same_dim_p(const FFTW_Image<T, FFTW_T>& Other) const;
    
  // ----------------------------------------
  // accessors

  size_t get_size() const { return array_size; }

  READER(size_t, x_dim);
  READER(size_t, y_dim);
  READER(bool, allocated_p);
  READWRITER(size_t, nb_threads);
  READWRITER(bool, threads_initialized_p);

  void report() const { timer_fourier -> stop(); cout << "for calls# :" << transform_calls; cout.flush();}


  // ----------------------------------------
  // fourier shit
  void set_fourier_forward() {fourier_perform(true); }
  void set_fourier_backward() {fourier_perform(false); }

  // cas particulier où on effectue spécifiquement une transformée avant real to complex
  // nécessaire d'utiliser import_R_r2c_from
  inline void set_fourier_real_forward(T* opt_source = NULL);


  // division des valeurs par la taille de l'image
  inline virtual void set_fourier_normalize();

  // ces deux fonctions sacrifient R pour y écrire:
  inline void compute_module_R();  // le module
  inline void compute_log_module_R();  // log(1 + module) => échelle log10
  inline void quantify_R(unsigned char max_val, float refactor = 1.0f);


  // par défaut, l'import et l'export gèrent le circshift centré.
  // désactivable par paramètre optionnel à false
  inline void import_from(T* data_R, T* data_I, bool shift = true); 
  inline virtual void export_to(T* data_R, T* data_I, bool shift = true ); 

  // variantes de import_from pour seulement copier R ou I 
  // dans une structure complexe et pour une transformation complexe vers complexe
  // (utile en hors-axe, évite un circshift)
  // NB: la partie omise est mise à zero!
  template <class U> inline void import_R_from(U* data_R, bool shift = true); // si transtypage requis
  inline virtual void import_R_from(T* data_R, bool shift = true);
  inline void import_I_from(T* data_I, bool shift = true);
 
  // l'export se fait dans tous les cas à partir d'une structure complexe
  template <class U> inline void export_R_to(U* data_R, bool shift = true ); // si transtypage requis
  inline void export_R_to(T* data_R, bool shift = true ); 
  inline void export_I_to(T* data_I, bool shift = true ) { MSG_ASSERT(false, "pas fait");}

  // import avec circshift non-centré (circshift2) 
  // complexe vers complexe
  // pour débog. essentiellement (hologrames recalés)
  inline void import_excentric(T* data_R, T* data_I, size_t peak_x, size_t peak_y);
  

  // dans le cas d'une transformation avant réel vers complexe, l'import est spécifique
  // il n'y a plus que le shift à gérer. Et si on n'en veut pas, ben c'est pas encore géré
  template <class U> inline void import_R_r2c_from(U* data_R);
  inline void import_R_r2c_from(T* data_R);
  


protected:
  char* wisdom_filename;

  inline void compute_wisdom_filename(char* result, bool is_forward);
  inline void compute_wisdom_r2c_fwd_filename(char* result);
  inline void write_radix(char *radix); 

  inline bool import_wisdom(const char* wisdom_filename);
  inline void account_wisdom(FILE* wisdom_file);


public:
  void set_wisdom_dir(const char* dir) { wisdom_dir = dir; }

  inline void create_wisdom(bool is_forward);
  inline void create_wisdom_r2c_fwd();
protected:
  inline void makeWisdom(bool is_forward);
  inline void make_r2c_wisdom();
};



#include "./FFTW_Image.cc"

// =============================================================================
// =============================================================================


#endif //  __FFTW_IMAGE2__




// =============================================================================
// example
// =============================================================================


/*
  int
  main(void)
  {

  FFTW_Image<double, fftw_complex> ImgD(512, 512);
  ImgD.allocate();

  FFTW_Image<float, fftwf_complex> ImgF(512, 512);
  ImgF.allocate();

  //FFTW_Image<double, fftwf_complex> ImgErr(512, 512);
  //ImgErr.allocate();

  ImgD.set_fourier_forward();
  }
*/



// =============================================================================
// usage
// =============================================================================

// on peut inclure le .h autant de fois que nécessaire, et le .cc doit être inclus 
// UNE SEULE FOIS PAR PROJET.
