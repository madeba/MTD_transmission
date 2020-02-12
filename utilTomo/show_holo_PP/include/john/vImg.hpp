#ifndef __VECTRA_IMAGE__ 
#define __VECTRA_IMAGE__


#ifndef OMP_THREADS
#define OMP_THREADS 3
#endif

#include "macros.h"
#include "vectra.h"
#include "AIR_Volume.h"

#include <pgmcode.h>

#include <iostream>
using namespace std;

#include "cv.h"  
#include "highgui.h"


// ********************************************************************************

/*

  Bibliothèque de base pour gérer le chargement, la sauvegarde et des
  opérations basiques sur les images 2D en NDG.

  A ne pas confondre avec cvDisplayArray, qui se focalise uniquement sur
  l'affichage en fenêtres, les évènements et en particulier la
  quantification à partir de données flottantes et doubles même pas normalisées


  On utilise comme container cv::Mat, simple et fonctionnel et universel. 

  Le Mat peut être entièrement géré par OpenCV, ou constitué à partir d'un tableau préalloué (e.g SSE) qui sera préservé dans tous les cas

  Une instance de classe contient nécessairement un objet Mat alloué, qui peut être retourné pour traitement/aff. opencv

  Il est possible de désallouer une image pour économiser la mémoire, mais c'est une opération irréversible

  Chargement et sauvegarde sont confiés par défaut à opencv, qui gère plein de filtres assez efficacement.
  Cas particulier du PGM 8 bits: 2x plus rapide que sur OpenCV, et permet un refresh image sans réallocations

*/

// ********************************************************************************




// ********************************************************************************


template <typename T>
class vImg{


  // ==================================================
private: 
  

  // ----------------------------------------
  // core 
  
  size_t x_dim;
  size_t y_dim;

  size_t i_array_size; 
  size_t i_nb_threads; 

  // cv::Mat est le conteneur de données principal
  cv::Mat i_cv_img;   // en cas de reload opencv, Mat est détruit!
  T* i_data;

  // une instance est construite allouée. Mais on autorise la désallocation pour économie, irréversible
  bool i_allocated_p;  
  // a-t-on alloué de manière externe (évt. SSE) ou bien en interne-opencv? nécessaire pour reload
  // nb: opencv ne désalloue pas un tableau quand il a construit un Mat sur ce tableau
  bool i_allocated_manual_p;
 
  // PGM spécifique
  bool i_startx_set_p;
  long i_startx;

  // ==================================================
public:

  
  // ----------------------------------------
  // constructeurs 

  // une instance est nécessairement allouée, en interne ou en externe.
  // le constructeur cv::Mat est utilisé pour allocations internes
  // on accepte les tableaux de données si allocation externe 

  inline vImg(); // interdit: une instance est créée allouée
  inline vImg(size_t x_dim, size_t y_dim, size_t nb_threads = 3);
  inline vImg(T* data, size_t x_dim, size_t y_dim, size_t nb_threads = 3);
  inline vImg(string filename, size_t nb_threads = 3);
  ///vImg(const vImg<T>& I);
  inline ~vImg();

private:
  inline void init(size_t x_dim, size_t y_dim, size_t nb_threads);

  // ----------------------------------------
  // memory 
public:
  inline void unallocate(); // liberté de désallouer, mais instance inutilisable.

private:
  inline void allocate();
  // méthodes surchargées selon type
  inline void allocate_mat(); // alloue un Mat
  inline void read_mat(string filename);
  inline void generate_mat(T* data_passed); // crée un Mat autour du tableau, préservé

  
  // ----------------------------------------
  // accessors

public:

  inline bool same_dim_p(const vImg<T>& img) const;
  READER(size_t, x_dim);
  READER(size_t, y_dim);
  NREADER(bool, i_allocated_p, get_allocated_p);
  NREADER(size_t, i_nb_threads, get_nb_threads);
  NWRITER(size_t, i_nb_threads, set_nb_threads);

 
  NREADER(T*, i_data, get_data); 
  // attention, le Mat est un objet temporaire, qui doit être MAJ en cas de reload_opencv
  //cv::Mat get_cv_Mat() const { return i_cv_img; }
  NREADER(cv::Mat, i_cv_img, get_cv_Mat);


  // ----------------------------------------
  // Memory I/O

  // recopie d'objets du même type

  inline void copy_from(T* data);
  // edge designe le coin supérieur gauche de la fenêtre de découpe dans l'image alien source 
  inline void copy_from_bigger(T* alien_data, size_t alien_w, size_t alien_h, size_t edge_x, size_t edge_y);
  inline void copy_from_bigger(const vImg<T>& bigger, size_t edge_x, size_t edge_y);
  // copie centrée
  void copy_from_bigger(const vImg<T>& bigger) 
  { int edge_x = (bigger.get_x_dim() - x_dim) / 2; 
    int edge_y = (bigger.get_y_dim() - y_dim) / 2; 
    copy_from_bigger( bigger, edge_x, edge_y );
  }
    

    
  // ----------------------------------------
  // File I/O

    
  // générique: l'image est sauvée en RAW quelque soit le type de données
  // visualisable avec cvVol, vtkVol, voire imageJ
  inline void save_AIR(string air_file_name);


  // met à jour l'instance courante à partir d'un nouveau fichier. si check_size, vérifie que même taille
  // l'implémentation opencv_imread oblige à allouer et désallouer, bien que testé sans perte temps
  inline void reload(string filename, bool check_size = false);
  // puisque notre bibliothèque PGM est 2x plus rapide qu'OpenCV et ne CONSe pas, on propose des méthodes
  inline void reload_pgm(string filename_pgm, bool checks = false);
  inline void save_pgm(string filename_pgm);
  void save(string filename)
  { cv::imwrite(filename, i_cv_img  ); }
  
  // EXPERIMENTAL
  // mettre à jour l'instance courante à partir d'un nouveau pointeur. Mat pas affecté
  void update_pointer_single( T* data ) 
  { 
    i_data = data; 
    i_cv_img.data = (unsigned char*) i_data; //.data est uc
  }
  
  
  
  // ----------------------------------------
  // circshift

  // centré (cf FFTW_Image)
  inline void circshift_to(vImg<T>& dst);
  inline void circshift_exc_to(vImg<T>& dst, size_t centreX, size_t centreY);


  // ----------------------------------------
  // image processing

  inline void fill(T value);
  inline void clear();  // fills zeros
  inline void tuckey(float alpha); // remplit image d'un masque de hamming / fenêtre de Tuckey
  // alpha < 2 (typiquement 0.20. 0 = pas d'apodisation);

  inline void module(vImg<T> real, vImg<T> imag, bool no_check = false);
  inline void module(T* real, T* imag); // no size checks

  inline void phase(vImg<T> real, vImg<T> imag, bool no_check = false);
  inline void phase(T* real, T* imag); // no size checks
  
  inline void draw_rectangle(size_t ul_x, size_t ul_y, size_t lr_x, size_t lr_y, T color);

  // ----------------------------------------


  // ----------------------------------------
  // fonctions auxilliaires 
  // ----------------------------------------

private:  

  inline void extract_subImage(T* src, T* dst, size_t src_dimx, size_t src_dimy, size_t edge_x, size_t edge_y, size_t dst_dimx, size_t dst_dimy);
  
  inline void reload_fast(string path); 
  inline void reload_full(string path); 
  inline void reload_pgm_fast(string path); 
  inline void reload_pgm_full(string path); 

  inline T norm(T r, T i);

  inline T phase(T r, T i); 

  inline T quantification_factor();
};



// =============================================================================
// implémentation
// =============================================================================

#include "vImg.cc"



#endif // __VECTRA_IMAGE__








/*

  8UC[1-4]
  8SC[1-4]
  16UC[1-4]
  16SC[1-4]
  32SC[1-4]
  32FC[1-4]
  64FC[1-4] 

  For popular image encodings, CvBridge will optionally do color or pixel depth conversions as necessary. To use this feature, specify the encoding to be one of the following strings:

  mono8: CV_8UC1, grayscale image

  mono16: CV_16UC1, 16-bit grayscale image

  bgr8: CV_8UC3, color image with blue-green-red color order

  rgb8: CV_8UC3, color image with red-green-blue color order

  bgra8: CV_8UC4, BGR color image with an alpha channel

  rgba8: CV_8UC4, RGB color image with an alpha channel 

*/
