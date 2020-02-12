// -*- mode: C++; coding: utf-8-unix ; fill-column: 80;

#ifndef _AIR_VOLUME_
#define _AIR_VOLUME_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */

   
/* PLEASE SEE AIR_Volume.txt for documentation */
   
/* ------------------------------------------------------------------------- */


#define HDR_SUFF ".hdr"
#define IMG_SUFF ".img"
#include "AIR/AIRmain.h" 
#include "AIR/src/HEADER.h" 


#include <avPoint.h>
#include <avVoxel.h>
#include <avVector.h>
#include <avSegment.h>
#include "macros.h"
#include "Filename.h"


#include <omp.h>


/* ------------------------------------------------------------------------- */


/* typedef struct hdr header; */


using namespace std; 


template <typename T> class AIR_Volume;
template <typename T> ostream& operator<< (ostream &o, const AIR_Volume<T> &v);


template <class T>
class AIR_Volume{

   


  // ===========================================================================
  // LEGACY, temporaire, à virer
 public:
   
  /* type local! */
  typedef struct {T* array; int size;} sized_array; // probablement à virer


  T* get_data_linear() { return get_data(); }
  template <class U>
    bool same_order_p(const AIR_Volume<U> &V) const
    {return true;}
  

  // ===========================================================================
 protected:

  
  // Volume number of BPP 
  size_t bits; // this is inferred by template-specialisation

  // Volume dimensions
  size_t x_dim;
  size_t y_dim;
  size_t z_dim;

  void update_voxels()
  { nb_voxels_layer = x_dim * y_dim; nb_voxels = nb_voxels_layer * z_dim; }

  // nombre total de voxels
  size_t nb_voxels; // calculé à lecture
  size_t nb_voxels_layer;

  // Volume voxel resolution
  double x_size;
  double y_size;
  double z_size;


  // Volume identifier codes
  int analyze_code;
  int AIR_code;

  // Object representing names of files used by current volume
  Filename VOL_files;

  // Volume comment
  string comment;
  
  /* for storage only */
  // Internal header structure used to communicate with .hdr files
  AIR_Key_info VOL_hdr; // foreign type
  FILE* fp_hdr;
  FILE* fp_img;
    
  T pen_color; // for drawings
  T default_color;  // color given when reading out ouf bounds

  
  // other internal state variables 
  bool allocated, allocated_external, files_attached, definitions_set, data_loaded;
  bool debug_mode;
  bool datacheck_mode;
  bool secure_mode; // if true, *some* operations will silently adopt
		    // different behaviour to overcome R/W out of
		    // volume bounds
  bool shifted;     // accounts for circular shifts performed. necessary for further Fourier classes
  bool unmatched_bits_p; // passe à vrai SSI le fichier donné pour source n'a pas le même nbre de bits
  // que le type donné en template. Un avertissement est envoyé par get_definitions_from et une erreur par update_data ou write_files;

  // ----------------------------------------------------------------------
  // data-related variables

  
  // true: data is allocated a single array and in natural order (C-order, Row-major order).
  // intended for volume processing with 1-pixel width window or for fourier algorithms.
  // T* data_linear pointer is used.
  // false: (default state)
  T *data_linear;
  // only one array allocated at once
  
  // image résidente pour traitements
  T* data_layer;

  // minima and maxima intensity values in current volume
  T glmin; 
  T glmax;

  // Misc
  char** c_type_names; // names of C type in file 
  // probablement à effacer


  // ----------------------------------------------------------------------
  /* private member functions */


  // initializes AIR volume values
  void AIR_init();


  // doubles a volume value - 1.
  T double_m1(T value);


  // associates ANALYZE datatype associated with bpp number.
  int bits_to_datatype(int bits); // déduit le code AIR à partir du nombre de bits (fixé par le format de fichier)
  bool check_bits(int bits);  // ne sait pas à quoi ca sert
  // le champ bits est lu à partir du .hdr d'un fichier existant
  // quand on CREE ex-nihilo un volume, il faut spécifier le nombre de bits 
  // cette fonction se charge de le faire pour des types connus 

 
  inline size_t infer_bits() const;
  // si inconnu, renvoie 0. Ca fait échouer la création du .hdr
 
  /*** WARING */
  // infer_bits est redéfinie pour chaque type supporté, et la version par défaut renvoie vers un code d'erreur
  // CEPENDANT, si ces fonctions ne sont pas inline, elles sont toutes compilées à chaque instanciation d'une instance
  // deux instanciations => bam, erreur au linkage
  
 public:

  ostream& display(ostream &o) const;
  
  // vérifie que les deux volumes ont les mêmes dimensions x/y/z
  template <class U>
  inline bool same_dim_p(const AIR_Volume<U> &V) const;
  // vérifie que les deux volumes ont les mêmes tailles de pixels
  template <class U>
  inline bool same_size_p(const AIR_Volume<U> &V) const;
  // détermine si volume est de type cubique
  inline bool cubic_p() const;
  // détermine si volume est isotrope
  inline bool isotropic_p() const;



  /* ------------------------------------------------------------------------- */
  /* ------------------------------------------------------------------------- */
 public:
  
  
  /* ------------------------------------------------------------------------- */
  /* constructors and destructors */
  /* ------------------------------------------------------------------------- */
  AIR_Volume();
  // input: radix(img-name|hdr-name) for which the couple radix.img/radix.hdr exists on disk
  // you MUST specify if data is supposed to be allocated in linear array mode 
  // or slice by slice, since data is loaded
  // default set to false to allow compatibility with older slice-by-slice programs
  AIR_Volume(const string &filename); 
  //, bool data_linear_mode_p = false);
  AIR_Volume(size_t x_dim, size_t y_dim, size_t z_dim)
    {AIR_init(); definitions_set = TRUE; this -> x_dim = x_dim; this -> y_dim = y_dim; this -> z_dim = z_dim; 
      update_voxels();}
  AIR_Volume(size_t x_dim, size_t y_dim, size_t z_dim, double voxel_size)
    {AIR_init(); definitions_set = TRUE; this -> x_dim = x_dim; this -> y_dim = y_dim; this -> z_dim = z_dim;
      this -> x_size = this -> y_size = this -> z_size = voxel_size;
      update_voxels();}
  AIR_Volume(size_t x_dim, size_t y_dim, size_t z_dim, double x_size, double y_size, double z_size)
    {AIR_init(); definitions_set = TRUE; this -> x_dim = x_dim; this -> y_dim = y_dim; this -> z_dim = z_dim;
      this -> x_size = x_size; this -> y_size = y_size; this -> z_size = z_size;
      update_voxels();}
    
      
  ~AIR_Volume();  


  /* ------------------------------------------------------------------------- */
  /* basic methods */
  /* ------------------------------------------------------------------------- */


  // Allocate/des. memory for AIR volume data 
  virtual void allocate();
  void unallocate();
  

  // Update AIR Volume definitions from attached hdr file.
  void update_definitions();
  // Update AIR Volume definitions from another AIR Volume.
  template <class U>
    void get_definitions_from(const AIR_Volume<U> &v);
  // Find minimum and maximum values the volume is featuring (if relevant).
  void set_glminmax(); //TODO
  // (over)write hdr file on disk from current Volume definitions.
  void write_hdr() const; 
  // Update AIR Volume data from attached img file.
  void update_data();
  // (over)write img file on disk from current volume data values.
  void write_data() const;
  // Change the files attached to AIR Volume 
  // (useful for saving current volume on other files)
  void change_files(const string &filename); 
  void write_files() const { write_data(); write_hdr(); }

  
  /* friend external fonctions */
  // '<<' overloading
  //friend ostream& operator<< <> (ostream &o, const AIR_Volume<T> &v); // nvcc cale

  
  // use on small matrices only
  ostream& show_values(ostream &o) const;


  /* ------------------------------------------------------------------------- */
  /* standard class accessors */
  /* ------------------------------------------------------------------------- */
  /*   xyz dim/size, glmin/man should be read-only */

  // !!!!! il faudrait interdire ces 3 là une fois l'allocation effectuée !!!!!!!!
  READWRITER(size_t, x_dim);      // int get_x_dim() created
  READWRITER(size_t, y_dim); 
  READWRITER(size_t, z_dim);

  READWRITER(double, x_size);
  READWRITER(double, y_size);
  READWRITER(double, z_size);
  READWRITER(T, pen_color);
  READWRITER(string, comment);  // string get_comment(), 
                                // void set_comment(string&) created
  READWRITER(int, bits);        
  // AIR old: values: 1/8/16/32 for bool/unsigned char/short int/int
  // AIR new: values: 8/16/32/64 for unsigned char/unsigned short int/float/double

  READWRITER(T, glmin);
  READWRITER(T, glmax);
  READER(Filename, VOL_files);
  READER(size_t, nb_voxels); 
  READER(bool, allocated);

  // more direct access to file names 
  const char* get_hdr_name_c() const { return VOL_files.get_hdr_name_c(); } 
  const char* get_img_name_c() const { return VOL_files.get_img_name_c(); } 
  // Give pointer to data, thus allowing input/output
  //T*** access_data() {return data_random;};

  void get_dims( size_t &dimX, size_t &dimY, size_t &dimZ)
  { dimX = x_dim; dimY = y_dim; dimZ = z_dim; }
  
  // --------------------------------------------------
  // deprecated
  READWRITER(bool, debug_mode);
  READWRITER(bool, secure_mode);
  READWRITER(bool, datacheck_mode);  

  // return the name of the C-language type contained by the file (depending on read bits value)
  // not used so far
  const char* get_c_type(); 

      
  /* ------------------------------------------------------------------------- */
  /* volume data accessors */
  /* ------------------------------------------------------------------------- */
  
  
  T& data(int idx) {return data_linear[idx];}

  

  T* get_data() const {
    MSG_ASSERT(allocated, "volume not allocated yet");
    return this -> data_linear;
  }
  

  // Ces fonctions permettent de construire une instance de AIR_Volume à partir 
  // de données déjà allouées. Le tableau est partagé et pas recopié.
  // L'instance courante ne doit pas avoir vu ses données alloués.
  // Avec des pointeurs, l'utilisateur est tenu de vérifier que les tableaux ont la même taille
  void set_data(T* data_linear);
  // Les données sont ici partagées avec un volume alloué
  void set_data(AIR_Volume<T> &V);

  // Permet d'échanger les données de deux volumes identiques et alloués. Contrôle que identiques.
  void permute_data(AIR_Volume<T> &V);
  
  // tout volume se retrouve alloué et accepte ce pointeur sans verification
  // unsafe mais nécessaire pour permutation
  void force_data(T* data_linear);
  
  /* ------------------------------------------------------------------------- */
  /* slice extractors */
  /* ------------------------------------------------------------------------- */

  // output is a pointer to an allocated image.

  void copy_XY_slice(size_t z, T* dest_slice) const;
  void copy_XZ_slice(size_t given_x, T* dest_slice) const;
  void copy_YZ_slice(size_t given_y, T* dest_slice) const;

  void copy_XY_MIP(T* dest_slice) const
  {copy_XY_IP(dest_slice, false);}
  void copy_XY_mIP(T* dest_slice) const
  {copy_XY_IP(dest_slice, true);}
  void copy_XZ_MIP(T* dest_slice) 
  {copy_XZ_IP(dest_slice, false);}
  void copy_XZ_mIP(T* dest_slice) 
  {copy_XZ_IP(dest_slice, true);}


 protected:
  virtual void tool_min_slice(T* mins_slice, T* values_slice) const;
  virtual void tool_max_slice(T* mins_slice, T* values_slice) const;
  void copy_XY_IP(T* dest_slice, bool do_min) const;
  void copy_XZ_IP(T* dest_slice, bool do_min); 
 public:

  /*   -------------------------------------------------- */
  // Legacy


  // copy nonzero voxels onto given volume
  // omp-optimized
  void copy_upon(AIR_Volume<T> &output_vol); //const

  // ATTENTION, ça n'a pas été testé ni optimisé à fond!!!
  // xy -> xy (défaut)
  void copy_slice(size_t in_slice, AIR_Volume<T> &output_vol, size_t out_slice) const;
  // xz -> xy
  void copy_yz_slice(size_t in_slice, AIR_Volume<T> &output_vol, size_t out_slice) const;  


  /* ------------------------------------------------------------------------- */
  /* operators */
  /* ------------------------------------------------------------------------- */


  // pour tout opérateur utilisé entre deux instances de volumes, par défaut:
  // * on ne tolèrera que le cas où les deux volumes sont chargés et identiques 
  //   en dimension et organisation des données
  // * on n'admet pas de transtypage, réservé à copy_from

  // opérateur de recopie
  // ici, on admet qu'aucun volume ne soit chargé (mais pas l'un XOR l'autre)
  void operator= (const AIR_Volume<T> &V);

  // ici, on admet le transtypage par défaut, à charge de l'utilisateur de vérifier 
  // que le cast s'effectue correctement
  template <class U>
    void operator+= (const AIR_Volume<U> &V);
  template <class U>
    void operator*= (const AIR_Volume<U> &V);

  void operator*= (double factor);
  void operator+= (T increment);  
  T* operator[] (size_t slice_number)
  { return (T *)(data_linear + (nb_voxels_layer * slice_number)); }

  /* ------------------------------------------------------------------------- */
  /* data operations (whole volume processed) */
  /* ------------------------------------------------------------------------- */

  
  // effectue une permutation circulaire des 8 semicubes afin de transposer 
  // le 0 situé en bas de l'image au centre
  // n'a de réel intérêt qu'en vue d'une TF3D (cuda cuFFT, fftw)
  void circshift_to(AIR_Volume<T> &V_dest) const;

  // inverse le flag *shifted* de la destination (faux à création).
  // on est hélas obligé de faire cette opération avec deux volumes alloués.
  
  // dans le cas où on charge un fichier que l'on sait shifté
  void declare_shifted()
  {shifted = true;}

  // swaps byte ordering (endianness) 
  void swap_endian();
  
  bool in_bounds(size_t x, size_t y, size_t z) const
  {return BOUNDEDP(x, 0, x_dim) && BOUNDEDP(y, 0, y_dim) && BOUNDEDP(z, 0, z_dim);}  

  // the 3 following functions are a little delicate to implement using map

  // Si data[i] > T / data[i] = value / sinon 0
  void binarize_data(T threshold, T value = 1);

  // retrieves (sets min & max to) minimum and maximum values of current volume
  void get_min_max(T &min, T &max);

  // data values coerced to provided interval with boundary adjustment
  void coerce_data(T min, T max);

  // data values coerced by cutoff
  void truncate_data(T min, T max);

  // fill volume with given value
  void fill(T value);
  void double_values();
  void fill_zero();

  // multiply volume elements by scalar
  void scalar_mult(T value);

  // sélection des valeurs par masque binaire. 
  // valeur != 0 dans le masque => conserve voxel correspondant dans le volume
  // sinon, mise à 0
  template <class U>
    void select_values(const AIR_Volume<U> &bin_mask);
  
  // map* a 1-arg function to current matrix
  // M.map(&sqr) => elevate M to square
  void map(T (*map_fun)(T));
  
  // applies* a 2-args function to current matrix and given matrix
  // M.apply(+, N) stores the result of M+N in M
  void map2_star(T (*map_fun)(T, T), const AIR_Volume<T> &V); 

  
  // retrieves sum of values
  T sum() const;

  // retrieves number of non-zero voxels
  size_t nzero_count() const;


  /* ------------------------------------------------------------------------- */
  /* data import/export */
  /* ------------------------------------------------------------------------- */
  
  // fills data values from data coming from another class instance, same type or not
  // checks: both instances must have data allocated, data of same dimensions, data ordered the same manner
  // if T != U, data is simply casted
  template <class U>
    void copy_from(const AIR_Volume<U> &V_Source);
  
 // those versions allow for importing external data in linear or slice/by/slice mode
  // user must ensure it matches the current context (size and data)
  template <typename U>
    void copy_from(const U* linear_array);
  // memcpy
  void copy_from(const T* linear_array);

  template <typename U>
    void copy_from(const U*** random_array); //vtk or not
  
  // import data from a smaller volume
  // imported data is centered
  // to this point, we only implemented linear allocation mode
  template <class U>
    void copy_from_smaller(const AIR_Volume<U> &V_Source);

  // conversely,
  template <class U>
    void copy_from_bigger(const AIR_Volume<U> &V_Source);
  
  // import data from bigger volume : X and Y are the same, 
  // Z source > Z and so upper and lower layers are squeezed 
  void copy_from_squeezing_layers(const AIR_Volume<T> &V_Source);

  // both volumes have same size, but we only copy a centered subpart of src into dst
  // edge: size of centered square
  template <class U>
    void copy_centered(const AIR_Volume<U> &V_Source, size_t edge);

  // the same, but applied destructively on current volume (borders set to 0);
  void crop_centered(size_t edge_size);


  /* TODO */
  // copied data is centered
  // to this point, we only implemented linear allocation mode
  // template <class U>
  /*   void copy_to_smaller(AIR_Volume<U> &V_Dest) const; */

  


  /* ------------------------------------------------------------------------- */
  /* Points section */
  /* ------------------------------------------------------------------------- */


  /* ------------------------------------------------------------------------- */
  // Section valable uniquement pour les volumes alloués en linéaire

  // basic data access, no control
  T get(size_t x, size_t y, size_t z) const;
  void set(size_t x, size_t y, size_t z, T value);

  // full secure access
  void write_Point_linear(size_t x, size_t y, size_t z, T color);
  void write_Point_linear(size_t x, size_t y, size_t z) {write_Point_linear(x, y, z, pen_color);}
 

  /* ------------------------------------------------------------------------- */
  /* Section morphologie mathématique LINEAIRE */
  /* ------------------------------------------------------------------------- */

 private:
  AIR_Volume<T> *V_morpho_cache;
  inline void fill_square(T* data, size_t pos, T color);
  
 public:
  
  // un appel à dilatation (fermeture= requiert l'allocation d'un volume cache persistant. Cet appel le libère
  void morpho_clear_cache(); 


  // calcule la somme des 6-voisins du voxel donné
  // nécessite volume binarisé 0-1
  size_t morpho_l_neighbors_6(size_t x, size_t y, size_t z);
  size_t morpho_l_neighbors_26(size_t x, size_t y, size_t z);

  // cette fonction effectue une pseudo-erosion
  // les voxels dont la connexité est < au seuil fourni sont effacés
  // i.e on y applique la couleur "noire" spécifiée
  void morpho_l_erode_6_under(size_t seuil_connex, T erode_color);
  void morpho_l_erode_26_under(size_t seuil_connex, T erode_color);
  // cette fonction alloue et libère un vol° temporaire

  // erosion tout court
  // seuls sont conservés les voxels dont la connexité est (respectivement) 6 / 26
  void morpho_l_erode_6() { morpho_l_erode_6_under(6, T(0)); }
  void morpho_l_erode_26() { morpho_l_erode_26_under(26, T(0)); }

  // dilatation 
  // chaque voxel allumé est remplacé par l'élément structurant choisi (6: une croix 3D, 26: un cube)
  // alloue un volume cache jusqu'à appel de morpho_clear_cache
  void morpho_l_dilate_6(T dilated_color);
  void morpho_l_dilate_26(T dilated_color);

  // ouverture: erosion suivie d'une dilatation
  void morpho_l_open_6(size_t erode_seuil_connex, T color);  
  void morpho_l_open_26(size_t erode_seuil_connex, T color);
  
  // fermeture: dilatation suivie d'une erosion
  void morpho_l_close_6(size_t erode_seuil_connex, T color); 
  void morpho_l_close_26(size_t erode_seuil_connex, T color);  

  // private
  void morpho_l_paint_6(size_t x, size_t y, size_t z);
  void morpho_l_paint_26(size_t x, size_t y, size_t z);


  /* ------------------------------------------------------------------------- */
  /* misc */

  // remplit le volume avec une aposidation (blanc au milieu, noir sur les bords) 
  // permet de ne pas analyser les bords si on fait la TF d'un volume spatial -> fréquentiel
  // il suffit ensuite de multiplier le volume image par le volume tukey
  void tukey_window(float alpha);
  // alpha < 2 (typiquement 0.20. 0 = pas d'apodisation);

  
  void flip();

};

template <typename T>
void circshift_do(T* src_array, T* dst_array, size_t x_dim, size_t y_dim, size_t z_dim);


#include "AIR_Volume.cc"




#endif /* _AIR_VOLUME_ */

