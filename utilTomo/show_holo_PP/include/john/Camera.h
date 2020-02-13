#ifndef __CAMERA__
#define __CAMERA__


#include <assert.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "msleep.h"
#include "macros.h"

#include <PvDeviceFinderWnd.h>
#include <PvDevice.h>
#include <PvBuffer.h>
#include <PvStream.h>
//#include <PvStreamRaw.h>
#include <PvBufferWriter.h>
#include <PvString.h>
#include <PvSystem.h>

#include <pgmcode.h>
#include "cv.h"
#include "highgui.h"


#include <boost/filesystem.hpp>
using namespace boost::filesystem;


#include "types_Camera.h"

using namespace std;



// =============================================================================
/* =============================================================================

   La classe Camera effectue un binding basique de l'API Pleora.
   Gère acquisition d'images en mode:
   * flux bufferisé (aussi vite que possible mais vitesse potentiellement variable)
   * synchrone (images prises une par une, au moment précis où l'API le demande. Gestion des attentes dans le programme)


   Le format par défaut de sauvegarde est désormais le PGM, la caméra étant en n&b. Le BMP est encore supporté.

   La caméra possède des méthodes genre read_ImgX et set_ImgX.
   * read_ImgX va lire dans la caméra connectée la valeur déjà configurée avec GEVplayer ou autre
   * set_ImgX effectue deux actions: d'une, il stocke la valeur demandée par l'utilisateur, qui conditionne l'allocation de buffers locaux, et de deux, il écrit cette valeur dans la mémoire de la caméra.

   La classe supporte un mode simulation permettant de faire fonctionner la classe à partir de données aléatoires ou bien lues dans images lues en boucle dans un répertoire.

   cf Camera.txt pour doc étendue et exemples.

   */
// =============================================================================
// =============================================================================




// =============================================================================
class Camera{

  // ----------------------------------------
 protected:
  // ----------------------------------------


  // ----------------------------------------
  // mode d'émulation pour fonctionner sans connexion
  //
  bool i_SIMULATE_p;  // simulation disque et non lecture caméra
  bool i_simulation_buffered_p; // la simulation lit-elle les fichiers ou bien une pile mémoire?
  unsigned char *i_simulation_stack; // pile mémoire créée au début à partir des fichiers si exigé.

  // devrait être privé
  size_t i_simu_images_count;
  size_t i_simu_image_current;
  vector<path> i_simu_path_vector;
  PvBuffer *i_simu_Buffer; // données brutes par i_simu_Buffer -> GetDataPointer();



  PvDeviceInfo* i_DeviceInfo;
  PvDevice i_Device;

  PvStream lStream;
  PvBuffer *lBuffers, *lLastBuffer;
  PvBuffer lBuffer;
  size_t lBufferCount;
  bool camera_connected, acquisition_snap_ready, acquisition_stream_ready, log_images, realtime_report, retries_report, pgm_buffer_allocated;


  /// IL REVIENT AU PROGRAMMEUR DE FIXER UN IT_EXPOSURE et un DELAY_EXPOSURE MUTUELLEMENT COHERENTS
  // IT_EXPOSURE : float, pour communiquer avec l'api Pleora de la camera
  // delay_exposure : int (deci-millisecondes) : pour gérer les attentes lors de la prise de vue
  // e.g 12 => 1.2ms

  // temps d'attente liés à l'exposition du CMOS et au readout de la caméra
  size_t delay_exposure, delay_readout;  // unité: déci-milliseconde, entier
  // distinct du temps d'exposition réel, qui doit être nécessairement inférieur au temps d'attente
  // défini par défaut à une valeur assez large pour éviter le stall. cf crunch_readout() pour accélérer

  // non vérifié
  float it_Exposure; // temps d'attente exact lu/ecrit par API Pleora


  // log images
  char *lLastBlockID, *lLastStamp;
  FILE* LOGFILE;
  PvGenCommand *lResetTimestamp;

  // realtime report
  PvGenParameterArray *lStreamParams;
  PvGenInteger *lCount;
  PvGenFloat *lFrameRate;
  PvGenFloat *lBandwidth;

  // commandes camera et variables état caméra
  PvGenInteger *lTLLocked, *lPayloadSize;
  PvGenCommand *lStart, *lStop, *lReset, *lTemps;
  PvGenEnum *lAcquisitionMode;

  // commandes permettant de changer les dimensions des images prises
  // par défaut, on prend celle déjà fixées par GEVplayer.
  // Mais si l'utilisateur prend la peine de modifier ces valeurs, on les envoie à la caméra
  PvGenInteger *ic_DimX, *ic_DimY, *ic_TempSensorLocal, *ic_TempSensorRemote, *ic_TempADCLocal;
  PvGenFloat *ic_Exposure;
  PvGenEnum *ic_PixDepth, *ic_Gain, *ic_LinLog;
  bool ib_DimX_set, ib_DimY_set, ib_PixDepth_set, ib_Exposure_set, ib_Gain_set, ib_Linlog_set;

  // paramètres effectifs
  size_t id_DimX, id_DimY;
  enum PixDepth ie_PixDepth; // 8/10/12 (bits) "Mono8" 10 12
  enum GainValue ie_Gain; // x1/x2/x4
  enum LinlogValue ie_Linlog; // compression: off/low/normal/high
  int itc_sensorLocal, itc_sensorRemote, itc_ADCLocal;

  const char* i_IpMac;

  // realtime report and log
  char lDoodle[];// = "|\\-|-/";
  int lDoodleIndex;
  int64_t lImageCountVal;
  double lFrameRateVal;
  double lBandwidthVal;

  PvResult lAcquisitionResult, lOperationResult, lSaveResult;

  PvBufferWriter lBufferWriter;
  unsigned char *lPgmBuffer;

  // info qui découle de i_Device, mais comme i_Device est conservé par adresse, on a besoin de stocker l'IP en dur.
  char* device_ip;

  size_t i_LAYER_MEM_SIZE;

  // reset internal variables
  void reset()
  {
    camera_connected = acquisition_snap_ready = acquisition_stream_ready = log_images = realtime_report = retries_report = pgm_buffer_allocated = i_simulation_buffered_p = false;
    ib_DimX_set = ib_DimY_set = ib_PixDepth_set = ib_Exposure_set = ib_Gain_set = ib_Linlog_set = false;
  }


  // ----------------------------------------
 public:
  // ----------------------------------------


  Camera();
  ~Camera();


  //------------------------------
  // active le mode simulation pour effectuer des essais sans caméra.
  //
  // les temporisations sont fidèlement reproduites, de même que les transferts mémoire et disque.
  // si un chemin correct est fourni vers des images png, elles vont servir de source de données (attention au coût de traitement)

  // charge les images une par une
  void set_simulation_mode(string pgm_dir_path="");
  // crée une pile d'images mémoire: long au démarrage mais timings plus précis (XOR unbuffered)
  void set_simulation_buffer_mode(string pgm_dir_path="");
  bool get_simulation_mode() const;

  // parfois nécessaire
  NREADER(size_t, i_simu_images_count, get_simu_images_count);


  //------------------------------
  // Soit on fait select_camera + connect, soit connect_firstCamera

  // choisir une caméra par menu QT parmi les devices connectées en ethernet
  bool select_camera();

  // réaliser la connexion ethernet directe avec la caméra sélectionnée
  bool connect();

  // connexion sans sélection, sur la première caméra disponible
  bool connect_firstCamera();

  // fermer la connexion
  void disconnect();



  // MAJ des valeurs des temperatures.
  void read_Temperatures();

  // lecture des temperatures (mise a jour effectuee à connexion et sur appel de read_Temperatures
  READER( int, itc_sensorLocal );
  READER( int, itc_sensorRemote );
  READER( int, itc_ADCLocal );
  READER( bool, camera_connected );
  // get_itc_sensorLocal(),  get_itc_sensorRemote(), get_itc_ADCLocal(), get_camera_connected


  // --------------------------------------------------
  // variables d'état

  // comportement par défaut: on ne change pas les valeurs en mémoire dans la caméra
  // sauf si explicité dans le programme par un set_TOTO
  // * si caméra connectée, va changer la valeur en ligne
  // * autrement, se contente de la stocker dans la classe et la changera à connexion de la caméra

  // 1312-896 by 32 multiples
  void set_ImgX(size_t width);
  // 1082-1
  void set_ImgY(size_t height);
  // _8b / _10b / _12b
  void set_ImgDepth(enum PixDepth bits);
  // 0.01-335 ms
  void set_Exposure(float t_ms);
  // x1 / x2 / x4
  void set_Gain(enum GainValue gain);
  // off / low / normal / high
  void set_LinlogMode(enum LinlogValue linlog);


  // va lire en ligne les valeurs contenues dans la camera,
  // de manière déconnectée des valeurs contenues dans la classe
  size_t read_ImgX() const;
  size_t read_ImgY() const;
  enum PixDepth read_ImgDepth() const;
  enum GainValue read_Gain() const;
  enum LinlogValue read_LinlogMode() const;
  float read_Exposure() const;


  // à n'utiliser qu'en urgence, caméra connectée
  void factory_reset();


  // ============================================================
  // mode flux bufférisé: les images sont prises aussi vite que possible
  // insérées dans un buffer et streamées
  // pas de synchro possible a priori
  // ============================================================


  // initialiser un flux d'image pour une prise de vue en flux et en asynchrone (thread séparé)
  // pas certain qu'elle puisse être apellée avant connection, à tester.
  // en général, on essaie d'allouer au moins 16 buffers si l'on peut: plus il y en a, plus vite on va
  size_t initiate_stream(size_t wished_bufCount);

  // marque le début de la phase d'acquisition en flux, mais n'acquiert pas d'images
  void start_stream()
  {
    assert(acquisition_stream_ready);
    lStart -> Execute();
  }

  // get one image from current open stream
  size_t get_image();

  // mettre fin à l'acquisition en flux
  void stop_stream() { lStop -> Execute(); }

  // après start et stop: libérer le flux
  void terminate_stream();


  // ============================================================
  // mode snapshot: une image est prise au moment indiqué
  // on est obligé de créer un stream, mais il ne contient qu'une image
  //
  // prépare une série de prises de vue uniques en mode synchrone
  // exposure correspond à l'attente observée par le programme pour exposition du CMOS
  // cela devrait correspondre avec la valeur fixée par un appel de set_Exposure
  // / valeur d'exposition déjà enregistrée dans la caméra
  // ============================================================


  // --------- SNAPSHOT CFG  -----------

  // helper functions: configures snapshot acquisition "properly"
  // fixe le readout correctement pour des cas limités, merde autrement
  void configure_snapshot();//size_t x_dim, size_t y_dim, size_t exposure_ms, enum PixDepth bits);

  // cette version nous laisse définir le readout au feeling
  void configure_snapshot(size_t readout_ms);//size_t x_dim, size_t y_dim, size_t exposure_ms, enum PixDepth bits, size_t readout_ms);

  // expérimental, unsafe (factor de 0à1)
  void crunch_readout(float factor);


 private:
  void initiate_snapshot(size_t exposure_decimillisecs);


  // temps de readout respecté en mode snap (deci-milliseconds)
  READWRITER(size_t, delay_readout);
  // get_delay_readout(), set_delay_readout(const size_t &val)

 public:


  // --------- SNAPSHOT BASICS -----------

  //* début de prise d'une image

  // start taking snapshot: halts during exposure time set at initialization
  void start_snap();

  // end of snapshot (TODO: should/may halt during readout time)
  void stop_snap();

  // retrieve snapshot from camera and copies it into lLastBuffer
  void get_snapshot();
  // if simulation activated and configured for a data source, fills lLastBuffer with provided image file / stack image
  void get_snapshot_simu();

  // effectuer ici appel de fonction pour sauver l'image dans lLastBuffer -> GetDataPointer();
  // save_pgm8image (raw, bmp), ou save_ramimage

  void release_image();

  //* fin de prise d'une image

  // à la fin de toutes les prises
  void terminate_snapshots();



  // ============================================================
  // sauvegarde d'images : commun aux modes FLUX et SNAP

  // sauvegarde la dernière image lue par la caméra
  // image contenue dans lLastBuffer -> GetDataPointer();
  // sur le disque au format BMP/PGM/etc
  // ou bien à l'adresse mémoire indiquée

  // échoue si on le tente avant le readout: dans ce cas, il faut releaser l'image et tout recommencer


  bool save_ramimage(unsigned char* data); // id_DimX * id_DimY
  void save_ramimage_nr(unsigned char* data); // id_DimX * id_DimY

  // code natif Pleora -> raw, bmp disque
  bool save_bmpimage(const char* filename); // natif Pleora
  bool save_rawimage(const char* filename); // natif Pleora

  // autres formats disque: sauvegarde sur tampon interne par save_ramimage, puis conversion par lib externe
  bool save_pgm8image(const char* filename);  // méthode PGMlib

  // court-circuit:
  bool save_buffer_pgm8image(const cv::Mat *frame ,const char* filename);
  bool save_buffer_pgm8image(const IplImage *frame ,const char* filename);

  // transfert dans une image OpenCV en mémoire (via save_ramimage)
  bool save_ram_cv8image(IplImage *dst_image);
  bool save_ram_cv8image(cv::Mat *dst_image);
  // soit pour display frame, soit pour sauvegarde (cf Camera.txt).


  // --------- HELPERS DE BASE EN MODE SNAPSHOT  -----------
  // implémentent le cycle complet de prise d'une image snapshot

  // NB: il faudrait faire une classe fille pour faire ça propre!!!!


  /* ------------------------------ */
  // 1) sauvegarde dans un fichier PGM sur (ram)disque, simu gérée
  void take_save_snapshot(const char* filename_pgm);

 private:
  void take_save_snapshot_std(const char* filename_pgm);
  void take_save_snapshot_simu(const char* filename_pgm);
 public:


  /* ------------------------------ */
  // 2) sauvegarde dans un buffer OpenCV affichable/sauvable, simu gérée
  void take_cv_snapshot(IplImage *dst_image);
  void take_cv_snapshot(cv::Mat *dst_image);
  void take_ram_snapshot(unsigned char *raw_image);

  // l'image peut être allouée par:
  IplImage* allocate_cv8image();
  cv::Mat* allocate_cv8image_Mat();

 private:
  void take_cv_snapshot_std(IplImage *dst_image);
  void take_cv_snapshot_std(cv::Mat *dst_image);
  void take_ram_snapshot_std(unsigned char* raw_image);
  void take_cv_snapshot_simu(IplImage *dst_image);
  void take_cv_snapshot_simu(cv::Mat *dst_image);
  void take_ram_snapshot_simu(unsigned char* raw_image);
 public:


  /* ------------------------------ */
  // 3) diverses versions accélérées de take_save_snapshot_std.
  // en versions pgm ou bmp. stabilité à voir.
  void fast_bmp(const char* filename_bmp);
  // same as before, but interleaving calls to spare time. readout is passed by
  void furious_bmp(const char* filename_bmp);
  void suicide_bmp(const char* filename_bmp);
  void suicide_pgm(const char* filename_pgm);


};



#endif /* __CAMERA__ */


/* *************************************************************************** */
/* *************************************************************************** */


