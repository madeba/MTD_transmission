#ifndef __VECTRA_UTILS
#define __VECTRA_UTILS

#include <sys/stat.h>
#include <stdio.h>


// namespace for hosting general-purpose utility functions


namespace vectra {

  // teste si le fichier existe ET est lisible
  bool file_exists_p(const char* filename);

  // renvoie la taille d'un fichier
  unsigned int file_size(FILE* fp);

  // teste si le répertoire existe (lisible ou non)
  bool dir_exists_p(char* filename);
  
  // while (! vectra::kbhit);
  bool kbhit(void);
  
  // vectra::makedir(argv[2], 0777);
  bool makedir(const char *path, mode_t mode);

};


#endif // __VECTRA_UTILS
