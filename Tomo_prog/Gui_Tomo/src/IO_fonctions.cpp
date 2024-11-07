#include "IO_fonctions.h"
#include <iostream>
//#include <cstring>
//#include <string>

#include <sys/stat.h>
///utliser la structure stat
using namespace std;
bool IsPathOk(const std::string &directory_name)
{
  struct stat buffer;
  return (stat (directory_name.c_str(), &buffer) == 0);
}
