#ifndef __MANIP_PARSER__
#define __MANIP_PARSER__


#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <exception>
#include <map>

#include "manipTomo.h"

using namespace std;



void
manip_parser_init(std::map<string, string> map_paths, std::map<string, bool> map_bools, \
		  std::map<string, size_t> map_dims);


void
manip_parse_file(const char* cfg_file, \
		 std::map<string, string> &map_paths, std::map<string, bool> &map_bools, \
		 std::map<string, size_t> &map_dims, std::map<string, float> &map_phys, \
		 ExpType &experiment, enum ScanType &scanning);


void
manip_parse_argts(size_t argc, char** argv, \
		  std::map<string, string> &map_paths, std::map<string, bool> &map_bools, \
		  std::map<string, size_t> &map_dims, std::map<string, float> &map_phys, \
		  ExpType &experiment, enum ScanType &scanning);



#endif //  __MANIP_PARSER__
