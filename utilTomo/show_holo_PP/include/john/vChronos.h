#ifndef __VECTRA_CHRONOS__
#define __VECTRA_CHRONOS__

// This class attemps to provide user-level functionnalities for timing programs using vChrono


#include "vChrono.h"
#include "macros.h"
#include <string.h>
using namespace std;

#define boostClock boost::chrono::system_clock
#define _MSG_LEN 500

class vChronos {

 private:
  
  char* i_message;
  bool i_verbose_p;
  bool i_paused_p;
  double i_elapsed_time;

  boostClock::time_point it_start;

  void init(char* message, bool verbose_p)
  {
    size_t len = strlen(message);
    ASSERT(len < _MSG_LEN);
    ARRAY_ALLOC(i_message, _MSG_LEN, char);

    strcpy(i_message, message);
    i_verbose_p = verbose_p;

    i_elapsed_time = 0;
    i_paused_p = false;
  }

 public:

  vChronos(const char* message, bool verbose_p)
    { init( (char*) message, verbose_p); }

  vChronos(const char* message)
    { init( (char*) message, true); }

  
  ~vChronos()
    { }
  
  void start()
  {
    if (i_verbose_p)
      { cout << "\n::: start " << i_message; cout.flush(); }
    it_start = boostClock::now(); 
  }

  void stop()
  { 
    if (! i_paused_p)
      i_elapsed_time += ( ( boostClock::now() - it_start ).count() *	\
			  ((double) boostClock::period::num / boostClock::period::den) );
    
    if (i_verbose_p)
      { 
	cout << "\n::: stop " << i_message << " in " << i_elapsed_time << "(s)"; 
	cout.flush();
      }
  }

  void pause()
  {
    if (! i_paused_p)
      {
	i_elapsed_time += ( ( boostClock::now() - it_start ).count() *	\
			    ((double) boostClock::period::num / boostClock::period::den) );
	i_paused_p = true;    
      }
  }

  void resume()
  {
    it_start = boostClock::now(); 	
    i_paused_p = false;
  }

  void clear()
  {
    i_elapsed_time = 0;
  }

  void change_message(char* new_message)
  {
    size_t len = strlen(new_message);
    ASSERT(len < _MSG_LEN);
    strcpy(i_message, new_message);
  }

  double get_elapsed_time()
  {return i_elapsed_time;}
  double get_elapsed_time_s()
  {return i_elapsed_time;}
  double get_elapsed_time_ms()
  {return 1000 * i_elapsed_time;}
  

  void set_verbose()
  {i_verbose_p = true;}
  void set_mute()
  {i_verbose_p = false;}
  
};


#endif //__VECTRA_CHRONOS__


//==============================================================================
// example
//==============================================================================

/*
#include <boost/chrono.hpp>
#include <iostream>

#include "vChronos.h"


using namespace std;


int main()
{
  vChronos toto("coucou", true);

  toto.start();

  std::cout << "Type the Enter key: ";
  std::cin.get();

  toto.pause();

  std::cout << "Type the Enter key: ";
  std::cin.get();
  
  toto.resume();

  std::cout << "Type the Enter key: ";
  std::cin.get();


  toto.stop();
  

  return 0;
  }
*/
