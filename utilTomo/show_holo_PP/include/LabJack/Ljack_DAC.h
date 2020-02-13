#ifndef __LJACK_DAC__
#define __LJACK_DAC__

#include "macros.h"
#include "msleep.h"

#include "Ljack.h"
// provisoire
#include "Ljack_DAC_plugin.h"



using namespace std;




class Ljack_DAC {

 private:

  u3TdacCalibrationInfo caliInfo;
  Ljack *host_labjack;
  HANDLE host_handle;
  // true if connected to free slot and calibrated OK
  bool connected;
  Ljack_DAC_plugin DAC_internal;

  
  // outputs a (coerced to) proper voltage on current exit
  void set_DAC_output(float voltage);


  // simulation
  bool simulate_p;
  size_t latency_microsecs; // utile pour simuler l'attente du labjack


 public:

  
  Ljack_DAC(Ljack &host_labjack)
    {
      connected = false;
      this -> host_labjack = &host_labjack;
      host_handle = host_labjack.get_handle();
      
      simulate_p = host_labjack.simulate_p; 
      if (! simulate_p)
	{
	  MSG_ASSERT(configIO_DAC(host_handle) == 0, \
		     "échec config DAC");
	}
      else
	latency_microsecs = host_labjack.latency_microsecs;
    }

  // connect current DAC expansion to specified slot of the labjack host
  // returns false if slot occupied or initialization failed
  bool connect(enum LjackSlot freeslot);

  // requires established connection
  // voltage will be coerced to [-10; +10v]
  void set_A_output(float voltage);
  void set_B_output(float voltage);

};


#endif // __LJACK_DAC__
