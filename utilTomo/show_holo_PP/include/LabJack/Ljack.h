#ifndef __LJACK_U3__
#define __LJACK_U3__

// define labjack U3 content
#include "u3.h"
#include <iostream>

// define labjack DAC plugin returned when DAC is slotted in
//#include "Ljack_DAC.h"

#define __SIMU__ false
#define LJ_SLOTS 4

using namespace std;



enum LjackSlot {
  FIO7_6 = 0,
  FIO5_4 = 1,
  AIN2_3 = 2,
  AIN0_1 = 3
};
//int calibration_indexes[4] = {6, 4, 2, 0};
//int cominit_indexes[4] = {7, 5, 3, 1};


class Ljack{


 private:  

  HANDLE hDevice;
  bool *occupied_ports;
  // Valeurs possibles: 7 (ports FIO7/6), 5 (FIO5/4), 3 (AIN2/3), 1 (AIN0/1)
  
 public:

  // simulation
  bool simulate_p;
  size_t latency_microsecs; // utile pour simuler l'attente du labjack


 public:

  Ljack(bool simulate = false);
  ~Ljack();
  
  void display(ostream &o) const;

  HANDLE get_handle() const
  { return hDevice; }
  
  // is given slot free?
  bool slotFree_p(enum LjackSlot slotID) const
  { return !occupied_ports[slotID]; }

  // marks free slot as occupied
  bool slotBook(enum LjackSlot slotID)
  {
    if (occupied_ports[slotID]) return false;
    return (occupied_ports[slotID] = true);
  }

  void set_simulate(size_t latency_microsecs = 1200)
  { 
    simulate_p = true;
    this -> latency_microsecs = latency_microsecs;
  }

};







#endif // LJack





/*
enum cali_index {
  FIO7_6 = 6,
  FIO5_4 = 4,
  AIN2_3 = 2,
  AIN0_1 = 0
};

enum com_index {
  FIO7_6 = 7,
  FIO5_4 = 5,
  AIN2_3 = 3,
  AIN0_1 = 1
};
*/
