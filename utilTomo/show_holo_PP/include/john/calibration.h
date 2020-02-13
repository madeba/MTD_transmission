// Code Hui pour autocalibration

#ifndef __HUI_AUTOCALIBRATION__
#define __HUI_AUTOCALIBRATION__


#include "Camera.h"
#include "Ljack.h"
#include "Ljack_DAC.h"



bool manip_calibration(Camera* GEV_cam, Ljack_DAC* DAC_phase, char* out_dir_calib, unsigned int nb_iters, float modphase_increment, unsigned int delay_phase_ms, unsigned int delay_camera_ms);


double calibration(const int premier_plan, const int Num_Angle_final, const int SautAngle, char *Chemin, size_t img_dimx, size_t img_dimy, size_t window_dimx, size_t window_dimy, size_t corner_x, size_t corner_y);


#endif // __HUI_AUTOCALIBRATION__
