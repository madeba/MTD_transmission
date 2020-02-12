#ifndef __CAMERA_TYPES
#define __CAMERA_TYPES


#define __CAM_DEPTHCOUNT 3
#define __CAM_GAINCOUNT  3
#define __CAM_LINLOGCOUNT  4

enum PixDepth {
  _8b = 0,
  _10b = 1,
  _12b = 2
};

const char i_PixDepthNames[__CAM_DEPTHCOUNT][7] = { "Mono8", "Mono10", "Mono12" }; 


enum GainValue {
  x1 = 0,
  x2 = 1,
  x4 = 2
};

const char i_GainNames[__CAM_GAINCOUNT][7] = { "Gain1x", "Gain2x", "Gain4x" }; 


enum LinlogValue {
  off = 0,
  low = 1,
  normal = 2,
  high = 3
};

const char i_LinlogNames[__CAM_LINLOGCOUNT][18] = { "Off", "LowCompression", "NormalCompression", "HighCompression" }; 



#endif  //  __CAMERA_TYPES
