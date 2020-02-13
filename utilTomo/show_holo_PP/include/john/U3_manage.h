#ifndef __U3MANAGE__
#define __U3MANAGE__


#include "u3.h"
#include <unistd.h>

int configIO_DAC(HANDLE hDevice);
int checkI2CErrorcode(uint8 errorcode);


#endif /* __U3MANAGE__ */
