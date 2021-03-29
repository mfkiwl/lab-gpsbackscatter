#ifndef   ___USER__DAC__H__
#define  ___USER__DAC__H__

#include "user_type.h"

typedef enum {
    CHANNEL_A = 0,
    CHANNEL_B = 1,
    CHANNEL_C = 2,
    CHANNEL_D = 3,

    RNG_RANGE_1X = 0,
    RNG_RANGE_2X = 1,

}ENUM_;


void user_dacInit(u8 RNG_range);
void user_dacTest(u8 RNG_range);




#endif
