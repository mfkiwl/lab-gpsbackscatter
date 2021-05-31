#ifndef   ___USER__DAC__H__
#define  ___USER__DAC__H__

#include "user_type.h"

#define TRUE 1
#define FALSE 0

#define VOL_REF_MV 2500 // (2.5V)
#define CMD_LOAD_A 0X10
#define CMD_LOAD_B 0X20
#define CMD_BUFF_A 0X00
#define CMD_BUFF_B 0X04

#define VOL_CHA_MV 150
#define VOL_CHB_MV 180

#define SW_PN_MODE  // FUNC DEFINE

#define PN_CODE     (0X74) // 0111_0100
#define PN_CODE_LEN (7)    // 
#define PN_MASK_MSB (0X40) // 7 BIT



typedef enum {
    CHANNEL_A  = 0,
    CHANNEL_B  = 1,

    MODE_SW    = 0,
    MODE_HOLD  = 1,
    MODE_SW_PN = 2,

}ENUM_;


void user_dacInit(u8 RNG_range);
void user_dacTest(void);




#endif
