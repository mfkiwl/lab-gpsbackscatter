#include "stdio.h"
#include "main.h"
#include "stm32l0xx_hal.h"
#include <math.h>
#include "user_dac.h"
#include "stm32l0xx_hal_tim.h"


#define SYNC_SET  {HAL_GPIO_WritePin(SYNC_GPIO_Port, SYNC_Pin, GPIO_PIN_SET);}
#define SYNC_CLR  {HAL_GPIO_WritePin(SYNC_GPIO_Port, SYNC_Pin, GPIO_PIN_RESET);}
#define SCLK_SET  {HAL_GPIO_WritePin(SCLK_GPIO_Port, SCLK_Pin, GPIO_PIN_SET);}
#define SCLK_CLR  {HAL_GPIO_WritePin(SCLK_GPIO_Port, SCLK_Pin, GPIO_PIN_RESET);}
#define DATA_SET  {HAL_GPIO_WritePin(DATA_GPIO_Port, DATA_Pin, GPIO_PIN_SET);}
#define DATA_CLR  {HAL_GPIO_WritePin(DATA_GPIO_Port, DATA_Pin, GPIO_PIN_RESET);}

//#define LDAC_SET  {HAL_GPIO_WritePin(LDAC_GPIO_Port, LDAC_Pin, GPIO_PIN_SET);}
//#define LDAC_CLR  {HAL_GPIO_WritePin(LDAC_GPIO_Port, LDAC_Pin, GPIO_PIN_RESET);}
//#define LOAD_SET  {HAL_GPIO_WritePin(LOAD_GPIO_Port, LOAD_Pin, GPIO_PIN_SET);}
//#define LOAD_CLR  {HAL_GPIO_WritePin(LOAD_GPIO_Port, LOAD_Pin, GPIO_PIN_RESET);}

#define CHANNEL_NUM 2

//volatile float volRef_mV[CHANNEL_NUM] = {2048.f, 2048.f, 2048.f, 2048.f};
volatile float volSet_mV[CHANNEL_NUM] = {0.f, 0.f};
// set = real + offset
volatile float volOffset_mV[CHANNEL_NUM] = {18.f, 12.5f};

volatile u8 cmdDac[CHANNEL_NUM] = {CMD_BUFF_A, CMD_LOAD_A|CMD_LOAD_B|CMD_BUFF_B};
volatile u16 volDac[CHANNEL_NUM] = {0};

volatile u8 flag_power = TRUE;



#define VOL2CODE(vol, ref) (u16)round(65536.f * vol /  ref)



void user_delay(u16 t)
{
    while(t--){}
}



void user_dacSetVol(float voltage_mV, u8 channel)
{
    u16 code = 0;
    if(channel > 1)
    {
        // wrong input
        return;
    }

    if(voltage_mV > 2*VOL_REF_MV)
    {
        voltage_mV = 2*VOL_REF_MV;
    }
    else if (voltage_mV < -2*VOL_REF_MV)
    {
        voltage_mV = -2*VOL_REF_MV;
    }

    // vol code
    volSet_mV[channel] = voltage_mV;
    voltage_mV -= volOffset_mV[channel];
    // Vo = 2(2Vi - Vref)
    voltage_mV = (voltage_mV + 2*VOL_REF_MV) / 4;
    code = VOL2CODE(voltage_mV, VOL_REF_MV);
    volDac[channel] = code;

    // cmd
    // default setting
}


// low pulse
void user_dacRefresh(void)
{

}

void user_dacSend(u8 channel)
{
    u32 data = 0;
    data = (cmdDac[channel] << 16) | volDac[channel];

    SYNC_SET;
    

    user_delay(5);
    SYNC_CLR;
    SCLK_SET;
    u8 i = 0;
    for(i=0;i<24;i++)
    {
        SCLK_SET;
        if(data & 0x0800000) // 24 bit MSB
        {
            DATA_SET;
        }
        else
        {
            DATA_CLR;
        }
        user_delay(5);
        data <<= 1;
        SCLK_CLR;
        user_delay(5);
    }
    SYNC_SET;

}




void user_dacInit(u8 RNG_range)
{






}

void user_dacTest()
{
    
    user_dacSetVol(0, CHANNEL_A);
    user_dacSetVol(0, CHANNEL_B);


}


// 500 ms
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
    static u16 cnt_sw = 0;
    static u8 i = 0;

    // switch freq
    // T / 2
    #define Tsw_2  1000  // 2 sec / 2 = 1000 ms
    if(++cnt_sw == (Tsw_2 / 500))
    {
        u8 mode = HAL_GPIO_ReadPin(MODE_GPIO_Port, MODE_Pin)? MODE_HOLD : MODE_SW;
        if(++i % 2 == 0)
        {
            // sw               
            HAL_GPIO_TogglePin(LED_GPIO_Port, LED_Pin);
            user_dacSend(CHANNEL_A);
            user_dacSend(CHANNEL_B);
        }
        else
        {
            // mid sw
            if((flag_power && mode == MODE_SW) || mode == MODE_HOLD)
            {
                user_dacSetVol(VOL_CHA_MV, CHANNEL_A);
                user_dacSetVol(VOL_CHB_MV, CHANNEL_B);
                flag_power = FALSE;
            }
            else
            {
                user_dacSetVol(0, CHANNEL_A);
                user_dacSetVol(0, CHANNEL_B);
                flag_power = TRUE;
            }
        }
        cnt_sw = 0;
    }
    

}


