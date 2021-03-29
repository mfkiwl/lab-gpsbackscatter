#include "stdio.h"
#include "main.h"
#include "stm32l0xx_hal.h"
#include <math.h>
#include "user_dac.h"
#include "stm32l0xx_hal_tim.h"


#define SCLK_SET  {HAL_GPIO_WritePin(SCLK_GPIO_Port, SCLK_Pin, GPIO_PIN_SET);}
#define SCLK_CLR  {HAL_GPIO_WritePin(SCLK_GPIO_Port, SCLK_Pin, GPIO_PIN_RESET);}
#define DATA_SET  {HAL_GPIO_WritePin(DATA_GPIO_Port, DATA_Pin, GPIO_PIN_SET);}
#define DATA_CLR  {HAL_GPIO_WritePin(DATA_GPIO_Port, DATA_Pin, GPIO_PIN_RESET);}
#define LDAC_SET  {HAL_GPIO_WritePin(LDAC_GPIO_Port, LDAC_Pin, GPIO_PIN_SET);}
#define LDAC_CLR  {HAL_GPIO_WritePin(LDAC_GPIO_Port, LDAC_Pin, GPIO_PIN_RESET);}
#define LOAD_SET  {HAL_GPIO_WritePin(LOAD_GPIO_Port, LOAD_Pin, GPIO_PIN_SET);}
#define LOAD_CLR  {HAL_GPIO_WritePin(LOAD_GPIO_Port, LOAD_Pin, GPIO_PIN_RESET);}


volatile u8 RNG_range = RNG_RANGE_1X;
volatile float volRef_mV[4] = {2048.f, 2048.f, 2048.f, 2048.f};
volatile float volSet_mV[4] = {0.f, 0.f, 0.f, 0.f};
volatile u8 volDac[4] = {0};


#define VOL2CODE(vol, ref) (u8)round(256.f * vol / ((1 + RNG_range) * ref))

void user_dacSetVol(float voltage_mV, u8 channel)
{
    if(channel > 4 || voltage_mV < 0 || voltage_mV > (volRef_mV[channel] * RNG_range))
    {
        // wrong input
        return;
    }
    volSet_mV[channel] = voltage_mV;
    volDac[channel] = VOL2CODE(voltage_mV, volRef_mV[channel]);
    

}


// low pulse
void user_dacRefresh(void)
{
    LDAC_SET;
    HAL_Delay(1);
    LDAC_SET;
    HAL_Delay(1);
    LDAC_SET;
    HAL_Delay(1);
}

void user_dacSend(u8 channel)
{
    u16 data = 0;
    data = (data | ((channel << 1) | RNG_range)) << 8;
    data |= volDac[channel];

    SCLK_SET;
    LDAC_SET;    
    LOAD_SET;
    HAL_Delay(1);
    u8 i = 0;
    for(i=0;i<11;i++)
    {
        SCLK_SET;
        if(data & 0x70)
        {
            DATA_SET;
        }
        else
        {
            DATA_CLR;
        }
        HAL_Delay(1);
        data <<= 1;
        SCLK_CLR;
        HAL_Delay(1);
    }
    
    LOAD_CLR;
    HAL_Delay(1);
    LOAD_SET;

}




void user_dacInit(u8 RNG_range)
{






}

void user_dacTest(u8 RNG_range)
{
    
    




}


// 500 ms
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
    static u16 cnt_sw = 0;
    static u8 i = 0;

    // switch freq
    // T / 2
    #define Tsw_2  1000  // 2000 ms / 2 = 1000 ms
    if(++cnt_sw == (Tsw_2 / 500))
    {
//        HAL_GPIO_TogglePin(SCLK_GPIO_Port, SCLK_Pin);
        if(++i % 2 == 0)
        {
            // sw               
//            user_dacRefresh();
        }
        else
        {
            // mid sw
//            user_dacSend(CHANNEL_A);
//            user_dacSend(CHANNEL_B);
//            user_dacSend(CHANNEL_C);
//            user_dacSend(CHANNEL_D);
        }
        cnt_sw = 0;
    }
    

}


