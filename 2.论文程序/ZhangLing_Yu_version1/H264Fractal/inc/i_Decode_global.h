/*
*i_decode_global.h
*/
#ifndef _I_DECODE_GLOBAL_H
#define _I_DECODE_GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include "defines_dec.h"
#include "i_codec_types.h"

#define W1 2841 /* 2048*sqrt(2)*cos(1*pi/16) */
#define W2 2676 /* 2048*sqrt(2)*cos(2*pi/16) */
#define W3 2408 /* 2048*sqrt(2)*cos(3*pi/16) */
#define W5 1609 /* 2048*sqrt(2)*cos(5*pi/16) */
#define W6 1108 /* 2048*sqrt(2)*cos(6*pi/16) */
#define W7 565  /* 2048*sqrt(2)*cos(7*pi/16) */

//macro definition
#define WIDTHBYTES(i)    ((i+31)/32*4)//??????????

//define return value of function
#define FUNC_OK 0
#define FUNC_MEMORY_ERROR 1
#define FUNC_FILE_ERROR 2
#define FUNC_FORMAT_ERROR 3

//#define bound(a)   ((a) < 0 ? 0 : ((a)>255? 255 : a))

extern double aanscalefactor[8]; 
extern BYTE zigzag[64];
extern BYTE std_luminance_qt[64];
extern BYTE std_chrominance_qt[64];
extern BYTE std_dc_luminance_nrcodes[17];
extern BYTE std_dc_luminance_values[12];
extern BYTE std_dc_chrominance_nrcodes[17];
extern BYTE std_dc_chrominance_values[12];
extern BYTE std_ac_luminance_nrcodes[17];
extern BYTE std_ac_luminance_values[162];
extern BYTE std_ac_chrominance_nrcodes[17];
extern BYTE std_ac_chrominance_values[162];	

//////////////////////////////////////////////////
//Jpeg functions

int  InitTag();
void init_Table();

int  Decode();

int  DecodeBlock(BYTE ,BYTE);

int  HufBlock(BYTE dchufindex,BYTE achufindex);
int  DecodeElement();

void IQtIZzBlock();
BYTE ReadByte();
void Initialize_Fast_IDCT();
void Fast_IDCT(int * block);
void idctrow(int * blk);
void idctcol(int * blk);
//////////////////////////////////////////////////

int ImgWidth , ImgHeight;
BYTE i_quality;
BYTE *plane;

extern unsigned char i_quality;
extern int imageHeight;
extern int imageWidth;
extern FILE** fp_out_region_rec;
extern FILE*fp_out_all_rec,h264_out;
extern int num_regions;
unsigned char   *lp;
BYTE   qt_table[2][64];
short   comp_num;

BYTE      YDcIndex,YAcIndex,UVDcIndex,UVAcIndex;
BYTE   HufTabIndex;
short      *YQtTable,*UQtTable,*VQtTable;
BYTE   And[9];
short      code_pos_table[4][16],code_len_table[4][16];
unsigned short code_value_table[4][256];
unsigned short huf_max_value[4][16],huf_min_value[4][16];
short   BitPos,CurByte;
short   rrun,vvalue;

short   BlockBuffer[64];
short   coef;
int   interval;

int    Y[4*64],U[4*64],V[4*64];
DWORD      sizei,sizej;
short    restart;
static  long iclip[1024];
static  long *iclp;
void i_Frm_Decoder();

#endif
