/*
 *jglobals.h
 */
#ifndef _I_GLOBALS_H_
#define _I_GLOBALS_H_

#include <stdio.h>
#include <string.h>
#include "i_codec_types.h"
#include "i_defines.h"

typedef struct { 
	BYTE length;
	WORD value;}
 bitstring;//3bytes

/*
 *All extern variable
 */
extern double aanscalefactor[8]; //dct算法的系数

///////////////////////////huffman表参数
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
extern WORD mask[16];

///////////////////////////

struct DQTINFOTYPE { //量化表

		 BYTE QTYinfo;// = 0:  bit 0..3: number of QT = 0 (table for Y)
				  //       bit 4..7: precision of QT, 0 = 8 bit
		 BYTE Ytable[64];
		 BYTE QTCbinfo; // = 1 (quantization table for Cb,Cr}
		 BYTE Cbtable[64];
};//8bytes

struct DHTINFOTYPE { //huffman表

		 BYTE HTYDCinfo; // bit 0..3: number of HT (0..3), for Y =0
			  //bit 4  :type of HT, 0 = DC table,1 = AC table
				//bit 5..7: not used, must be 0
		 BYTE YDC_nrcodes[16]; //at index i = nr of codes with length i
		 BYTE YDC_values[12];
		 BYTE HTYACinfo; // = 0x10
		 BYTE YAC_nrcodes[16];
		 BYTE YAC_values[162];//we'll use the standard Huffman tables
		 BYTE HTCbDCinfo; // = 1
		 BYTE CbDC_nrcodes[16];
		 BYTE CbDC_values[12];
		 BYTE HTCbACinfo; //  = 0x11
		 BYTE CbAC_nrcodes[16];
		 BYTE CbAC_values[162];
};//16bytes

		 
typedef struct jpeg_encode{
	//variable

	int  Width;
	int  Height;
	int  quality;//[0,255], good to bad
	
	bitstring YDC_HT[12];    //HT表示huffman table
	bitstring CbDC_HT[12];
	bitstring YAC_HT[256];
	bitstring CbAC_HT[256];

	BYTE *category_alloc;//(BYTE *)malloc(65535*sizeof(BYTE));
	BYTE *category; 		//Here we'll keep the category of the numbers in range: -32767..32767
	bitstring *bitcode_alloc;//(bitstring *)malloc(65535*sizeof(bitstring));
    bitstring *bitcode; 	// their bitcoded representation

    float fdtbl_Y[64];
    float fdtbl_Cb[64]; 	//the same with the fdtbl_Cr[64]

    SWORD DU_DCT[64]; 	// Current DU (after DCT and quantization) which we'll zigzag
    SWORD DU[64]; 		//zigzag reordered DU which will be Huffman coded

	BYTE *Y_Buffer;//(DWORD *)(malloc(Width*Height*4));
	BYTE *Cb_Buffer;//(DWORD *)(malloc(Width*Height*4));
	BYTE *Cr_Buffer;//(DWORD *)(malloc(Width*Height*4));

	BYTE *Y_Buffer_plane;//(DWORD *)(malloc(Width*Height*4));
	BYTE *Cb_Buffer_plane;//(DWORD *)(malloc(Width*Height*4));
	BYTE *Cr_Buffer_plane;//(DWORD *)(malloc(Width*Height*4));

    SBYTE  Y_block[64];
    SBYTE  Cb_block[64];
    SBYTE  Cr_block[64];

	SBYTE  Y_block_plane[64];
    SBYTE  Cb_block_plane[64];
    SBYTE  Cr_block_plane[64];

	BYTE bytenew;
	SBYTE bytepos;
	//structure

	struct DQTINFOTYPE  DQTinfo;
	struct DHTINFOTYPE  DHTinfo; 

	//input structure
	//output file
	BYTE **fp_i_stream;//(BYTE *)(malloc(Width*Height*3));
	WORD *i_length;
}I_FRM_ENC_GLOBAL;//143bytes,即36DWORDS

void init_I_all(I_FRM_ENC_GLOBAL*);
void free_I_all(I_FRM_ENC_GLOBAL*);
BYTE i_Frm_Encoder(I_FRM_ENC_GLOBAL*);
#endif //_JGLOBALS_H_

