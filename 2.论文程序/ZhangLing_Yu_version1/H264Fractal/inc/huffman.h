#ifndef __HUFFMAN2__H_
#define __HUFFMAN2__H_

#include <math.h>
#include "defines_enc.h"
#include "global.h"
#define MAX_CL 12 
#define MASK   (1<<MAX_CL)-1

typedef struct Huffman_code
{
  int code;
  unsigned char bits;
}HCODE;
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// long **xy_p,**s_p,**o_p,**redual_ple,**redual_pru,**redual_ple_u,**redual_ple_v,**redual_pru_u,**redual_pru_v;// ÁÙÊ±±äÁ¿

extern int search_range;

HCODE ***Huff_ptr;
long   **xy_p, /* probability for CPM and NCIM */
       **s_p,
       **o_p;

void CreateHuffmanCodeBook (HCODE*,int,long *);
void HuffmanEncoder	(unsigned char *,int, int*, HCODE*, long*, int);

void GetInverseCodeBook	(int, HCODE*, int*);
void HuffmanDecoder (int*,HCODE*, int*, unsigned char *, long, int);

void huffman_Code();
#endif
