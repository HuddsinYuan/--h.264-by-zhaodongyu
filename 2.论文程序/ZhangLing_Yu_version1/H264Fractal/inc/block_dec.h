#ifndef _BLOCK_DEC_H_
#define _BLOCK_DEC_H_

#include "global.h"

extern long **xy_p,**s_p,**o_p;

void decode_one_macroblock(int,TRANS_NODE**,int);
void decode_block_rect(int b_x,int b_y,int Curblock_8,TRANS_NODE trans,int con,int mode,int reg,int depth);
void decode_block_8(int b_x,int b_y,int Curblock_8,TRANS_NODE trans,int con,int reg);
void decode_block_4(int b_x,int b_y,TRANS_NODE trans,int con,int reg);

#endif