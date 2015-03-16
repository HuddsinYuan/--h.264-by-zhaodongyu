#ifndef _BLOCK_ENC_H_
#define _BLOCK_ENC_H_

#include "global.h"

extern long **xy_p,**s_p,**o_p;
int S_PP;
unsigned int S_P;
void encode_one_macroblock(int CurMb,TRANS_NODE **trans,int con);
void p_encode_one_macroblock(int CurMb,TRANS_NODE **trans,int con);

void changeReferenceFrame(char);
double full_search(int,int,int,int,int,TRANS_NODE *);
double full_search_R(int,int,int,int,int,TRANS_NODE *);
double new_hexagon_block_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans);
double hexagon_block_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans);
double UMHEXIntegerPelBlockMotionSearch  (int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans);
int bound_chk(int m, int n, int cx, int cy,int,int,int);
void encode_block_8(int b_x,int b_y,int Curblock_8,TRANS_NODE *trans,int con);
void encode_block_4(int b_x,int b_y,TRANS_NODE *trans,int con);
int encode_block_rect(int b_x,int b_y,int Curblock_8,TRANS_NODE *trans,int con,int mode,int depth);
#endif