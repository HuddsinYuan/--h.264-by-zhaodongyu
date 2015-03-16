#ifndef _COMPUTE_H_
#define _COMPUTE_H_
#include "global.h"

void compute_range_Sum();
void compute_domain_Sum();
double compute_rdSum(int,int,int,int,int,int,int);
double compute_rms(int,int,int,int,double *,double*,int,int,int);
void classify(int block_x,int block_y,int m,int n,int obj,int,int,int con);
void partitionBlock(int,int,int,int,unsigned char,TRANS_NODE **,int);

#endif




















