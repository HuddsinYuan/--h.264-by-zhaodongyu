#ifndef _CODE_H_
#define _CODE_H_
#include "global.h"

void Configure(char **);

void write_CompInfo(compressionInfo );
void Memoalloc();
void freeTrans(TRANS_NODE **trans);
void Memofree();
long pack(int option, long data, int bits, FILE * output);

void start_oneframe(char);/*start frame,给全局变量赋值*/
int encode_oneframe (int);
void decode_Oneframe();

void process_oneframe();
void dct_luma_fractal_uv();

void write_Istream();
void write_Codestream();

#endif