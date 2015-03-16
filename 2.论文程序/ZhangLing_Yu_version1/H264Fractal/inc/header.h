
/*!
 *************************************************************************************
 * \file header.h
 *
 * \brief
 *    Prototypes for header.c
 *************************************************************************************
 */

#ifndef _HEADER_H_
#define _HEADER_H_

int SliceHeader();
int Partition_BC_Header();

int  writeERPS(SyntaxElement *sym, DataPartition *partition);
// int  SequenceHeader(FILE *outf);
void write_terminating_bit (short);
/////////////////////////////////////>


int FirstPartOfSliceHeader();
int RestOfSliceHeader();

void dec_ref_pic_marking_dec(Bitstream *currStream);
void decode_poc(struct img_par *img);
// int  dumppoc   (struct img_par *img);



/////////////////////////////////////<
#endif

