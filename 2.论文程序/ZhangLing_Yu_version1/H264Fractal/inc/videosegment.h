#ifndef _VIDEOSEGMENT_H
#define _VIDEOSEGMENT_H
// #include <stdio.h>
// #include <math.h>

#include "parsetcommon.h"

#define TRUE1 1
#define  FALSE1 0




#define b(i) ((i%2)==0?(i/2):((i-1)/2))
#include "global.h"
void grayerosion(unsigned char *imputimage, int *cake, int cakeheight, int cakewidth, unsigned char *outputimage);//�Ҷȸ�ʴ�ӳ���
void grayinflation(unsigned char *imputimage, int *cake, int cakeheight, int cakewidth, unsigned char *outputimage);
void median(unsigned char *inputbmp,unsigned char *outputbmp);//��ֵ�˲��ӳ���
  
void VideoSegment(FILE*fp,unsigned char*binarybuf,int tol);
#endif