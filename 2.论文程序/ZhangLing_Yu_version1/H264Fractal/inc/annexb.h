
/*!
 **************************************************************************************
 * \file
 *    annexb.h
 * \brief
 *    Byte stream operations support
 *    This code reflects JVT version xxx
 *  \date 7 December 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */


/////////////////////////////////////>
#ifndef _ANNEXB_H_
#define _ANNEXB_H_

/////////////////////////////////////<
#include <stdio.h>
#include "nalucommon.h"

/////////////////////////////////////>
extern int IsFirstByteStreamNALU;
extern int LastAccessUnitExists;
extern int NALUCount;
/////////////////////////////////////<

int WriteAnnexbNALU (NALU_t *n);
void CloseAnnexbFile();
void OpenAnnexbFile (char *Filename);

/////////////////////////////////////>
int  GetAnnexbNALU (NALU_t *nalu);
void OpenBitstreamFile (char *fn);
void CloseBitstreamFile();
void CheckZeroByteNonVCL(NALU_t *nalu, int * ret);
void CheckZeroByteVCL(NALU_t *nalu, int * ret);

#endif
/////////////////////////////////////<