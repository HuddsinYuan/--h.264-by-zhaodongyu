
/*!
 ***************************************************************************
 *
 * \file fmo.h
 *
 * \brief
 *    Support for Flexible Macroblock Ordering
 *
 * \date
 *    16 June 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifndef _FMO_H_
#define _FMO_H_

#define MAXSLICEGROUPIDS 8

int FmoInit(ImageParameters * img, pic_parameter_set_rbsp_t * pps, seq_parameter_set_rbsp_t * sps);
void FmoUninit ();
int FmoFinit (seq_parameter_set_rbsp_t * sps);
int FmoMB2SliceGroup (int mb);
int FmoGetFirstMBOfSliceGroup (int SliceGroupID);
int FmoGetFirstMacroblockInSlice (int SliceGroup);
int FmoGetNextMBNr (int CurrentMbNr);
int FmoGetLastCodedMBOfSliceGroup (int SliceGroupID);
int FmoStartPicture ();
int FmoEndPicture();
int FmoSliceGroupCompletelyCoded(int SliceGroupID);
void FmoSetLastMacroblockInSlice (int mb);

int FmoGetPreviousMBNr (int CurrentMbNr);

extern int *MBAmap; 
/////////////////////////////////////>


int FmoInit_dec (pic_parameter_set_rbsp_t* pps, seq_parameter_set_rbsp_t* sps);
int FmoFinit_dec ();

int FmoGetNumberOfSliceGroup();
int FmoGetLastMBOfPicture();
int FmoGetLastMBInSliceGroup(int SliceGroup);
int FmoGetSliceGroupId (int mb);
int FmoGetNextMBNr_dec (int CurrentMbNr);



/////////////////////////////////////<

#endif
