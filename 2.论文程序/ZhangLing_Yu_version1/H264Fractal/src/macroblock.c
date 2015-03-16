/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////


/*!
 *************************************************************************************
 * \file macroblock.c
 *
 * \brief
 *    Process one macroblock
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Lang鴜               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Ragip Kurceren                  <ragip.kurceren@nokia.com>
 *************************************************************************************
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "elements.h"
#include "macroblock.h"
#include "refbuf.h"
#include "fmo.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"
#include "ratectl.h"              // head file for rate control
#include "cabac.h"
#include "i_global.h"
#include "i_defines.h"

#include "biaridecod.h"
#include "global.h"
#include "mbuffer.h"


#include "stdio.h"
#include "defines_enc.h"

//Rate control
int predict_error,dq;
extern int DELTA_QP,DELTA_QP2;
extern int QP,QP2;



/////////////////////////////////////>
#if TRACE
#define TRACE_STRING(s) strncpy(currSE.tracestring, s, TRACESTRING_SIZE)
#else
#define TRACE_STRING(s) // do nothing
#endif

extern ColocatedParams *Co_located;



/////////////////////////////////////<


/*!
************************************************************************
* \brief
*    updates the coordinates for the next macroblock to be processed
*
* \param mb_addr
*    MB address in scan order
************************************************************************
*/

void set_MB_parameters (int mb_addr)
{
	img->current_mb_nr = mb_addr;
	
	get_mb_block_pos(mb_addr, &img->mb_x, &img->mb_y);//宏块的光栅扫描坐标转换成二维x,y坐标,赋给img->mb_x，img->mb_y（以宏块为单位）
	
	img->block_x = img->mb_x << 2;//当前宏块的x坐标（以8*8块为单位）
	img->block_y = img->mb_y << 2;//当前宏块的y坐标（以8*8块为单位）
	
	img->pix_x   = img->block_x << 2;//当前宏块的x坐标（以像素为单位）
	img->pix_y   = img->block_y << 2;//当前宏块的y坐标（以8*8块为单位）
	
	img->opix_x   = img->pix_x;//当前宏块的x坐标（以像素为单位）
	
	if (img->MbaffFrameFlag)//no
	{
		if (img->mb_data[mb_addr].mb_field)
		{
			
			imgY_org  = (mb_addr % 2) ? imgY_org_bot  : imgY_org_top;
			imgUV_org = (mb_addr % 2) ? imgUV_org_bot : imgUV_org_top;
			img->opix_y   = (img->mb_y >> 1 ) << 4;
		}
		else
		{
			imgY_org  = imgY_org_frm;
			imgUV_org = imgUV_org_frm;
			img->opix_y   = img->block_y << 2;
		}
	}
	else
	{
		img->opix_y   = img->block_y << 2;////当前宏块的y坐标（以1/4像素为单位）
	}
	
	img->pix_c_x = img->pix_x >> 1;//色度宏块的x坐标（以像素为单位）
	img->pix_c_y = img->pix_y >> 1;//色度宏块的y坐标（以像素为单位）
	
	img->opix_c_x = img->opix_x >> 1;//色度宏块的x坐标（以像素为单位）
	img->opix_c_y = img->opix_y >> 1;//色度宏块的x坐标（以1/4像素为单位）
	
	//  printf ("set_MB_parameters: mb %d,  mb_x %d,  mb_y %d\n", mb_addr, img->mb_x, img->mb_y);
}


int clip1a(int a)
{
	return ((a)>255?255:((a)<0?0:(a)));
}

/*!
 ************************************************************************
 * \brief
 *    updates the coordinates and statistics parameter for the
 *    next macroblock
 ************************************************************************
 */
void proceed2nextMacroblock()
{
#if TRACE
  int i;
  int use_bitstream_backing = (input->slice_mode == FIXED_RATE || input->slice_mode == CALLBACK);
#endif
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int*        bitCount = currMB->bitcounter;

#if TRACE

  if (p_trace)
  {
    fprintf(p_trace, "\n*********** Pic: %i (I/P) MB: %i Slice: %i **********\n\n", frame_no, img->current_mb_nr, img->current_slice_nr);
    if(use_bitstream_backing)
      fprintf(p_trace, "\n*********** Pic: %i (I/P) MB: %i Slice: %i **********\n\n", frame_no, img->current_mb_nr, img->current_slice_nr);
   // Write out the tracestring for each symbol
    for (i=0; i<currMB->currSEnr; i++)
      trace2out(&(img->MB_SyntaxElements[i]));
  }
#endif

  // Update the statistics
  stat->bit_use_mb_type    [img->type]  += bitCount[BITS_MB_MODE];
  stat->bit_use_coeffY     [img->type]  += bitCount[BITS_COEFF_Y_MB] ;
  stat->tmp_bit_use_cbp    [img->type]  += bitCount[BITS_CBP_MB];
  stat->bit_use_coeffC     [img->type]  += bitCount[BITS_COEFF_UV_MB];
  stat->bit_use_delta_quant[img->type]  += bitCount[BITS_DELTA_QUANT_MB];

  ++stat->mode_use[img->type][currMB->mb_type];
  stat->bit_use_mode[img->type][currMB->mb_type]+= bitCount[BITS_INTER_MB];

  // Statistics
  if ((img->type == P_SLICE)||(img->type==SP_SLICE) )
  {
    ++stat->quant0;
    stat->quant1 += currMB->qp;      // to find average quant for inter frames
  }
}
/*!
************************************************************************
* \brief
*    initializes the current macroblock
************************************************************************
*/
void start_macroblock(int mb_addr, int mb_field)
{
	int i,j,k,l;
	int use_bitstream_backing = (input->slice_mode == FIXED_RATE || input->slice_mode == CALLBACK1);//0
	Macroblock *currMB = &img->mb_data[mb_addr];//当前宏块数居
	Slice *curr_slice = img->currentSlice;//当前条带
	DataPartition *dataPart;
	Bitstream *currStream;
	EncodingEnvironmentPtr eep;
	
	currMB->mb_field = mb_field;//0
	
	enc_picture->mb_field[mb_addr] = mb_field;//0
	
	set_MB_parameters (mb_addr);
	
	if(use_bitstream_backing)//no
	{
		// Keep the current state of the bitstreams
		if(!img->cod_counter)
			for (i=0; i<curr_slice->max_part_nr; i++)
			{
				dataPart = &(curr_slice->partArr[i]);
				currStream = dataPart->bitstream;
				currStream->stored_bits_to_go   = currStream->bits_to_go;
				currStream->stored_byte_pos   = currStream->byte_pos;
				currStream->stored_byte_buf   = currStream->byte_buf;
				
				if (input->symbol_mode ==CABAC)
				{
					eep = &(dataPart->ee_cabac);
					eep->ElowS            = eep->Elow;
					eep->ErangeS           = eep->Erange;
					eep->EbufferS         = eep->Ebuffer;
					eep->Ebits_to_goS     = eep->Ebits_to_go;
					eep->Ebits_to_followS = eep->Ebits_to_follow;
					eep->EcodestrmS       = eep->Ecodestrm;
					eep->Ecodestrm_lenS   = eep->Ecodestrm_len;
					eep->CS               = eep->C;
					eep->BS               = eep->B;
					eep->ES               = eep->E;
				}
			}
	}
	
	// Save the slice number of this macroblock. When the macroblock below
	// is coded it will use this to decide if prediction for above is possible
	currMB->slice_nr = img->current_slice_nr;//0
	
	// Initialize delta qp change from last macroblock. Feature may be used for future rate control
	// Rate control
	currMB->qpsp       = img->qpsp;//0
	if(input->RCEnable)//no
	{
		if (img->current_mb_nr==0)
		{
			currMB->prev_qp = img->qp;
			currMB->prev_delta_qp = 0;
		}
		else
		{    
			currMB->prev_qp = img->mb_data[img->current_mb_nr-1].qp;
			currMB->prev_delta_qp = img->mb_data[img->current_mb_nr-1].delta_qp;
		}
		/*frame layer rate control*/
		if(input->basicunit==img->Frame_Total_Number_MB)
		{
			currMB->delta_qp = 0;
			currMB->qp       = img->qp;
		}
		/*basic unit layer rate control*/
		else
		{
			/*each I or B frame has only one QP*/
			if((img->type==I_SLICE)||(img->type==B_SLICE))
			{
				currMB->delta_qp = 0;
				currMB->qp       = img->qp;
			}
			else if(img->type==P_SLICE)
			{
				
				if (!img->write_macroblock) //write macroblock
				{
					if (!currMB->mb_field)  //frame macroblock
					{
						if (img->current_mb_nr == 0) //first macroblock
						{
							// Initialize delta qp change from last macroblock. Feature may be used for future rate control
							currMB->delta_qp = 0;
							currMB->qp       = img->qp;
							DELTA_QP = DELTA_QP2 = currMB->delta_qp;
							QP = QP2 = currMB->qp;
						}
						else
						{
							if (!((input->MbInterlace) && img->bot_MB)) //top macroblock
							{
								if (img->mb_data[img->current_mb_nr-1].prev_cbp == 1)
								{
									currMB->delta_qp = 0;
									currMB->qp       = img->qp;
								}
								else
								{
									currMB->qp = img->mb_data[img->current_mb_nr-1].prev_qp;
									currMB->delta_qp = currMB->qp - img->mb_data[img->current_mb_nr-1].qp;
									img->qp = currMB->qp;
								}
								DELTA_QP = DELTA_QP2 = currMB->delta_qp;
								QP = QP2 = currMB->qp;
							}
							else //bottom macroblock
							{
								// Initialize delta qp change from last macroblock. Feature may be used for future rate control
								currMB->delta_qp = 0;
								currMB->qp       = img->qp;       // needed in loop filter (even if constant QP is used)
							}
						}
					}
					else  //field macroblock
					{
						if (!img->bot_MB) //top macroblock 
						{
							currMB->delta_qp = DELTA_QP2;
							currMB->qp   = img->qp    = QP2;
						}
						else//bottom macroblock
						{
							currMB->qp = img->qp;
							currMB->delta_qp = 0;
						}
						
					}
					
				}
				else 
				{
					if (!img->bot_MB) //write top macroblock
					{
						if (img->write_macroblock_frame)
						{
							currMB->delta_qp = DELTA_QP;
							img->qp = currMB->qp = QP;
						}
						else
						{
							currMB->delta_qp = DELTA_QP2;
							img->qp = currMB->qp = QP2;
						}
					}
					else //write bottom macroblock
					{
						currMB->delta_qp = 0;
						currMB->qp = img->qp;
					}
				}
				
				/*compute the quantization parameter for each basic unit of P frame*/
				
				if(!((input->MbInterlace)&&img->bot_MB))
				{
					if(!currMB->mb_field)
					{
						
						if((img->NumberofCodedMacroBlocks>0)\
							&&(img->NumberofCodedMacroBlocks%img->BasicUnit==0))
						{
							
							/*frame coding*/
							if(active_sps->frame_mbs_only_flag)
							{
								updateRCModel();
								img->BasicUnitQP=updateQuantizationParameter(img->TopFieldFlag);
							}
							/*adaptive field/frame coding*/
							else if((input->PicInterlace==ADAPTIVE_CODING)&&(!input->MbInterlace)&&(img->IFLAG==0))
							{
								updateRCModel();
								img->BasicUnitQP=updateQuantizationParameter(img->TopFieldFlag);
							}
							/*field coding*/
							else if((input->PicInterlace==FIELD_CODING)&&(!input->MbInterlace)&&(img->IFLAG==0))
							{
								updateRCModel();
								img->BasicUnitQP=updateQuantizationParameter(img->TopFieldFlag);
							}
							/*mb adaptive f/f coding, field coding*/
							else if((input->MbInterlace)&&(img->IFLAG==0)&&(img->FieldControl==1))
							{
								updateRCModel();
								img->BasicUnitQP=updateQuantizationParameter(img->TopFieldFlag);
							}
							/*mb adaptive f/f coding, frame coding*/
							else if((input->MbInterlace)&&(img->IFLAG==0)&&(img->FieldControl==0))
							{
								updateRCModel();
								img->BasicUnitQP=updateQuantizationParameter(img->TopFieldFlag);
							} 
							
							
						}
						
						
						if(img->current_mb_nr==0)
							img->BasicUnitQP=img->qp;
						
						currMB->predict_qp=img->BasicUnitQP;
						
						if(currMB->predict_qp>currMB->qp+25)
							currMB->predict_qp=currMB->qp+25;
						else if(currMB->predict_qp<currMB->qp-26)
							currMB->predict_qp=currMB->qp-26; 
						
						currMB->prev_qp = currMB->predict_qp;
						
						dq = currMB->delta_qp + currMB->predict_qp-currMB->qp;
						if(dq < -26) 
						{
							dq = -26;
							predict_error = dq-currMB->delta_qp;
							img->qp = img->qp+predict_error;
							currMB->delta_qp = -26;
						}
						else if(dq > 25)
						{
							dq = 25;
							predict_error = dq - currMB->delta_qp;
							img->qp = img->qp + predict_error;
							currMB->delta_qp = 25;
						}
						else
						{
							currMB->delta_qp = dq;
							predict_error=currMB->predict_qp-currMB->qp;
							img->qp = currMB->predict_qp;
						}
						currMB->qp =  img->qp;
						if (input->MbInterlace)
						{
							DELTA_QP = DELTA_QP2 = currMB->delta_qp;
							QP = QP2     = currMB->qp;
							DELTA_QP2 = currMB->delta_qp;
						}
						currMB->predict_error=predict_error;
					}
					else
						predict_error=currMB->predict_error;
				}
				else
					currMB->prev_qp=img->qp;
       }
    }   
  }
  else//yes
  {
	  Slice* currSlice = img->currentSlice;
	  
	  int prev_mb = FmoGetPreviousMBNr(img->current_mb_nr);//-1
	  if (prev_mb>-1)
	  {
		  currMB->prev_qp = img->mb_data[prev_mb].qp;
		  currMB->prev_delta_qp = img->mb_data[prev_mb].delta_qp;
	  }
	  else//1
	  {
		  currMB->prev_qp = currSlice->qp;//输入的QP
		  currMB->prev_delta_qp = 0;
	  }
	  
	  currMB->qp       = currSlice->qp ;//输入的QP
	  
	  currMB->delta_qp = currMB->qp - currMB->prev_qp;//0
	  DELTA_QP = DELTA_QP2 = currMB->delta_qp;//0
	  QP = QP2 = currMB->qp;//输入的QP
	  
  }
  // Initialize counter for MB symbols 对宏块编码符号计数器进行初始化
  currMB->currSEnr=0;
  
  // loop filter parameter环路滤波参数
  if (active_pps->deblocking_filter_control_present_flag)//no
  {
	  currMB->LFDisableIdc    = img->LFDisableIdc;
	  currMB->LFAlphaC0Offset = img->LFAlphaC0Offset;
	  currMB->LFBetaOffset    = img->LFBetaOffset;
  }
  else//yes
  {
	  currMB->LFDisableIdc    = 0;
	  currMB->LFAlphaC0Offset = 0;
	  currMB->LFBetaOffset    = 0;
  }
  
  // If MB is next to a slice boundary, mark neighboring blocks unavailable for prediction
  CheckAvailabilityOfNeighbors();//标记相邻块的可用性
  
  if (input->symbol_mode == CABAC)
	  CheckAvailabilityOfNeighborsCABAC();
  
  // Reset vectors before doing motion search in motion_search(). //运动估计之前将MV初始化为(0,0)，l包含list0,list1  
  for (l=0; l<2; l++)//遍历list0,list1
  {
	  for (j=0; j < BLOCK_MULTIPLE; j++)//BLOCK_MULTIPLE=4,所以我猜测是以4x4块为基本单位初始化的
		  for (i=0; i < BLOCK_MULTIPLE; i++)
			  for (k=0; k < 2; k++)//代表x分量和y分量
				  enc_picture->mv[l][img->block_x+i][img->block_y+j][k]=0;
  }
  
  //initialize reference index 初始化listo,list1的参考索引为-1，我猜测是以4x4块为基本单位
  for (j=0; j < BLOCK_MULTIPLE; j++)
  {
	  for (i=0; i < BLOCK_MULTIPLE; i++)
		  for (l=0; l<2; l++)//遍历list0,list1
		  {
			  enc_picture->ref_idx[ l ][ img->block_x+i ][ img->block_y+j ] = -1;
			  enc_picture->ref_pic_id[ l ][ img->block_x+i ][ img->block_y+j ] = -1;
		  }
  }
  
  // Reset syntax element entries in MB struct 重置宏块的输入句法元素
  currMB->mb_type   = 0;
  currMB->cbp_blk   = 0;
  currMB->cbp       = 0;
  
  for (l=0; l < 2; l++) //初始化list0,list1的MVD为0，我猜测是以4x4块为基本单位
	  for (j=0; j < BLOCK_MULTIPLE; j++)
		  for (i=0; i < BLOCK_MULTIPLE; i++)
			  for (k=0; k < 2; k++)
				  currMB->mvd[l][j][i][k] = 0;
			  
  currMB->cbp_bits   = 0;
  currMB->c_ipred_mode = DC_PRED_8; //预测模式初始化为DC_PRED_8
			  
  for (i=0; i < (BLOCK_MULTIPLE*BLOCK_MULTIPLE); i++)
    currMB->intra_pred_modes[i] = DC_PRED;
			  
  //initialize the whole MB as INTRA coded
  //Blocks ar set to notINTRA in write_one_macroblock
  if (input->UseConstrainedIntraPred)//no
  {
	  img->intra_block[img->current_mb_nr] = 1;
  }
  // store filtering parameters for this MB; For now, we are using the //滤波初始化
  // same offset throughout the sequence
  currMB->lf_disable = img->LFDisableIdc;//0
  currMB->lf_alpha_c0_offset = img->LFAlphaC0Offset;//0
  currMB->lf_beta_offset = img->LFBetaOffset;//0
			  
			  
  // Initialize bitcounters for this macroblock
  if(img->current_mb_nr == 0)//yes // No slice header to account for
  {
	  currMB->bitcounter[BITS_HEADER] = 0;
  }
  else if (currMB->slice_nr == img->mb_data[img->current_mb_nr-1].slice_nr) // current MB belongs to the
  // same slice as the last MB
  {
	  currMB->bitcounter[BITS_HEADER] = 0;//BITS_HEADER=0
  }
			  
  currMB->bitcounter[BITS_MB_MODE] = 0;//BITS_MB_MODE=2
  currMB->bitcounter[BITS_COEFF_Y_MB] = 0;//BITS_COEFF_Y_MB=5
  currMB->bitcounter[BITS_INTER_MB] = 0;//BITS_INTER_MB=3
  currMB->bitcounter[BITS_CBP_MB] = 0;//BITS_CBP_MB=4
  currMB->bitcounter[BITS_DELTA_QUANT_MB] = 0;//BITS_DELTA_QUANT_MB=7
  currMB->bitcounter[BITS_COEFF_UV_MB] = 0;//BITS_COEFF_UV_MB=6
			  
#ifdef _FAST_FULL_ME_
  if(!input->FMEnable)
    ResetFastFullIntegerSearch ();
#endif
}

/*!
 ************************************************************************
 * \brief
 *    terminates processing of the current macroblock depending
 *    on the chosen slice mode
 ************************************************************************
 */
void terminate_macroblock(Boolean *end_of_slice, Boolean *recode_macroblock)
{
  int i;
  Slice *currSlice = img->currentSlice;
  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int *partMap = assignSE2partition[input->partition_mode];
  DataPartition *dataPart;
  Bitstream *currStream;
  int rlc_bits=0;
  EncodingEnvironmentPtr eep;
  int use_bitstream_backing = (input->slice_mode == FIXED_RATE || input->slice_mode == CALLBACK1);
  int new_slice;
  static int skip = FALSE1;

	 
  // if previous mb in the same slice group has different slice number as the current, it's the
  // the start of new slice
  new_slice=0;
  if ( (img->current_mb_nr==0) || (FmoGetPreviousMBNr(img->current_mb_nr)<0) )
    new_slice=1;
  else if( img->mb_data[FmoGetPreviousMBNr(img->current_mb_nr)].slice_nr != img->current_slice_nr )
    new_slice=1;
	
  
  *recode_macroblock=FALSE1;

  switch(input->slice_mode)
  {
  case NO_SLICES:
    currSlice->num_mb++;
    *recode_macroblock = FALSE1;
    if ((currSlice->num_mb) == img->total_number_mb) // maximum number of MBs reached
      *end_of_slice = TRUE1;

    // if it's end of current slice group, slice ends too
    *end_of_slice |= (img->current_mb_nr == FmoGetLastCodedMBOfSliceGroup (FmoMB2SliceGroup (img->current_mb_nr)));
		
    break;
  case FIXED_MB:
    // For slice mode one, check if a new slice boundary follows
    currSlice->num_mb++;
    *recode_macroblock = FALSE1;
    //! Check end-of-slice group condition first
    *end_of_slice = (img->current_mb_nr == FmoGetLastCodedMBOfSliceGroup (FmoMB2SliceGroup (img->current_mb_nr)));
    //! Now check maximum # of MBs in slice
    *end_of_slice |= (currSlice->num_mb >= input->slice_argument);
		
    break;

    // For slice modes two and three, check if coding of this macroblock
    // resulted in too many bits for this slice. If so, indicate slice
    // boundary before this macroblock and code the macroblock again
  case FIXED_RATE:
     // in case of skip MBs check if there is a slice boundary
     // only for UVLC (img->cod_counter is always 0 in case of CABAC)
     if(img->cod_counter)
     {
       // write out the skip MBs to know how many bits we need for the RLC
       currSE->value1 = img->cod_counter;
       currSE->value2 = 0;
       currSE->mapping = ue_linfo;
       currSE->type = SE_MBTYPE;
       dataPart = &(currSlice->partArr[partMap[currSE->type]]);

       dataPart->writeSyntaxElement(  currSE, dataPart);
       rlc_bits=currSE->len;

       currStream = dataPart->bitstream;
       // save the bitstream as it would be if we write the skip MBs
       currStream->bits_to_go_skip  = currStream->bits_to_go;
       currStream->byte_pos_skip    = currStream->byte_pos;
       currStream->byte_buf_skip    = currStream->byte_buf;
       // restore the bitstream
       currStream->bits_to_go = currStream->stored_bits_to_go;
       currStream->byte_pos = currStream->stored_byte_pos;
       currStream->byte_buf = currStream->stored_byte_buf;
       skip = TRUE1;
     }
     //! Check if the last coded macroblock fits into the size of the slice
     //! But only if this is not the first macroblock of this slice
     if (!new_slice)
     {
       if(slice_too_big(rlc_bits))
       {
         *recode_macroblock = TRUE1;
         *end_of_slice = TRUE1;
       }
       else if(!img->cod_counter)
         skip = FALSE1;
     }
     // maximum number of MBs
		 
     // check if current slice group is finished
     if ((*recode_macroblock == FALSE1) && (img->current_mb_nr == FmoGetLastCodedMBOfSliceGroup (FmoMB2SliceGroup (img->current_mb_nr)))) 
     {
       *end_of_slice = TRUE1;
       if(!img->cod_counter)
         skip = FALSE1;
     }
   
     //! (first MB OR first MB in a slice) AND bigger that maximum size of slice
     if (new_slice && slice_too_big(rlc_bits))
     {
       *end_of_slice = TRUE1;
       if(!img->cod_counter)
         skip = FALSE1;
     }
     if (!*recode_macroblock)
       currSlice->num_mb++;
     break;

  case  CALLBACK1:
    if (img->current_mb_nr > 0 && !new_slice)
    {
      if (currSlice->slice_too_big(rlc_bits))
      {
        *recode_macroblock = TRUE1;
        *end_of_slice = TRUE1;
      }
    }
		
    if ( (*recode_macroblock == FALSE1) && (img->current_mb_nr == FmoGetLastCodedMBOfSliceGroup (FmoMB2SliceGroup (img->current_mb_nr)))) 
      *end_of_slice = TRUE1;
    break;

  default:
    printf(errortext, ET_SIZE, "Slice Mode %d not supported", input->slice_mode);
    error(errortext, 600);
  }

  if(*recode_macroblock == TRUE1)
  {
    // Restore everything
    for (i=0; i<currSlice->max_part_nr; i++)
    {
      dataPart = &(currSlice->partArr[i]);
      currStream = dataPart->bitstream;
      currStream->bits_to_go = currStream->stored_bits_to_go;
      currStream->byte_pos  = currStream->stored_byte_pos;
      currStream->byte_buf  = currStream->stored_byte_buf;
      if (input->symbol_mode == CABAC)
      {
        eep = &(dataPart->ee_cabac);
        eep->Elow            = eep->ElowS;
        eep->Erange           = eep->ErangeS;
        eep->Ebuffer         = eep->EbufferS;
        eep->Ebits_to_go     = eep->Ebits_to_goS;
        eep->Ebits_to_follow = eep->Ebits_to_followS;
        eep->Ecodestrm       = eep->EcodestrmS;
        eep->Ecodestrm_len   = eep->Ecodestrm_lenS;
        eep->C               = eep->CS;
        eep->B               = eep->BS;
        eep->E               = eep->ES;       
      }
    }
  }

  if(*end_of_slice == TRUE1  && skip == TRUE1) //! TO 4.11.2001 Skip MBs at the end of this slice
  { 
    //! only for Slice Mode 2 or 3
    // If we still have to write the skip, let's do it!
    if(img->cod_counter && *recode_macroblock == TRUE1) //! MB that did not fit in this slice
    { 
      // If recoding is true and we have had skip, 
      // we have to reduce the counter in case of recoding
      img->cod_counter--;
      if(img->cod_counter)
      {
        currSE->value1 = img->cod_counter;
        currSE->value2 = 0;
        currSE->mapping = ue_linfo;
        currSE->type = SE_MBTYPE;
        dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        dataPart->writeSyntaxElement(  currSE, dataPart);
        rlc_bits=currSE->len;
        currMB->bitcounter[BITS_MB_MODE]+=rlc_bits;
        img->cod_counter = 0;
      }
    }
    else //! MB that did not fit in this slice anymore is not a Skip MB
    {
      dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);       
      currStream = dataPart->bitstream;
        // update the bitstream
      currStream->bits_to_go = currStream->bits_to_go_skip;
      currStream->byte_pos  = currStream->byte_pos_skip;
      currStream->byte_buf  = currStream->byte_buf_skip;

      // update the statistics
      img->cod_counter = 0;
      skip = FALSE1;
    }
  }
  
  //! TO 4.11.2001 Skip MBs at the end of this slice for Slice Mode 0 or 1
  if(*end_of_slice == TRUE1 && img->cod_counter && !use_bitstream_backing)
  {
    currSE->value1 = img->cod_counter;
    currSE->value2 = 0;
    currSE->mapping = ue_linfo;
    currSE->type = SE_MBTYPE;
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    dataPart->writeSyntaxElement(  currSE, dataPart);
    rlc_bits=currSE->len;
    currMB->bitcounter[BITS_MB_MODE]+=rlc_bits;
    img->cod_counter = 0;
  }
}
/*!
 *****************************************************************************
 *
 * \brief 
 *    For Slice Mode 2: Checks if one partition of one slice exceeds the 
 *    allowed size
 * 
 * \return
 *    FALSE if all Partitions of this slice are smaller than the allowed size
 *    TRUE is at least one Partition exceeds the limit
 *
 * \par Side effects
 *    none
 *
 * \date
 *    4 November 2001
 *
 * \author
 *    Tobias Oelbaum      drehvial@gmx.net
 *****************************************************************************/
 
 int slice_too_big(int rlc_bits)
 {
   Slice *currSlice = img->currentSlice;
   DataPartition *dataPart;
   Bitstream *currStream;
   EncodingEnvironmentPtr eep;
   int i;
   int size_in_bytes;
  
   //! UVLC
   if (input->symbol_mode == UVLC)
   {
     for (i=0; i<currSlice->max_part_nr; i++)
     {
       dataPart = &(currSlice->partArr[i]);
       currStream = dataPart->bitstream;
       size_in_bytes = currStream->byte_pos /*- currStream->tmp_byte_pos*/;

       if (currStream->bits_to_go < 8)
         size_in_bytes++;
       if (currStream->bits_to_go < rlc_bits)
         size_in_bytes++;
       if(size_in_bytes > input->slice_argument)
         return TRUE1;
     }
   }
    
   //! CABAC
   if (input->symbol_mode ==CABAC)
   {
     for (i=0; i<currSlice->max_part_nr; i++)
     {
        dataPart= &(currSlice->partArr[i]);
        eep = &(dataPart->ee_cabac);
      
       if( arienco_bits_written(eep) > (input->slice_argument*8))
          return TRUE1;
     }
   }
   return FALSE1;
 }


/*!
 ************************************************************************
 * \brief
 *    Predict one component of a 4x4 Luma block
 ************************************************************************
 */
void
OneComponentLumaPrediction4x4 (int*   mpred,      //  --> array of prediction values (row by row)
                               int    pic_pix_x,  // <--  absolute horizontal coordinate of 4x4 block
                               int    pic_pix_y,  // <--  absolute vertical   coordinate of 4x4 block
                               int*   mv,         // <--  motion vector
                               int    ref,        // <--  reference frame 
                               StorablePicture **list)
{
  pel_t** ref_pic;
  int     pix_add = 4;
  int     j0      = (pic_pix_y << 2) + mv[1], j1=j0+pix_add, j2=j1+pix_add, j3=j2+pix_add;
  int     i0      = (pic_pix_x << 2) + mv[0], i1=i0+pix_add, i2=i1+pix_add, i3=i2+pix_add;
  
  pel_t (*get_pel) (pel_t**, int, int, int, int) = UMVPelY_14;

  int img_width =list[ref]->size_x;
  int img_height=list[ref]->size_y;

  ref_pic   = list[ref]->imgY_ups;
  
  *mpred++ = get_pel (ref_pic, j0, i0, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j0, i1, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j0, i2, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j0, i3, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j1, i0, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j1, i1, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j1, i2, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j1, i3, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j2, i0, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j2, i1, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j2, i2, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j2, i3, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j3, i0, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j3, i1, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j3, i2, img_height, img_width);
  *mpred++ = get_pel (ref_pic, j3, i3, img_height, img_width);

}


/*!
 ************************************************************************
 * \brief
 *    Predict one 4x4 Luma block
 ************************************************************************
 */
void LumaPrediction4x4 (int  block_x,    // <--  relative horizontal block coordinate of 4x4 block
                        int  block_y,    // <--  relative vertical   block coordinate of 4x4 block
                        int  p_dir,      // <--  prediction direction (0=forward, 1=backward, 2=bidir)
                        int  fw_mode,    // <--  forward  prediction mode (1-7, 0=DIRECT if bw_mode=0)
                        int  bw_mode,    // <--  backward prediction mode (1-7, 0=DIRECT if fw_mode=0)
                        int  fw_ref_idx, // <--  reference frame for forward prediction (-1: Intra4x4 pred. with fw_mode)
                        int  bw_ref_idx  )    
{
  static int fw_pred[16];
  static int bw_pred[16];

  int  i, j;
  int  block_x4  = block_x+4;
  int  block_y4  = block_y+4;
  int  pic_opix_x = img->opix_x + block_x;
  int  pic_opix_y = img->opix_y + block_y;
  int  bx        = block_x >> 2;
  int  by        = block_y >> 2;
  int* fpred     = fw_pred;
  int* bpred     = bw_pred;
//  int  direct    = (fw_mode == 0 && bw_mode == 0 && (img->type == B_SLICE));
//  int  skipped   = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE));

//  int  apply_weights = ( (input->WeightedPrediction && (img->type == P_SLICE || img->type == SP_SLICE)) ||
//                         (input->WeightedBiprediction && (img->type ==B_SLICE)));  
  int  apply_weights = ( (active_pps->weighted_pred_flag && (img->type== P_SLICE || img->type == SP_SLICE)) ||
                         (active_pps->weighted_bipred_idc && (img->type== B_SLICE)));  

  
  int  list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;

  if ((p_dir==0)||(p_dir==2))
  {
    OneComponentLumaPrediction4x4 (fw_pred, pic_opix_x, pic_opix_y, img->all_mv[bx][by][LIST_0][fw_ref_idx][fw_mode], fw_ref_idx, listX[0+list_offset]);   
  }

  if ((p_dir==1)||(p_dir==2))
  { 
    OneComponentLumaPrediction4x4 (bw_pred, pic_opix_x, pic_opix_y, img->all_mv[bx][by][LIST_1][bw_ref_idx][bw_mode], bw_ref_idx, listX[1+list_offset]);   
  }

  if (apply_weights)
  {

    if (p_dir==2)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][0] * *fpred++ + 
                                    wbp_weight[1][fw_ref_idx][bw_ref_idx][0] * *bpred++ + 
                                    2*wp_luma_round) >> (luma_log_weight_denom + 1)) + 
                                    ((wp_offset[0][fw_ref_idx][0] + wp_offset[1][bw_ref_idx][0] + 1)>>1)); 
    }
    else if (p_dir==0)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a(((wp_weight[0][fw_ref_idx][0] * *fpred++  + wp_luma_round) >> luma_log_weight_denom) +
          + wp_offset[0][fw_ref_idx][0] );
    }
    else // p_dir==1
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a(((wp_weight[1][bw_ref_idx][0] * *bpred++  + wp_luma_round) >> luma_log_weight_denom) +
          wp_offset[1][bw_ref_idx][0] );
    }


  }
  else
  {
    if (p_dir==2)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
    }
    else if (p_dir==0)
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *fpred++;
    }
    else // p_dir==1
    {
      for   (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *bpred++;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Residual Coding of an 8x8 Luma block (not for intra)
 ************************************************************************
 */
int                                       //  ==> coefficient cost
LumaResidualCoding8x8 (int  *cbp,         //  --> cbp (updated according to processed 8x8 luminance block)8*8是否有残差
                       int  *cbp_blk,     //  --> block cbp (updated according to processed 8x8 luminance block)4*4是否有残差
                       int  block8x8,     // <--  block number of 8x8 block
                       int  p_dir,        // <--  prediction direction
                       int  fw_mode,      // <--  forward  prediction mode (1-7, 0=DIRECT)
                       int  bw_mode,      // <--  backward prediction mode (1-7, 0=DIRECT)
                       int  fw_refframe,  // <--  reference frame for forward prediction
                       int  bw_refframe   // <--  reference frame for backward prediction
                       )
{
  int    block_y, block_x, pic_pix_y, pic_pix_x, i, j, nonzero, cbp_blk_mask;
  int    coeff_cost = 0;
  int    mb_y       = (block8x8 / 2) << 3;
  int    mb_x       = (block8x8 % 2) << 3;
  int    cbp_mask   = 1 << block8x8;
  int    bxx, byy;                   // indexing curr_blk
  int    scrFlag = 0;                // 0=noSCR, 1=strongSCR, 2=jmSCR
  int    skipped    = (fw_mode == 0 && bw_mode == 0 && (img->type != B_SLICE));

  if (img->type==B_SLICE)
    scrFlag = 1;

  //===== loop over 4x4 blocks =====
  for (byy=0, block_y=mb_y; block_y<mb_y+8; byy+=4, block_y+=4)
  {
    pic_pix_y = img->opix_y + block_y;

    for (bxx=0, block_x=mb_x; block_x<mb_x+8; bxx+=4, block_x+=4)
    {
      pic_pix_x = img->opix_x + block_x;

      cbp_blk_mask = (block_x>>2) + block_y;

      //===== prediction of 4x4 block =====没有帧内的预测 暂时删除
      LumaPrediction4x4 (block_x, block_y, p_dir, fw_mode, bw_mode, fw_refframe, bw_refframe);

      //===== get displaced frame difference ======                
      for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        img->m7[i][j] = imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[i+block_x][j+block_y];
      }

      //===== DCT, Quantization, inverse Quantization, IDCT, Reconstruction =====      
      if (img->NoResidueDirect != 1 && !skipped  )
      {
        //===== DCT, Quantization, inverse Quantization, IDCT, Reconstruction =====
        if (img->type!=SP_SLICE)  nonzero = dct_luma   (block_x, block_y, &coeff_cost, 0);
        else                      nonzero = dct_luma_sp(block_x, block_y, &coeff_cost);
        if (nonzero)
        {
          (*cbp_blk) |= 1 << cbp_blk_mask;  // one bit for every 4x4 block
          (*cbp)     |= cbp_mask;           // one bit for the 4x4 blocks of an 8x8 block
        }
      }
    }
  }
//   /*
//   The purpose of the action below is to prevent that single or 'expensive' coefficients are coded.
//   With 4x4 transform there is larger chance that a single coefficient in a 8x8 or 16x16 block may be nonzero.
//   A single small (level=1) coefficient in a 8x8 block will cost: 3 or more bits for the coefficient,
//   4 bits for EOBs for the 4x4 blocks,possibly also more bits for CBP.  Hence the total 'cost' of that single
//   coefficient will typically be 10-12 bits which in a RD consideration is too much to justify the distortion improvement.
//   The action below is to watch such 'single' coefficients and set the reconstructed block equal to the prediction according
//   to a given criterium.  The action is taken only for inter luma blocks.
// 
//   Notice that this is a pure encoder issue and hence does not have any implication on the standard.
//   coeff_cost is a parameter set in dct_luma() and accumulated for each 8x8 block.  If level=1 for a coefficient,
//   coeff_cost is increased by a number depending on RUN for that coefficient.The numbers are (see also dct_luma()): 3,2,2,1,1,1,0,0,...
//   when RUN equals 0,1,2,3,4,5,6, etc.
//   If level >1 coeff_cost is increased by 9 (or any number above 3). The threshold is set to 3. This means for example:
//   1: If there is one coefficient with (RUN,level)=(0,1) in a 8x8 block this coefficient is discarded.
//   2: If there are two coefficients with (RUN,level)=(1,1) and (4,1) the coefficients are also discarded
//   sum_cnt_nonz is the accumulation of coeff_cost over a whole macro block.  If sum_cnt_nonz is 5 or less for the whole MB,
//   all nonzero coefficients are discarded for the MB and the reconstructed block is set equal to the prediction.
//   */

  if (img->NoResidueDirect != 1 && !skipped && coeff_cost <= _LUMA_COEFF_COST_)
  {
    coeff_cost  = 0;
    (*cbp)     &=  (63 - cbp_mask);
    (*cbp_blk) &= ~(51 << (4*block8x8-2*(block8x8%2)));

    for (i=mb_x; i<mb_x+8; i++)
    for (j=mb_y; j<mb_y+8; j++)
    {
      enc_picture->imgY[img->pix_y+j][img->pix_x+i] = img->mpr[i][j];
    }
    if (img->type==SP_SLICE)
    {
      for (i=mb_x; i < mb_x+BLOCK_SIZE*2; i+=BLOCK_SIZE)
        for (j=mb_y; j < mb_y+BLOCK_SIZE*2; j+=BLOCK_SIZE)
          copyblock_sp(i,j);
    }
  }

  return coeff_cost;
}

int                                        //  ==> coefficient cost
LumaResidualCoding8x8_fract (int  *cbp,         //  --> cbp (updated according to processed 8x8 luminance block)8*8是否有残差
							 int  *cbp_blk,     //  --> block cbp (updated according to processed 8x8 luminance block)4*4是否有残差
							 int  block8x8//,     // <--  block number of 8x8 block
							 //                        int  p_dir,        // <--  prediction direction
							 //                        int  fw_mode,      // <--  forward  prediction mode (1-7, 0=DIRECT)
							 //                        int  bw_mode,      // <--  backward prediction mode (1-7, 0=DIRECT)
							 //                        int  fw_refframe//,  // <--  reference frame for forward prediction
							 //                        int  bw_refframe   // <--  reference frame for backward prediction
							 )
{
	
	int    block_y, block_x, pic_pix_y, pic_pix_x, i, j, nonzero, cbp_blk_mask;
	int    coeff_cost = 0;
	int    mb_y       = (block8x8 / 2) << 3;//0,8 //当前要处理的8*8块在当前宏块中的坐标（以像素为单位）
	int    mb_x       = (block8x8 % 2) << 3;//0,8 //当前要处理的8*8块在当前宏块中的坐标（以像素为单位）
	int    cbp_mask   = 1 << block8x8;//1,2,4,8
	int    bxx, byy;                   // indexing curr_blk
	

// 	printf("mb_x=%d,mb_y=%d",mb_x,mb_y);

	//===== loop over 4x4 blocks =====将一个8*8块分成4个4*4块，分别处理
	for (block_y=mb_y; block_y<mb_y+8; block_y+=4)
	{
		for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
		{
			pic_pix_y = img->pix_y + block_y;//当前4*4块在图像中的像素坐标
			pic_pix_x = img->pix_x + block_x;
			
			cbp_blk_mask = (block_x>>2) + block_y;//什么意思？？？
// 			printf("\nblock_y=%2d,block_x=%2d",block_y,block_x);
			//       //===== prediction of 4x4 block =====
			//       LumaPrediction4x4 (block_x, block_y, p_dir, fw_mode, bw_mode, fw_refframe, bw_refframe);
			
			//===== get displaced frame difference ====== 
//            printf("\nrec(org)");
			for (j=0; j<4; j++)
			{				
// 				printf("\n");
				for (i=0; i<4; i++)
				{
					img->mpr[i+block_x][j+block_y] = imgY_rec[pic_pix_y+j][pic_pix_x+i];//重建值赋给了img->mpr
					img->m7[i][j] = imgY_org[pic_pix_y+j][pic_pix_x+i] - imgY_rec[pic_pix_y+j][pic_pix_x+i]; //计算一个4*4块的残差，原始值减去重建值，赋给了img->m7[i][j];
// 					printf("%3d(-%3d)=%3d ",imgY_rec[pic_pix_y+j][pic_pix_x+i],imgY_org[pic_pix_y+j][pic_pix_x+i],img->m7[i][j]);
				}
			}	

			//=====残差 DCT, 量化, 反量化, IDCT =====      
			if (img->NoResidueDirect != 1  )//yes
			{
			
				if (img->type!=SP_SLICE)  nonzero = dct_luma   (block_x, block_y, &coeff_cost, 0);//yes
				else                      nonzero = dct_luma_sp(block_x, block_y, &coeff_cost);
				if (nonzero)//有非零系数
				{
					(*cbp_blk) |= 1 << cbp_blk_mask;  // one bit for every 4x4 block
					(*cbp)     |= cbp_mask;           // one bit for the 4x4 blocks of an 8x8 block
				}
			}
		}
	}
	
	
	if (img->NoResidueDirect != 1 && coeff_cost <= _LUMA_COEFF_COST_)//_LUMA_COEFF_COST_=4// && !skipped
	{
		coeff_cost  = 0;
		(*cbp)     &=  (63 - cbp_mask);
		(*cbp_blk) &= ~(51 << (4*block8x8-2*(block8x8%2)));
		
		for (i=mb_x; i<mb_x+8; i++)
			for (j=mb_y; j<mb_y+8; j++)
			{
				enc_picture->imgY[img->pix_y+j][img->pix_x+i] = img->mpr[i][j];//img->mpr[i][j]是什么意思呢？？？
			}
			//     if (img->type==SP_SLICE)
			//     {
			//       for (i=mb_x; i < mb_x+BLOCK_SIZE*2; i+=BLOCK_SIZE)
			//         for (j=mb_y; j < mb_y+BLOCK_SIZE*2; j+=BLOCK_SIZE)
			//           copyblock_sp(i,j);
			//     }
	}
	
	return coeff_cost;
	
}


/*!
 ************************************************************************
 * \brief
 *    Set mode parameters and reference frames for an 8x8 block
 ************************************************************************
 */
void
SetModesAndRefframe (int b8, int* p_dir, int* fw_mode, int* bw_mode, int* fw_ref, int* bw_ref)
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];
  int         j      = 2*(b8/2);
  int         i      = 2*(b8%2);

  *fw_mode = *bw_mode = *fw_ref = *bw_ref = -1;

  *p_dir  = currMB->b8pdir[b8];

  if (img->type!=B_SLICE)
  {
//     *fw_ref = enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j];
    *bw_ref = 0;
    *fw_mode  = currMB->b8mode[b8];
    *bw_mode  = 0;
  }
  else
  {
    if (currMB->b8pdir[b8]==-1)
    {
      *fw_ref   = -1;
      *bw_ref   = -1;
      *fw_mode  =  0;
      *bw_mode  =  0;
    }
    else if (currMB->b8pdir[b8]==0)
    {
//       *fw_ref   = enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j];
      *bw_ref   = 0;
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = 0;
    }
    else if (currMB->b8pdir[b8]==1)
    {
      *fw_ref   = 0;
//       *bw_ref   = enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j];
      *fw_mode  = 0;
      *bw_mode  = currMB->b8mode[b8];
    }
    else
    {
//       *fw_ref   = enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j];
//       *bw_ref   = enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j];
      *fw_mode  = currMB->b8mode[b8];
      *bw_mode  = currMB->b8mode[b8];
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Residual Coding of a Luma macroblock (not for intra)
 ************************************************************************
 */
void
LumaResidualCoding ()
{
  int i,j,block8x8,b8_x,b8_y;
  int p_dir, fw_mode, bw_mode, refframe;
  int sum_cnt_nonz;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  currMB->cbp     = 0 ;
  currMB->cbp_blk = 0 ;
  sum_cnt_nonz    = 0 ;

  for (block8x8=0; block8x8<4; block8x8++)
  {
    int bw_ref;
    SetModesAndRefframe (block8x8, &p_dir, &fw_mode, &bw_mode, &refframe, &bw_ref);

    sum_cnt_nonz += LumaResidualCoding8x8 (&(currMB->cbp), &(currMB->cbp_blk), block8x8,
                                           p_dir, fw_mode, bw_mode, refframe, bw_ref);
  }

  if (sum_cnt_nonz <= 5 )
  {
     currMB->cbp     &= 0xfffff0 ;
     currMB->cbp_blk &= 0xff0000 ;
     for (i=0; i < MB_BLOCK_SIZE; i++)
     {
       for (j=0; j < MB_BLOCK_SIZE; j++)
       {
         enc_picture->imgY[img->pix_y+j][img->pix_x+i]=img->mpr[i][j];
       }
     }
     if (img->type==SP_SLICE)
     {
       for(block8x8=0;block8x8<4;block8x8++)
       {
         b8_x=(block8x8&1)<<3;
         b8_y=(block8x8&2)<<2;
         for (i=0;i<8;i+=4)
           for (j=0;j<8;j+=4)
             copyblock_sp(b8_x+i,b8_y+j);
       }
     }
   }
}

void LumaResidualCoding_fract(int block16_x,int block16_y)
{
	int i,j,block8x8,b8_x,b8_y;
	int p_dir, fw_mode, bw_mode, refframe;
	int sum_cnt_nonz;
	Macroblock *currMB = &img->mb_data[img->current_mb_nr];//当前宏块数据
	
	currMB->cbp     = 0 ;
	currMB->cbp_blk = 0 ;
	sum_cnt_nonz    = 0 ;
	
	img->pix_x=block16_x*16;//当前宏块的x坐标（以像素为单位）
	img->pix_y=block16_y*16;//当前宏块的y坐标（以像素为单位）
	
	for (block8x8=0; block8x8<4; block8x8++)//以8*8为基本单位进行处理，所以要循环4次
	{
// 		printf("\n\nblock8=%d\n",block8x8);
		sum_cnt_nonz += LumaResidualCoding8x8_fract (&(currMB->cbp), &(currMB->cbp_blk), block8x8);//currMB->cbp_blk是什么意思？？？？
	}
	
	if (sum_cnt_nonz <= 5 )//添加这个之后PSNR减少不到0.1，这段的用处是什么？？？
	{
		currMB->cbp     &= 0xfffff0 ;
		currMB->cbp_blk &= 0xff0000 ;
		for (i=0; i < MB_BLOCK_SIZE; i++)
		{
			for (j=0; j < MB_BLOCK_SIZE; j++)
			{
				enc_picture->imgY[img->pix_y+j][img->pix_x+i]=img->mpr[i][j];//img->mpr[i][j]是什么意思？？？
			}
		}
		if (img->type==SP_SLICE)
		{
			for(block8x8=0;block8x8<4;block8x8++)
			{
				b8_x=(block8x8&1)<<3;
				b8_y=(block8x8&2)<<2;
				for (i=0;i<8;i+=4)
					for (j=0;j<8;j+=4)
						copyblock_sp(b8_x+i,b8_y+j);
			}
		}
	}
	
}

/*!
 ************************************************************************
 * \brief
 *    Predict one component of a chroma 4x4 block
 ************************************************************************
 */
void
OneComponentChromaPrediction4x4 (int*      mpred,      //!< array to store prediction values
                                 int       block_c_x,  //!< horizontal pixel coordinate of 4x4 block
                                 int       block_c_y,  //!< vertical   pixel coordinate of 4x4 block
                                 int****** mv,         //!< motion vector array
                                 int       list_idx,   //!< reference picture list
                                 int       ref,        //!< reference index
                                 int       blocktype,  //!< block type
                                 int       uv)         //!< chroma component
{
  int     i, j, ii, jj, ii0, jj0, ii1, jj1, if0, if1, jf0, jf1;
  int*    mvb;
  pel_t** refimage;
  int     f1        = 8 , f2=f1-1, f3=f1*f1, f4=f3>>1;
  int     s1        = 3;
  int     list_offset;
  int     max_y_cr;

  StorablePicture **list;

  int curr_mb_field = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field));
  // find out the correct list offsets
  if (curr_mb_field)
  {
    if(img->current_mb_nr%2)
      list_offset = 4; // top field mb
    else
      list_offset = 2; // bottom field mb
    max_y_cr = img->height_cr/2-1;
  }
  else
  {
    list_offset = 0;  // no mb aff or frame mb
    max_y_cr = img->height_cr-1;
  }


  list      = listX[list_idx + list_offset];

  refimage  = list[ref]->imgUV[uv];
  
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    //ref_idx   = enc_picture->ref_idx[list_idx][(img->pix_c_x + block_c_x)>>1][(img->pix_c_y + block_c_y)>>1];  

    mvb  = mv [(i+block_c_x)>>1][(j+block_c_y)>>1][list_idx][ref][blocktype];

    ii   = ((i+block_c_x+img->opix_c_x)<<s1) + mvb[0];
    jj   = ((j+block_c_y+img->opix_c_y)<<s1) + mvb[1];
    
    jj  += list[ref]->chroma_vector_adjustment;

    ii0  = max (0, min (img->width_cr -1, ii>>s1     ));
    jj0  = max (0, min (max_y_cr, jj>>s1     ));
    ii1  = max (0, min (img->width_cr -1, (ii>>s1)+1));
    jj1  = max (0, min (max_y_cr, (jj>>s1)+1));

    if1  = (ii&f2);  if0 = f1-if1;
    jf1  = (jj&f2);  jf0 = f1-jf1;

    *mpred++ = (if0 * jf0 * refimage[jj0][ii0] +
                if1 * jf0 * refimage[jj0][ii1] +
                if0 * jf1 * refimage[jj1][ii0] +
                if1 * jf1 * refimage[jj1][ii1] + f4) >> 6;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Predict an intra chroma 4x4 block
 ************************************************************************
 */

void IntraChromaPrediction4x4 (int uv,       // <-- colour component
                               int block_x,  // <-- relative horizontal block coordinate of 4x4 block
                               int block_y)  // <-- relative vertical   block coordinate of 4x4 block
{
  int mode = img->mb_data[img->current_mb_nr].c_ipred_mode;
  int i, j;

  //===== prediction =====
  for (j=block_y; j<block_y+4; j++)
  for (i=block_x; i<block_x+4; i++)
  {
    img->mpr[i][j] = img->mprr_c[uv][mode][i][j];
  }
}

/*!
 ************************************************************************
 * \brief
 *    Predict one chroma 4x4 block
 ************************************************************************
 */
void
ChromaPrediction4x4 (int  uv,           // <-- colour component
                     int  block_x,      // <-- relative horizontal block coordinate of 4x4 block
                     int  block_y,      // <-- relative vertical   block coordinate of 4x4 block
                     int  p_dir,        // <-- prediction direction
                     int  fw_mode,      // <-- forward  prediction mode (1-7, 0=DIRECT if bw_mode=0)
                     int  bw_mode,      // <-- backward prediction mode (1-7, 0=DIRECT if fw_mode=0)
                     int  fw_ref_idx,   // <-- reference frame for forward prediction (if (<0) -> intra prediction)
                     int  bw_ref_idx)   // <-- reference frame for backward prediction 
{
  static int fw_pred[16];
  static int bw_pred[16];

  int  i, j;
  int  block_x4   = block_x+4;
  int  block_y4   = block_y+4;
  int* fpred      = fw_pred;
  int* bpred      = bw_pred;
  int****** mv_array = img->all_mv;

  //int apply_weights = ( (input->WeightedPrediction && (img->type == P_SLICE||img->type == SP_SLICE)) ||
//                     (input->WeightedBiprediction && (img->type == B_SLICE)));
  int  apply_weights = ( (active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
                         (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));  


  //===== INTRA PREDICTION =====
  if (p_dir==-1)
  {
    IntraChromaPrediction4x4 (uv, block_x, block_y);
    return;
  }
  
  //===== INTER PREDICTION =====
  if ((p_dir==0) || (p_dir==2))
  {
    OneComponentChromaPrediction4x4 (fw_pred, block_x, block_y, mv_array, LIST_0, fw_ref_idx, fw_mode, uv);
  }
  if ((p_dir==1) || (p_dir==2))
  {
    OneComponentChromaPrediction4x4 (bw_pred, block_x, block_y, mv_array, LIST_1, bw_ref_idx, bw_mode, uv);
  }

  if (apply_weights)
  {
    if (p_dir==2)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
            img->mpr[i][j] =  clip1a(((wbp_weight[0][fw_ref_idx][bw_ref_idx][uv+1] * *fpred++ + wbp_weight[1][fw_ref_idx][bw_ref_idx][uv+1] * *bpred++ 
                  + 2*wp_chroma_round) >> (chroma_log_weight_denom + 1)) + ((wp_offset[0][fw_ref_idx][uv+1] + wp_offset[1][bw_ref_idx][uv+1] + 1)>>1) );
    }
    else if (p_dir==0)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
           img->mpr[i][j] = clip1a(((wp_weight[0][fw_ref_idx][uv+1] * *fpred++ + wp_chroma_round) >> chroma_log_weight_denom) +  wp_offset[0][fw_ref_idx][uv+1]);
    }
    else // (p_dir==1)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = clip1a(((wp_weight[1][bw_ref_idx][uv+1] * *bpred++ + wp_chroma_round) >> chroma_log_weight_denom) + wp_offset[1][bw_ref_idx][uv+1]);
    }
  }
  else
  {
    if (p_dir==2)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  
          img->mpr[i][j] = (*fpred++ + *bpred++ + 1) / 2; 
    }
    else if (p_dir==0)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *fpred++;
    }
    else // (p_dir==1)
    {
      for (j=block_y; j<block_y4; j++)
        for (i=block_x; i<block_x4; i++)  img->mpr[i][j] = *bpred++;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Chroma residual coding for an macroblock
 ************************************************************************
 */
void ChromaResidualCoding (int* cr_cbp)
{
  int   uv, block8, block_y, block_x, j, i;
  int   p_dir, fw_mode, bw_mode, refframe;
  int   skipped = (img->mb_data[img->current_mb_nr].mb_type == 0 && (img->type == P_SLICE || img->type == SP_SLICE));
  int   bw_ref;

  for (*cr_cbp=0, uv=0; uv<2; uv++)
  {
    //===== prediction of chrominance blocks ===d==
    block8 = 0;
    for (block_y=0; block_y<8; block_y+=4)
    for (block_x=0; block_x<8; block_x+=4, block8++)
    {
      SetModesAndRefframe (block8, &p_dir, &fw_mode, &bw_mode, &refframe, &bw_ref);

      ChromaPrediction4x4 (uv, block_x, block_y, p_dir, fw_mode, bw_mode, refframe, bw_ref);
    }

        // ==== set chroma residue to zero for skip Mode in SP frames 
    if (img->NoResidueDirect)
    {
      for (j=0; j<8; j++)
        for (i=0; i<8; i++)
        {
          enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i] = img->mpr[i][j];
        }
    }
    else
      if (skipped && img->type==SP_SLICE)
      {
        for (j=0; j<8; j++)
          for (i=0; i<8; i++)
          {
            img->m7[i][j] = 0;
          }
      }
      else
        if (skipped)
        {
          for (j=0; j<8; j++)
            for (i=0; i<8; i++)
            {
              enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i] = img->mpr[i][j];
            }
        }
        else
        {
          for (j=0; j<8; j++)
            for (i=0; i<8; i++)
            {
              img->m7[i][j] = imgUV_org[uv][img->opix_c_y+j][img->opix_c_x+i] - img->mpr[i][j];
            }
        }

    //===== DCT, Quantization, inverse Quantization, IDCT, and Reconstruction =====
    //===== Call function for skip mode in SP frames to properly process frame ====
    
    if (skipped && img->type==SP_SLICE)
    {
      *cr_cbp=dct_chroma_sp(uv,*cr_cbp);
    }
    else
    {
      if (!img->NoResidueDirect && !skipped)
      {
        if (img->type!=SP_SLICE || IS_INTRA (&img->mb_data[img->current_mb_nr]))
          *cr_cbp=dct_chroma   (uv,*cr_cbp);
        else
          *cr_cbp=dct_chroma_sp(uv,*cr_cbp);
      }
    }
  }
  
  //===== update currMB->cbp =====
  img->mb_data[img->current_mb_nr].cbp += ((*cr_cbp)<<4);  
}






/*!
************************************************************************
* \brief
*    Chroma residual coding for an macroblock
************************************************************************
*/
void ChromaResidualCoding_fract(int* cr_cbp)
{
	int   uv, block8, block_y, block_x, j, i;
	int   p_dir, fw_mode, bw_mode, refframe;
	int   skipped = (img->mb_data[img->current_mb_nr].mb_type == 0 && (img->type == P_SLICE || img->type == SP_SLICE));
	int   bw_ref;
	
	for (*cr_cbp=0, uv=0; uv<2; uv++)
	{
		//===== prediction of chrominance blocks ===d==
// 		block8 = 0;
// 		for (block_y=0; block_y<8; block_y+=4)
// 			for (block_x=0; block_x<8; block_x+=4, block8++)
// 			{
// 				SetModesAndRefframe (block8, &p_dir, &fw_mode, &bw_mode, &refframe, &bw_ref);
// 				
// 				ChromaPrediction4x4 (uv, block_x, block_y, p_dir, fw_mode, bw_mode, refframe, bw_ref);
// 			}
			
			// ==== set chroma residue to zero for skip Mode in SP frames 
			if (img->NoResidueDirect)
			{
				for (j=0; j<8; j++)
					for (i=0; i<8; i++)
					{
						enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i] = img->mpr[i][j];
					}
			}
			else
			if (skipped && img->type==SP_SLICE)
			{
				for (j=0; j<8; j++)
					for (i=0; i<8; i++)
					{
						img->m7[i][j] = 0;
					}
			}
			else
			if (skipped)
			{
				for (j=0; j<8; j++)
					for (i=0; i<8; i++)
					{
						enc_picture->imgUV[uv][img->pix_c_y+j][img->pix_c_x+i] = img->mpr[i][j];
					}
			}
			else
			{
// 				printf("imgUV_org(imgUV_org-imgUV_rec)");
				for (j=0; j<8; j++)
				{
// 					printf("\n");
					for (i=0; i<8; i++)
					{

						img->mpr[i][j]= imgUV_rec[uv][img->opix_c_y+j][img->opix_c_x+i];
						img->m7[i][j] = imgUV_org[uv][img->opix_c_y+j][img->opix_c_x+i] - img->mpr[i][j];
// 					printf(" %3d(%3d) ",imgUV_org[uv][img->opix_c_y+j][img->opix_c_x+i],img->m7[i][j]);
					}
				}
// 				printf("\n");
			}
			
			//===== DCT, Quantization, inverse Quantization, IDCT, and Reconstruction =====
			//===== Call function for skip mode in SP frames to properly process frame ====
				
			if (skipped && img->type==SP_SLICE)
			{
				*cr_cbp=dct_chroma_sp(uv,*cr_cbp);
			}
			else
			{
				if (!img->NoResidueDirect && !skipped)
				{
					if (img->type!=SP_SLICE || IS_INTRA (&img->mb_data[img->current_mb_nr]))
						*cr_cbp=dct_chroma   (uv,*cr_cbp);
					else
						*cr_cbp=dct_chroma_sp(uv,*cr_cbp);

// 					printf("\nenc_picture->imgUV(enc_picture->imgUV-imgUV_org)");
					for (j=0; j<8; j++)
					{
// 						printf("\n");
						for (i=0; i<8; i++)
						{
// 							printf(" %3d(%3d) ",enc_picture->imgUV[uv][img->opix_c_y+j][img->opix_c_x+i],imgUV_org[uv][img->opix_c_y+j][img->opix_c_x+i] - imgUV_rec[uv][img->opix_c_y+j][img->opix_c_x+i]);
						}
					}
// 					printf("\n");

				}
			}
	}
	
	//===== update currMB->cbp =====
	img->mb_data[img->current_mb_nr].cbp += ((*cr_cbp)<<4);  
}






/*!
 ************************************************************************
 * \brief
 *    Predict an intra chroma 8x8 block
 ************************************************************************
 */
void IntraChromaPrediction8x8 (int *mb_up, int *mb_left, int*mb_up_left)
{

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int     s, s0, s1, s2, s3, i, j, k;
  pel_t** image;
  int     block_x, block_y;
  int     mb_nr             = img->current_mb_nr;
  int     mb_available_up;
  int     mb_available_left[2];
  int     mb_available_up_left;
  int     ih,iv;
  int     ib,ic,iaa;
  int     uv;
  int     hline[8], vline[9];
  int     mode;
  int     best_mode = DC_PRED_8;     //0    //just an initilaization here, should always be overwritten
  int     cost;
  int     min_cost;
  int     diff[16];
  PixelPos up;       //!< pixel position p(0,-1)
  PixelPos left[9];  //!< pixel positions p(-1, -1..8)


  for (i=0;i<9;i++)
  {
    getNeighbour(mb_nr, -1 ,  i-1 , 0, &left[i]);
  }
  
  getNeighbour(mb_nr, 0     ,  -1 , 0, &up);


  mb_available_up                             = up.available;//0
  mb_available_up_left                        = left[0].available;//0
  mb_available_left[0] = mb_available_left[1] = left[1].available;//0

  if(input->UseConstrainedIntraPred)
  {
    mb_available_up      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, mb_available_left[0]=1; i<5;i++)
      mb_available_left[0]  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    for (i=5, mb_available_left[1]=1; i<9;i++)
      mb_available_left[1]  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    mb_available_up_left = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  if (mb_up)
    *mb_up = mb_available_up;
  if (mb_left)
    *mb_left = mb_available_left[0] && mb_available_left[1];
  if( mb_up_left )
    *mb_up_left = mb_available_up_left;

  // compute all chroma intra prediction modes for both U and V
  for (uv=0; uv<2; uv++)
  {
    image = enc_picture->imgUV[uv];

    // DC prediction
    for (block_y=0; block_y<8; block_y+=4)
    for (block_x=0; block_x<8; block_x+=4)
    {
      s=128;
      s0=s1=s2=s3=0;
      //===== get prediction value =====
      switch ((block_y>>1) + (block_x>>2))
      {
      case 0:  //===== TOP LEFT =====
        if      (mb_available_up)       for (i=0;i<4;i++)  s0 += image[up.pos_y][up.pos_x + i];
        if      (mb_available_left[0])  for (i=1;i<5;i++)  s2 += image[left[i].pos_y][left[i].pos_x];
        if      (mb_available_up && mb_available_left[0])  s  = (s0+s2+4) >> 3;
        else if (mb_available_up)                          s  = (s0   +2) >> 2;
        else if (mb_available_left[0])                     s  = (s2   +2) >> 2;
        break;
      case 1: //===== TOP RIGHT =====
        if      (mb_available_up)       for (i=4;i<8;i++)  s1 += image[up.pos_y][up.pos_x + i];
        else if (mb_available_left[0])  for (i=1;i<5;i++)  s2 += image[left[i].pos_y][left[i].pos_x];
        if      (mb_available_up)                          s  = (s1   +2) >> 2;
        else if (mb_available_left[0])                     s  = (s2   +2) >> 2;
        break;
      case 2: //===== BOTTOM LEFT =====
        if      (mb_available_left[1])  for (i=5;i<9;i++)  s3 += image[left[i].pos_y][left[i].pos_x];
        else if (mb_available_up)       for (i=0;i<4;i++)  s0 += image[up.pos_y][up.pos_x + i];
        if      (mb_available_left[1])                     s  = (s3   +2) >> 2;
        else if (mb_available_up)                          s  = (s0   +2) >> 2;
        break;
      case 3: //===== BOTTOM RIGHT =====
        if      (mb_available_up)       for (i=4;i<8;i++)  s1 += image[up.pos_y][up.pos_x + i];
        if      (mb_available_left[1])  for (i=5;i<9;i++)  s3 += image[left[i].pos_y][left[i].pos_x];
        if      (mb_available_up && mb_available_left[1])  s  = (s1+s3+4) >> 3;
        else if (mb_available_up)                          s  = (s1   +2) >> 2;
        else if (mb_available_left[1])                     s  = (s3   +2) >> 2;
        break;
      }


      //===== prediction =====
      for (j=block_y; j<block_y+4; j++)
      for (i=block_x; i<block_x+4; i++)
      {
        img->mprr_c[uv][DC_PRED_8][i][j] = s;
      }
    }

    // vertical prediction
    if (mb_available_up)
    {
      for (i=0; i<8; i++)
        hline[i] = image[up.pos_y][up.pos_x + i];
      for (i=0; i<8; i++)
      for (j=0; j<8; j++)
        img->mprr_c[uv][VERT_PRED_8][i][j] = hline[i];
    }

    // horizontal prediction 
    if (mb_available_left[0] && mb_available_left[1])
    {
      for (i=1; i<9; i++)
        vline[i] = image[left[i].pos_y][left[i].pos_x];
      for (i=0; i<8; i++)
      for (j=0; j<8; j++)
        img->mprr_c[uv][HOR_PRED_8][i][j] = vline[j+1]; 
    }

    // plane prediction 
    if (mb_available_left[0] && mb_available_left[1] && mb_available_up && mb_available_up_left)
    {
      ih = 4*(hline[7] - image[left[0].pos_y][left[0].pos_x]);
      iv = 4*(vline[7+1] - image[left[0].pos_y][left[0].pos_x]);
      for (i=1;i<4;i++)
      {
        ih += i*(hline[3+i] - hline[3-i]);
        iv += i*(vline[3+i+1] - vline[3-i+1]);
      }
      ib=(17*ih+16)>>5;
      ic=(17*iv+16)>>5;

      iaa=16*(hline[7]+vline[7+1]);
      for (j=0; j<8; j++)
      for (i=0; i<8; i++)
        img->mprr_c[uv][PLANE_8][i][j]=max(0,min(255,(iaa+(i-3)*ib +(j-3)*ic + 16)/32));// store plane prediction
    }
  }

  if (!input->rdopt) // the rd-opt part does not work correctly (see  i_encode_one_macroblock)
  {                       // since ipredmodes could be overwritten => encoder-decoder-mismatches
    // pick lowest cost prediction mode
    min_cost = 1<<20;
    for (i=0;i<8;i++)
    {
      getNeighbour(mb_nr, 0 ,  i , 0, &left[i]);
    }
    for (mode=DC_PRED_8; mode<=PLANE_8; mode++)
    {
      if ((mode==VERT_PRED_8 && !mb_available_up) ||
          (mode==HOR_PRED_8 && (!mb_available_left[0] || !mb_available_left[1])) ||
          (mode==PLANE_8 && (!mb_available_left[0] || !mb_available_left[1] || !mb_available_up || !mb_available_up_left)))
        continue;

      cost = 0;
      for (uv=0; uv<2; uv++)
      {
        image = imgUV_org[uv];
        for (block_y=0; block_y<8; block_y+=4)
        for (block_x=0; block_x<8; block_x+=4)
        {
          for (k=0,j=block_y; j<block_y+4; j++)
          for (i=block_x; i<block_x+4; i++,k++)
          {
            diff[k] = image[left[j].pos_y][left[j].pos_x+i] - img->mprr_c[uv][mode][i][j];
          }
          cost += SATD(diff, input->hadamard);
        }
      }
      if (cost < min_cost)
      {
        best_mode = mode;
        min_cost = cost;
      }
    }

    currMB->c_ipred_mode = best_mode;
  }
 
}


/*!
 ************************************************************************
 * \brief
 *    Check if all reference frames for a macroblock are zero
 ************************************************************************
 */
int
ZeroRef (Macroblock* currMB)
{
  int i,j;

  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    if (enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]!=0)
    {
        return 0;
    }
  }
  return 1;
}

/*!
 ************************************************************************
 * \brief
 *    Converts macroblock type to coding value
 ************************************************************************
 */
int
MBType2Value (Macroblock* currMB)
{
  static const int dir1offset[3]    =  { 1,  2, 3};
  static const int dir2offset[3][3] = {{ 0,  4,  8},   // 1. block forward
                                       { 6,  2, 10},   // 1. block backward
                                       {12, 14, 16}};  // 1. block bi-directional

  int mbtype, pdir0, pdir1;

  if (img->type!=B_SLICE)
  {
    if      (currMB->mb_type==I4MB)     return (img->type==I_SLICE ? 0 : 6);
    else if (currMB->mb_type==I16MB)    return (img->type==I_SLICE ? 0 : 6) + img->i16offset;
    else if (currMB->mb_type==P8x8)
    {
      if (input->symbol_mode==UVLC && ZeroRef (currMB))  return 5;
      else                                               return 4;
    }
    else                                return currMB->mb_type;
  }
  else
  {
    mbtype = currMB->mb_type;
    pdir0  = currMB->b8pdir[0];
    pdir1  = currMB->b8pdir[3];

    if      (mbtype==0)       return 0;
    else if (mbtype==I4MB)    return 23;
    else if (mbtype==I16MB)   return 23 + img->i16offset;
    else if (mbtype==P8x8)    return 22;
    else if (mbtype==1)       return dir1offset[pdir0];
    else if (mbtype==2)       return 4 + dir2offset[pdir0][pdir1];
    else                      return 5 + dir2offset[pdir0][pdir1];
  }
}

/*!
 ************************************************************************
 * \brief
 *    Writes intra prediction modes for an 8x8 block
 ************************************************************************
 */
int writeIntra4x4Modes(int only_this_block)
{
  int i,j,bs_x,bs_y,ii,jj;
  int block8x8;
  int rate;
  int ipred_array[16],cont_array[16],ipred_number;
  Macroblock    *currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  int           *bitCount   = currMB->bitcounter;
  Slice         *currSlice  = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap    = assignSE2partition[input->partition_mode];

  ipred_number=0;
  for(block8x8=0;block8x8<4;block8x8++)
  {
    if( currMB->b8mode[block8x8]==IBLOCK && (only_this_block<0||only_this_block==block8x8) )
    {
      bs_x=bs_y=4;
      ii=(bs_x>>2); // bug fix for solaris. mwi 
      jj=(bs_y>>2); // bug fix for solaris. mwi
      
      for(j=0;j<2;j+=jj)
      {
        for(i=0;i<2;i+=ii)
        {
          ipred_array[ipred_number]=currMB->intra_pred_modes[(block8x8<<2)|(j<<1)|i];
          cont_array[ipred_number]=(block8x8<<2)+(j<<1)+i;
          ipred_number++;
        }
      }
    }
  }
  rate=0;

  for(i=0;i<ipred_number;i++)
  {
    currMB->IntraChromaPredModeFlag = 1;
    currSE->context = cont_array[i];
    currSE->value1  = ipred_array[i];
    currSE->value2  = 0;

#if TRACE
    printf(currSE->tracestring, TRACESTRING_SIZE, "Intra mode     = %3d %d",currSE->value1,currSE->context);
#endif

//     /*--- set symbol type and function pointers ---*/
    if (input->symbol_mode != UVLC)    currSE->writing = writeIntraPredMode_CABAC;       
    currSE->type = SE_INTRAPREDMODE;

    /*--- choose data partition ---*/
    dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);   
    /*--- encode and update rate ---*/
    if (input->symbol_mode == UVLC)    writeSyntaxElement_Intra4x4PredictionMode(currSE, dataPart);
    else                               dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB]+=currSE->len;
    rate += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }

  return rate;
}

/*!
 ************************************************************************
 * \brief
 *    Converts 8x8 block tyoe to coding value
 ************************************************************************
 */
int
B8Mode2Value (int b8mode, int b8pdir)
{
  static const int b8start[8] = {0,0,0,0, 1, 4, 5, 10};
  static const int b8inc  [8] = {0,0,0,0, 1, 2, 2, 1};

  if (img->type!=B_SLICE)
  {
    return (b8mode-4);
  }
  else
  {
    return b8start[b8mode] + b8inc[b8mode] * b8pdir;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Codes macroblock header//熵编码宏块头（宏块类型）
 ************************************************************************
 */
int writeMBHeader(int rdopt)  // GB CHROMA !!!!!!!!
{
  int             i,j;
  int             mb_nr     = img->current_mb_nr;
  Macroblock*     currMB    = &img->mb_data[mb_nr];//当前宏块
  Macroblock*     prevMB    = mb_nr ? (&img->mb_data[mb_nr-1]) : NULL;//前一宏块
  SyntaxElement  *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount  = currMB->bitcounter;
  Slice*          currSlice = img->currentSlice;
  DataPartition*  dataPart;
  const int*      partMap   = assignSE2partition[input->partition_mode];//=assignSE2partition_NoDP 对编码东西根据重要性分区ＡＢＣ
  int             no_bits   = 0;
  int             skip      = currMB->mb_type ? 0:((img->type == B_SLICE) ? !currMB->cbp:1);//0????不一样啊
  int             mb_type;
  int			  prevMbSkipped = 0;
  int			  mb_field_tmp;
  Macroblock	 *topMB = NULL;
  
  int             WriteFrameFieldMBInHeader = 0;
//    printf("\nwriteMBHeader");

  if (img->MbaffFrameFlag)//no
  {
    if (0==(mb_nr%2))
    {
      WriteFrameFieldMBInHeader = 1; // top field

      prevMbSkipped = 0;
    }
    else
    {
      if (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1))
      {
        WriteFrameFieldMBInHeader = 1; // bottom, if top was skipped
      }

      topMB= &img->mb_data[img->current_mb_nr-1];
      if(!(img->type == B_SLICE))
        prevMbSkipped = (topMB->mb_type == 0);
      else 
        prevMbSkipped = (topMB->mb_type == 0 && topMB->cbp == 0);
    }
  }
  currMB->IntraChromaPredModeFlag = IS_INTRA(currMB);//0

  // choose the appropriate data partition
  dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);//
  
  //=====  BITS FOR MACROBLOCK MODE =====
  if(img->type == I_SLICE)//GB//no
  {
    // write mb_aff ---------------------------熵编码 帧场编码模式
    if(img->MbaffFrameFlag && !skip) // check for copy mode, Krit            //no
    {
      if(WriteFrameFieldMBInHeader)
      {
        currSE->value1 =  currMB->mb_field;
        currSE->value2 = 0;
        currSE->type   =  SE_MBTYPE;
        
        if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
        else                            currSE->writing = writeFieldModeInfo_CABAC;
          
#if TRACE
        printf(currSE->tracestring, TRACESTRING_SIZE, "Field mode = %3d",currMB->mb_field);
#endif		
        if( input->symbol_mode==UVLC)
        {
          currSE->bitpattern = (currMB->mb_field ? 1 : 0);
          currSE->len = 1;
          writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
        }
        else
        {
          dataPart->writeSyntaxElement(currSE, dataPart);
        }

        bitCount[BITS_MB_MODE] += currSE->len;
//  		printf("\n1、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
    
    // write mb_type//-----------------熵编码宏块类型
    currSE->value1  = MBType2Value (currMB);//23//0
    currSE->value2  = 0;
    currSE->type    = SE_MBTYPE;//2
    
    if (input->symbol_mode == UVLC)  currSE->mapping = ue_linfo;
    else                             currSE->writing = writeMB_typeInfo_CABAC;
    
    dataPart->writeSyntaxElement( currSE, dataPart);//语法元素写入码流
#if TRACE
    printf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
#endif
    bitCount[BITS_MB_MODE] += currSE->len;//11//0
//  	printf("\n2、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
   no_bits                += currSE->len;//11//0
    currSE++;
    currMB->currSEnr++;//1
  }
  else if (input->symbol_mode == CABAC)//CABAC熵编码
  {
    
    if (img->MbaffFrameFlag && (img->current_mb_nr%2 == 0||prevMbSkipped))
    {
      mb_field_tmp = currMB->mb_field;
      currMB->mb_field = field_flag_inference();
      CheckAvailabilityOfNeighborsCABAC();
      currMB->mb_field = mb_field_tmp;
    }
    
    // write mb_skip_flag
    mb_type         = MBType2Value (currMB);
    currSE->value1  = mb_type;
    currSE->value2  = currMB->cbp;
    currSE->type    = SE_MBTYPE;
    currSE->writing = writeMB_skip_flagInfo_CABAC;
    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    if (img->type == B_SLICE)  printf(currSE->tracestring, TRACESTRING_SIZE, "B_MB skipflag(%2d,%2d) = %3d",img->mb_x, img->mb_y, (mb_type!=0 ||currMB->cbp!=0));
    else                       printf(currSE->tracestring, TRACESTRING_SIZE, "MB skipflag(%2d,%2d,%d) = %3d",img->mb_x, img->mb_y, currSE->context,(mb_type!=0));
#endif
    bitCount[BITS_MB_MODE] += currSE->len;
//  	printf("\n3、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
   no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
    
    CheckAvailabilityOfNeighborsCABAC();
			
    // write mb_aff
    if(img->MbaffFrameFlag && !skip) // check for copy mode, Krit
    {
      if(WriteFrameFieldMBInHeader)
      {
        currSE->value1 = currMB->mb_field;
        currSE->value2 = 0;
        currSE->type   =  SE_MBTYPE;
        
        if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
        else                            currSE->writing = writeFieldModeInfo_CABAC;
        
        if( input->symbol_mode==UVLC)
        {
          currSE->bitpattern = (currMB->mb_field ? 1 : 0);
          currSE->len = 1;
          writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
        }
        else
        {
          dataPart->writeSyntaxElement(currSE, dataPart);
        }
#if TRACE
        printf(currSE->tracestring, TRACESTRING_SIZE, "Field mode = %3d",currMB->mb_field);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
//  		printf("\n4、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
      
    // write mb_type
    if (currMB->mb_type != 0 || ((img->type == B_SLICE) && currMB->cbp != 0))
    {
      currSE->value1  = mb_type;
      currSE->value2  = 0;
      currSE->type    = SE_MBTYPE;
      currSE->writing = writeMB_typeInfo_CABAC;
      dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
      if (img->type == B_SLICE)  printf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
      else                     printf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;
//  	  printf("\n5、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
      no_bits                += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
  }

  else if (currMB->mb_type != 0 || ((img->type == B_SLICE) && currMB->cbp != 0))//yes
  {
    //===== Run Length Coding: Non-Skipped macorblock =====
    currSE->value1  = img->cod_counter;
    currSE->value2  = 0;
    currSE->mapping = ue_linfo;
    currSE->type    = SE_MBTYPE;
    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    printf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
    bitCount[BITS_MB_MODE] += currSE->len;
//  	printf("\n6、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
    
    // Reset cod counter
    img->cod_counter = 0;
    
    // write mb_aff
    if(img->MbaffFrameFlag && !skip) //no  // check for copy mode, Krit
    {
      if(WriteFrameFieldMBInHeader)
      {
        currSE->value1 = currMB->mb_field;
        currSE->type   =  SE_MBTYPE;
        currSE->mapping = ue_linfo;
        
        //dataPart->writeSyntaxElement(currSE, dataPart);
        currSE->bitpattern = (currMB->mb_field ? 1 : 0);
        currSE->len = 1;
        writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);

#if TRACE
        printf(currSE->tracestring, TRACESTRING_SIZE, "Field mode = %3d",currMB->mb_field);
#endif
        bitCount[BITS_MB_MODE] += currSE->len;
//  		printf("\n7、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
        no_bits                += currSE->len;
        currSE++;
        currMB->currSEnr++;
      }
    }
    // Put out mb mode
    currSE->value1  = MBType2Value (currMB);//0
    
    if (img->type != B_SLICE)
    {
      currSE->value1--;
    }
    currSE->mapping = ue_linfo;
    currSE->type    = SE_MBTYPE;
    currSE->value2  = 0;
    
    dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
    if (img->type == B_SLICE)   printf(currSE->tracestring, TRACESTRING_SIZE, "B_MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->mb_type);
    else                      printf(currSE->tracestring, TRACESTRING_SIZE,   "MB mode(%2d,%2d) = %3d",img->mb_x, img->mb_y,currMB->mb_type);
#endif
    bitCount[BITS_MB_MODE] += currSE->len;
//  	printf("\n8、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
    no_bits                += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }
  else
  {
    //Run Length Coding: Skipped macroblock
    img->cod_counter++;

    // CAVLC
    for (j=0; j < 6; j++)
      for (i=0; i < 4; i++)
        img->nz_coeff [img->current_mb_nr][i][j]=0;


    if(img->current_mb_nr == img->total_number_mb)
    {
      // Put out run
      currSE->value1  = img->cod_counter;
      currSE->value2  = 0;
      currSE->mapping = ue_linfo;
      currSE->type    = SE_MBTYPE;
      
      dataPart->writeSyntaxElement( currSE, dataPart);
#if TRACE
      printf(currSE->tracestring, TRACESTRING_SIZE, "MB runlength = %3d",img->cod_counter);
#endif
      bitCount[BITS_MB_MODE] += currSE->len;//11
//  	  printf("\n9、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
      no_bits                += currSE->len;//11
      currSE++;
      currMB->currSEnr++;
      
      // Reset cod counter
      img->cod_counter = 0;
    }
  }
  
  //===== BITS FOR 8x8 SUB-PARTITION MODES =====
  if (IS_P8x8 (currMB))//8×8划分子块
  {
    dataPart = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    
    for (i=0; i<4; i++)
    {
      if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
      else                            currSE->writing = writeB8_typeInfo_CABAC;

      currSE->value1  = B8Mode2Value (currMB->b8mode[i], currMB->b8pdir[i]);
      currSE->value2  = 0;
      currSE->type    = SE_MBTYPE;
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      printf(currSE->tracestring, TRACESTRING_SIZE, "8x8 mode/pdir(%2d) = %3d/%d",
        i,currMB->b8mode[i],currMB->b8pdir[i]);
#endif
      bitCount[BITS_MB_MODE]+= currSE->len;
//  	  printf("\n10、currSE->len=%3d,bitcount=%3d",currSE->len,bitCount[BITS_MB_MODE]);
      no_bits               += currSE->len;
      currSE++;
      currMB->currSEnr++;
    }
  }

 //===== BITS FOR INTRA PREDICTION MODES ====
  no_bits += writeIntra4x4Modes(-1);//11//51
  //===== BITS FOR CHROMA INTRA PREDICTION MODE ====
  if (currMB->IntraChromaPredModeFlag)//1
    no_bits += writeChromaIntraPredMode();//12//53
  else if(!rdopt) //GB CHROMA !!!!!
    currMB->c_ipred_mode = DC_PRED_8; //setting c_ipred_mode to default is not the right place here
                                      //resetting in rdopt.c (but where ??)
                                      //with cabac and bframes maybe it could crash without this default
                                      //since cabac needs the right neighborhood for the later MBs

  return no_bits;
}

void write_terminating_bit (short bit)
{
	DataPartition*          dataPart;
	const int*              partMap   = assignSE2partition[input->partition_mode];
	EncodingEnvironmentPtr  eep_dp;
	
	//--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
	dataPart = &(img->currentSlice->partArr[partMap[SE_MBTYPE]]);
	dataPart->bitstream->write_flag = 1;
	eep_dp                          = &(dataPart->ee_cabac);
	
	biari_encode_symbol_final(eep_dp, bit); 
#if TRACE
	fprintf (p_trace, "      CABAC terminating bit = %d\n",bit);
#endif
	
}


/*!
 ************************************************************************
 * \brief
 *    Write chroma intra prediction mode.
 ************************************************************************
 */
int writeChromaIntraPredMode()
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  int*            bitCount  = currMB->bitcounter;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             rate      = 0;
  DataPartition*  dataPart;

  //===== BITS FOR CHROMA INTRA PREDICTION MODES
  if (input->symbol_mode==UVLC)   currSE->mapping = ue_linfo;
  else                            currSE->writing = writeCIPredMode_CABAC;

  currSE->value1 = currMB->c_ipred_mode;
  currSE->value2 = 0;
  currSE->type = SE_INTRAPREDMODE;
  dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);

  dataPart->writeSyntaxElement (currSE, dataPart);
  bitCount[BITS_COEFF_UV_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  printf(currSE->tracestring, TRACESTRING_SIZE, "Chroma intra pred mode");
#endif
  currSE++;
  currMB->currSEnr++;

  return rate;
}





/*!
************************************************************************
* \brief
*    Set global last_dquant according to macroblock delta qp 将当前宏块的delta qp赋给全局变量last_dquant 
************************************************************************
*/

extern int last_dquant;

void set_last_dquant()
{
	Macroblock *currMB = &img->mb_data[img->current_mb_nr];
	
	if ((IS_INTERMV (currMB)  || IS_INTRA (currMB)  ) ||
		((img->type==B_SLICE)     && currMB->cbp != 0)  )
	{
		// non-skip
		last_dquant = currMB->delta_qp;
	}
	else
	{
		// skip
		last_dquant = 0;
	}
}


/*!
************************************************************************
* \brief
*    Passes the chosen syntax elements to the NAL 熵编码句法元素，写入码流
************************************************************************
*/
void write_one_macroblock (int eos_bit)
{
	Macroblock* currMB   = &img->mb_data[img->current_mb_nr];//当前宏块
	int*        bitCount = currMB->bitcounter;//比特计数
	int i,j;
	
	extern int cabac_encoding;
	
	//===== init and update number of intra macroblocks =====
	if (img->current_mb_nr==0)//如果是第一个宏块
		intras=0;
	//if ((img->type==P_SLICE || img->type==SP_SLICE || (img->type==B_SLICE && img->nal_reference_idc>0)) && IS_INTRA(currMB))
	if (IS_INTRA(currMB))
		intras++;
	
	//如果当前宏块不是片中的第一个宏块--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
	if (input->symbol_mode==CABAC && img->current_mb_nr!=img->currentSlice->start_mb_nr && eos_bit)
	{
		write_terminating_bit (0);//熵编码非片结束符号0
	}
	
	cabac_encoding = 1;
	
	//--- write header ---
	writeMBHeader (0); //8//熵编码宏块头信息
	
	//  Do nothing more if copy and inter mode
	if ((IS_INTERMV (currMB)  || IS_INTRA (currMB)  ) ||
		((img->type==B_SLICE)     && currMB->cbp != 0)  )//yes
	{
		writeMotionInfo2NAL  ();//熵编码运动信息
		writeCBPandLumaCoeff ();//熵编码当前宏块的CBP，与亮度系数。
		writeChromaCoeff     ();//熵编码色度系数
	}
	else//no
	{ 
		for (j=0; j < 6; j++)
			for (i=0; i < 4; i++)
				img->nz_coeff [img->current_mb_nr][i][j]=0;  // CAVLC
	}
	
	set_last_dquant();//将当前宏块的delta qp赋给全局变量last_dquant 
	
	//--- constrain intra prediction ---
	if(input->UseConstrainedIntraPred && (img->type==P_SLICE || img->type==B_SLICE))//no
	{
		if( !IS_NEWINTRA( currMB ) && currMB->mb_type!=I4MB )
		{
			img->intra_block[img->current_mb_nr] = 0;
		}
	}
	
	//--- set total bit-counter ---
	bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE] + bitCount[BITS_COEFF_Y_MB]     + bitCount[BITS_INTER_MB]
		+ bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
//  	printf("\n%4d、%3d、%3d、%3d、%3d、%3d、%3d",bitCount[BITS_TOTAL_MB],bitCount[BITS_MB_MODE],bitCount[BITS_COEFF_Y_MB],bitCount[BITS_INTER_MB]
//  		,bitCount[BITS_CBP_MB],bitCount[BITS_DELTA_QUANT_MB],bitCount[BITS_COEFF_UV_MB]);
	//Rate control
	img->NumberofMBHeaderBits=bitCount[BITS_MB_MODE]   + bitCount[BITS_INTER_MB]
		+ bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB];//宏块头的比特数
	img->NumberofMBTextureBits= bitCount[BITS_COEFF_Y_MB]+ bitCount[BITS_COEFF_UV_MB];//残差系数的比特数
	img->NumberofTextureBits +=img->NumberofMBTextureBits;//累加所有宏块的残差系数比特数
	img->NumberofHeaderBits +=img->NumberofMBHeaderBits;//累加所有宏块的宏块头比特数
	/*basic unit layer rate control*/
	if(img->BasicUnit<img->Frame_Total_Number_MB)//no
	{
		img->NumberofBasicUnitHeaderBits +=img->NumberofMBHeaderBits;
		img->NumberofBasicUnitTextureBits +=img->NumberofMBTextureBits;
	}
	/*record the total number of MBs*/
	img->NumberofCodedMacroBlocks++;//已编码的宏块数目
	
	stat->bit_slice += bitCount[BITS_TOTAL_MB];
//  	printf("\nstat->bit_slice=%4d",stat->bit_slice);
	
	cabac_encoding = 0;//把CABAC编码标记cabac_encoding置为0
}


/*!
************************************************************************
* \brief
*    Passes the chosen syntax elements to the NAL把句法元素传递到NAL
************************************************************************
*/
void write_one_macroblock_fract (int eos_bit)//---zl---分形参数的编码---
{
	Macroblock* currMB   = &img->mb_data[img->current_mb_nr];//当前宏块数据
	int*        bitCount = currMB->bitcounter;
	int i,j;
	int mvbit;
	extern int cabac_encoding;
// 	mvbit=0;
	//===== init and update number of intra macroblocks =====
	if (img->current_mb_nr==0)//一帧中第一个宏块
		intras=0;
	//if ((img->type==P_SLICE || img->type==SP_SLICE || (img->type==B_SLICE && img->nal_reference_idc>0)) && IS_INTRA(currMB))
	if (IS_INTRA(currMB))//no
		intras++;
	
	//如果当前宏块不是片中的第一个宏块--- write non-slice termination symbol if the macroblock is not the first one in its slice ---
	if (input->symbol_mode==CABAC && img->current_mb_nr!=img->currentSlice->start_mb_nr && eos_bit)
	{
		write_terminating_bit (0);//对非片结束符号0进行编码
	}
	
	cabac_encoding = 1;
	
	//--- write header ---
	writeMBHeader (0); //写宏块头信息（宏块类型）
	
	//  Do nothing more if copy and inter mode
	if ((IS_INTERMV (currMB)  || IS_INTRA (currMB)  ) ||//---zl---分形参数的编码---
		((img->type==B_SLICE)     && currMB->cbp != 0)  )//yes
	{
// 		mvbit+=
		writeXYAndSOInfo2NAL_fract(1);//熵编码x,y和scale，offset
		writeCBPandLumaCoeff ();//熵编码当前宏块的块编码模式，与亮度残差系数。
		writeChromaCoeff     ();//熵编码色度系数
	}
	else
	{ 
		for (j=0; j < 6; j++)
			for (i=0; i < 4; i++)
				img->nz_coeff [img->current_mb_nr][i][j]=0;  // CAVLC
	}
	
	set_last_dquant();//什么意思???
	
	//--- constrain intra prediction ---
	if(input->UseConstrainedIntraPred && (img->type==P_SLICE || img->type==B_SLICE))//no
	{
		if( !IS_NEWINTRA( currMB ) && currMB->mb_type!=I4MB )
		{
			img->intra_block[img->current_mb_nr] = 0;
		}
	}
	
	//--- set total bit-counter ---计算编码的比特数
	bitCount[BITS_TOTAL_MB] = bitCount[BITS_MB_MODE] + bitCount[BITS_COEFF_Y_MB]     + bitCount[BITS_INTER_MB]
		+ bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB] + bitCount[BITS_COEFF_UV_MB];
	//  	printf("\n%4d、%3d、%3d、%3d、%3d、%3d、%3d",bitCount[BITS_TOTAL_MB],bitCount[BITS_MB_MODE],bitCount[BITS_COEFF_Y_MB],bitCount[BITS_INTER_MB]
	//  		,bitCount[BITS_CBP_MB],bitCount[BITS_DELTA_QUANT_MB],bitCount[BITS_COEFF_UV_MB]);
	//Rate control
	img->NumberofMBHeaderBits=bitCount[BITS_MB_MODE]   + bitCount[BITS_INTER_MB]
		+ bitCount[BITS_CBP_MB]  + bitCount[BITS_DELTA_QUANT_MB];//当前宏块
	img->NumberofMBTextureBits= bitCount[BITS_COEFF_Y_MB]+ bitCount[BITS_COEFF_UV_MB];//当前宏块
	img->NumberofTextureBits +=img->NumberofMBTextureBits;//一帧所有宏块
	img->NumberofHeaderBits +=img->NumberofMBHeaderBits;//一帧所有宏块
	/*basic unit layer rate control*/
	if(img->BasicUnit<img->Frame_Total_Number_MB)//no
	{
		img->NumberofBasicUnitHeaderBits +=img->NumberofMBHeaderBits;
		img->NumberofBasicUnitTextureBits +=img->NumberofMBTextureBits;
	}
	/*record the total number of MBs*/
	img->NumberofCodedMacroBlocks++;//所有已经编码完的宏块个数，包括之前几帧的宏块
	
	stat->bit_slice += bitCount[BITS_TOTAL_MB];
	//  	printf("\nstat->bit_slice=%4d",stat->bit_slice);
	
	cabac_encoding = 0;//把CABAC编码标记cabac_encoding置为0
// 	return mvbit;
}

void Scale_Offset_code(int value/*,int current_mb_nr*/)
{
	int rate;
	DataPartition* dataPart;

// 	img->current_mb_nr=current_mb_nr;
	Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
	SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
	Slice*         currSlice  = img->currentSlice;
	int*           bitCount   = currMB->bitcounter;
	const int*     partMap    = assignSE2partition[input->partition_mode];
	dataPart = &(currSlice->partArr[partMap[currSE->type]]);
	
	currSE->value1 = value;
	currSE->value2 = 0;
	currSE->type   = SE_REFFRAME;
	if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
	{
		currSE->mapping = se_linfo;//有符号0阶指数哥伦布编码
// 		dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
		dataPart->writeSyntaxElement (currSE, dataPart);//语法元素写入码流
	}
// 				else//删除了内容，下同。
// 	bitCount[BITS_INTER_MB] += currSE->len;
// 	rate                    += currSE->len;
// 	currSE++;  
// 	currMB->currSEnr++;

}

int writeScaleAndOffset_new(int con)
{
	int k, j, i, refframe,m,n,avg_8=0;
	int block16_x_n,block16_y_n,block16_x,block16_y;
	int pre_block16_x,pre_block16_y;
	DataPartition* dataPart;
	int Y,U,V,yuv;
	int            rate       = 0;
	int            curr_mvd,curr_offset,curr_scale;
	Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
	SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
	Slice*         currSlice  = img->currentSlice;
	int*           bitCount   = currMB->bitcounter;
	const int*     partMap    = assignSE2partition[input->partition_mode];
	// 	int             no_bits   = 0;
	int list_idx;
	dataPart = &(currSlice->partArr[partMap[currSE->type]]);
	Y=U=V=0;
	yuv=4;
	if (2==con)//把UV分量存在Y分量的后面了
	{
		U=img->width/4;//88
		// 		currMB->mb_type=1;
// 		yuv=2;
	} 
	if (3==con)
	{
		U=img->width/4;//88
		V=img->height/8;//36
		// 		currMB->mb_type=1;
// 		yuv=2;
	}//else
	block16_x_n=img->width/16/min(con,2);//图像宽度/16
	if ((img->current_mb_nr)!=0)//不是当前帧中第一个宏块
	{
		pre_block16_x=(img->current_mb_nr-1)%block16_x_n;//当前宏块前面一个宏块的x坐标
		pre_block16_y=(img->current_mb_nr-1)/block16_x_n;//当前宏块前面一个宏块的y坐标
	}
	block16_x=img->current_mb_nr%block16_x_n;
	block16_y=img->current_mb_nr/block16_x_n;
	if ((img->current_frame%input->I_frame)==1)//(264里面的参考帧的编码)如果是紧接着I帧后的一个P帧，对此P的scale和offset的预测利用前一块的值
	{
		//是以MB来进行的，先编码offset,和scale一起了吧
		if (img->current_mb_nr==0)//如果是第一个块，第一块减去(10+6*yuv)在编码。
		{
			if (currMB->mb_type==1)//16*16
			{
				curr_scale=enc_picture->ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv];
				curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv]-(10+6*yuv/min(2,con));
				Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
				Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 				currSE->value1 = curr_offset;
// 				currSE->value2 = curr_scale;
// 				currSE->type   = SE_REFFRAME;
// 				if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 				{
// 					currSE->mapping = se_linfo;
// 					dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
// 					dataPart->writeSyntaxElement (currSE, dataPart);
// 				}
// 				// 				else//删除了内容，下同。
// 				bitCount[BITS_INTER_MB] += currSE->len;
// 				rate                    += currSE->len;
// 				currSE++;  
// 				currMB->currSEnr++;
			}else//是第一块，但是不是16*16块,可以改进，但是用处不大吧？？？？？？？
			{
				for (j=0;j<yuv;j++)
				{
					for (i=0;i<yuv;i++)
					{
						curr_scale     = enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j];
						curr_offset    = enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j];
						Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
						Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 						currSE->value1 = curr_offset;
// 						currSE->value2 = curr_scale;
// 						currSE->type   = SE_REFFRAME;
// 						if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 						{
// 							currSE->mapping = ue_linfo;
// 							dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
// 							dataPart->writeSyntaxElement (currSE, dataPart);
// 						}
// 						// 						else
// 						bitCount[BITS_INTER_MB] += currSE->len;
// 						rate                    += currSE->len;
// 						currSE++;  
// 						currMB->currSEnr++;
						
					}	
				}
			}
		}else//是紧接着I帧后的一个P帧;不是第一块的情况
		{
			//////////////////改进改进//////
			if (1 == currMB->mb_type)//16*16
			{//现在是16*16，不管之前是什么就直接减去这个块的第一个值？以后再说吧///////////改进改进//////////待改进地方
				curr_scale     = enc_picture->ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv][V+pre_block16_y*yuv];//[0][0];
				curr_offset    = enc_picture->ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv][V+pre_block16_y*yuv];//[0][0];
				Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
				Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

				
// 				currSE->value1 = curr_offset;
// 				currSE->value2 = curr_scale;
// 				currSE->type   = SE_REFFRAME;
// 				if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，只考虑有一个参考了。
// 				{
// 					currSE->mapping = ue_linfo;
// 					dataPart->writeSyntaxElement (currSE, dataPart);
// 				}
// 				// 				else
// 				bitCount[BITS_INTER_MB] += currSE->len;
// 				rate                    += currSE->len;
// 				currSE++;  
// 				currMB->currSEnr++;
			}else  if (2 == currMB->mb_type)//16*8的情况
			{
				for (j=0,i=0;i<yuv;i+yuv/2)
				{
					curr_scale     = enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+i][V+pre_block16_y*yuv+j];//[0][0];//
					curr_offset    = enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+i][V+pre_block16_y*yuv+j];//[0][0];//
					Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
					Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

					
// 					currSE->value1 = curr_offset;
// 					currSE->value2 = curr_scale;
// 					currSE->type   = SE_REFFRAME;
// 					if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，只考虑有一个参考了。
// 					{
// 						currSE->mapping = ue_linfo;
// 						dataPart->writeSyntaxElement (currSE, dataPart);
// 					}
// 					// 						else
// 					bitCount[BITS_INTER_MB] += currSE->len;
// 					rate                    += currSE->len;
// 					currSE++;  
// 					currMB->currSEnr++;
				}
				// 				}
			}else  if (3==currMB->mb_type)//8*16
			{
				for (j=0,i=0;j<yuv;j+yuv/2)
				{
					curr_scale     = enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+i][V+pre_block16_y*yuv+j];//[0][0];//
					curr_offset    = enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+i][V+pre_block16_y*yuv+j];//[0][0];//
					Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
					Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 					
// 					currSE->value1 = curr_offset;
// 					currSE->value2 = curr_scale;
// 					currSE->type   = SE_REFFRAME;
// 					if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 					{
// 						currSE->mapping = ue_linfo;
// 						dataPart->writeSyntaxElement (currSE, dataPart);
// 					}
// 					// 						else
// 					bitCount[BITS_INTER_MB] += currSE->len;
// 					rate                    += currSE->len;
// 					currSE++;  
// 					currMB->currSEnr++;
				}
				
			}else if(8 == currMB->mb_type)//8*8
			{
				for(m=0;m<2;m++)//列,2个
					for(n=0;n<2;n++)//行，2个
					{
						k=2*m+n;
						if (4==currMB->b8mode[k])//8*8,一共四个要编码
						{
							curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n][V+block16_y*yuv+m]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+n][V+pre_block16_y*yuv+m];//[0][0];//
							curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n][V+block16_y*yuv+m]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+n][V+pre_block16_y*yuv+m];//[0][0];//avg_8;//
							Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
							Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 							
// 							currSE->value1 = curr_offset;
// 							currSE->value2 = curr_scale;
// 							currSE->type   = SE_REFFRAME;
// 							if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 							{
// 								currSE->mapping = ue_linfo;
// 								dataPart->writeSyntaxElement (currSE, dataPart);
// 							}
// 							// 						else
// 							bitCount[BITS_INTER_MB] += currSE->len;
// 							rate                    += currSE->len;
// 							currSE++;  
// 							currMB->currSEnr++;
						}//8*8
						if (5==currMB->b8mode[k])//8*4
						{
							for (i=0;i<2;i++)
							{
								curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+n*2+i][V+pre_block16_y*yuv+m*2];//[0][0];//
								curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+n*2+i][V+pre_block16_y*yuv+m*2];//[0][0];//-avg_8;
								Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
								Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

								
// 								currSE->value1 = curr_offset;
// 								currSE->value2 = curr_scale;
// 								currSE->type   = SE_REFFRAME;
// 								if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 								{
// 									currSE->mapping = ue_linfo;
// 									dataPart->writeSyntaxElement (currSE, dataPart);
// 								}
// 								// 						else
// 								bitCount[BITS_INTER_MB] += currSE->len;
// 								rate                    += currSE->len;
// 								currSE++;  
// 								currMB->currSEnr++;
							}
				
						}//8*4
						if (6==currMB->b8mode[k])//4*8
						{
							for (j=0;j<2;j++)
							{
								curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2+j]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+n*2][V+pre_block16_y*yuv+m*2+j];//[0][0];//
								curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2+j]/*-avg_8;*/-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+n*2][V+pre_block16_y*yuv+m*2+j];//[0][0];//
								Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
								Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

								
// 								currSE->value1 = curr_offset;
// 								currSE->value2 = curr_scale;
// 								currSE->type   = SE_REFFRAME;
// 								if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 								{
// 									currSE->mapping = ue_linfo;
// 									dataPart->writeSyntaxElement (currSE, dataPart);
// 								}
// 								// 						else
// 								bitCount[BITS_INTER_MB] += currSE->len;
// 								rate                    += currSE->len;
// 								currSE++;  
// 								currMB->currSEnr++;
							}
						}
						if (7==currMB->b8mode[k])//4*4
						{
							for (j=0;j<2;j++)
							{
								for (i=0;i<2;i++)
								{
// 									curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n][V+block16_y*yuv+m]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+n][V+pre_block16_y*yuv+m];//[0][0];//+i+j
// 									curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n][V+block16_y*yuv+m]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+n][V+pre_block16_y*yuv+m];//[0][0];//+i+j-avg_8;//
									
									curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2+j]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+n*2][V+pre_block16_y*yuv+m*2];//[0][0];//+i+j
									curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2+j]-/*avg_8;//*/enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+n*2][V+pre_block16_y*yuv+m*2];//[0][0];//+i+j
									Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
									Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

									
// 									currSE->value1 = curr_offset;
// 									currSE->value2 = curr_scale;
// 									currSE->type   = SE_REFFRAME;
// 									if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 									{
// 										currSE->mapping = ue_linfo;
// 										dataPart->writeSyntaxElement (currSE, dataPart);
// 									}
// 									// 						else
// 									bitCount[BITS_INTER_MB] += currSE->len;
// 									rate                    += currSE->len;
// 									currSE++;  
// 									currMB->currSEnr++;
								}
							}
						}//4*4
					}//全体8*8	
			}//是不是16*16	
		}//是不是第一块
	}//是I帧后的第一帧
	else//不是I帧后的第一帧
	{//不用分是不是第一块了，反正是和前面一帧进行预测的~~~前后帧的块类型是不是一致。还要存储前一帧的块类型。
		if (1==currMB->mb_type)//16*16
		{
			curr_scale=enc_picture->ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv];
			curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv];
			Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
			Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

			
// 			currSE->value1 = curr_offset;
// 			currSE->value2 = curr_scale;
// 			currSE->type   = SE_REFFRAME;
// 			if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 			{
// 				currSE->mapping = ue_linfo;
// 				dataPart->writeSyntaxElement (currSE, dataPart);
// 			}
// 			// 				else
// 			bitCount[BITS_INTER_MB] += currSE->len;
// 			rate                    += currSE->len;
// 			currSE++;  
// 			currMB->currSEnr++;
		} 
		else if (2==currMB->mb_type)//16*8的情况
		{
			for (j=0,i=0;i<yuv;i+yuv/2)
			{
				curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j];//[0][0];//
				curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j];//[0][0];//
				Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
				Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 				currSE->value1 = curr_offset;
// 				currSE->value2 = curr_scale;
// 				currSE->type   = SE_REFFRAME;
// 				if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 				{
// 					currSE->mapping = ue_linfo;
// 					dataPart->writeSyntaxElement (currSE, dataPart);
// 				}
// 				// 						else
// 				bitCount[BITS_INTER_MB] += currSE->len;
// 				rate                    += currSE->len;
// 				currSE++;  
// 				currMB->currSEnr++;
			}
		}else if (3==currMB->mb_type)
		{
			for (j=0,i=0;j<yuv;j+yuv/2)
			{
				curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j];//[0][0];//
				curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j];//[0][0];//
				Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
				Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 				currSE->value1 = curr_offset;
// 				currSE->value2 = curr_scale;
// 				currSE->type   = SE_REFFRAME;
// 				if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 				{
// 					currSE->mapping = ue_linfo;
// 					dataPart->writeSyntaxElement (currSE, dataPart);
// 				}
// 				// 						else
// 				bitCount[BITS_INTER_MB] += currSE->len;
// 				rate                    += currSE->len;
// 				currSE++;  
// 				currMB->currSEnr++;
			}
			
		}else if(8 == currMB->mb_type)//8*8
		{
			// 			for(m=0;m<4;m++)//列
			// 			{	
			// 				for(n=0;n<4;n++)//行
			// 				{
			// 					avg_8+=enc_picture->ref_pic_id[1][U+block16_x*yuv+n][V+block16_y*yuv+m];
			// 				}
			// 			}
			// 			avg_8=avg_8/16;
			for(m=0;m<2;m++)//列,2个
			{	
				for(n=0;n<2;n++)//行，2个
				{
					k=2*m+n;
					if (4==currMB->b8mode[k])//8*8,一共四个要编码
					{
						curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n][V+block16_y*yuv+m]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+n][V+block16_y*yuv+m];//[0][0];//
						curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n][V+block16_y*yuv+m]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+n][V+block16_y*yuv+m];//[0][0];//
						Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
						Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 						currSE->value1 = curr_offset;
// 						currSE->value2 = curr_scale;
// 						currSE->type   = SE_REFFRAME;
// 						if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 						{
// 							currSE->mapping = ue_linfo;
// 							dataPart->writeSyntaxElement (currSE, dataPart);
// 						}
// 						// 						else
// 						bitCount[BITS_INTER_MB] += currSE->len;
// 						rate                    += currSE->len;
// 						currSE++;  
// 						currMB->currSEnr++;
					}//8*8
					if (5==currMB->b8mode[k])//8*4
					{
						for (i=0;i<2;i++)
						{
							curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2];//[0][0];//
							curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2];//[0][0];//
							Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
							Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 							currSE->value1 = curr_offset;
// 							currSE->value2 = curr_scale;
// 							currSE->type   = SE_REFFRAME;
// 							if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 							{
// 								currSE->mapping = ue_linfo;
// 								dataPart->writeSyntaxElement (currSE, dataPart);
// 							}
// 							// 						else
// 							bitCount[BITS_INTER_MB] += currSE->len;
// 							rate                    += currSE->len;
// 							currSE++;  
// 							currMB->currSEnr++;
						}
						
					}//8*4
					if (6==currMB->b8mode[k])//4*8
					{
						for (j=0;j<2;j++)
						{
							curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2+j]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2+j];//[0][0];//
							curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2+j]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2+j];//[0][0];//
							Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
							Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 							currSE->value1 = curr_offset;
// 							currSE->value2 = curr_scale;
// 							currSE->type   = SE_REFFRAME;
// 							if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 							{
// 								currSE->mapping = ue_linfo;
// 								dataPart->writeSyntaxElement (currSE, dataPart);
// 							}
// 							// 						else
// 							bitCount[BITS_INTER_MB] += currSE->len;
// 							rate                    += currSE->len;
// 							currSE++;  
// 							currMB->currSEnr++;
						}
					}
					if (7==currMB->b8mode[k])//4*4
					{
						for (j=0;j<2;j++)
						{
							for (i=0;i<2;i++)
							{
								curr_scale =enc_picture->ref_pic_id[0][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2+j]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2];//[0][0];//+i+j
								curr_offset =enc_picture->ref_pic_id[1][U+block16_x*yuv+n*2+i][V+block16_y*yuv+m*2+j]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+n*2][V+block16_y*yuv+m*2];//[0][0];//+i+j
								Scale_Offset_code(curr_scale/*,int current_mb_nr*/);
								Scale_Offset_code(curr_offset/*,int current_mb_nr*/);

// 								currSE->value1 = curr_offset;
// 								currSE->value2 = curr_scale;
// 								currSE->type   = SE_REFFRAME;
// 								if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
// 								{
// 									currSE->mapping = ue_linfo;
// 									dataPart->writeSyntaxElement (currSE, dataPart);
// 								}
// 								// 						else
// 								bitCount[BITS_INTER_MB] += currSE->len;
// 								rate                    += currSE->len;
// 								currSE++;  
// 								currMB->currSEnr++;
							}
						}
					}//4*4
				}//行全体8*8
			}//列全体8*8
		}//mb_type种类
	}//是不是I帧后第一个P
	return rate;
}


// /*!
//  ************************************************************************
//  * \brief
//  *    Sets context for reference frame parameter
//  ************************************************************************
//  */
// int BType2CtxRef (int btype)
// {
//   if (btype<4)   return 0;
//   else           return 1;
// }
// 

int writeScaleAndOffset(int con)
{
	int k, j, i, refframe;
	int block16_x_n,block16_y_n,block16_x,block16_y;
	int pre_block16_x,pre_block16_y;
	DataPartition* dataPart;
	int Y,U,V,yuv;
	int            rate       = 0;
	int            curr_mvd,curr_offset,curr_scale;
	Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
	SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
	Slice*         currSlice  = img->currentSlice;
	int*           bitCount   = currMB->bitcounter;
	const int*     partMap    = assignSE2partition[input->partition_mode];
	// 	int             no_bits   = 0;
	int list_idx;
	dataPart = &(currSlice->partArr[partMap[currSE->type]]);
	Y=U=V=0;
	yuv=4;
	if (2==con)//把UV分量存在Y分量的后面了
	{
		U=img->width/4;//88
		yuv=2;
	} 
	if (3==con)
	{
		U=img->width/4;//88
		V=img->height/8;//36
		yuv=2;
	}//else
	block16_x_n=img->width/16;
	if ((img->current_mb_nr)!=0)
	{
		pre_block16_x=(img->current_mb_nr-1)%block16_x_n;//21
		pre_block16_y=(img->current_mb_nr-1)/block16_x_n;//0
	}
	block16_x=img->current_mb_nr%block16_x_n;
	block16_y=img->current_mb_nr/block16_x_n;
	if ((img->current_frame%input->I_frame)==1)//(264里面的参考帧的编码)如果是紧接着I帧后的一个P帧，对此P的scale和offset的预测利用前一块的值
	{
		//是以MB来进行的，先编码offset,和scale一起了吧
		if (img->current_mb_nr==0)//如果是第一个块，第一块减去35在编码。
		{
			if (currMB->mb_type==1)//16*16
			{
				curr_offset=enc_picture->ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv];
				curr_scale =enc_picture->ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv]-(10+6*yuv);
				currSE->value1 = curr_offset;
				currSE->value2 = curr_scale;
				currSE->type   = SE_REFFRAME;
				if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
				{
					currSE->mapping = ue_linfo;
					dataPart->writeSyntaxElement (currSE, dataPart);
				}
// 				else//删除了内容，下同。
				bitCount[BITS_INTER_MB] += currSE->len;
				rate                    += currSE->len;
				currSE++;  
				currMB->currSEnr++;
			}else//不是16*16块
			{
				for (j=0;j<yuv;j++)
				{
					for (i=0;i<yuv;i++)
					{
						curr_offset=enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j];
						curr_scale =enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j];
						currSE->value1 = curr_offset;
						currSE->value2 = curr_scale;
						currSE->type   = SE_REFFRAME;
						if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
						{
							currSE->mapping = ue_linfo;
							dataPart->writeSyntaxElement (currSE, dataPart);
						}
// 						else
						bitCount[BITS_INTER_MB] += currSE->len;
						rate                    += currSE->len;
						currSE++;  
						currMB->currSEnr++;
						
					}	
				}
			}
		}else//不是第一块的情况，是I帧后紧跟的P帧
		{
//////////////////改进改进//////
			if (currMB->mb_type==1)//16*16
			{//现在是16*16，不管之前是什么就直接减去这个块的第一个值？以后再说吧///////////改进改进////////////////这样效果好点，少了将近400bit吧，以后可以算一个平均值。待改进地方
				curr_offset=enc_picture->ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv][V+pre_block16_y*yuv];//[0][0];
	            curr_scale =enc_picture->ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv][V+pre_block16_y*yuv];//[0][0];
				currSE->value1 = curr_offset;
				currSE->value2 = curr_scale;
				currSE->type   = SE_REFFRAME;
				if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
				{
					currSE->mapping = ue_linfo;
					dataPart->writeSyntaxElement (currSE, dataPart);
				}
// 				else
				bitCount[BITS_INTER_MB] += currSE->len;
				rate                    += currSE->len;
				currSE++;  
				currMB->currSEnr++;
			}else//不是16*16的情况
			{
				for (j=0;j<yuv;j++)
				{
					for (i=0;i<yuv;i++)
					{
						curr_offset=enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->ref_pic_id[0][U+pre_block16_x*yuv+i][V+pre_block16_y*yuv+j];//[0][0];//
						curr_scale =enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->ref_pic_id[1][U+pre_block16_x*yuv+i][V+pre_block16_y*yuv+j];//[0][0];//
						currSE->value1 = curr_offset;
						currSE->value2 = curr_scale;
						currSE->type   = SE_REFFRAME;
						if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
						{
							currSE->mapping = ue_linfo;
							dataPart->writeSyntaxElement (currSE, dataPart);
						}
// 						else
						bitCount[BITS_INTER_MB] += currSE->len;
						rate                    += currSE->len;
						currSE++;  
						currMB->currSEnr++;
					}
				}
			}//是不是16*16	
		}//是不是第一块
	}//是I帧后的第一帧
	else//不是I帧后的第一帧
	{//不用分是不是第一块了，反正是和前面一帧进行预测的~~~前后帧的块类型是不是一致。还要存储前一帧的块类型。
		if (currMB->mb_type==1)//16*16
		{
			curr_offset=enc_picture->ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv][V+block16_y*yuv];
			curr_scale =enc_picture->ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv][V+block16_y*yuv];
			currSE->value1 = curr_offset;
			currSE->value2 = curr_scale;
			currSE->type   = SE_REFFRAME;
			if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
			{
				currSE->mapping = ue_linfo;
				dataPart->writeSyntaxElement (currSE, dataPart);
			}
			// 				else
			bitCount[BITS_INTER_MB] += currSE->len;
			rate                    += currSE->len;
			currSE++;  
			currMB->currSEnr++;
		} 
		else
		{
			for (j=0;j<yuv;j++)
			{
				for (i=0;i<yuv;i++)
				{
					curr_offset=enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->pre_ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j];
					curr_scale =enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]-enc_picture->pre_ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j];
					currSE->value1 = curr_offset;
					currSE->value2 = curr_scale;
					currSE->type   = SE_REFFRAME;
					if (input->symbol_mode == UVLC)//和264有差别，去掉了前后参考帧的差别吧，值考虑有一个参考了。
					{
						currSE->mapping = ue_linfo;
						dataPart->writeSyntaxElement (currSE, dataPart);
					}
					// 						else
					bitCount[BITS_INTER_MB] += currSE->len;
					rate                    += currSE->len;
					currSE++;  
					currMB->currSEnr++;
				}
			}
		}
	}
	return rate;
}
int   writeXAndY (int con)
{
	int k, j, i, refframe;
	int block16_x_n,block16_y_n,block16_x,block16_y;
	int block16_x_pre,block16_y_pre;
	int Y,U,V;
	int yuv;
	
	DataPartition* dataPart;
	
	int            rate       = 0;
	int            curr_mvd,curr_offset;
	Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
	SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
	Slice*         currSlice  = img->currentSlice;
	int*           bitCount   = currMB->bitcounter;
	const int*     partMap    = assignSE2partition[input->partition_mode];
	// 	int             no_bits   = 0;
	int   bframe          = (img->type==B_SLICE);
	int   step_h0         = (input->blc_size[ IS_P8x8(currMB) ? 4 : currMB->mb_type ][0] >> 2);
	int   step_v0         = (input->blc_size[ IS_P8x8(currMB) ? 4 : currMB->mb_type ][1] >> 2);
	int   list_idx;
	block16_x_n=img->width/16;
	block16_y_n=img->height/16;
	block16_x=img->current_mb_nr%block16_x_n;
	block16_y=img->current_mb_nr/block16_x_n;
	block16_x_pre=(img->current_mb_nr-1)%block16_x_n;
	block16_y_pre=(img->current_mb_nr-1)/block16_x_n;
	Y=U=V=0;
	yuv=4;
	if (2==con)//把UV分量存在Y分量的后面了
	{
		// 		con=1;
		U=img->width/4;//88
		yuv=2;
	} 
	if (3==con)
	{
		// 		con=1;
		U=img->width/4;//88
		V=img->height/8;//36
		yuv=2;
	}//else
	list_idx=0;//我自己设定的，因为把分形的参考帧都设置为0了,后面也不会进入到含有这个变量的循环里
	if (IS_INTERMV (currMB))///??????有什么用途？？？？
	{
		if (currMB->mb_type==1)//16*16
		{
			for (k=0;k<2;k++)//这样就是包括X和Y了
			{
				curr_mvd=enc_picture->mv[0][U+block16_x*yuv][V+block16_y*yuv][k];
				
				currSE->value1 = curr_mvd;
				currSE->value2 = 0;
				currSE->type   = SE_MVD;
				if (input->symbol_mode == UVLC)
				{
					currSE->mapping = se_linfo;//有符号0阶指数哥伦布编码
				}
				else
				{
					img->subblock_x = i; // position used for context determination
					img->subblock_y = j; // position used for context determination
					currSE->value2  = 2*k+list_idx; // identifies the component and the direction; only used for context determination
					currSE->writing = writeMVD_CABAC;
				}  
				dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
				dataPart->writeSyntaxElement (currSE, dataPart);//语法元素写入码流
				bitCount[BITS_INTER_MB] += currSE->len;
				rate                    += currSE->len;
				currSE++;  
				currMB->currSEnr++;
			}
		}else
		{
			for (j=0;j<yuv;j++)
			{
				for (i=0;i<yuv;i++)
				{	
					for (k=0;k<2;k++)
					{
						curr_mvd=enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][k];
						
						currSE->value1 = curr_mvd;
						currSE->value2 = 0;
						currSE->type   = SE_MVD;
						if (input->symbol_mode == UVLC)
						{
							currSE->mapping = se_linfo;
						}
						else
						{
							img->subblock_x = i; // position used for context determination
							img->subblock_y = j; // position used for context determination
							currSE->value2  = 2*k+list_idx; // identifies the component and the direction; only used for context determination
							currSE->writing = writeMVD_CABAC;
						}  
						dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
						dataPart->writeSyntaxElement (currSE, dataPart);
						bitCount[BITS_INTER_MB] += currSE->len;
						rate                    += currSE->len;
						currSE++;  
						currMB->currSEnr++;
					}
				}
			}
		}
	}
    return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Codes the reference frame
 ************************************************************************
 */
int writeReferenceFrame (int mode, int i, int j, int fwd_flag, int  ref)
{
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  int*            bitCount  = currMB->bitcounter;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             rate      = 0;
  DataPartition*  dataPart;
  int             list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
  int             num_ref   = ( fwd_flag ? listXsize[LIST_0+list_offset]: listXsize[LIST_1+list_offset]);
  int             flag_mode = 0;

  if( num_ref == 1 )
  {
    return 0;
  }
  if ( num_ref == 2 )
  {
    flag_mode = 1;
  }

  currSE->value1 = ref;
  currSE->value2  = 0;
  currSE->type   = SE_REFFRAME;

  dataPart = &(currSlice->partArr[partMap[currSE->type]]);
  if (input->symbol_mode == UVLC)
  {
    if( flag_mode )
    {
      currSE->bitpattern = 1 - currSE->value1;
      currSE->len = 1;
      writeSyntaxElement2Buf_Fixed(currSE, dataPart->bitstream);
    }
    else
    {
      currSE->mapping = ue_linfo;
      dataPart->writeSyntaxElement (currSE, dataPart);
    }
  }
  else
  {
    currSE->context = BType2CtxRef (mode);
    img->subblock_x = i; // position used for context determination
    img->subblock_y = j; // position used for context determination
    currSE->writing = writeRefFrame_CABAC;
    currSE->value2 = (fwd_flag)? LIST_0:LIST_1;
    dataPart->writeSyntaxElement (currSE, dataPart);
  }

  bitCount[BITS_INTER_MB] += currSE->len;
  rate                    += currSE->len;
#if TRACE
  if (fwd_flag)
  {
    printf(currSE->tracestring, TRACESTRING_SIZE, "ref_idx_l0 = %d", currSE->value1);
  }
  else
  {
    printf(currSE->tracestring, TRACESTRING_SIZE, "ref_idx_l1 = %d", currSE->value1);
  }
#endif
  currSE++;
  currMB->currSEnr++;

  return rate;
}



/*!
 ************************************************************************
 * \brief
 *    Writes motion vectors of an 8x8 block
 ************************************************************************
 */
int writeMotionVector8x8 (int  i0,
                          int  j0,
                          int  i1,
                          int  j1,
                          int  refframe,
                          int  list_idx,
                          int  mv_mode)
{
  int            i, j, k, l, m;
  int            curr_mvd;
  DataPartition* dataPart;

  int            rate       = 0;
  int            step_h     = input->blc_size[mv_mode][0] >> 2;
  int            step_v     = input->blc_size[mv_mode][1] >> 2;
  Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
  SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*         currSlice  = img->currentSlice;
  int*           bitCount   = currMB->bitcounter;
  const int*     partMap    = assignSE2partition[input->partition_mode];
  int            refindex   = refframe;

  int******      all_mv     = img->all_mv;
  int******      pred_mv    = img->pred_mv;


  for (j=j0; j<j1; j+=step_v)
  for (i=i0; i<i1; i+=step_h)
  {
    for (k=0; k<2; k++) 
    {
      curr_mvd = all_mv[i][j][list_idx][refindex][mv_mode][k] - pred_mv[i][j][list_idx][refindex][mv_mode][k];


      //--- store (oversampled) mvd ---
      for (l=0; l < step_v; l++)
        for (m=0; m < step_h; m++)
          currMB->mvd[list_idx][j+l][i+m][k] = curr_mvd;

      currSE->value1 = curr_mvd;
      currSE->value2 = 0;
      currSE->type   = SE_MVD;
      if (input->symbol_mode == UVLC)
      {
        currSE->mapping = se_linfo;
      }
      else
      {
        img->subblock_x = i; // position used for context determination
        img->subblock_y = j; // position used for context determination
        currSE->value2  = 2*k+list_idx; // identifies the component and the direction; only used for context determination
        currSE->writing = writeMVD_CABAC;
      }  
      dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
      dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
      if (!list_idx)
      {
        printf(currSE->tracestring, TRACESTRING_SIZE, "mvd_l0 (%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][list_idx][refindex][mv_mode][k], pred_mv[i][j][list_idx][refindex][mv_mode][k]);
      }
      else
      {
        
        printf(currSE->tracestring, TRACESTRING_SIZE, "mvd_l1 (%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][list_idx][refindex][mv_mode][k], pred_mv[i][j][list_idx][refindex][mv_mode][k]);
      }

#endif
      bitCount[BITS_INTER_MB] += currSE->len;
      rate                    += currSE->len;
      currSE++;  
      currMB->currSEnr++;
    }
  }

  return rate;
}

/*!
************************************************************************
* \brief
*    Writes motion vectors of an 8x8 block
************************************************************************
*/
int writeMotionVector8x8_fract (int  i0,
                          int  j0,
                          int  i1,
                          int  j1,
//                           int  refframe,
//                           int  list_idx,
                          int  mv_mode)
{
	int            i, j, k, l, m;
	int            curr_mvd;
	DataPartition* dataPart;
	
	int            rate       = 0;
	int            step_h     = input->blc_size[mv_mode][0] >> 2;
	int            step_v     = input->blc_size[mv_mode][1] >> 2;
	Macroblock*    currMB     = &img->mb_data[img->current_mb_nr];
	SyntaxElement* currSE     = &img->MB_SyntaxElements[currMB->currSEnr];
	Slice*         currSlice  = img->currentSlice;
	int*           bitCount   = currMB->bitcounter;
	const int*     partMap    = assignSE2partition[input->partition_mode];
// 	int            refindex   = refframe;
	
	int******      all_mv     = img->all_mv;
	int******      pred_mv    = img->pred_mv;
	
	
	for (j=j0; j<j1; j+=step_v)
		for (i=i0; i<i1; i+=step_h)
		{
			for (k=0; k<2; k++) 
			{
// 				curr_mvd = all_mv[i][j][list_idx][refindex][mv_mode][k] - pred_mv[i][j][list_idx][refindex][mv_mode][k];
				
				
				//--- store (oversampled) mvd ---
				for (l=0; l < step_v; l++)
					for (m=0; m < step_h; m++)
						currMB->mvd[0][j+l][i+m][k] = curr_mvd;
// 						currMB->mvd[list_idx][j+l][i+m][k] = curr_mvd;
					currSE->value1 = curr_mvd;
					currSE->value2 = 0;
					currSE->type   = SE_MVD;
					if (input->symbol_mode == UVLC)
					{
						currSE->mapping = se_linfo;
					}
					else
					{
// 						img->subblock_x = i; // position used for context determination
// 						img->subblock_y = j; // position used for context determination
// 						currSE->value2  = 2*k+list_idx; // identifies the component and the direction; only used for context determination
						currSE->writing = writeMVD_CABAC;
					}  
					dataPart = &(currSlice->partArr[partMap[SE_MVD]]);
					dataPart->writeSyntaxElement (currSE, dataPart);
#if TRACE
					if (!list_idx)
					{
						printf(currSE->tracestring, TRACESTRING_SIZE, "mvd_l0 (%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][list_idx][refindex][mv_mode][k], pred_mv[i][j][list_idx][refindex][mv_mode][k]);
					}
					else
					{
						
						printf(currSE->tracestring, TRACESTRING_SIZE, "mvd_l1 (%d) = %3d  (org_mv %3d pred_mv %3d)",k, curr_mvd, all_mv[i][j][list_idx][refindex][mv_mode][k], pred_mv[i][j][list_idx][refindex][mv_mode][k]);
					}
					
#endif
					bitCount[BITS_INTER_MB] += currSE->len;
					rate                    += currSE->len;
					currSE++;  
					currMB->currSEnr++;
			}
		}
		
		return rate;
}

/*!
 ************************************************************************
 * \brief
 *    Writes motion info
 ************************************************************************
 */
int writeMotionInfo2NAL ()
{
  int k, j0, i0, refframe;

  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int             no_bits   = 0;
  int   bframe          = (img->type==B_SLICE);
  int   step_h0         = (input->blc_size[ IS_P8x8(currMB) ? 4 : currMB->mb_type ][0] >> 2);
  int   step_v0         = (input->blc_size[ IS_P8x8(currMB) ? 4 : currMB->mb_type ][1] >> 2);

  //=== If multiple ref. frames, write reference frame for the MB ===
  if (IS_INTERMV (currMB))
  {
    // if UVLC is turned on, a 8x8 macroblock with all ref=0 in a P-frame is signalled in macroblock mode
    if (!IS_P8x8 (currMB) || !ZeroRef (currMB) || input->symbol_mode==CABAC || bframe)//不应该进来
    {
      for (j0=0; j0<4; j0+=step_v0)
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=j0+(i0/2);

        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
        {
          no_bits += writeReferenceFrame (currMB->b8mode[k], i0, j0, 1, enc_picture->ref_idx[LIST_0][img->block_x+i0][img->block_y+j0]);
        }
      }
      for (j0=0; j0<4; j0+=step_v0)
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=j0+(i0/2);
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
        {
          no_bits += writeReferenceFrame (currMB->b8mode[k], i0, j0, 0, enc_picture->ref_idx[LIST_1][img->block_x+i0][img->block_y+j0]);
        }
      }
    }
  }

  //===== write forward motion vectors =====
  if (IS_INTERMV (currMB))///??????有什么用途？？？？
  {
    for (j0=0; j0<4; j0+=step_v0)
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=j0+(i0/2);
      if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has forward vector
      {
        refframe  = enc_picture->ref_idx[LIST_0][img->block_x+i0][img->block_y+j0];
        no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, LIST_0, currMB->b8mode[k]);
      }
    }
  }


  //===== write backward motion vectors =====
  if (IS_INTERMV (currMB) && bframe)
  {
    for (j0=0; j0<4; j0+=step_v0)
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=j0+(i0/2);
      if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)//has backward vector
      {
        refframe  = enc_picture->ref_idx[LIST_1][img->block_x+i0][img->block_y+j0];
        no_bits  += writeMotionVector8x8 (i0, j0, i0+step_h0, j0+step_v0, refframe, LIST_1, currMB->b8mode[k]);
      }
    }
  }
  return no_bits;
}

/*!
************************************************************************
* \brief
*    Writes motion info
************************************************************************
*/
void writeXYAndSOInfo2NAL_fract (int con)//这里con没有什么用
{

	//===== write ref(scale and offset) =====
// 	writeScaleAndOffset(1);//从这里改为YUV都写入的。
// 	if (img->current_mb_nr %2 == 0 && img->current_mb_nr !=0)
// 	{
// 	
// 		img->current_mb_nr /=2;
// 		writeScaleAndOffset(2);//从这里改为YUV都写入的。
// 		writeScaleAndOffset(3);//从这里改为YUV都写入的。
// 		img->current_mb_nr *=2;
// 	}

    writeScaleAndOffset_new(1);
	if (img->current_mb_nr %2 == 0 && img->current_mb_nr !=0)
	{
		img->current_mb_nr =img->current_mb_nr/2-1;
        writeScaleAndOffset_new(2);
        writeScaleAndOffset_new(3);
		img->current_mb_nr =(img->current_mb_nr+1)*2;
	}


	//=====     write mv（x and y ）    =====
	writeXAndY(1);
	if (img->current_mb_nr %2 == 0 && img->current_mb_nr !=0)
	{
		img->current_mb_nr =img->current_mb_nr/2-1;
		
//		img->current_mb_nr /=2;
		writeXAndY(2);
		writeXAndY(3);
		img->current_mb_nr =(img->current_mb_nr+1)*2;
	}


}

/*!
 ************************************************************************
 * \brief
 *    Writes chrominance coefficients
 ************************************************************************
 */
int writeChromaCoeff ()
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int*            bitCount  = currMB->bitcounter;
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             cbp       = currMB->cbp;
  DataPartition*  dataPart;

  int   level, run;
  int   i, j, k, uv, mb_x, mb_y, i1, ii, j1, jj;
  int   b8, b4, param;
  int*  ACLevel;
  int*  ACRun;
  int*  DCLevel;
  int*  DCRun;


  //=====
  //=====   D C - C O E F F I C I E N T S
  //=====
  if (cbp > 15)  // check if any chroma bits in coded block pattern is set
  {
    for (uv=0; uv < 2; uv++)
    {

      if (input->symbol_mode == UVLC)
      {
        param = uv;
        rate += writeCoeff4x4_CAVLC (CHROMA_DC, 0, 0, param);//6
          // CAVLC
      }
      else
      {

        DCLevel = img->cofDC[uv+1][0];
        DCRun   = img->cofDC[uv+1][1];

        level=1;
        for (k=0; k < 5 && level != 0; ++k)
        {
          level = currSE->value1 = DCLevel[k]; // level
          run   = currSE->value2 = DCRun  [k]; // run

          if (input->symbol_mode == UVLC)   currSE->mapping = levrun_linfo_c2x2;
          else                              currSE->writing = writeRunLevel_CABAC;

          currSE->context     = CHROMA_DC;//6
          currSE->type        = (IS_INTRA(currMB) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER);
          img->is_intra_block =  IS_INTRA(currMB);
          img->is_v_block     = uv;
    
          // choose the appropriate data partition
          dataPart = &(currSlice->partArr[partMap[currSE->type]]);    
          dataPart->writeSyntaxElement (currSE, dataPart);
          bitCount[BITS_COEFF_UV_MB] += currSE->len;
          rate                       += currSE->len;
#if TRACE
          printf(currSE->tracestring, TRACESTRING_SIZE, "2x2 DC Chroma %2d: level =%3d run =%2d",k, level, run);
#endif
          // proceed to next SE 
          currSE++;  
          currMB->currSEnr++;
        }
      }
    }
  }


  //=====
  //=====   A C - C O E F F I C I E N T S
  //=====
  uv=-1;   
  if (cbp >> 4 == 2) // check if chroma bits in coded block pattern = 10b
  {  
    for (mb_y=4; mb_y < 6; mb_y += 2)
    for (mb_x=0; mb_x < 4; mb_x += 2)
    {
      for (j=mb_y; j < mb_y+2; j++)
      {
        jj=j/2;
        j1=j-4;
        for (i=mb_x; i < mb_x+2; i++)
        {
          b8      = 4 + i/2;
          b4      = 2*(j/5)+ (i%2);

          if (input->symbol_mode == UVLC)
          {
            param = i << 4 | j;
            rate += writeCoeff4x4_CAVLC (CHROMA_AC, b8, b4, param);//7
            // CAVLC
          }
          else
          {

            ACLevel = img->cofAC[b8][b4][0];
            ACRun   = img->cofAC[b8][b4][1];

            ii=i/2;
            i1=i%2;
            level=1;
            uv++;

            img->subblock_y = b4/2;
            img->subblock_x = b4%2;

            for (k=0; k < 16 && level != 0; k++)
            {
              level = currSE->value1 = ACLevel[k]; // level
              run   = currSE->value2 = ACRun  [k]; // run

              if (input->symbol_mode == UVLC)   currSE->mapping = levrun_linfo_inter;
              else                              currSE->writing = writeRunLevel_CABAC;
            
              currSE->context     = CHROMA_AC;
              currSE->type        = (IS_INTRA(currMB) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER);
              img->is_intra_block =  IS_INTRA(currMB);
              img->is_v_block     = (uv>=4);

              // choose the appropriate data partition
              dataPart = &(currSlice->partArr[partMap[currSE->type]]); 
              dataPart->writeSyntaxElement (currSE, dataPart);
              bitCount[BITS_COEFF_UV_MB] += currSE->len;
              rate                       += currSE->len;
#if TRACE
              printf(currSE->tracestring, TRACESTRING_SIZE, "AC Chroma %2d: level =%3d run =%2d",k, level, run);
#endif

              // proceed to next SE 
              currSE++;  
              currMB->currSEnr++;
            }
          }
        }
      }
    }
  }

  return rate;
}

/*!
 ************************************************************************
 * \brief
 *    Writes Luma coeff of an 4x4 block
 ************************************************************************
 */
int writeLumaCoeff4x4_CABAC (int b8, int b4, int intra4x4mode)
{
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int*            bitCount  = currMB->bitcounter;
  DataPartition*  dataPart;

  int   level, run;
  int   k;
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];

  img->subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); // horiz. position for coeff_count context
  img->subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); // vert.  position for coeff_count context

  level=1; // get inside loop
  for(k=0; k<=16 && level !=0; k++)
  {
    level = currSE->value1 = ACLevel[k]; // level
    run   = currSE->value2 = ACRun  [k]; // run
      
    currSE->writing = writeRunLevel_CABAC;

    currSE->context     = LUMA_4x4;
    currSE->type        = (k==0 ? (intra4x4mode?SE_LUM_DC_INTRA:SE_LUM_DC_INTER) : (intra4x4mode?SE_LUM_AC_INTRA:SE_LUM_AC_INTER));
    img->is_intra_block = intra4x4mode;

    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);          
    dataPart->writeSyntaxElement (currSE, dataPart);
    bitCount[BITS_COEFF_Y_MB] += currSE->len;
    rate                      += currSE->len;
#if TRACE
    printf(currSE->tracestring, TRACESTRING_SIZE, "Luma sng(%2d) level =%3d run =%2d", k, level,run);
#endif
    /* proceed to next SE */
    currSE++;  
    currMB->currSEnr++;
  }

  return rate;
}
/*!
 ************************************************************************
 * \brief
 *    Writes Luma Coeff of an 8x8 block
 ************************************************************************
 */
int writeLumaCoeff8x8 (int block8x8, int intra4x4mode)
{
  int  block4x4, rate = 0;

  for (block4x4=0; block4x4<4; block4x4++)
  {
    if (input->symbol_mode == UVLC )
      rate += writeCoeff4x4_CAVLC (LUMA, block8x8, block4x4, 0);// CAVLC
    else
      rate += writeLumaCoeff4x4_CABAC (block8x8, block4x4, intra4x4mode);
  }

  return rate;
}

/*!
 ************************************************************************
 * \brief
 *    Writes CBP, DQUANT, and Luma Coefficients of an macroblock
 ************************************************************************
 */
int writeCBPandLumaCoeff ()
{
  int             mb_x, mb_y, i, j, k;
  int             level, run;
  int             rate      = 0;
  Macroblock*     currMB    = &img->mb_data[img->current_mb_nr];
  int*            bitCount  = currMB->bitcounter;
  SyntaxElement*  currSE    = &img->MB_SyntaxElements[currMB->currSEnr];//将句法元素结构体与img对应起来
  Slice*          currSlice = img->currentSlice;
  const int*      partMap   = assignSE2partition[input->partition_mode];
  int             cbp       = currMB->cbp;
  DataPartition*  dataPart;

  int   b8, b4;
  int*  DCLevel = img->cofDC[0][0];
  int*  DCRun   = img->cofDC[0][1];
  int*  ACLevel;
  int*  ACRun;

  if (!IS_NEWINTRA (currMB))//如果不是I16MB
  {
    //=====   C B P   =====
    //---------------------
    currSE->value1 = cbp;
    
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB)
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_intra;
      currSE->type = SE_CBP_INTRA;
    }
    else
    {
      if (input->symbol_mode == UVLC)  currSE->mapping = cbp_linfo_inter;
      currSE->type = SE_CBP_INTER;
    }
    if (input->symbol_mode == CABAC)   currSE->writing = writeCBP_CABAC;
                      
    // choose the appropriate data partition
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    
    dataPart->writeSyntaxElement(currSE, dataPart);
    bitCount[BITS_CBP_MB] += currSE->len;
    rate                  += currSE->len;
#if TRACE
    printf(currSE->tracestring, TRACESTRING_SIZE, "CBP (%2d,%2d) = %3d",img->mb_x, img->mb_y, cbp);
#endif
    // proceed to next SE
    currSE++;  
    currMB->currSEnr++;
  }

  //=====   DQUANT   =====//编码语法元素deltaqp
  //----------------------
  if (cbp!=0 || IS_NEWINTRA (currMB))//如果CBP≠0或者是I16MB
  {
    currSE->value1 = currMB->delta_qp;//deltaqp

    if (input->symbol_mode==UVLC)   currSE->mapping = se_linfo;//UVLC就等同于有符号0阶指数哥伦布编码
    else                            currSE->writing = writeDquant_CABAC;

    if (IS_INTER (currMB))  currSE->type = SE_DELTA_QUANT_INTER;
    else                    currSE->type = SE_DELTA_QUANT_INTRA;


    // choose the appropriate data partition
    dataPart = &(img->currentSlice->partArr[partMap[currSE->type]]);
    dataPart->writeSyntaxElement(  currSE, dataPart);
    bitCount[BITS_DELTA_QUANT_MB] += currSE->len;
    rate                          += currSE->len;
#if TRACE
    printf(currSE->tracestring, TRACESTRING_SIZE, "Delta QP (%2d,%2d) = %3d",img->mb_x, img->mb_y, currMB->delta_qp);
#endif
    // proceed to next SE
    currSE++;
    currMB->currSEnr++;
  }

  for (j=0; j < 6; j++)
    for (i=0; i < 4; i++)
      img->nz_coeff [img->current_mb_nr][i][j]=0;  // CAVLC

  if (!IS_NEWINTRA (currMB))//如果不是I16MB
  {
    //=====  L U M I N A N C E   =====
    //--------------------------------
    for (i=0; i<4; i++)  if (cbp & (1<<i))
    {
        rate += writeLumaCoeff8x8 (i, (currMB->b8mode[i]==IBLOCK));//首先编码器是把一个宏块分成4个block8x8，
		//然后再将每个block8x8分成4个block4x4，对每个block4x4调用writeCoeff4x4_CAVLC（这里我们以主要档次为例）
    }
  }
  else//如果是I16MB
  {
    //=====  L U M I N A N C E   f o r   1 6 x 1 6   =====
    //----------------------------------------------------
    // DC coeffs
    if (input->symbol_mode == UVLC)//熵编码模式=UVLC
    {
      rate += writeCoeff4x4_CAVLC (LUMA_INTRA16x16DC, 0, 0, 0);  // CAVLC熵编码4×4块的残差DC系数
    }
    else//熵编码模式=CABAC
    {
      level=1; // get inside loop
      for (k=0; k<=16 && level!=0; k++)
      {
        level = currSE->value1 = DCLevel[k]; // level
        run   = currSE->value2 = DCRun  [k]; // run

        if (input->symbol_mode == UVLC)
        {
          currSE->mapping = levrun_linfo_inter;
        }
        else
        {
          currSE->writing = writeRunLevel_CABAC;
        }

        currSE->context     = LUMA_16DC;
        currSE->type        = SE_LUM_DC_INTRA;   // element is of type DC
        img->is_intra_block = 1;

        // choose the appropriate data partition
        dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    
        dataPart->writeSyntaxElement (currSE, dataPart);
        bitCount[BITS_COEFF_Y_MB] += currSE->len;
        rate                      += currSE->len;
#if TRACE
        printf(currSE->tracestring, TRACESTRING_SIZE, "DC luma 16x16 sng(%2d) level =%3d run =%2d", k, level, run);
#endif
        // proceed to next SE
        currSE++;
        currMB->currSEnr++;
      }
    }

    // AC coeffs
    if (cbp & 15)
    {
      for (mb_y=0; mb_y < 4; mb_y += 2)
      for (mb_x=0; mb_x < 4; mb_x += 2)
      for (j=mb_y; j < mb_y+2; j++)
      for (i=mb_x; i < mb_x+2; i++)
      {
        b8      = 2*(j/2) + (i/2);
        b4      = 2*(j%2) + (i%2);
        if (input->symbol_mode == UVLC)
        {
          rate += writeCoeff4x4_CAVLC (LUMA_INTRA16x16AC, b8, b4, 0);  // CAVLC
        }
        else
        {
          ACLevel = img->cofAC[b8][b4][0];
          ACRun   = img->cofAC[b8][b4][1];

          img->subblock_y = j;
          img->subblock_x = i;

          level=1; // get inside loop
          for (k=0;k<16 && level !=0;k++)
          {
            level = currSE->value1 = ACLevel[k]; // level
            run   = currSE->value2 = ACRun  [k]; // run

            if (input->symbol_mode == UVLC)
            {
              currSE->mapping = levrun_linfo_inter;
            }
            else
            {
              currSE->writing = writeRunLevel_CABAC;
            }
            currSE->context     = LUMA_16AC;
            currSE->type        = SE_LUM_AC_INTRA;   // element is of type AC
            img->is_intra_block = 1;

            // choose the appropriate data partition
           dataPart = &(currSlice->partArr[partMap[currSE->type]]);

            dataPart->writeSyntaxElement (currSE, dataPart);
            bitCount[BITS_COEFF_Y_MB] += currSE->len;
            rate                      += currSE->len;
#if TRACE
            printf(currSE->tracestring, TRACESTRING_SIZE, "AC luma 16x16 sng(%2d) level =%3d run =%2d", k, level, run);
#endif
            // proceed to next SE
            currSE++;
            currMB->currSEnr++;
          }
        }
      }
    }
  }

  return rate;
}

/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients //从邻居块预测非零系数的数目
 *    
 *    Luma Blocks
 ************************************************************************
 */
int predict_nnz(int i,int j)
{
  PixelPos pix;

  int pred_nnz = 0;
  int cnt      = 0;
  int mb_nr    = img->current_mb_nr;

  // left block//左相邻块
  getLuma4x4Neighbour(mb_nr, i, j, -1, 0, &pix);
/* to be inserted only for dp
  if (pix.available && img->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz = img->nz_coeff [pix.mb_addr ][pix.x][pix.y];
    cnt++;
  }

  // top block//上相邻块
  getLuma4x4Neighbour(mb_nr, i, j, 0, -1, &pix);
/* to be inserted only for dp
  if (pix.available && img->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz += img->nz_coeff [pix.mb_addr ][pix.x][pix.y];
    cnt++;
  }

  if (cnt==2)//左相邻块和上相邻块均可用
  {
    pred_nnz++;
    pred_nnz/=cnt; //取平均值
  }

  return pred_nnz;//返回预测的非零系数的个数
}

/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients 
 *    
 *    Chroma Blocks   
 ************************************************************************
 */
int predict_nnz_chroma(int i,int j)
{
  PixelPos pix;

  int pred_nnz = 0;
  int cnt      =0;
  int mb_nr    = img->current_mb_nr;

  // left block
  getChroma4x4Neighbour(mb_nr, i%2, j-4, -1, 0, &pix);
/*  to be inserted only for dp
  if (pix.available && img->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz = img->nz_coeff [pix.mb_addr ][2 * (i/2) + pix.x][4 + pix.y];
    cnt++;
  }
  
  // top block
  getChroma4x4Neighbour(mb_nr, i%2, j-4, 0, -1, &pix);
/*  to be inserted only for dp
  if (pix.available && img->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz += img->nz_coeff [pix.mb_addr ][2 * (i/2) + pix.x][4 + pix.y];
    cnt++;
  }

  if (cnt==2)
  {
    pred_nnz++;
    pred_nnz/=cnt; 
  }

  return pred_nnz;
}

/*!
 ************************************************************************
 * \brief
 *    Writes coeff of an 4x4 block (CAVLC)//以CAVLC方式熵编码一个4×4块的系数
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */

int writeCoeff4x4_CAVLC (int reference, int b8, int b4, int param)
{
  int           no_bits    = 0;
  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  int           *bitCount  = currMB->bitcounter;
  Slice         *currSlice = img->currentSlice;
  DataPartition *dataPart;
  int           *partMap   = assignSE2partition[input->partition_mode];

  int k,level,run,vlcnum;
  int numcoeff, lastcoeff, numtrailingones; 
  int numones, totzeros, zerosleft, numcoef;
  int numcoeff_vlc;
  int code, level_two_or_higher;
  int dptype = 0, bitcounttype = 0;
  int nnz, max_coeff_num = 0, cdc=0, cac=0;
  int subblock_x, subblock_y;
  char type[15];

  int incVlc[] = {0,3,6,12,24,48,32768};  // maximum vlc = 6


  int*  pLevel = NULL;
  int*  pRun = NULL;

  switch (reference)//1
  {
  case LUMA://0
    max_coeff_num = 16;
    bitcounttype = BITS_COEFF_Y_MB;

    pLevel = img->cofAC[b8][b4][0];
    pRun   = img->cofAC[b8][b4][1];

    sprintf(type, "%s", "Luma");
    if (IS_INTRA (currMB))
    {
      dptype = SE_LUM_AC_INTRA;
    }
    else
    {
      dptype = SE_LUM_AC_INTER;
    }
    break;
  case LUMA_INTRA16x16DC://1
    max_coeff_num = 16;
    bitcounttype = BITS_COEFF_Y_MB;//5

    pLevel = img->cofDC[0][0];//存放Level的值
    pRun   = img->cofDC[0][1];//存放Run的值

    sprintf(type, "%s", "Lum16DC");
    dptype = SE_LUM_DC_INTRA;//7
    break;
  case LUMA_INTRA16x16AC://2
    max_coeff_num = 15;
    bitcounttype = BITS_COEFF_Y_MB;

    pLevel = img->cofAC[b8][b4][0];
    pRun   = img->cofAC[b8][b4][1];

    sprintf(type, "%s", "Lum16AC");
    dptype = SE_LUM_AC_INTRA;
    break;

  case CHROMA_DC://6
    max_coeff_num = 4;
    bitcounttype = BITS_COEFF_UV_MB;
    cdc = 1;//色度的cdc=1,亮度的cdc=0

    pLevel = img->cofDC[param+1][0];
    pRun   = img->cofDC[param+1][1];

    sprintf(type, "%s", "ChrDC");
    if (IS_INTRA (currMB))
    {
      dptype = SE_CHR_DC_INTRA;
    }
    else
    {
      dptype = SE_CHR_DC_INTER;
    }
    break;
  case CHROMA_AC://7
    max_coeff_num = 15;
    bitcounttype = BITS_COEFF_UV_MB;
    cac = 1;

    pLevel = img->cofAC[b8][b4][0];
    pRun   = img->cofAC[b8][b4][1];

    sprintf(type, "%s", "ChrAC");
    if (IS_INTRA (currMB))
    {
      dptype = SE_CHR_AC_INTRA;
    }
    else
    {
      dptype = SE_CHR_AC_INTER;
    }
    break;
  default:
//     error("writeCoeff4x4_CAVLC: Invalid block type", 600);
    break;
  }

  dataPart = &(currSlice->partArr[partMap[dptype]]);

  numcoeff = 0;
  numtrailingones = 0;
  numones = 0;
  lastcoeff = 0;
  totzeros = 0;
  level = 1;

  for(k = 0; (k <= cdc?4:16)&& level !=0; k++)
  {
    level = pLevel[k]; // level
    run   = pRun[k];   // run

    if (level)
    {
      if (run)
        totzeros += run;
      if (abs(level) == 1)
      {
        numtrailingones ++;
        numones ++;
        if (numtrailingones > 3)
        {
          numtrailingones = 3; /* clip to 3 */
        }
      }
      else
      {
        numtrailingones = 0;
      }
      numcoeff ++;
      lastcoeff = k;
    }
  }

   if (!cdc)//!0
  {
    if (!cac)//!0
    {
      // luma
      subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); 
        // horiz. position for coeff_count context
      subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); 
        // vert.  position for coeff_count context
      nnz = predict_nnz(subblock_x,subblock_y);//预测的非零系数的个数
    }
    else
    {
      // chroma AC
      subblock_x = param >> 4;
      subblock_y = param & 15;
      nnz = predict_nnz_chroma(subblock_x,subblock_y);
    }

    img->nz_coeff [img->current_mb_nr ][subblock_x][subblock_y] = numcoeff;//真实非零系数的个数


    if (nnz < 2)//预测的非零系数的个数<2
    {
      numcoeff_vlc = 0;
    }
    else if (nnz < 4)
    {
      numcoeff_vlc = 1;
    }
    else if (nnz < 8)
    {
      numcoeff_vlc = 2;
    }
    else 
    {
      numcoeff_vlc = 3;
    }


  }
  else
  {
    // chroma DC (has its own VLC)
    // numcoeff_vlc not relevant
    numcoeff_vlc = 0;

    subblock_x = param;
    subblock_y = param;
  }

  currSE->type  = dptype;   

  currSE->value1 = numcoeff;//非零系数的个数
  currSE->value2 = numtrailingones;//拖尾系数的个数
  currSE->len = numcoeff_vlc; /* use len to pass vlcnum */

#if TRACE
  printf(currSE->tracestring, 
    TRACESTRING_SIZE, "%s # c & tr.1s(%d,%d) vlc=%d #c=%d #t1=%d",
    type, subblock_x, subblock_y, numcoeff_vlc, numcoeff, numtrailingones);
#endif

  if (!cdc)//1
    writeSyntaxElement_NumCoeffTrailingOnes(currSE, dataPart);//熵编码亮度分量的拖尾系数
  else
    writeSyntaxElement_NumCoeffTrailingOnesChromaDC(currSE, dataPart);

  bitCount[bitcounttype]+=currSE->len;
  no_bits               +=currSE->len;

  // proceed to next SE
  currSE++;
  currMB->currSEnr++;


  if (!numcoeff)
    return no_bits;

  if (numcoeff)
  {
    code = 0;
    for (k = lastcoeff; k > lastcoeff-numtrailingones; k--)
    {
      level = pLevel[k]; // level
      if (abs(level) > 1)
      {
        printf("ERROR: level > 1\n");
        exit(-1);
      }
      code <<= 1;
      if (level < 0)
      {
        code |= 0x1;
      }
    }

    if (numtrailingones)
    {
      currSE->type  = dptype;   

      currSE->value2 = numtrailingones;
      currSE->value1 = code;

#if TRACE
      printf(currSE->tracestring, 
        TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)", 
        type, subblock_x, subblock_y);
#endif

      writeSyntaxElement_VLC (currSE, dataPart);//生成VLC码字写入比特流
      bitCount[bitcounttype]+=currSE->len;
      no_bits               +=currSE->len;

      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }

    // encode levels//编码levels
    level_two_or_higher = 1;
    if (numcoeff > 3 && numtrailingones == 3)
      level_two_or_higher = 0;

    if (numcoeff > 10 && numtrailingones < 3)
      vlcnum = 1;
    else
      vlcnum = 0;

    for (k = lastcoeff - numtrailingones; k >= 0; k--)
    {
      level = pLevel[k]; // level

      currSE->value1 = level;
      currSE->type  = dptype;   
      
#if TRACE
      printf(currSE->tracestring, 
        TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d lev=%3d",
        type, subblock_x, subblock_y, k, vlcnum, level);
#endif
      
      if (level_two_or_higher)
      {
        if (currSE->value1 > 0)
          currSE->value1 --;
        else
          currSE->value1 ++;
        level_two_or_higher = 0;
      }

      //    encode level
      if (vlcnum == 0)
        writeSyntaxElement_Level_VLC1(currSE, dataPart);//熵编码Level，写入比特流
      else
        writeSyntaxElement_Level_VLCN(currSE, vlcnum, dataPart);

      // update VLC table
      if (abs(level)>incVlc[vlcnum])
        vlcnum++;

      if (k == lastcoeff - numtrailingones && abs(level)>3)
        vlcnum = 2;

      bitCount[bitcounttype]+=currSE->len;
      no_bits               +=currSE->len;

      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }


    // encode total zeroes
    if (numcoeff < max_coeff_num)
    {

      currSE->type  = dptype;   
      currSE->value1 = totzeros;

      vlcnum = numcoeff-1;

      currSE->len = vlcnum;

#if TRACE
      printf(currSE->tracestring, 
        TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d totzeros=%3d",
          type, subblock_x, subblock_y, vlcnum, totzeros);
#endif
      if (!cdc)
        writeSyntaxElement_TotalZeros(currSE, dataPart);//熵编码TotalZeros，并写入码流
      else
        writeSyntaxElement_TotalZerosChromaDC(currSE, dataPart);

      bitCount[bitcounttype]+=currSE->len;
      no_bits               +=currSE->len;

      // proceed to next SE
      currSE++;
      currMB->currSEnr++;
    }

    // encode run before each coefficient
    zerosleft = totzeros;
    numcoef = numcoeff;
    for (k = lastcoeff; k >= 0; k--)
    {
      run = pRun[k]; // run

      currSE->value1 = run;
      currSE->type  = dptype;   

      // for last coeff, run is remaining totzeros
      // when zerosleft is zero, remaining coeffs have 0 run
      if (numcoeff <= 1 || !zerosleft)
        break;

      if (numcoef > 1 && zerosleft) 
      {

        vlcnum = zerosleft - 1;
        if (vlcnum > RUNBEFORE_NUM-1)
          vlcnum = RUNBEFORE_NUM-1;

        currSE->len = vlcnum;

#if TRACE
        printf(currSE->tracestring, 
          TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d run=%2d",
            type, subblock_x, subblock_y, k, vlcnum, run);
#endif

        writeSyntaxElement_Run(currSE, dataPart);//熵编码Run，并写入码流

        bitCount[bitcounttype]+=currSE->len;
        no_bits               +=currSE->len;

        zerosleft -= run;
        numcoef --;

        // proceed to next SE
		currSE++;
        currMB->currSEnr++;
      }
    }
  }

  return no_bits;
}






/*!
 ************************************************************************
 * \brief
 *    Find best 16x16 based intra mode
 *
 * \par Input:
 *    Image parameters, pointer to best 16x16 intra mode
 *
 * \par Output:
 *    best 16x16 based SAD
 ************************************************************************/
int find_sad_16x16(int *intra_mode)
{
  int current_intra_sad_2,best_intra_sad2;
  int M1[16][16],M0[4][4][4][4],M3[4],M4[4][4];

  int i,j,k;
  int ii,jj;
  int mb_nr = img->current_mb_nr;
  
  PixelPos up;          //!< pixel position p(0,-1)
  PixelPos left[17];    //!< pixel positions p(-1, -1..15)

  int up_avail, left_avail, left_up_avail;

  for (i=0;i<17;i++)
  {
    getNeighbour(mb_nr, -1 ,  i-1 , 1, &left[i]);
  }
  
  getNeighbour(mb_nr, 0     ,  -1 , 1, &up);

  if (!(input->UseConstrainedIntraPred))
  {
    up_avail   = up.available;
    left_avail = left[1].available;
    left_up_avail = left[0].available;
  }
  else
  {
    up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
    for (i=1, left_avail=1; i<17;i++)
      left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
    left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
  }

  best_intra_sad2=MAX_VALUE;

  for (k=0;k<4;k++)
  {
    //check if there are neighbours to predict from
    if ((k==0 && !up_avail) || (k==1 && !left_avail) || (k==3 && (!left_avail || !up_avail || !left_up_avail)))
    {
      ; // edge, do nothing
    }
    else
    {
      for (j=0;j<16;j++)
      {
        for (i=0;i<16;i++)
        {
          M1[i][j]=imgY_org[img->opix_y+j][img->opix_x+i]-img->mprr_2[k][j][i];
          M0[i%4][i/4][j%4][j/4]=M1[i][j];
        }
      }
      current_intra_sad_2=0;              // no SAD start handicap here
      for (jj=0;jj<4;jj++)
      {
        for (ii=0;ii<4;ii++)
        {
          for (j=0;j<4;j++)
          {
            M3[0]=M0[0][ii][j][jj]+M0[3][ii][j][jj];
            M3[1]=M0[1][ii][j][jj]+M0[2][ii][j][jj];
            M3[2]=M0[1][ii][j][jj]-M0[2][ii][j][jj];
            M3[3]=M0[0][ii][j][jj]-M0[3][ii][j][jj];

            M0[0][ii][j][jj]=M3[0]+M3[1];
            M0[2][ii][j][jj]=M3[0]-M3[1];
            M0[1][ii][j][jj]=M3[2]+M3[3];
            M0[3][ii][j][jj]=M3[3]-M3[2];
          }

          for (i=0;i<4;i++)
          {
            M3[0]=M0[i][ii][0][jj]+M0[i][ii][3][jj];
            M3[1]=M0[i][ii][1][jj]+M0[i][ii][2][jj];
            M3[2]=M0[i][ii][1][jj]-M0[i][ii][2][jj];
            M3[3]=M0[i][ii][0][jj]-M0[i][ii][3][jj];

            M0[i][ii][0][jj]=M3[0]+M3[1];
            M0[i][ii][2][jj]=M3[0]-M3[1];
            M0[i][ii][1][jj]=M3[2]+M3[3];
            M0[i][ii][3][jj]=M3[3]-M3[2];
            for (j=0;j<4;j++)
              if ((i+j)!=0)
                current_intra_sad_2 += abs(M0[i][ii][j][jj]);
          }
        }
      }

      for (j=0;j<4;j++)
        for (i=0;i<4;i++)
          M4[i][j]=M0[0][i][0][j]/4;

        // Hadamard of DC koeff
        for (j=0;j<4;j++)
        {
          M3[0]=M4[0][j]+M4[3][j];
          M3[1]=M4[1][j]+M4[2][j];
          M3[2]=M4[1][j]-M4[2][j];
          M3[3]=M4[0][j]-M4[3][j];

          M4[0][j]=M3[0]+M3[1];
          M4[2][j]=M3[0]-M3[1];
          M4[1][j]=M3[2]+M3[3];
          M4[3][j]=M3[3]-M3[2];
        }

        for (i=0;i<4;i++)
        {
          M3[0]=M4[i][0]+M4[i][3];
          M3[1]=M4[i][1]+M4[i][2];
          M3[2]=M4[i][1]-M4[i][2];
          M3[3]=M4[i][0]-M4[i][3];

          M4[i][0]=M3[0]+M3[1];
          M4[i][2]=M3[0]-M3[1];
          M4[i][1]=M3[2]+M3[3];
          M4[i][3]=M3[3]-M3[2];

          for (j=0;j<4;j++)
            current_intra_sad_2 += abs(M4[i][j]);
        }
        if(current_intra_sad_2 < best_intra_sad2)
        {
          best_intra_sad2=current_intra_sad_2;
          *intra_mode = k; // update best intra mode

        }
    }
  }
  best_intra_sad2 = best_intra_sad2/2;

  return best_intra_sad2;

}
/////////////////////////////////////>

static void SetMotionVectorPredictor (struct img_par  *img,
                                      int             *pmv_x,
                                      int             *pmv_y,
                                      int             ref_frame,
                                      int             list,
                                      int             ***refPic,
                                      int             ****tmp_mv,
                                      int             block_x,
                                      int             block_y,
                                      int             blockshape_x,
                                      int             blockshape_y);


/*!
 ************************************************************************
 * \brief
 *    initializes the current macroblock
 ************************************************************************
 */
void start_macroblock_dec(struct img_par *img,struct inp_par_dec *inp, int CurrentMBInScanOrder)
{
	int i,j,k,l;
	Macroblock *currMB;   // intialization code deleted, see below, StW
	
	assert (img->current_mb_nr>=0 && img->current_mb_nr<img->PicSizeInMbs);
	
	currMB = &img->mb_data[img->current_mb_nr];
	
	/* Update coordinates of the current macroblock */
	if (img->MbaffFrameFlag)
	{
		img->mb_x = (img->current_mb_nr)%((2*img->width)/MB_BLOCK_SIZE);
		img->mb_y = 2*((img->current_mb_nr)/((2*img->width)/MB_BLOCK_SIZE));
		
		if (img->mb_x % 2)
		{
			img->mb_y++;
		}
		
		img->mb_x /= 2;
	}
	else
	{
		img->mb_x = (img->current_mb_nr)%(img->width/MB_BLOCK_SIZE);
		img->mb_y = (img->current_mb_nr)/(img->width/MB_BLOCK_SIZE);
	}
	
	/* Define vertical positions */
	img->block_y = img->mb_y * BLOCK_SIZE;      /* luma block position */
	img->pix_y   = img->mb_y * MB_BLOCK_SIZE;   /* luma macroblock position */
	img->pix_c_y = img->mb_y * MB_BLOCK_SIZE/2; /* chroma macroblock position */
	
	/* Define horizontal positions */
	img->block_x = img->mb_x * BLOCK_SIZE;      /* luma block position */
	img->pix_x   = img->mb_x * MB_BLOCK_SIZE;   /* luma pixel position */
	img->pix_c_x = img->mb_x * MB_BLOCK_SIZE/2; /* chroma pixel position */
	
	// Save the slice number of this macroblock. When the macroblock below
	// is coded it will use this to decide if prediction for above is possible
	currMB->slice_nr = img->current_slice_nr;
	
	if (img->current_slice_nr >= MAX_NUM_SLICES)
	{
		error ("maximum number of supported slices exceeded, please recompile with increased value for MAX_NUM_SLICES", 200);
	}
	
	dec_picture->slice_id[img->mb_x][img->mb_y] = img->current_slice_nr;
	if (img->current_slice_nr > dec_picture->max_slice_id)
	{
		dec_picture->max_slice_id=img->current_slice_nr;
	}
	
	CheckAvailabilityOfNeighbors(img);
	
	// Reset syntax element entries in MB struct
	currMB->qp          = img->qp ;
	currMB->mb_type     = 0;
	currMB->delta_quant = 0;
	currMB->cbp         = 0;
	currMB->cbp_blk     = 0;
	currMB->c_ipred_mode= DC_PRED_8; //GB
	
	for (l=0; l < 2; l++)
		for (j=0; j < BLOCK_MULTIPLE; j++)
			for (i=0; i < BLOCK_MULTIPLE; i++)
				for (k=0; k < 2; k++)
					currMB->mvd[l][j][i][k] = 0;
				
				currMB->cbp_bits   = 0;
				
				// initialize img->m7 for ABT
				for (j=0; j<MB_BLOCK_SIZE; j++)
					for (i=0; i<MB_BLOCK_SIZE; i++)
						img->m7[i][j] = 0;
					
					// store filtering parameters for this MB 
					currMB->LFDisableIdc = img->currentSlice->LFDisableIdc;
					currMB->LFAlphaC0Offset = img->currentSlice->LFAlphaC0Offset;
					currMB->LFBetaOffset = img->currentSlice->LFBetaOffset;
					
}

/*!
************************************************************************
* \brief
*    set coordinates of the next macroblock
*    check end_of_slice condition 
************************************************************************
*/
int exit_macroblock(struct img_par *img,struct inp_par_dec *inp,int eos_bit)
{
	
	//! The if() statement below resembles the original code, which tested 
	//! img->current_mb_nr == img->PicSizeInMbs.  Both is, of course, nonsense
	//! In an error prone environment, one can only be sure to have a new
	//! picture by checking the tr of the next slice header!
	
	// printf ("exit_macroblock: FmoGetLastMBOfPicture %d, img->current_mb_nr %d\n", FmoGetLastMBOfPicture(), img->current_mb_nr);
	img->num_dec_mb++;
	
	if (img->num_dec_mb == img->PicSizeInMbs)
		//  if (img->current_mb_nr == FmoGetLastMBOfPicture(currSlice->structure))
	{
		//thb
		/*
		if (currSlice->next_header != EOS)
		currSlice->next_header = SOP;
		*/
		//the
		assert (nal_startcode_follows (img, inp, eos_bit) == TRUE1);
		return TRUE1;
	}
	// ask for last mb in the slice  UVLC
	else
	{
		// printf ("exit_macroblock: Slice %d old MB %d, now using MB %d\n", img->current_slice_nr, img->current_mb_nr, FmoGetNextMBNr (img->current_mb_nr));
		
		img->current_mb_nr = FmoGetNextMBNr_dec (img->current_mb_nr);
		
		if (img->current_mb_nr == -1)     // End of Slice group, MUST be end of slice
		{
			assert (nal_startcode_follows (img, inp, eos_bit) == TRUE1);
			return TRUE1;
		}
		
		if(nal_startcode_follows(img, inp, eos_bit) == FALSE1) 
			return FALSE1;
		
		if(img->type == I_SLICE  || img->type == SI_SLICE || active_pps->entropy_coding_mode_flag == CABAC)
			return TRUE1;
		if(img->cod_counter<=0)
			return TRUE1;
		return FALSE1;
	}
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for P-Frames
 ************************************************************************
 */
void interpret_mb_mode_P(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int         mbmode = currMB->mb_type;

#define ZERO_P8x8     (mbmode==5)
#define MODE_IS_P8x8  (mbmode==4 || mbmode==5)
#define MODE_IS_I4x4  (mbmode==6)
#define I16OFFSET     (mbmode-7)
#define MODE_IS_IPCM  (mbmode==31)

  if(mbmode <4)
  {
    currMB->mb_type = mbmode;
    for (i=0;i<4;i++)
    {
      currMB->b8mode[i]   = mbmode;
      currMB->b8pdir[i]   = 0;
    }
  }
  else if(MODE_IS_P8x8)
  {
    currMB->mb_type = P8x8;
    img->allrefzero = ZERO_P8x8;
  }
  else if(MODE_IS_I4x4)
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++)
    {
      currMB->b8mode[i] = IBLOCK;
      currMB->b8pdir[i] = -1;
    }
  }
  else if(MODE_IS_IPCM)
  {
    currMB->mb_type=IPCM;
    
    for (i=0;i<4;i++) 
    {
      currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; 
    }
    currMB->cbp= -1;
    currMB->i16mode = 0;
  }
  else
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(I16OFFSET)>>2];
    currMB->i16mode = (I16OFFSET) & 0x03;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for I-Frames
 ************************************************************************
 */
void interpret_mb_mode_I(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int         mbmode   = currMB->mb_type;

  if (mbmode==0)
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
  }
  else if(mbmode==25)
  {
    currMB->mb_type=IPCM;

    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= -1;
    currMB->i16mode = 0;

  }
  else
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(mbmode-1)>>2];
    currMB->i16mode = (mbmode-1) & 0x03;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for B-Frames
 ************************************************************************
 */
void interpret_mb_mode_B(struct img_par *img)
{
  static const int offset2pdir16x16[12]   = {0, 0, 1, 2, 0,0,0,0,0,0,0,0};
  static const int offset2pdir16x8[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},{1,0},
                                             {0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2},{0,0}};
  static const int offset2pdir8x16[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},
                                             {1,0},{0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2}};

  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int i, mbmode;
  int mbtype  = currMB->mb_type;
  int *b8mode = currMB->b8mode;
  int *b8pdir = currMB->b8pdir;

  //--- set mbtype, b8type, and b8pdir ---
  if (mbtype==0)       // direct
  {
      mbmode=0;       for(i=0;i<4;i++) {b8mode[i]=0;          b8pdir[i]=2; }
  }
  else if (mbtype==23) // intra4x4
  {
    mbmode=I4MB;    for(i=0;i<4;i++) {b8mode[i]=IBLOCK;     b8pdir[i]=-1; }
  }
  else if ((mbtype>23) && (mbtype<48) ) // intra16x16
  {
    mbmode=I16MB;   for(i=0;i<4;i++) {b8mode[i]=0;          b8pdir[i]=-1; }
    currMB->cbp     = ICBPTAB[(mbtype-24)>>2];
    currMB->i16mode = (mbtype-24) & 0x03;
  }
  else if (mbtype==22) // 8x8(+split)
  {
    mbmode=P8x8;       // b8mode and pdir is transmitted in additional codewords
  }
  else if (mbtype<4)   // 16x16
  {
    mbmode=1;       for(i=0;i<4;i++) {b8mode[i]=1;          b8pdir[i]=offset2pdir16x16[mbtype]; }
  }
  else if(mbtype==48)
  {
    mbmode=IPCM;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= -1;
    currMB->i16mode = 0;
  }

  else if (mbtype%2==0) // 16x8
  {
    mbmode=2;       for(i=0;i<4;i++) {b8mode[i]=2;          b8pdir[i]=offset2pdir16x8 [mbtype][i/2]; }
  }
  else
  {
    mbmode=3;       for(i=0;i<4;i++) {b8mode[i]=3;          b8pdir[i]=offset2pdir8x16 [mbtype][i%2]; }
  }
  currMB->mb_type = mbmode;
}
/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for SI-Frames
 ************************************************************************
 */
void interpret_mb_mode_SI(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int         mbmode   = currMB->mb_type;

  if (mbmode==0)
  {
    currMB->mb_type = SI4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
    img->siblock[img->mb_x][img->mb_y]=1;
  }
  else if (mbmode==1)
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
  }
  else if(mbmode==26)
  {
    currMB->mb_type=IPCM;

    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= -1;
    currMB->i16mode = 0;

  }

  else
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(mbmode-1)>>2];
    currMB->i16mode = (mbmode-2) & 0x03;
  }
}



/*!
************************************************************************
* \brief
*    init macroblock I and P frames
************************************************************************
*/
void init_macroblock(struct img_par *img)
{
	int i,j;
	
	for (i=0;i<BLOCK_SIZE;i++)
	{                           // reset vectors and pred. modes
		for(j=0;j<BLOCK_SIZE;j++)
		{
			dec_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0]=0;
			dec_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1]=0;
			dec_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0]=0;
			dec_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1]=0;
			
			img->ipredmode[img->block_x+i][img->block_y+j] = DC_PRED;
		}
	}
	
	for (j=0; j<BLOCK_SIZE; j++)
		for (i=0; i<BLOCK_SIZE; i++)
		{
			dec_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = -1;
			dec_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
			dec_picture->ref_pic_id[LIST_0][img->block_x+i][img->block_y+j] = INT64_MIN;
			dec_picture->ref_pic_id[LIST_1][img->block_x+i][img->block_y+j] = INT64_MIN;
		}
}


/*!
************************************************************************
* \brief
*    Sets mode for 8x8 block
************************************************************************
*/
void SetB8Mode (struct img_par* img, Macroblock* currMB, int value, int i)
{
	static const int p_v2b8 [ 5] = {4, 5, 6, 7, IBLOCK};
	static const int p_v2pd [ 5] = {0, 0, 0, 0, -1};
	static const int b_v2b8 [14] = {0, 4, 4, 4, 5, 6, 5, 6, 5, 6, 7, 7, 7, IBLOCK};
	static const int b_v2pd [14] = {2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1, 2, -1};
	
	if (img->type==B_SLICE)
	{
		currMB->b8mode[i]   = b_v2b8[value];
		currMB->b8pdir[i]   = b_v2pd[value];
		
	}
	else
	{
		currMB->b8mode[i]   = p_v2b8[value];
		currMB->b8pdir[i]   = p_v2pd[value];
	}
	
}


void reset_coeffs()
{
	int i, j, iii, jjj;
	
	// reset luma coeffs
	for (i=0;i<BLOCK_SIZE;i++)
	{ 
		for (j=0;j<BLOCK_SIZE;j++)
			for(iii=0;iii<BLOCK_SIZE;iii++)
				for(jjj=0;jjj<BLOCK_SIZE;jjj++)
					img->cof[i][j][iii][jjj]=0;
	}
	
	// reset chroma coeffs
	for (j=4;j<6;j++)
	{ 
		for (i=0;i<4;i++)
			for (iii=0;iii<4;iii++)
				for (jjj=0;jjj<4;jjj++)
					img->cof[i][j][iii][jjj]=0;
	}
	
	// CAVLC
	for (i=0; i < 4; i++)
		for (j=0; j < 6; j++)
			img->nz_coeff[img->current_mb_nr][i][j]=0;  
		
}

/*!
************************************************************************
* \brief
*    Get the syntax elements from the NAL
************************************************************************
*/
int read_one_macroblock(struct img_par *img,struct inp_par_dec *inp)//img在此已经有值了，看前一个函数呗。
{
	int i;
	
	SyntaxElement currSE;
	Macroblock *currMB = &img->mb_data[img->current_mb_nr];
	
	Slice *currSlice = img->currentSlice;
	DataPartition *dP;
	int *partMap = assignSE2partition[currSlice->dp_mode];
	Macroblock *topMB = NULL;
	int  prevMbSkipped = 0;
	int  img_block_y;
	int  check_bottom, read_bottom, read_top;
	
	if (img->MbaffFrameFlag)
	{
		if (img->current_mb_nr%2)
		{
			topMB= &img->mb_data[img->current_mb_nr-1];
			if(!(img->type == B_SLICE))
				prevMbSkipped = (topMB->mb_type == 0);
			else 
				prevMbSkipped = (topMB->mb_type == 0 && topMB->cbp == 0);
		}
		else 
			prevMbSkipped = 0;
	}
	
	if (img->current_mb_nr%2 == 0)
		currMB->mb_field = 0;
	else
		currMB->mb_field = img->mb_data[img->current_mb_nr-1].mb_field;
	
	
	currMB->qp = img->qp ;
	
	currSE.type = SE_MBTYPE;
	
	//  read MB mode *****************************************************************
	dP = &(currSlice->partArr[partMap[currSE.type]]);
	
	if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping_dec = linfo_ue_dec;
	
	if(img->type == I_SLICE || img->type == SI_SLICE)
	{
		// read MB aff
		if (img->MbaffFrameFlag && img->current_mb_nr%2==0)
		{
			TRACE_STRING("mb_field_decoding_flag");
			if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
			{
				currSE.len = 1;
				readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
			}
			else
			{
				currSE.reading = readFieldModeInfo_CABAC;
				dP->readSyntaxElement(&currSE,img,inp,dP);
			}
			currMB->mb_field = currSE.value1;
		}
		if(active_pps->entropy_coding_mode_flag  == CABAC)
			CheckAvailabilityOfNeighborsCABAC();
		
		//  read MB type
		TRACE_STRING("mb_type");
		currSE.reading = readMB_typeInfo_CABAC;
		dP->readSyntaxElement(&currSE,img,inp,dP);
		
		currMB->mb_type = currSE.value1;
		if(!dP->bitstream->ei_flag)
			currMB->ei_flag = 0;
	}
	
	// non I/SI-slice CABAC
	else if (active_pps->entropy_coding_mode_flag == CABAC)
	{
		// read MB skipflag
		if (img->MbaffFrameFlag && (img->current_mb_nr%2 == 0||prevMbSkipped))
			field_flag_inference();
		
		CheckAvailabilityOfNeighborsCABAC();
		TRACE_STRING("mb_skip_flag");
		currSE.reading = readMB_skip_flagInfo_CABAC;
		dP->readSyntaxElement(&currSE,img,inp,dP);
		currMB->mb_type = currSE.value1;
		
		if (img->type==B_SLICE)
			currMB->cbp = currSE.value2;
		if(!dP->bitstream->ei_flag)
			currMB->ei_flag = 0;
		
		if ((img->type==B_SLICE) && currSE.value1==0 && currSE.value2==0)
			img->cod_counter=0;
		
		// read MB aff
		if (img->MbaffFrameFlag) 
		{
			check_bottom=read_bottom=read_top=0;
			if (img->current_mb_nr%2==0)
			{
				check_bottom =  (img->type!=B_SLICE)? 
					(currMB->mb_type == 0):
				(currMB->mb_type == 0 && currMB->cbp == 0);
				read_top = !check_bottom;
			}
			else
				read_bottom = (img->type!=B_SLICE)? 
				(topMB->mb_type == 0 && currMB->mb_type != 0) :
			((topMB->mb_type == 0 && topMB->cbp == 0) && (currMB->mb_type != 0 || currMB->cbp != 0));
			
			if (read_bottom || read_top)
			{
				TRACE_STRING("mb_field_decoding_flag");
				currSE.reading = readFieldModeInfo_CABAC;
				dP->readSyntaxElement(&currSE,img,inp,dP);
				currMB->mb_field = currSE.value1;
			}
			if (check_bottom)
				check_next_mb_and_get_field_mode_CABAC(&currSE,img,inp,dP);
			
		}
		if(active_pps->entropy_coding_mode_flag  == CABAC)
			CheckAvailabilityOfNeighborsCABAC();
		
		// read MB type
		if (currMB->mb_type != 0 )
		{
			currSE.reading = readMB_typeInfo_CABAC;
			TRACE_STRING("mb_type");
			dP->readSyntaxElement(&currSE,img,inp,dP);
			currMB->mb_type = currSE.value1;
			if(!dP->bitstream->ei_flag)
				currMB->ei_flag = 0;
		}
	}
	// VLC Non-Intra
	else
	{
		if(img->cod_counter == -1)
		{
			TRACE_STRING("mb_skip_run");
			dP->readSyntaxElement(&currSE,img,inp,dP);
			img->cod_counter = currSE.value1;
		}
		if (img->cod_counter==0)
		{
			// read MB aff
			if ((img->MbaffFrameFlag) && ((img->current_mb_nr%2==0) || ((img->current_mb_nr%2) && prevMbSkipped)))
			{
				TRACE_STRING("mb_field_decoding_flag");
				currSE.len = 1;
				readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
				currMB->mb_field = currSE.value1;
			}
			
			// read MB type
			TRACE_STRING("mb_type");
			dP->readSyntaxElement(&currSE,img,inp,dP);
			if(img->type == P_SLICE || img->type == SP_SLICE)
				currSE.value1++;
			currMB->mb_type = currSE.value1;
			if(!dP->bitstream->ei_flag)
				currMB->ei_flag = 0;
			img->cod_counter--;
		} 
		else
		{
			img->cod_counter--;
			currMB->mb_type = 0;
			currMB->ei_flag = 0;
			
			// read field flag of bottom block
			if(img->MbaffFrameFlag)
			{
				if(img->cod_counter == 0 && (img->current_mb_nr%2 == 0))
				{
					TRACE_STRING("mb_field_decoding_flag (of coded bottom mb)");
					currSE.len = 1;
					readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
					dP->bitstream->frame_bitoffset--;
					currMB->mb_field = currSE.value1;
				}
				else if(img->cod_counter > 0 && (img->current_mb_nr%2 == 0))
				{
					// check left macroblock pair first
					if (mb_is_available(img->current_mb_nr-2, img->current_mb_nr)&&((img->current_mb_nr%(img->PicWidthInMbs*2))!=0))
					{
						currMB->mb_field = img->mb_data[img->current_mb_nr-2].mb_field;
					}
					else
					{
						// check top macroblock pair
						if (mb_is_available(img->current_mb_nr-2*img->PicWidthInMbs, img->current_mb_nr))
						{
							currMB->mb_field = img->mb_data[img->current_mb_nr-2*img->PicWidthInMbs].mb_field;
						}
						else
							currMB->mb_field = 0;
					}
				}
			}
		}
	}
	
	dec_picture->mb_field[img->current_mb_nr] = currMB->mb_field;
	
	img->siblock[img->mb_x][img->mb_y]=0;
	
	if ((img->type==P_SLICE ))    // inter frame
		interpret_mb_mode_P(img);
	else if (img->type==I_SLICE)                                  // intra frame
		interpret_mb_mode_I(img);
	else if ((img->type==B_SLICE))       // B frame
		interpret_mb_mode_B(img);
	else if ((img->type==SP_SLICE))     // SP frame
		interpret_mb_mode_P(img);
	else if (img->type==SI_SLICE)     // SI frame
		interpret_mb_mode_SI(img);
	
	if(img->MbaffFrameFlag)
	{
		if(currMB->mb_field)
		{
			img->num_ref_idx_l0_active <<=1;
			img->num_ref_idx_l1_active <<=1;
		}
	}
	
	//====== READ 8x8 SUB-PARTITION MODES (modes of 8x8 blocks) and Intra VBST block modes ======
	if (IS_P8x8 (currMB))
	{
		currSE.type    = SE_MBTYPE;
		dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
		
		for (i=0; i<4; i++)
		{
			if (active_pps->entropy_coding_mode_flag ==UVLC || dP->bitstream->ei_flag) currSE.mapping_dec = linfo_ue_dec;
			else                                                  currSE.reading = readB8_typeInfo_CABAC;
			
			TRACE_STRING("sub_mb_type");
			dP->readSyntaxElement (&currSE, img, inp, dP);
			SetB8Mode (img, currMB, currSE.value1, i);
		}
	}
	
	if(active_pps->constrained_intra_pred_flag && (img->type==P_SLICE|| img->type==B_SLICE))        // inter frame
	{
		if( !IS_INTRA(currMB) )
		{
			img->intra_block[img->current_mb_nr] = 0;
		}
	}
	
	//! TO for Error Concelament
	//! If we have an INTRA Macroblock and we lost the partition
	//! which contains the intra coefficients Copy MB would be better 
	//! than just a grey block.
	//! Seems to be a bit at the wrong place to do this right here, but for this case 
	//! up to now there is no other way.
	dP = &(currSlice->partArr[partMap[SE_CBP_INTRA]]);
	if(IS_INTRA (currMB) && dP->bitstream->ei_flag && img->number)
	{
		currMB->mb_type = 0;
		currMB->ei_flag = 1;
		for (i=0;i<4;i++) {currMB->b8mode[i]=currMB->b8pdir[i]=0; }
	}
	dP = &(currSlice->partArr[partMap[currSE.type]]);
	//! End TO
	
	
	//--- init macroblock data ---
	init_macroblock       (img);
	
	if (IS_DIRECT (currMB) && img->cod_counter >= 0)
	{
		currMB->cbp = 0;
		reset_coeffs();
		
		if (active_pps->entropy_coding_mode_flag ==CABAC)
			img->cod_counter=-1;
		
		return DECODE_MB;
	}
	
// 	if (IS_COPY (currMB))   //keep last macroblock   //ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ有问题
// 	{
// 		int i, j, k, pmv[2];
// 		int zeroMotionAbove;
// 		int zeroMotionLeft;
// 		PixelPos mb_a, mb_b;
// 		int      a_mv_y = 0;
// 		int      a_ref_idx = 0;
// 		int      b_mv_y = 0;
// 		int      b_ref_idx = 0;
// 		int      list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
// 		
// 		getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_a);
// 		getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_b);
// 		
// 		if (mb_a.available)
// 		{
// 			a_mv_y    = dec_picture->mv[LIST_0][mb_a.pos_x][mb_a.pos_y][1];
// 			a_ref_idx = dec_picture->ref_idx[LIST_0][mb_a.pos_x][mb_a.pos_y];
// 			
// 			if (currMB->mb_field && !img->mb_data[mb_a.mb_addr].mb_field)
// 			{
// 				a_mv_y    /=2;
// 				a_ref_idx *=2;
// 			}
// 			if (!currMB->mb_field && img->mb_data[mb_a.mb_addr].mb_field)
// 			{
// 				a_mv_y    *=2;
// 				a_ref_idx >>=1;
// 			}
// 		}
// 		
// 		if (mb_b.available)
// 		{
// 			b_mv_y    = dec_picture->mv[LIST_0][mb_b.pos_x][mb_b.pos_y][1];
// 			b_ref_idx = dec_picture->ref_idx[LIST_0][mb_b.pos_x][mb_b.pos_y];
// 			
// 			if (currMB->mb_field && !img->mb_data[mb_b.mb_addr].mb_field)
// 			{
// 				b_mv_y    /=2;
// 				b_ref_idx *=2;
// 			}
// 			if (!currMB->mb_field && img->mb_data[mb_b.mb_addr].mb_field)
// 			{
// 				b_mv_y    *=2;
// 				b_ref_idx >>=1;
// 			}
// 		}
// 		
// 		zeroMotionLeft  = !mb_a.available ? 1 : a_ref_idx==0 && dec_picture->mv[LIST_0][mb_a.pos_x][mb_a.pos_y][0]==0 && a_mv_y==0 ? 1 : 0;
// 		zeroMotionAbove = !mb_b.available ? 1 : b_ref_idx==0 && dec_picture->mv[LIST_0][mb_b.pos_x][mb_b.pos_y][0]==0 && b_mv_y==0 ? 1 : 0;
// 		
// 		currMB->cbp = 0;
// 		reset_coeffs();
// 		
// 		img_block_y   = img->block_y;
// 		
// 		if (zeroMotionAbove || zeroMotionLeft)
// 		{
// 			for(i=0;i<BLOCK_SIZE;i++)
// 				for(j=0;j<BLOCK_SIZE;j++)
// 					for (k=0;k<2;k++)
// 						dec_picture->mv[LIST_0][img->block_x+i][img->block_y+j][k] = 0;
// 		}
// 		else
// 		{
// 			SetMotionVectorPredictor (img, pmv, pmv+1, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
// 			
// 			for(i=0;i<BLOCK_SIZE;i++)
// 				for(j=0;j<BLOCK_SIZE;j++)
// 					for (k=0;k<2;k++)
// 					{
// 						dec_picture->mv[LIST_0][img->block_x+i][img_block_y+j][k] = pmv[k];
// 					}
// 		}
// 		
// 		for(i=0;i<BLOCK_SIZE;i++)
// 			for(j=0;j<BLOCK_SIZE;j++)
// 			{
// 				dec_picture->ref_idx[LIST_0][img->block_x+i][img_block_y+j] = 0;
// 				dec_picture->ref_pic_id[LIST_0][img->block_x+i][img_block_y+j] = dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_0 + list_offset][dec_picture->ref_idx[LIST_0][img->block_x+i][img_block_y+j]];
// 			}
// 			
// 		return DECODE_MB;
// 	}

	if(currMB->mb_type!=IPCM)
	{
		
		// intra prediction modes for a macroblock 4x4 **********************************************
		read_ipred_modes(img,inp);//查找img的赋值。
		
		// read inter frame vector data *********************************************************
		if (IS_INTERMV (currMB))
		{
			readMotionInfoFromNAL (img, inp);
		}
		// read CBP and Coeffs  ***************************************************************
		readCBPandCoeffsFromNAL (img,inp);
	}
	else
	{
		//read pcm_alignment_zero_bit and pcm_byte[i] 
		
		// here dP is assigned with the same dP as SE_MBTYPE, because IPCM syntax is in the 
		// same category as MBTYPE
		dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
		readIPCMcoeffsFromNAL(img,inp,dP);
	}
	return DECODE_MB;//1
}



/*!
************************************************************************
* \brief
*    Read IPCM pcm_alignment_zero_bit and pcm_byte[i] from stream to img->cof
*    (for IPCM CABAC and IPCM CAVLC  28/11/2003)
*
* \author
*    Dong Wang <Dong.Wang@bristol.ac.uk>
************************************************************************
*/

void readIPCMcoeffsFromNAL(struct img_par *img, struct inp_par_dec *inp, struct datapartition *dP)
{
	SyntaxElement currSE;
	int i,j;
	
	//For CABAC, we don't need to read bits to let stream byte aligned
	//  because we have variable for integer bytes position
	if(active_pps->entropy_coding_mode_flag  == CABAC)
	{
		//read luma and chroma IPCM coefficients
		currSE.len=8;
		
		for(i=0;i<16;i++)
			for(j=0;j<16;j++)
			{
				readIPCMBytes_CABAC(&currSE, dP->bitstream);
				img->cof[i/4][j/4][i%4][j%4]=currSE.value1;
			}
			
			for(i=0;i<8;i++)
				for(j=0;j<8;j++)
				{
					readIPCMBytes_CABAC(&currSE, dP->bitstream);
					img->cof[i/4][j/4+4][i%4][j%4]=currSE.value1;
				}
				
				for(i=0;i<8;i++)
					for(j=0;j<8;j++)
					{
						readIPCMBytes_CABAC(&currSE, dP->bitstream);
						img->cof[i/4+2][j/4+4][i%4][j%4]=currSE.value1;
					}
					
					//If the decoded MB is IPCM MB, decoding engine is initialized
					
					// here the decoding engine is directly initialized without checking End of Slice
					// The reason is that, whether current MB is the last MB in slice or not, there is
					// at least one 'end of slice' syntax after this MB. So when fetching bytes in this 
					// initialisation process, we can guarantee there is bits available in bitstream. 
// 					init_decoding_engine_IPCM(img);ZZZZZZZZZZZZZZZZZZZZZZZZZZ检查是否经过再决定是否添加函数
					
	}
	else
	{ 
		
		//read bits to let stream byte aligned
		
		
		if((dP->bitstream->frame_bitoffset)%8!=0)
		{
			
			currSE.len=8-(dP->bitstream->frame_bitoffset)%8;
			readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
			
		}
		
		//read luma and chroma IPCM coefficients
		currSE.len=8;
		
		for(i=0;i<16;i++)
			for(j=0;j<16;j++)
			{
				readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
				img->cof[i/4][j/4][i%4][j%4]=currSE.value1;
			}
			
			for(i=0;i<8;i++)
				for(j=0;j<8;j++)
				{
					readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
					img->cof[i/4][j/4+4][i%4][j%4]=currSE.value1;
				}
				
				for(i=0;i<8;i++)
					for(j=0;j<8;j++)
					{
						readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
						img->cof[i/4+2][j/4+4][i%4][j%4]=currSE.value1;
					}
	}
}



void read_ipred_modes(struct img_par *img,struct inp_par_dec *inp)
{
	int b8,i,j,bi,bj,bx,by,dec;
	SyntaxElement currSE;
	Slice *currSlice;
	DataPartition *dP;
	int *partMap;
	Macroblock *currMB;
	int ts, ls;
	int mostProbableIntraPredMode;
	int upIntraPredMode;
	int leftIntraPredMode;
	int IntraChromaPredModeFlag;
	
	PixelPos left_block;
	PixelPos top_block;
	
	currMB = &img->mb_data[img->current_mb_nr];
	
	IntraChromaPredModeFlag = IS_INTRA(currMB);//1
	
	currSlice = img->currentSlice;
	partMap = assignSE2partition[currSlice->dp_mode];
	
	currSE.type = SE_INTRAPREDMODE;
	
	TRACE_STRING("intra4x4_pred_mode");
	dP = &(currSlice->partArr[partMap[currSE.type]]);//注意dp的赋值。
	
	if (!(active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)) 
		currSE.reading = readIntraPredMode_CABAC;
	
	for(b8=0;b8<4;b8++)  //loop 8x8 blocks
	{
		if( currMB->b8mode[b8]==IBLOCK )
		{
			IntraChromaPredModeFlag = 1;
			
			for(j=0;j<2;j++)  //loop subblocks
				for(i=0;i<2;i++)
				{
					//get from stream
					if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
						readSyntaxElement_Intra4x4PredictionMode_dec(&currSE,img,inp,dP);
					else 
					{
						currSE.context=(b8<<2)+(j<<1)+i;
						dP->readSyntaxElement(&currSE,img,inp,dP);//注意dp，dep_dp=&(dp->de_cabac)
					}
					
					bx = ((b8&1)<<1) + i;
					by = (b8&2)      + j;
					
					getLuma4x4Neighbour(img->current_mb_nr, bx, by, -1,  0, &left_block);
					getLuma4x4Neighbour(img->current_mb_nr, bx, by,  0, -1, &top_block);
					
					//get from array and decode
					bi = img->block_x + bx;
					bj = img->block_y + by;
					
					if (active_pps->constrained_intra_pred_flag)
					{
						left_block.available = left_block.available ? img->intra_block[left_block.mb_addr] : 0;
						top_block.available  = top_block.available  ? img->intra_block[top_block.mb_addr]  : 0;
					}
					
					// !! KS: not sure if the follwing is still correct...
					ts=ls=0;   // Check to see if the neighboring block is SI
					if (IS_OLDINTRA(currMB) && img->type == SI_SLICE)           // need support for MBINTLC1
					{
						if (left_block.available)
							if (img->siblock [left_block.pos_x][left_block.pos_y])
								ls=1;
							
							if (top_block.available)
								if (img->siblock [top_block.pos_x][top_block.pos_y])
									ts=1;
					}
					
					upIntraPredMode            = (top_block.available  &&(ts == 0)) ? img->ipredmode[top_block.pos_x ][top_block.pos_y ] : -1;
					leftIntraPredMode          = (left_block.available &&(ls == 0)) ? img->ipredmode[left_block.pos_x][left_block.pos_y] : -1;
					
					mostProbableIntraPredMode  = (upIntraPredMode < 0 || leftIntraPredMode < 0) ? DC_PRED : upIntraPredMode < leftIntraPredMode ? upIntraPredMode : leftIntraPredMode;
					
					dec = (currSE.value1 == -1) ? mostProbableIntraPredMode : currSE.value1 + (currSE.value1 >= mostProbableIntraPredMode);//1
					//  =        1                          2                       1             1                          2
// 					printf("%d=(%d==-1)？%d:%d+(%d>=%d)",dec,currSE.value1,mostProbableIntraPredMode,currSE.value1,currSE.value1, mostProbableIntraPredMode);
					//set
					img->ipredmode[bi][bj]=dec; //1
// 					printf("img->ipredmode[%d][%d]=%d\n",bi,bj,img->ipredmode[bi][bj]);
				}
		}
	}
	
	if (IntraChromaPredModeFlag)
	{
		currSE.type = SE_INTRAPREDMODE;
		TRACE_STRING("intra_chroma_pred_mode");
		dP = &(currSlice->partArr[partMap[currSE.type]]);
		
		if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping_dec = linfo_ue_dec;
		else                                                    currSE.reading = readCIPredMode_CABAC;
		
		dP->readSyntaxElement(&currSE,img,inp,dP);
		currMB->c_ipred_mode = currSE.value1;
		
		if (currMB->c_ipred_mode < DC_PRED_8 || currMB->c_ipred_mode > PLANE_8)
		{
			error("illegal chroma intra pred mode!\n", 600);
		}
	}
}



/*!
************************************************************************
* \brief
*    Set motion vector predictor
************************************************************************
*/
static void SetMotionVectorPredictor (struct img_par  *img,
                                      int             *pmv_x,
                                      int             *pmv_y,
                                      int             ref_frame,
                                      int             list,
                                      int             ***refPic,
                                      int             ****tmp_mv,
                                      int             block_x,
                                      int             block_y,
                                      int             blockshape_x,
                                      int             blockshape_y)
{
	int mb_x                 = BLOCK_SIZE*block_x;
	int mb_y                 = BLOCK_SIZE*block_y;
	int mb_nr                = img->current_mb_nr;
	
	int mv_a, mv_b, mv_c, pred_vec=0;
	int mvPredType, rFrameL, rFrameU, rFrameUR;
	int hv;
	
	
	PixelPos block_a, block_b, block_c, block_d;
	
	getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1,  0, &block_a);
	getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
	getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
	getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);
	
	if (mb_y > 0)
	{
		if (mb_x < 8)  // first column of 8x8 blocks
		{
			if (mb_y==8)
			{
				if (blockshape_x == 16)      block_c.available  = 0;
				else                         block_c.available &= 1;
			}
			else
			{
				if (mb_x+blockshape_x != 8)  block_c.available &= 1;
				else                         block_c.available  = 0;
			}
		}
		else
		{
			if (mb_x+blockshape_x != 16)   block_c.available &= 1;
			else                           block_c.available  = 0;
		}
	}
	
	if (!block_c.available)
	{
		block_c=block_d;
	}
	
	mvPredType = MVPRED_MEDIAN;
	
	if (!img->MbaffFrameFlag)
	{
		rFrameL    = block_a.available    ? refPic[list][block_a.pos_x][block_a.pos_y] : -1;
		rFrameU    = block_b.available    ? refPic[list][block_b.pos_x][block_b.pos_y] : -1;
		rFrameUR   = block_c.available    ? refPic[list][block_c.pos_x][block_c.pos_y] : -1;
	}
	else
	{
		if (img->mb_data[img->current_mb_nr].mb_field)
		{
			rFrameL    = block_a.available    ? 
				img->mb_data[block_a.mb_addr].mb_field ? 
				refPic[list][block_a.pos_x][block_a.pos_y]:
			refPic[list][block_a.pos_x][block_a.pos_y] * 2: 
			-1;
			rFrameU    = block_b.available    ? 
				img->mb_data[block_b.mb_addr].mb_field ? 
				refPic[list][block_b.pos_x][block_b.pos_y]:
			refPic[list][block_b.pos_x][block_b.pos_y] * 2: 
			-1;
			rFrameUR    = block_c.available    ? 
				img->mb_data[block_c.mb_addr].mb_field ? 
				refPic[list][block_c.pos_x][block_c.pos_y]:
			refPic[list][block_c.pos_x][block_c.pos_y] * 2: 
			-1;
		}
		else
		{
			rFrameL    = block_a.available    ? 
				img->mb_data[block_a.mb_addr].mb_field ? 
				refPic[list][block_a.pos_x][block_a.pos_y] >>1:
			refPic[list][block_a.pos_x][block_a.pos_y] : 
			-1;
			rFrameU    = block_b.available    ? 
				img->mb_data[block_b.mb_addr].mb_field ? 
				refPic[list][block_b.pos_x][block_b.pos_y] >>1:
			refPic[list][block_b.pos_x][block_b.pos_y] : 
			-1;
			rFrameUR    = block_c.available    ? 
				img->mb_data[block_c.mb_addr].mb_field ? 
				refPic[list][block_c.pos_x][block_c.pos_y] >>1:
			refPic[list][block_c.pos_x][block_c.pos_y] : 
			-1;
		}
	}
	
	
	/* Prediction if only one of the neighbors uses the reference frame
	* we are checking
	*/
	if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       mvPredType = MVPRED_L;
	else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  mvPredType = MVPRED_U;
	else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  mvPredType = MVPRED_UR;
	// Directional predictions 
	if(blockshape_x == 8 && blockshape_y == 16)
	{
		if(mb_x == 0)
		{
			if(rFrameL == ref_frame)
				mvPredType = MVPRED_L;
		}
		else
		{
			if( rFrameUR == ref_frame)
				mvPredType = MVPRED_UR;
		}
	}
	else if(blockshape_x == 16 && blockshape_y == 8)
	{
		if(mb_y == 0)
		{
			if(rFrameU == ref_frame)
				mvPredType = MVPRED_U;
		}
		else
		{
			if(rFrameL == ref_frame)
				mvPredType = MVPRED_L;
		}
	}
	
	for (hv=0; hv < 2; hv++)
	{
		if (!img->MbaffFrameFlag || hv==0)
		{
			mv_a = block_a.available  ? tmp_mv[list][block_a.pos_x][block_a.pos_y][hv] : 0;
			mv_b = block_b.available  ? tmp_mv[list][block_b.pos_x][block_b.pos_y][hv] : 0;
			mv_c = block_c.available  ? tmp_mv[list][block_c.pos_x][block_c.pos_y][hv] : 0;
		}
		else
		{
			if (img->mb_data[img->current_mb_nr].mb_field)
			{
				mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field?
					tmp_mv[list][block_a.pos_x][block_a.pos_y][hv]:
				tmp_mv[list][block_a.pos_x][block_a.pos_y][hv] / 2: 
				0;
				mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field?
					tmp_mv[list][block_b.pos_x][block_b.pos_y][hv]:
				tmp_mv[list][block_b.pos_x][block_b.pos_y][hv] / 2: 
				0;
				mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field?
					tmp_mv[list][block_c.pos_x][block_c.pos_y][hv]:
				tmp_mv[list][block_c.pos_x][block_c.pos_y][hv] / 2: 
				0;
			}
			else
			{
				mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field?
					tmp_mv[list][block_a.pos_x][block_a.pos_y][hv] * 2:
				tmp_mv[list][block_a.pos_x][block_a.pos_y][hv]: 
				0;
				mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field?
					tmp_mv[list][block_b.pos_x][block_b.pos_y][hv] * 2:
				tmp_mv[list][block_b.pos_x][block_b.pos_y][hv]: 
				0;
				mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field?
					tmp_mv[list][block_c.pos_x][block_c.pos_y][hv] * 2:
				tmp_mv[list][block_c.pos_x][block_c.pos_y][hv]: 
				0;
			}
		}
		
		
		switch (mvPredType)
		{
		case MVPRED_MEDIAN:
			if(!(block_b.available || block_c.available))
				pred_vec = mv_a;
			else
				pred_vec = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
			break;
		case MVPRED_L:
			pred_vec = mv_a;
			break;
		case MVPRED_U:
			pred_vec = mv_b;
			break;
		case MVPRED_UR:
			pred_vec = mv_c;
			break;
		default:
			break;
		}
		
		if (hv==0)  *pmv_x = pred_vec;
		else        *pmv_y = pred_vec;
		
	}
}


/*!
************************************************************************
* \brief
*    Set context for reference frames
************************************************************************
*/
int
BType2CtxRef (int btype)
{
	if (btype<4)  return 0;
	else          return 1;
}

/*!
************************************************************************
* \brief
*    Read motion info
************************************************************************
*/
void readMotionInfoFromNAL (struct img_par *img, struct inp_par_dec *inp)
{
	int i,j,k;
	int step_h,step_v;
	int curr_mvd;
	Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
	SyntaxElement currSE;
	Slice *currSlice    = img->currentSlice;
	DataPartition *dP;
	int *partMap        = assignSE2partition[currSlice->dp_mode];
	int bframe          = (img->type==B_SLICE);
	int partmode        = (IS_P8x8(currMB)?4:currMB->mb_type);
	int step_h0         = BLOCK_STEP [partmode][0];
	int step_v0         = BLOCK_STEP [partmode][1];
	
	int mv_mode, i0, j0, refframe;
	int pmv[2];
	int j4, i4, ii,jj;
	int vec;
	
	int mv_scale = 0;
	
	int flag_mode;
	
	int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
	
	byte **    moving_block;
	int ****   co_located_mv;
	int ***    co_located_ref_idx;
	int64 ***    co_located_ref_id;
	
	if ((img->MbaffFrameFlag)&&(currMB->mb_field))
	{
		if(img->current_mb_nr%2)
		{
			moving_block = Co_located->bottom_moving_block;
			co_located_mv = Co_located->bottom_mv;
			co_located_ref_idx = Co_located->bottom_ref_idx;
			co_located_ref_id = Co_located->bottom_ref_pic_id;
		}
		else
		{
			moving_block = Co_located->top_moving_block;
			co_located_mv = Co_located->top_mv;
			co_located_ref_idx = Co_located->top_ref_idx;
			co_located_ref_id = Co_located->top_ref_pic_id;
		}
	}
	else
	{
		moving_block = Co_located->moving_block;
		co_located_mv = Co_located->mv;
		co_located_ref_idx = Co_located->ref_idx;
		co_located_ref_id = Co_located->ref_pic_id;
	}
	
	if (bframe && IS_P8x8 (currMB))
	{
		if (img->direct_type)
		{
			int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
			int fw_rFrameL, fw_rFrameU, fw_rFrameUL, fw_rFrameUR;
			int bw_rFrameL, bw_rFrameU, bw_rFrameUL, bw_rFrameUR;
			
			PixelPos mb_left, mb_up, mb_upleft, mb_upright;
			
			int fw_rFrame,bw_rFrame;
			int pmvfw[2]={0,0},pmvbw[2]={0,0};
			
			
			getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1,  0, &mb_left);
			getLuma4x4Neighbour(img->current_mb_nr, 0, 0,  0, -1, &mb_up);
			getLuma4x4Neighbour(img->current_mb_nr, 0, 0, 16, -1, &mb_upright);
			getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1, -1, &mb_upleft);
			
			if (!img->MbaffFrameFlag)
			{
				fw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : -1;
				fw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
				fw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
				fw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
				
				bw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
				bw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
				bw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
				bw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
			}
			else
			{
				if (img->mb_data[img->current_mb_nr].mb_field)
				{
					fw_rFrameL = mb_left.available ? 
						img->mb_data[mb_left.mb_addr].mb_field  || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] < 0? 
						dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : 
					dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] * 2: -1;
					fw_rFrameU = mb_up.available ? 
						img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] < 0? 
						dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : 
					dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] * 2: -1;
					
					fw_rFrameUL = mb_upleft.available ? 
						img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0?
						dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : 
					dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      
					
					fw_rFrameUR = mb_upright.available ? 
						img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] < 0 ?
						dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : 
					dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] * 2: fw_rFrameUL;      
					
					bw_rFrameL = mb_left.available ? 
						img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : 
					dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] * 2: -1;
					
					bw_rFrameU = mb_up.available ? 
						img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : 
					dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] * 2: -1;
					
					bw_rFrameUL = mb_upleft.available ? 
						img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : 
					dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      
					
					bw_rFrameUR = mb_upright.available ? 
						img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : 
					dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] * 2: bw_rFrameUL;      
					
				}
				else
				{
					fw_rFrameL = mb_left.available ? 
						img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] >> 1 : 
					dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]: -1;
					
					fw_rFrameU = mb_up.available ? 
						img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] >> 1 :  
					dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
					
					fw_rFrameUL = mb_upleft.available ? 
						img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0 ?
						dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y]>> 1 : 
					dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      
					
					fw_rFrameUR = mb_upright.available ? 
						img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] >> 1 :  
					dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
					
					bw_rFrameL = mb_left.available ? 
						img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] < 0 ?
						dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] >> 1 :  
					dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
					bw_rFrameU = mb_up.available ? 
						img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] < 0 ?
						dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] >> 1 : 
					dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
					
					bw_rFrameUL = mb_upleft.available ? 
						img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y]  < 0 ?
						dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] >> 1 : 
					dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      
					
					bw_rFrameUR = mb_upright.available ? 
						img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] < 0 ?
						dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] >> 1: 
					dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
				}
			}
			
			fw_rFrame = (fw_rFrameL >= 0 && fw_rFrameU >= 0) ? min(fw_rFrameL,fw_rFrameU): max(fw_rFrameL,fw_rFrameU);
			fw_rFrame = (fw_rFrame >= 0 && fw_rFrameUR >= 0) ? min(fw_rFrame,fw_rFrameUR): max(fw_rFrame,fw_rFrameUR);
			
			bw_rFrame = (bw_rFrameL >= 0 && bw_rFrameU >= 0) ? min(bw_rFrameL,bw_rFrameU): max(bw_rFrameL,bw_rFrameU);
			bw_rFrame = (bw_rFrame >= 0 && bw_rFrameUR >= 0) ? min(bw_rFrame,bw_rFrameUR): max(bw_rFrame,bw_rFrameUR);
			
			
			if (fw_rFrame >=0)
				SetMotionVectorPredictor (img, pmvfw, pmvfw+1, fw_rFrame, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
			
			if (bw_rFrame >=0)
				SetMotionVectorPredictor (img, pmvbw, pmvbw+1, bw_rFrame, LIST_1, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
			
			
			for (i=0;i<4;i++)
			{
				if (currMB->b8mode[i] == 0)
					for(j=2*(i/2);j<2*(i/2)+2;j++)
						for(k=2*(i%2);k<2*(i%2)+2;k++)
						{
							int j6 = imgblock_y+j;
							j4 = img->block_y+j;
							i4 = img->block_x+k;
							
							
							if (fw_rFrame >= 0)
							{
								
								if  (!fw_rFrame  && ((!moving_block[i4][j6]) && (!listX[1+list_offset][0]->is_long_term)))
								{                    
									dec_picture->mv  [LIST_0][i4][j4][0] = 0;
									dec_picture->mv  [LIST_0][i4][j4][1] = 0;
									dec_picture->ref_idx[LIST_0][i4][j4] = 0;                    
								}
								else
								{
									
									dec_picture->mv  [LIST_0][i4][j4][0] = pmvfw[0];
									dec_picture->mv  [LIST_0][i4][j4][1] = pmvfw[1];
									dec_picture->ref_idx[LIST_0][i4][j4] = fw_rFrame;
								}
							}
							else
							{
								dec_picture->mv  [LIST_0][i4][j4][0] = 0;
								dec_picture->mv  [LIST_0][i4][j4][1] = 0;
								dec_picture->ref_idx[LIST_0][i4][j4] = -1;
							}
							if (bw_rFrame >= 0)
							{
								if  (bw_rFrame==0 && ((!moving_block[i4][j6])&& (!listX[1+list_offset][0]->is_long_term)))
								{
									dec_picture->mv  [LIST_1][i4][j4][0] = 0;
									dec_picture->mv  [LIST_1][i4][j4][1] = 0;
									dec_picture->ref_idx[LIST_1][i4][j4] = 0;
								}
								else
								{
									dec_picture->mv  [LIST_1][i4][j4][0] = pmvbw[0];
									dec_picture->mv  [LIST_1][i4][j4][1] = pmvbw[1];
									dec_picture->ref_idx[LIST_1][i4][j4] = bw_rFrame;
								}
							}
							else
							{
								dec_picture->mv  [LIST_1][i4][j4][0] = 0;
								dec_picture->mv  [LIST_1][i4][j4][1] = 0;
								dec_picture->ref_idx[LIST_1][i4][j4] = -1;                               
							}
							
							if (fw_rFrame <0 && bw_rFrame <0)
							{
								dec_picture->ref_idx[LIST_0][i4][j4] = 0;
								dec_picture->ref_idx[LIST_1][i4][j4] = 0;                  
							}
						}
			}
    }
    else
    {
		for (i=0;i<4;i++)
		{
			if (currMB->b8mode[i] == 0)
			{
				
				for(j=2*(i/2);j<2*(i/2)+2;j++)
				{
					for(k=2*(i%2);k<2*(i%2)+2;k++)
					{
						
						int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
						int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2 : img->block_y/2 : img->block_y;
						int refList = co_located_ref_idx[LIST_0 ][img->block_x+k][imgblock_y+j]== -1 ? LIST_1 : LIST_0;
						int ref_idx = co_located_ref_idx[refList][img->block_x + k][imgblock_y + j];
						int mapped_idx=-1, iref;                             
						
						
						
						if (ref_idx == -1)
						{
							dec_picture->ref_idx [LIST_0][img->block_x + k][img->block_y + j] = 0;
							dec_picture->ref_idx [LIST_1][img->block_x + k][img->block_y + j] = 0;                
						}
						else
						{
							for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
							{
								if (dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][img->block_x + k][imgblock_y + j])
								{
									mapped_idx=iref;
									break;
								}
								else //! invalid index. Default to zero even though this case should not happen
									mapped_idx=INVALIDINDEX;
							}
							if (INVALIDINDEX == mapped_idx)
							{
								error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
							}
							dec_picture->ref_idx [LIST_0][img->block_x + k][img->block_y + j] = mapped_idx;
							dec_picture->ref_idx [LIST_1][img->block_x + k][img->block_y + j] = 0;                
						}
					}
				}
			}
		}
    }
  } 
  
  //  If multiple ref. frames, read reference frame for the MB *********************************
  if(img->num_ref_idx_l0_active>1) 
  {
	  flag_mode = ( img->num_ref_idx_l0_active == 2 ? 1 : 0);
	  
	  currSE.type = SE_REFFRAME;
	  dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);
	  
	  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping_dec = linfo_ue_dec;
	  else                                                      currSE.reading = readRefFrame_CABAC;
	  
	  for (j0=0; j0<4; j0+=step_v0)
	  {
		  for (i0=0; i0<4; i0+=step_h0)
		  {
			  k=2*(j0/2)+(i0/2);
			  if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
			  {
				  TRACE_STRING("ref_idx_l0");
				  
				  img->subblock_x = i0;
				  img->subblock_y = j0;
				  
				  if (!IS_P8x8 (currMB) || bframe || (!bframe && !img->allrefzero))
				  {
					  currSE.context = BType2CtxRef (currMB->b8mode[k]);
					  if( (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) && flag_mode )
					  {
						  currSE.len = 1;
						  readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
						  currSE.value1 = 1 - currSE.value1;
					  }
					  else
					  {
						  currSE.value2 = LIST_0;
						  dP->readSyntaxElement (&currSE,img,inp,dP);
					  }
					  refframe = currSE.value1;
					  
				  }
				  else
				  {
					  refframe = 0;
				  }
				  
				  /*
				  if (bframe && refframe>img->buf_cycle)    // img->buf_cycle should be correct for field MBs now
				  {
				  set_ec_flag(SE_REFFRAME);
				  refframe = 1;
				  }
				  */
				  
				  for (j=j0; j<j0+step_v0;j++)
					  for (i=i0; i<i0+step_h0;i++)
					  {
						  dec_picture->ref_idx[LIST_0][img->block_x + i][img->block_y + j] = refframe;
					  }
					  
			  }
		  }
	  }
  }
  else
  {
	  for (j0=0; j0<4; j0+=step_v0)
	  {
		  for (i0=0; i0<4; i0+=step_h0)
		  {
			  k=2*(j0/2)+(i0/2);
			  if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
			  {
				  for (j=j0; j<j0+step_v0;j++)
					  for (i=i0; i<i0+step_h0;i++)
					  {
						  dec_picture->ref_idx[LIST_0][img->block_x + i][img->block_y + j] = 0;
					  }
			  }
		  }
	  }
  }
  
  //  If backward multiple ref. frames, read backward reference frame for the MB *********************************
  if(img->num_ref_idx_l1_active>1)
  {
	  flag_mode = ( img->num_ref_idx_l1_active == 2 ? 1 : 0);
	  
	  currSE.type = SE_REFFRAME;
	  dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);
	  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping_dec = linfo_ue_dec;
	  else                                                      currSE.reading = readRefFrame_CABAC;
	  
	  for (j0=0; j0<4; j0+=step_v0)
	  {
		  for (i0=0; i0<4; i0+=step_h0)
		  {
			  k=2*(j0/2)+(i0/2);
			  if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
			  {
				  TRACE_STRING("ref_idx_l1");
				  
				  img->subblock_x = i0;
				  img->subblock_y = j0;
				  
				  currSE.context = BType2CtxRef (currMB->b8mode[k]);
				  if( (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) && flag_mode )
				  {
					  currSE.len = 1;
					  readSyntaxElement_FLC_dec(&currSE, dP->bitstream);
					  currSE.value1 = 1-currSE.value1;
				  }
				  else
				  {
					  currSE.value2 = LIST_1;
					  dP->readSyntaxElement (&currSE,img,inp,dP);
				  }
				  refframe = currSE.value1;
				  
				  for (j=j0; j<j0+step_v0;j++)
				  {
					  for (i=i0; i<i0+step_h0;i++)
					  {
						  dec_picture->ref_idx[LIST_1][img->block_x + i][img->block_y + j] = refframe;
					  }
				  }
			  }
		  }
	  }
  }
  else
  {
	  for (j0=0; j0<4; j0+=step_v0)
	  {
		  for (i0=0; i0<4; i0+=step_h0)
		  {
			  k=2*(j0/2)+(i0/2);
			  if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
			  {
				  for (j=j0; j<j0+step_v0;j++)
					  for (i=i0; i<i0+step_h0;i++)
					  {
						  dec_picture->ref_idx[LIST_1][img->block_x + i][img->block_y + j] = 0;
					  }
			  }
		  }
	  }
  }
  
  //=====  READ FORWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP = &(currSlice->partArr[partMap[SE_MVD]]);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping_dec = linfo_se_dec;
  else                                                  currSE.reading = readMVD_CABAC;
  
  for (j0=0; j0<4; j0+=step_v0)
	  for (i0=0; i0<4; i0+=step_h0)
	  {
		  k=2*(j0/2)+(i0/2);
		  
		  if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && (currMB->b8mode[k] !=0))//has forward vector
		  {
			  mv_mode  = currMB->b8mode[k];
			  step_h   = BLOCK_STEP [mv_mode][0];
			  step_v   = BLOCK_STEP [mv_mode][1];
			  
			  refframe = dec_picture->ref_idx[LIST_0][img->block_x+i0][img->block_y+j0];
			  
			  for (j=j0; j<j0+step_v0; j+=step_v)
			  {
				  for (i=i0; i<i0+step_h0; i+=step_h)
				  {
					  j4 = img->block_y+j;
					  i4 = img->block_x+i;
					  
					  // first make mv-prediction
					  SetMotionVectorPredictor (img, pmv, pmv+1, refframe, LIST_0, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v);
					  
					  for (k=0; k < 2; k++) 
					  {
						  TRACE_STRING("mvd_l0");
						  
						  img->subblock_x = i; // position used for context determination
						  img->subblock_y = j; // position used for context determination
						  currSE.value2 = k<<1; // identifies the component; only used for context determination
						  dP->readSyntaxElement(&currSE,img,inp,dP);
						  curr_mvd = currSE.value1; 
						  
						  vec=curr_mvd+pmv[k];           /* find motion vector */
						  
						  for(ii=0;ii<step_h;ii++)
						  {
							  for(jj=0;jj<step_v;jj++)
							  {
								  dec_picture->mv  [LIST_0][i4+ii][j4+jj][k] = vec;
								  currMB->mvd      [LIST_0][j+jj] [i+ii] [k] = curr_mvd;
							  }
						  }
					  }
				  }
			  }
		  }
		  else if (currMB->b8mode[k=2*(j0/2)+(i0/2)]==0)      
		  {  
			  if (!img->direct_type)
			  {
				  int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
				  int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2 : img->block_y;
				  
				  int refList = (co_located_ref_idx[LIST_0][img->block_x+i0][imgblock_y+j0]== -1 ? LIST_1 : LIST_0);
				  int ref_idx =  co_located_ref_idx[refList][img->block_x+i0][imgblock_y+j0];          
				  
				  if (ref_idx==-1)
				  {
					  for (j=j0; j<j0+step_v0; j++)
						  for (i=i0; i<i0+step_h0; i++)
						  {            
							  dec_picture->ref_idx [LIST_1][img->block_x+i][img->block_y+j]=0;
							  dec_picture->ref_idx [LIST_0][img->block_x+i][img->block_y+j]=0; 
							  j4 = img->block_y+j;
							  i4 = img->block_x+i;            
							  for (ii=0; ii < 2; ii++) 
							  {                                    
								  dec_picture->mv [LIST_0][i4][j4][ii]=0;
								  dec_picture->mv [LIST_1][i4][j4][ii]=0;                  
							  }
						  }
				  }
				  else 
				  {        
					  int mapped_idx=-1, iref;                             
					  int j6;
					  
					  for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
					  {
						  
						  if (dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][img->block_x+i0][imgblock_y+j0])
						  {
							  mapped_idx=iref;
							  break;
						  }
						  else //! invalid index. Default to zero even though this case should not happen
							  mapped_idx=INVALIDINDEX;
					  }
					  
					  if (INVALIDINDEX == mapped_idx)
					  {
						  error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
					  }
					  
					  
					  for (j=j0; j<j0+step_v0; j++)
						  for (i=i0; i<i0+step_h0; i++)
						  {
							  {
								  mv_scale = img->mvscale[LIST_0 + list_offset][mapped_idx];
								  
								  dec_picture->ref_idx [LIST_0][img->block_x+i][img->block_y+j] = mapped_idx;
								  dec_picture->ref_idx [LIST_1][img->block_x+i][img->block_y+j] = 0;
								  
								  j4 = img->block_y+j;
								  j6 = imgblock_y+j;
								  i4 = img->block_x+i;
								  
								  for (ii=0; ii < 2; ii++) 
								  {              
									  //if (iTRp==0)
									  if (mv_scale == 9999 || listX[LIST_0+list_offset][mapped_idx]->is_long_term)
										  //                    if (mv_scale==9999 || Co_located->is_long_term)
									  {                      
										  dec_picture->mv  [LIST_0][i4][j4][ii]=co_located_mv[refList][i4][j6][ii];
										  dec_picture->mv  [LIST_1][i4][j4][ii]=0;
									  }
									  else
									  {
										  dec_picture->mv  [LIST_0][i4][j4][ii]=(mv_scale * co_located_mv[refList][i4][j6][ii] + 128 ) >> 8;
										  dec_picture->mv  [LIST_1][i4][j4][ii]=dec_picture->mv[LIST_0][i4][j4][ii] - co_located_mv[refList][i4][j6][ii];
									  }
								  }
							  } 
						  }
				  }  
			  } 
		  }
  }
  
  //=====  READ BACKWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP          = &(currSlice->partArr[partMap[SE_MVD]]);
  
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping_dec = linfo_se_dec;
  else                                                    currSE.reading = readMVD_CABAC;
  
  for (j0=0; j0<4; j0+=step_v0)
  {
	  for (i0=0; i0<4; i0+=step_h0)
	  {
		  k=2*(j0/2)+(i0/2);
		  if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && (currMB->b8mode[k]!=0))//has backward vector
		  {
			  mv_mode  = currMB->b8mode[k];
			  step_h   = BLOCK_STEP [mv_mode][0];
			  step_v   = BLOCK_STEP [mv_mode][1];
			  
			  refframe = dec_picture->ref_idx[LIST_1][img->block_x+i0][img->block_y+j0];
			  
			  for (j=j0; j<j0+step_v0; j+=step_v)
			  {
				  for (i=i0; i<i0+step_h0; i+=step_h)
				  {
					  j4 = img->block_y+j;
					  i4 = img->block_x+i;
					  
					  // first make mv-prediction
					  SetMotionVectorPredictor (img, pmv, pmv+1, refframe, LIST_1, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v);
					  
					  for (k=0; k < 2; k++) 
					  {
						  TRACE_STRING("mvd_l1");
						  
						  img->subblock_x = i; // position used for context determination
						  img->subblock_y = j; // position used for context determination
						  currSE.value2   = (k<<1) +1; // identifies the component; only used for context determination
						  dP->readSyntaxElement(&currSE,img,inp,dP);
						  curr_mvd = currSE.value1; 
						  
						  vec=curr_mvd+pmv[k];           /* find motion vector */
						  
						  for(ii=0;ii<step_h;ii++)
						  {
							  for(jj=0;jj<step_v;jj++)
							  {
								  dec_picture->mv  [LIST_1][i4+ii][j4+jj][k] = vec;
								  currMB->mvd      [LIST_1][j+jj] [i+ii] [k] = curr_mvd;
							  }
						  }
					  }
				  }
			  }
		  }
	  }
  }
  // record reference picture Ids for deblocking decisions
  
  for(i4=img->block_x;i4<(img->block_x+4);i4++)
	  for(j4=img->block_y;j4<(img->block_y+4);j4++)
	  {
		  if(dec_picture->ref_idx[LIST_0][i4][j4]>=0)
			  dec_picture->ref_pic_id[LIST_0][i4][j4] = dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_0 + list_offset][dec_picture->ref_idx[LIST_0][i4][j4]];
		  else
			  dec_picture->ref_pic_id[LIST_0][i4][j4] = INT64_MIN;
		  if(dec_picture->ref_idx[LIST_1][i4][j4]>=0)
			  dec_picture->ref_pic_id[LIST_1][i4][j4] = dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_1 + list_offset][dec_picture->ref_idx[LIST_1][i4][j4]];  
		  else
			  dec_picture->ref_pic_id[LIST_1][i4][j4] = INT64_MIN;  
	  }
}



/*!
************************************************************************
* \brief
*    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients 
*    
*    Luma Blocks
************************************************************************
*/
int predict_nnz_dec(struct img_par *img, int i,int j)
{
	PixelPos pix;
	
	int pred_nnz = 0;
	int cnt      = 0;
	int mb_nr    = img->current_mb_nr;
	
	// left block
	getLuma4x4Neighbour(mb_nr, i, j, -1, 0, &pix);
	/* to be inserted only for dp
	if (pix.available && active_pps->constrained_intra_pred_flag)
	{
    pix.available &= img->intra_block[pix.mb_addr];
	}
	*/  
	if (pix.available)
	{
		pred_nnz = img->nz_coeff [pix.mb_addr ][pix.x][pix.y];
		cnt++;
	}
	
	// top block
	getLuma4x4Neighbour(mb_nr, i, j, 0, -1, &pix);
	/* to be inserted only for dp
	if (pix.available && active_pps->constrained_intra_pred_flag)
	{
    pix.available &= img->intra_block[pix.mb_addr];
	}
	*/  
	if (pix.available)
	{
		pred_nnz += img->nz_coeff [pix.mb_addr ][pix.x][pix.y];
		cnt++;
	}
	
	if (cnt==2)
	{
		pred_nnz++;
		pred_nnz/=cnt; 
	}
	
	return pred_nnz;
}


/*!
************************************************************************
* \brief
*    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients 
*    
*    Chroma Blocks   
************************************************************************
*/
int predict_nnz_chroma_dec(struct img_par *img, int i,int j)
{
	PixelPos pix;
	
	int pred_nnz = 0;
	int cnt      =0;
	int mb_nr    = img->current_mb_nr;
	
	// left block
	getChroma4x4Neighbour(mb_nr, i%2, j-4, -1, 0, &pix);
	/*  to be inserted only for dp
	if (pix.available && active_pps->constrained_intra_pred_flag)
	{
    pix.available &= img->intra_block[pix.mb_addr];
	}
	*/  
	if (pix.available)
	{
		pred_nnz = img->nz_coeff [pix.mb_addr ][2 * (i/2) + pix.x][4 + pix.y];
		cnt++;
	}
	
	// top block
	getChroma4x4Neighbour(mb_nr, i%2, j-4, 0, -1, &pix);
	/*  to be inserted only for dp
	if (pix.available && active_pps->constrained_intra_pred_flag)
	{
    pix.available &= img->intra_block[pix.mb_addr];
	}
	*/  
	if (pix.available)
	{
		pred_nnz += img->nz_coeff [pix.mb_addr ][2 * (i/2) + pix.x][4 + pix.y];
		cnt++;
	}
	
	if (cnt==2)
	{
		pred_nnz++;
		pred_nnz/=cnt; 
	}
	
	return pred_nnz;
}


/*!
************************************************************************
* \brief
*    Reads coeff of an 4x4 block (CAVLC)
*
* \author
*    Karl Lillevold <karll@real.com>
*    contributions by James Au <james@ubvideo.com>
************************************************************************
*/


void readCoeff4x4_CAVLC (struct img_par *img,struct inp_par_dec *inp,
						 int reference, 
						 int i, int j, int levarr[16], int runarr[16],
						 int *number_coefficients)
{
	int mb_nr = img->current_mb_nr;
	Macroblock *currMB = &img->mb_data[mb_nr];
	SyntaxElement currSE;
	Slice *currSlice = img->currentSlice;
	DataPartition *dP;
	int *partMap = assignSE2partition[currSlice->dp_mode];
	
	
	int k, code, vlcnum;
	int numcoeff, numtrailingones, numcoeff_vlc;
	int level_two_or_higher;
	int numones, totzeros, level, cdc=0, cac=0;
	int zerosleft, ntr, dptype = 0;
	int max_coeff_num = 0, nnz;
	char type[15];
	int incVlc[] = {0,3,6,12,24,48,32768};    // maximum vlc = 6
	
	numcoeff = 0;
	
	switch (reference)
	{
	case LUMA:
		max_coeff_num = 16;
		sprintf(type, "%s", "Luma");
		if (IS_INTRA (currMB))
		{
			dptype = SE_LUM_AC_INTRA;
		}
		else
		{
			dptype = SE_LUM_AC_INTER;
		}
		break;
	case LUMA_INTRA16x16DC:
		max_coeff_num = 16;
		sprintf(type, "%s", "Lum16DC");
		dptype = SE_LUM_DC_INTRA;
		break;
	case LUMA_INTRA16x16AC:
		max_coeff_num = 15;
		sprintf(type, "%s", "Lum16AC");
		dptype = SE_LUM_AC_INTRA;
		break;
		
	case CHROMA_DC:
		max_coeff_num = 4;
		cdc = 1;
		
		sprintf(type, "%s", "ChrDC");
		if (IS_INTRA (currMB))
		{
			dptype = SE_CHR_DC_INTRA;
		}
		else
		{
			dptype = SE_CHR_DC_INTER;
		}
		break;
	case CHROMA_AC:
		max_coeff_num = 15;
		cac = 1;
		sprintf(type, "%s", "ChrAC");
		if (IS_INTRA (currMB))
		{
			dptype = SE_CHR_AC_INTRA;
		}
		else
		{
			dptype = SE_CHR_AC_INTER;
		}
		break;
	default:
		error ("readCoeff4x4_CAVLC: invalid block type", 600);
		break;
	}
	
	currSE.type = dptype;
	dP = &(currSlice->partArr[partMap[dptype]]);
	
	img->nz_coeff[img->current_mb_nr][i][j] = 0;
	
	
	if (!cdc)
	{
		// luma or chroma AC
		if (!cac)
		{
			nnz = predict_nnz_dec(img, i, j);
		}
		else
		{
			nnz = predict_nnz_chroma_dec(img, i, j);
		}
		
		if (nnz < 2)
		{
			numcoeff_vlc = 0;
		}
		else if (nnz < 4)
		{
			numcoeff_vlc = 1;
		}
		else if (nnz < 8)
		{
			numcoeff_vlc = 2;
		}
		else //
		{
			numcoeff_vlc = 3;
		}
		
		currSE.value1 = numcoeff_vlc;
		
		readSyntaxElement_NumCoeffTrailingOnes_dec(&currSE, dP, type);
		
		numcoeff =  currSE.value1;
		numtrailingones =  currSE.value2;
		
		img->nz_coeff[img->current_mb_nr][i][j] = numcoeff;
	}
	else
	{
		// chroma DC
		readSyntaxElement_NumCoeffTrailingOnesChromaDC_dec(&currSE, dP);
		
		numcoeff =  currSE.value1;
		numtrailingones =  currSE.value2;
	}
	
	
	for (k = 0; k < max_coeff_num; k++)
	{
		levarr[k] = 0;
		runarr[k] = 0;
	}
	
	numones = numtrailingones;
	*number_coefficients = numcoeff;
	
	if (numcoeff)
	{
		if (numtrailingones)
		{
			
			currSE.len = numtrailingones;
			
#if TRACE
			printf(currSE.tracestring, 
				TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)", type, i, j);
#endif
			
			readSyntaxElement_FLC_dec (&currSE, dP->bitstream);
			
			code = currSE.inf;
			ntr = numtrailingones;
			for (k = numcoeff-1; k > numcoeff-1-numtrailingones; k--)
			{
				ntr --;
				if ((code>>ntr)&1)
					levarr[k] = -1;
				else
					levarr[k] = 1;
			}
		}
		
		// decode levels
		level_two_or_higher = 1;
		if (numcoeff > 3 && numtrailingones == 3)
			level_two_or_higher = 0;
		
		if (numcoeff > 10 && numtrailingones < 3)
			vlcnum = 1;
		else
			vlcnum = 0;
		
		for (k = numcoeff - 1 - numtrailingones; k >= 0; k--)
		{
			
#if TRACE
			printf(currSE.tracestring, 
				TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d ", type,
				i, j, k, vlcnum);
#endif
			
			if (vlcnum == 0)
				readSyntaxElement_Level_VLC0_dec(&currSE, dP);
			else
				readSyntaxElement_Level_VLCN_dec(&currSE, vlcnum, dP);
			
			if (level_two_or_higher)
			{
				if (currSE.inf > 0)
					currSE.inf ++;
				else
					currSE.inf --;
				level_two_or_higher = 0;
			}
			
			level = levarr[k] = currSE.inf;
			if (abs(level) == 1)
				numones ++;
			
			// update VLC table
			if (abs(level)>incVlc[vlcnum])
				vlcnum++;
			
			if (k == numcoeff - 1 - numtrailingones && abs(level)>3)
				vlcnum = 2;
			
		}
		
		if (numcoeff < max_coeff_num)
		{
			// decode total run
			vlcnum = numcoeff-1;
			currSE.value1 = vlcnum;
			
#if TRACE
			printf(currSE.tracestring, 
				TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d ", type, i,j, vlcnum);
#endif
			if (cdc)
				readSyntaxElement_TotalZerosChromaDC_dec(&currSE, dP);
			else
				readSyntaxElement_TotalZeros_dec(&currSE, dP);
			
			totzeros = currSE.value1;
		}
		else
		{
			totzeros = 0;
		}
		
		// decode run before each coefficient
		zerosleft = totzeros;
		i = numcoeff-1;
		if (zerosleft > 0 && i > 0)
		{
			do 
			{
				// select VLC for runbefore
				vlcnum = zerosleft - 1;
				if (vlcnum > RUNBEFORE_NUM-1)
					vlcnum = RUNBEFORE_NUM-1;
				
				currSE.value1 = vlcnum;
#if TRACE
				printf(currSE.tracestring, 
					TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d ",
					type, i, j, i, vlcnum);
#endif
				
				readSyntaxElement_Run_dec(&currSE, dP);
				runarr[i] = currSE.value1;
				
				zerosleft -= runarr[i];
				i --;
			} while (zerosleft != 0 && i != 0);
		}
		runarr[i] = zerosleft;
		
  } // if numcoeff
}



/*!
************************************************************************
* \brief
*    Get coded block pattern and coefficients (run/level)
*    from the NAL
************************************************************************
*/
void readCBPandCoeffsFromNAL(struct img_par *img,struct inp_par_dec *inp)
{
	int i,j,k;
	int level;
	int mb_nr = img->current_mb_nr;
	int ii,jj;
	int i1,j1, m2,jg2;
	Macroblock *currMB = &img->mb_data[mb_nr];
	int cbp;
	SyntaxElement currSE;
	Slice *currSlice = img->currentSlice;
	DataPartition *dP;
	int *partMap = assignSE2partition[currSlice->dp_mode];
	int iii,jjj;
	int coef_ctr, i0, j0, b8;
	int ll;
	int block_x,block_y;
	int start_scan;
	int uv;
	int qp_uv;
	int run, len;
	int levarr[16], runarr[16], numcoeff;
	
	int qp_per    = (img->qp-MIN_QP)/6;
	int qp_rem    = (img->qp-MIN_QP)%6;
	int qp_per_uv = QP_SCALE_CR[img->qp-MIN_QP]/6;
	int qp_rem_uv = QP_SCALE_CR[img->qp-MIN_QP]%6;
	int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
	
	// QPI
	qp_uv = img->qp + active_pps->chroma_qp_index_offset;
	qp_uv = Clip3(0, 51, qp_uv);
	qp_per_uv = QP_SCALE_CR[qp_uv-MIN_QP]/6;
	qp_rem_uv = QP_SCALE_CR[qp_uv-MIN_QP]%6;
	
	// read CBP if not new intra mode
	if (!IS_NEWINTRA (currMB))
	{
		if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB )   currSE.type = SE_CBP_INTRA;
		else                        currSE.type = SE_CBP_INTER;
		
		dP = &(currSlice->partArr[partMap[currSE.type]]);
		
		if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
		{
			if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB)  currSE.mapping_dec = linfo_cbp_intra_dec;
			else                       currSE.mapping_dec = linfo_cbp_inter_dec;
		}
		else
		{
			currSE.reading = readCBP_CABAC;
		}
		
		TRACE_STRING("coded_block_pattern");
		
		dP->readSyntaxElement(&currSE,img,inp,dP);
		currMB->cbp = cbp = currSE.value1;
		// Delta quant only if nonzero coeffs
		if (cbp !=0)
		{
			if (IS_INTER (currMB))  currSE.type = SE_DELTA_QUANT_INTER;
			else                    currSE.type = SE_DELTA_QUANT_INTRA;
			
			dP = &(currSlice->partArr[partMap[currSE.type]]);
			
			if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
			{
				currSE.mapping_dec = linfo_se_dec;
			}
			else
                currSE.reading= readDquant_CABAC; //gabi
			
			TRACE_STRING("mb_qp_delta");
			
			dP->readSyntaxElement(&currSE,img,inp,dP);
			currMB->delta_quant = currSE.value1;
			img->qp= (img->qp-MIN_QP+currMB->delta_quant+(MAX_QP-MIN_QP+1))%(MAX_QP-MIN_QP+1)+MIN_QP;
		}
	}
	else
	{
		cbp = currMB->cbp;
	}
	
	for (i=0;i<BLOCK_SIZE;i++)
		for (j=0;j<BLOCK_SIZE;j++)
			for(iii=0;iii<BLOCK_SIZE;iii++)
				for(jjj=0;jjj<BLOCK_SIZE;jjj++)
					img->cof[i][j][iii][jjj]=0;// reset luma coeffs
				
				
				if (IS_NEWINTRA (currMB)) // read DC coeffs for new intra modes
				{
					currSE.type = SE_DELTA_QUANT_INTRA;
					
					dP = &(currSlice->partArr[partMap[currSE.type]]);
					
					if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
					{
						currSE.mapping_dec = linfo_se_dec;
					}
					else
					{
						currSE.reading= readDquant_CABAC;
					}
#if TRACE
					printf(currSE.tracestring, TRACESTRING_SIZE, "Delta quant ");
#endif
					dP->readSyntaxElement(&currSE,img,inp,dP);
					currMB->delta_quant = currSE.value1;
					img->qp= (img->qp-MIN_QP+currMB->delta_quant+(MAX_QP-MIN_QP+1))%(MAX_QP-MIN_QP+1)+MIN_QP;
					
					for (i=0;i<BLOCK_SIZE;i++)
						for (j=0;j<BLOCK_SIZE;j++)
							img->ipredmode[img->block_x+i][img->block_y+j]=DC_PRED;
						
						
						if (active_pps->entropy_coding_mode_flag == UVLC)
						{
							readCoeff4x4_CAVLC(img, inp, LUMA_INTRA16x16DC, 0, 0,
								levarr, runarr, &numcoeff);
							
							coef_ctr=-1;
							level = 1;                            // just to get inside the loop
							for(k = 0; k < numcoeff; k++)
							{
								if (levarr[k] != 0)                     // leave if len=1
								{
									coef_ctr=coef_ctr+runarr[k]+1;
									
									if ((img->structure == FRAME) && (!currMB->mb_field)) 
									{
										i0=SNGL_SCAN[coef_ctr][0];
										j0=SNGL_SCAN[coef_ctr][1];
									}
									else 
									{ // Alternate scan for field coding
										i0=FIELD_SCAN[coef_ctr][0];
										j0=FIELD_SCAN[coef_ctr][1];
									}
									
									img->cof[i0][j0][0][0]=levarr[k];// add new intra DC coeff
								}
							}
						}
						else
						{
							
							currSE.type = SE_LUM_DC_INTRA;
							dP = &(currSlice->partArr[partMap[currSE.type]]);
							
							currSE.context      = LUMA_16DC;
							currSE.type         = SE_LUM_DC_INTRA;
							img->is_intra_block = 1;
							
							if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
							{
								currSE.mapping_dec = linfo_levrun_inter_dec;
							}
							else
							{
								currSE.reading = readRunLevel_CABAC;
							}
							
							
							
							coef_ctr=-1;
							level = 1;                            // just to get inside the loop
							for(k=0;(k<17) && (level!=0);k++)
							{
#if TRACE
								printf(currSE.tracestring, TRACESTRING_SIZE, "DC luma 16x16 ");
#endif
								dP->readSyntaxElement(&currSE,img,inp,dP);
								level = currSE.value1;
								run   = currSE.value2;
								len   = currSE.len;
								
								if (level != 0)                     // leave if len=1
								{
									coef_ctr=coef_ctr+run+1;
									
									if ((img->structure == FRAME) && (!currMB->mb_field)) 
									{
										i0=SNGL_SCAN[coef_ctr][0];
										j0=SNGL_SCAN[coef_ctr][1];
									}
									else 
									{ // Alternate scan for field coding
										i0=FIELD_SCAN[coef_ctr][0];
										j0=FIELD_SCAN[coef_ctr][1];
									}
									
									img->cof[i0][j0][0][0]=level;// add new intra DC coeff
								}
							}
						}
						itrans_2(img);// transform new intra DC
  }
  
  qp_per    = (img->qp-MIN_QP)/6;
  qp_rem    = (img->qp-MIN_QP)%6;
  qp_uv = img->qp + active_pps->chroma_qp_index_offset;
  qp_uv = Clip3(0, 51, qp_uv);
  qp_per_uv = QP_SCALE_CR[qp_uv-MIN_QP]/6;
  qp_rem_uv = QP_SCALE_CR[qp_uv-MIN_QP]%6;
  currMB->qp = img->qp;
  
  // luma coefficients
  for (block_y=0; block_y < 4; block_y += 2) /* all modes */
  {
	  for (block_x=0; block_x < 4; block_x += 2)
	  {
		  
		  b8 = 2*(block_y/2) + block_x/2;
		  if (active_pps->entropy_coding_mode_flag == UVLC)
		  {
			  for (j=block_y; j < block_y+2; j++)
			  {
				  for (i=block_x; i < block_x+2; i++)
				  {
					  ii = block_x/2; jj = block_y/2;
					  b8 = 2*jj+ii;
					  
					  if (cbp & (1<<b8))  /* are there any coeff in current block at all */
					  {
						  if (IS_NEWINTRA(currMB))
						  {
							  readCoeff4x4_CAVLC(img, inp, LUMA_INTRA16x16AC, i, j,
								  levarr, runarr, &numcoeff);
							  
							  start_scan = 1;
						  }
						  else
						  {
							  readCoeff4x4_CAVLC(img, inp, LUMA, i, j,
								  levarr, runarr, &numcoeff);
							  start_scan = 0;
						  }
						  
						  coef_ctr = start_scan-1;
						  for (k = 0; k < numcoeff; k++)
						  {
							  if (levarr[k] != 0)
							  {
								  coef_ctr             += runarr[k]+1;
								  
								  if ((img->structure == FRAME) && (!currMB->mb_field)) 
								  {
									  i0=SNGL_SCAN[coef_ctr][0];
									  j0=SNGL_SCAN[coef_ctr][1];
								  }
								  else { // Alternate scan for field coding
									  i0=FIELD_SCAN[coef_ctr][0];
									  j0=FIELD_SCAN[coef_ctr][1];
								  }
								  currMB->cbp_blk      |= 1 << ((j<<2) + i) ;
								  img->cof[i][j][i0][j0]= levarr[k]*dequant_coef[qp_rem][i0][j0]<<qp_per;
							  }
						  }
					  }
					  else
					  {
						  img->nz_coeff[img->current_mb_nr][i][j] = 0;
					  }
				  }
			  }
		  } // VLC
		  else
		  { 
			  b8 = 2*(block_y/2) + block_x/2;
			  // CABAC 
			  for (j=block_y; j < block_y+2; j++)
			  {
				  //jj=j/2;
				  for (i=block_x; i < block_x+2; i++)
				  {
					  //ii = i/2;
					  //b8 = 2*jj+ii;
					  
					  if (IS_NEWINTRA (currMB))   start_scan = 1; // skip DC coeff
					  else                        start_scan = 0; // take all coeffs
					  
					  img->subblock_x = i; // position for coeff_count ctx
					  img->subblock_y = j; // position for coeff_count ctx
					  if (cbp & (1<<b8))  // are there any coeff in current block at all
					  {
						  coef_ctr = start_scan-1;
						  level    = 1;
						  for(k=start_scan;(k<17) && (level!=0);k++)
						  {
						  /*
						  * make distinction between INTRA and INTER coded
						  * luminance coefficients
							  */
							  currSE.context      = (IS_NEWINTRA(currMB) ? LUMA_16AC : LUMA_4x4);
							  currSE.type         = (IS_INTRA(currMB) ?
								  (k==0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) :
							  (k==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER));
							  img->is_intra_block = IS_INTRA(currMB);
							  
#if TRACE
							  sprintf(currSE.tracestring, "Luma sng ");
#endif
							  dP = &(currSlice->partArr[partMap[currSE.type]]);
							  
							  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)  currSE.mapping_dec = linfo_levrun_inter_dec;
							  else                                                     currSE.reading = readRunLevel_CABAC;
							  
							  dP->readSyntaxElement(&currSE,img,inp,dP);
							  level = currSE.value1;
							  run   = currSE.value2;
							  len   = currSE.len;
							  
							  if (level != 0)    /* leave if len=1 */
							  {
								  coef_ctr             += run+1;
								  
								  if ((img->structure == FRAME) && (!currMB->mb_field)) 
								  {
									  i0=SNGL_SCAN[coef_ctr][0];
									  j0=SNGL_SCAN[coef_ctr][1];
								  }
								  else { // Alternate scan for field coding
									  i0=FIELD_SCAN[coef_ctr][0];
									  j0=FIELD_SCAN[coef_ctr][1];
								  }
								  currMB->cbp_blk      |= 1 << ((j<<2) + i) ;
								  img->cof[i][j][i0][j0]= level*dequant_coef[qp_rem][i0][j0]<<qp_per;
							  }
						  }
					  }
				  }
			  }
		  } 
    }
  }
  
  for (j=4;j<6;j++) // reset all chroma coeffs before read
	  for (i=0;i<4;i++)
		  for (iii=0;iii<4;iii++)
			  for (jjj=0;jjj<4;jjj++)
				  img->cof[i][j][iii][jjj]=0;
			  
			  m2 =img->mb_x*2;
			  jg2=img->mb_y*2;
			  
			  // chroma 2x2 DC coeff
			  if(cbp>15)
			  {
				  for (ll=0;ll<3;ll+=2)
				  {
					  for (i=0;i<4;i++)
						  img->cofu[i]=0;
					  
					  
					  if (active_pps->entropy_coding_mode_flag == UVLC)
					  {
						  
						  readCoeff4x4_CAVLC(img, inp, CHROMA_DC, 0, 0,
							  levarr, runarr, &numcoeff);
						  coef_ctr=-1;
						  level=1;
						  for(k = 0; k < numcoeff; k++)
						  {
							  if (levarr[k] != 0)
							  {
								  currMB->cbp_blk |= 0xf0000 << (ll<<1) ;
								  coef_ctr=coef_ctr+runarr[k]+1;
								  img->cofu[coef_ctr]=levarr[k];
							  }
						  }
					  }
					  else
					  {
						  coef_ctr=-1;
						  level=1;
						  for(k=0;(k<5)&&(level!=0);k++)
						  {
							  currSE.context      = CHROMA_DC;
							  currSE.type         = (IS_INTRA(currMB) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER);
							  img->is_intra_block =  IS_INTRA(currMB);
							  img->is_v_block     = ll;
							  
#if TRACE
							  printf(currSE.tracestring, TRACESTRING_SIZE, " 2x2 DC Chroma ");
#endif
							  dP = &(currSlice->partArr[partMap[currSE.type]]);
							  
							  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
								  currSE.mapping_dec = linfo_levrun_c2x2_dec;
							  else
								  currSE.reading = readRunLevel_CABAC;
							  
							  dP->readSyntaxElement(&currSE,img,inp,dP);
							  level = currSE.value1;
							  run = currSE.value2;
							  len = currSE.len;
							  if (level != 0)
							  {
								  currMB->cbp_blk |= 0xf0000 << (ll<<1) ;
								  coef_ctr=coef_ctr+run+1;
								  // Bug: img->cofu has only 4 entries, hence coef_ctr MUST be <4 (which is
								  // caught by the assert().  If it is bigger than 4, it starts patching the
								  // img->predmode pointer, which leads to bugs later on.
								  //
								  // This assert() should be left in the code, because it captures a very likely
								  // bug early when testing in error prone environments (or when testing NAL
								  // functionality).
								  assert (coef_ctr < 4);
								  img->cofu[coef_ctr]=level;
							  }
						  }
					  }
					  
					  if (smb) // check to see if MB type is SPred or SIntra4x4 
					  {
						  img->cof[0+ll][4][0][0]=img->cofu[0];   img->cof[1+ll][4][0][0]=img->cofu[1];
						  img->cof[0+ll][5][0][0]=img->cofu[2];   img->cof[1+ll][5][0][0]=img->cofu[3];
					  }
					  else
					  { 
						  for (i=0;i<4;i++)
							  img->cofu[i]*=dequant_coef[qp_rem_uv][0][0]<<qp_per_uv;
						  img->cof[0+ll][4][0][0]=(img->cofu[0]+img->cofu[1]+img->cofu[2]+img->cofu[3])>>1;
						  img->cof[1+ll][4][0][0]=(img->cofu[0]-img->cofu[1]+img->cofu[2]-img->cofu[3])>>1;
						  img->cof[0+ll][5][0][0]=(img->cofu[0]+img->cofu[1]-img->cofu[2]-img->cofu[3])>>1;
						  img->cof[1+ll][5][0][0]=(img->cofu[0]-img->cofu[1]-img->cofu[2]+img->cofu[3])>>1;
					  }
				  }
			  }
			  
			  // chroma AC coeff, all zero fram start_scan
			  if (cbp<=31)
				  for (j=4; j < 6; j++)
					  for (i=0; i < 4; i++)
						  img->nz_coeff [img->current_mb_nr ][i][j]=0;
					  
					  
					  // chroma AC coeff, all zero fram start_scan
					  uv=-1;
					  if (cbp>31)
					  {
						  block_y=4;
						  for (block_x=0; block_x < 4; block_x += 2)
						  {
							  for (j=block_y; j < block_y+2; j++)
							  {
								  jj=j/2;
								  j1=j-4;
								  for (i=block_x; i < block_x+2; i++)
								  {
									  
									  ii=i/2;
									  i1=i%2;
									  
									  if (active_pps->entropy_coding_mode_flag == UVLC)
									  {
										  readCoeff4x4_CAVLC(img, inp, CHROMA_AC, i, j,
											  levarr, runarr, &numcoeff);
										  coef_ctr=0;
										  level=1;
										  uv++;
										  for(k = 0; k < numcoeff;k++)
										  {
											  if (levarr[k] != 0)
											  {
												  currMB->cbp_blk |= 1 << (16 + (j1<<1) + i1 + (block_x<<1) ) ;
												  coef_ctr=coef_ctr+runarr[k]+1;
												  
												  if ((img->structure == FRAME) && (!currMB->mb_field)) 
												  {
													  i0=SNGL_SCAN[coef_ctr][0];
													  j0=SNGL_SCAN[coef_ctr][1];
												  }
												  else { // Alternate scan for field coding
													  i0=FIELD_SCAN[coef_ctr][0];
													  j0=FIELD_SCAN[coef_ctr][1];
												  }
												  img->cof[i][j][i0][j0]=levarr[k]*dequant_coef[qp_rem_uv][i0][j0]<<qp_per_uv;
											  }
										  }
									  }
									  
									  else
									  {
										  coef_ctr=0;
										  level=1;
										  uv++;
										  
										  img->subblock_y = j/5;
										  img->subblock_x = i%2;
										  
										  for(k=0;(k<16)&&(level!=0);k++)
										  {
											  currSE.context      = CHROMA_AC;
											  currSE.type         = (IS_INTRA(currMB) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER);
											  img->is_intra_block =  IS_INTRA(currMB);
											  img->is_v_block     = (uv>=4);
											  
#if TRACE
											  printf(currSE.tracestring, TRACESTRING_SIZE, " AC Chroma ");
#endif
											  dP = &(currSlice->partArr[partMap[currSE.type]]);
											  
											  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
												  currSE.mapping_dec = linfo_levrun_inter_dec;
											  else
												  currSE.reading = readRunLevel_CABAC;
											  
											  dP->readSyntaxElement(&currSE,img,inp,dP);
											  level = currSE.value1;
											  run = currSE.value2;
											  len = currSE.len;
											  
											  if (level != 0)
											  {
												  currMB->cbp_blk |= 1 << (16 + (j1<<1) + i1 + (block_x<<1) ) ;
												  coef_ctr=coef_ctr+run+1;
												  
												  if ((img->structure == FRAME) && (!currMB->mb_field)) 
												  {
													  i0=SNGL_SCAN[coef_ctr][0];
													  j0=SNGL_SCAN[coef_ctr][1];
												  }
												  else { // Alternate scan for field coding
													  i0=FIELD_SCAN[coef_ctr][0];
													  j0=FIELD_SCAN[coef_ctr][1];
												  }
												  img->cof[i][j][i0][j0]=level*dequant_coef[qp_rem_uv][i0][j0]<<qp_per_uv;
											  }
										  }
									  }
								  }
							  }
						  }
					  }
}


/*!
************************************************************************
* \brief
*    Copy IPCM coefficients to decoded picture buffer and set parameters for this MB
*    (for IPCM CABAC and IPCM CAVLC  28/11/2003)
*
* \author
*    Dong Wang <Dong.Wang@bristol.ac.uk>
************************************************************************
*/

void decode_ipcm_mb(struct img_par *img)
{
	int i,j;
	
	Macroblock *currMb = &img->mb_data[img->current_mb_nr];    
	
	//Copy coefficents to decoded picture buffer
	//IPCM coefficents are stored in img->cof which is set in function readIPCMcoeffsFromNAL()
	
	for(i=0;i<16;i++)
		for(j=0;j<16;j++)
			dec_picture->imgY[img->pix_y+i][img->pix_x+j]=img->cof[i/4][j/4][i%4][j%4];
		
		for(i=0;i<8;i++)
			for(j=0;j<8;j++)
				dec_picture->imgUV[0][img->pix_c_y+i][img->pix_c_x+j]=img->cof[i/4][j/4+4][i%4][j%4];
			
			for(i=0;i<8;i++)
				for(j=0;j<8;j++)
					dec_picture->imgUV[1][img->pix_c_y+i][img->pix_c_x+j]=img->cof[i/4+2][j/4+4][i%4][j%4];
				
				
				
				//For Deblocking Filter  16/08/2003
				if (currMb->mb_type==IPCM)
					currMb->qp=0;
				
				//For CAVLC
				//Set the nz_coeff to 16. 
				//These parameters are to be used in CAVLC decoding of neighbour blocks
				for(i=0;i<4;i++)
					for (j=0;j<6;j++)
						img->nz_coeff[img->current_mb_nr][i][j]=16;
					
					
					//For CABAC decoding of MB skip flag 
					if (currMb->mb_type==IPCM)
						currMb->skip_flag=1;
					
					//for Loop filter CABAC
					if (currMb->mb_type==IPCM)
						currMb->cbp_blk=0xFFFF;
					
					//For CABAC decoding of Dquant
					last_dquant=0;
}

/*!
 ************************************************************************
 * \brief
 *    decode one macroblock
 ************************************************************************
 */

int decode_one_macroblock_dec(struct img_par *img,struct inp_par_dec *inp)
{
  int tmp_block[BLOCK_SIZE][BLOCK_SIZE];
  int tmp_blockbw[BLOCK_SIZE][BLOCK_SIZE];
  int i=0,j=0,k,l,ii=0,jj=0,i1=0,j1=0,j4=0,i4=0;
  int uv, hv;
  int vec1_x=0,vec1_y=0,vec2_x=0,vec2_y=0;
  int ioff,joff;
  int block8x8;   // needed for ABT

  int bw_pred=0, fw_pred=0, pred, ifx;
  int ii0,jj0,ii1,jj1,if1,jf1,if0,jf0;
  int mv_mul,f1,f2,f3,f4;

  const byte decode_block_scan[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};

  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int ref_idx, fw_refframe=-1, bw_refframe=-1, mv_mode, pred_dir, intra_prediction; // = currMB->ref_frame;
  int fw_ref_idx=-1, bw_ref_idx=-1;

  int  *** mv_array, ***fw_mv_array, ***bw_mv_array;

  int mv_scale;

  int mb_nr             = img->current_mb_nr;
  int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
  int list_offset;
  int max_y_cr;

  StorablePicture **list;

  int jf;

  int fw_rFrame=-1,bw_rFrame=-1;
  int pmvfw[2]={0,0},pmvbw[2]={0,0};              

  int direct_pdir=-1;

  int curr_mb_field = ((img->MbaffFrameFlag)&&(currMB->mb_field));
  
  byte **     moving_block;
  int ****     co_located_mv;
  int ***     co_located_ref_idx;
  int64 ***    co_located_ref_id;

  if(currMB->mb_type==IPCM)
  {
    //copy readed data into imgY and set parameters
    decode_ipcm_mb(img);
    return 0;
  }

//////////////////////////

  // find out the correct list offsets
  if (curr_mb_field)
  {
    if(mb_nr%2)
    {
      list_offset = 4; // top field mb
      moving_block = Co_located->bottom_moving_block;
      co_located_mv = Co_located->bottom_mv;
      co_located_ref_idx = Co_located->bottom_ref_idx;
      co_located_ref_id = Co_located->bottom_ref_pic_id;
    }
    else
    {
      list_offset = 2; // bottom field mb
      moving_block = Co_located->top_moving_block;
      co_located_mv = Co_located->top_mv;
      co_located_ref_idx = Co_located->top_ref_idx;
      co_located_ref_id = Co_located->top_ref_pic_id;
    }
    max_y_cr = dec_picture->size_y_cr/2-1;
  }
  else
  {
    list_offset = 0;  // no mb aff or frame mb
    moving_block = Co_located->moving_block;
    co_located_mv = Co_located->mv;
    co_located_ref_idx = Co_located->ref_idx;
    co_located_ref_id = Co_located->ref_pic_id;
    max_y_cr = dec_picture->size_y_cr-1;
  }



  if (!img->MbaffFrameFlag)
  {
    for (l=0+list_offset;l<(2+list_offset);l++)
    {
      for(k = 0; k < listXsize[l]; k++)
      {
        listX[l][k]->chroma_vector_adjustment= 0;
        if(img->structure == TOP_FIELD && img->structure != listX[l][k]->structure)
          listX[l][k]->chroma_vector_adjustment = -2;
        if(img->structure == BOTTOM_FIELD && img->structure != listX[l][k]->structure)
          listX[l][k]->chroma_vector_adjustment = 2;
      }
    }
  }
  else
  {
    if (curr_mb_field)
    {
      for (l=0+list_offset;l<(2+list_offset);l++)
      {
        for(k = 0; k < listXsize[l]; k++)
        {
          listX[l][k]->chroma_vector_adjustment= 0;
          if(img->current_mb_nr % 2 == 0 && listX[l][k]->structure == BOTTOM_FIELD)
            listX[l][k]->chroma_vector_adjustment = -2;
          if(img->current_mb_nr % 2 == 1 && listX[l][k]->structure == TOP_FIELD)
            listX[l][k]->chroma_vector_adjustment = 2;
        }
      }
    }
    else
    {
      for (l=0+list_offset;l<(2+list_offset);l++)
      {
        for(k = 0; k < listXsize[l]; k++)
        {
          listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
  }

  mv_mul=4;
  f1=8;
  f2=7;

  f3=f1*f1;
  f4=f3/2;

  // luma decoding **************************************************

  // get prediction for INTRA_MB_16x16
  if (IS_NEWINTRA (currMB))
  {
    intrapred_luma_16x16_dec(img, currMB->i16mode);
  }

  if (img->type==B_SLICE && img->direct_type && (IS_DIRECT (currMB) || 
    (IS_P8x8(currMB) && !(currMB->b8mode[0] && currMB->b8mode[1] && currMB->b8mode[2] && currMB->b8mode[3]))))
  {
    int fw_rFrameL, fw_rFrameU, fw_rFrameUL, fw_rFrameUR;
    int bw_rFrameL, bw_rFrameU, bw_rFrameUL, bw_rFrameUR;    
    
    PixelPos mb_left, mb_up, mb_upleft, mb_upright;              
    
    getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_left);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_up);
    getLuma4x4Neighbour(img->current_mb_nr,0,0,16, -1,&mb_upright);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, -1,-1,&mb_upleft);

    if (!img->MbaffFrameFlag)
    {
      fw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : -1;
      fw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
      fw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
      fw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
      
      bw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
      bw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
      bw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
      bw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
    }
    else
    {
      if (img->mb_data[img->current_mb_nr].mb_field)
      {
        fw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field  || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] < 0? 
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : 
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] * 2: -1;

        fw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] < 0? 
          dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : 
        dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] * 2: -1;

        fw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0?         
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : 
        dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      

        fw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] < 0?
          dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : 
        dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] * 2: fw_rFrameUL;      
        
        bw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] < 0? 
          dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] * 2: -1;

        bw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] < 0? 
          dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] * 2: -1;

        bw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] < 0?         
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      

        bw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] < 0?         
          dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] * 2: bw_rFrameUL;              
      }
      else
      {
        fw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]  < 0 ?
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] >> 1 : 
        dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]: -1;
        
        fw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] >> 1 :  
        dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
        
        fw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y]>> 1 : 
        dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      
        
        fw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] < 0 ? 
          dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] >> 1 :  
        dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
        
        bw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] >> 1 :  
        dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
        
        bw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] >> 1 : 
        dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
        
        bw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] >> 1 : 
        dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      
        
        bw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] >> 1: 
        dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
      }
    }
    

    fw_rFrame = (fw_rFrameL >= 0 && fw_rFrameU >= 0) ? min(fw_rFrameL,fw_rFrameU): max(fw_rFrameL,fw_rFrameU);
    fw_rFrame = (fw_rFrame >= 0 && fw_rFrameUR >= 0) ? min(fw_rFrame,fw_rFrameUR): max(fw_rFrame,fw_rFrameUR);
    
    bw_rFrame = (bw_rFrameL >= 0 && bw_rFrameU >= 0) ? min(bw_rFrameL,bw_rFrameU): max(bw_rFrameL,bw_rFrameU);
    bw_rFrame = (bw_rFrame >= 0 && bw_rFrameUR >= 0) ? min(bw_rFrame,bw_rFrameUR): max(bw_rFrame,bw_rFrameUR);        
    
    if (fw_rFrame >=0)
      SetMotionVectorPredictor (img, pmvfw, pmvfw+1, fw_rFrame, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
    
    if (bw_rFrame >=0)
      SetMotionVectorPredictor (img, pmvbw, pmvbw+1, bw_rFrame, LIST_1, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);

  }

  for (block8x8=0; block8x8<4; block8x8++)
  {
    for (k = block8x8*4; k < block8x8*4+4; k ++)
    {
      i = (decode_block_scan[k] & 3);
      j = ((decode_block_scan[k] >> 2) & 3);
      
      ioff=i*4;
      i4=img->block_x+i;
      
      joff=j*4;
      j4=img->block_y+j;

      mv_mode  = currMB->b8mode[2*(j/2)+(i/2)];
      pred_dir = currMB->b8pdir[2*(j/2)+(i/2)];
      
      assert (pred_dir<=2);

      // PREDICTION
      if (mv_mode==IBLOCK)
      {
        //===== INTRA PREDICTION =====
        if (intrapred(img,ioff,joff,i4,j4)==SEARCH_SYNC)  /* make 4x4 prediction block mpr from given prediction img->mb_mode */
          return SEARCH_SYNC;                   /* bit error */
      }
      else if (!IS_NEWINTRA (currMB))
      {
        if (pred_dir != 2)
        {
          //===== FORWARD/BACKWARD PREDICTION =====
          fw_refframe = ref_idx  = dec_picture->ref_idx[LIST_0 + pred_dir][i4][j4];
          mv_array = dec_picture->mv[LIST_0 + pred_dir];
          list     = listX[0+list_offset+ pred_dir];
          vec1_x = i4*4*mv_mul + mv_array[i4][j4][0];

          if (!curr_mb_field)
          {
            vec1_y = j4*4*mv_mul + mv_array[i4][j4][1];
          }
          else
          {
            if (mb_nr%2 == 0) 
              vec1_y = (img->block_y * 2 + joff) * mv_mul + mv_array[i4][j4][1];
            else
              vec1_y = ((img->block_y-4) * 2 + joff)* mv_mul + mv_array[i4][j4][1];
          }

          get_block (ref_idx, list, vec1_x, vec1_y, img, tmp_block);

          if (img->apply_weights)
          {
            if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
            {
              ref_idx >>=1;
            }

            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = Clip1(((img->wp_weight[pred_dir][ref_idx][0] *  tmp_block[ii][jj]+ img->wp_round_luma) >>img->luma_log2_weight_denom)  + img->wp_offset[pred_dir][fw_refframe>>curr_mb_field][0] );
          }
          else
          {
            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
          }
        }
        else
        {
          if (mv_mode != 0)
          {
            //===== BI-DIRECTIONAL PREDICTION =====
            fw_mv_array = dec_picture->mv[LIST_0];
            bw_mv_array = dec_picture->mv[LIST_1];
            
            fw_refframe = dec_picture->ref_idx[LIST_0][i4][j4];
            bw_refframe = dec_picture->ref_idx[LIST_1][i4][j4];
            fw_ref_idx = fw_refframe;
            bw_ref_idx = bw_refframe;
          }
          else
          {
            //===== DIRECT PREDICTION =====
            fw_mv_array = dec_picture->mv[LIST_0];
            bw_mv_array = dec_picture->mv[LIST_1];
            bw_refframe = 0;

            if (img->direct_type )
            {
              int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
              int j6=imgblock_y+j;


              if (fw_rFrame >=0)
              {
                if (!fw_rFrame  && ((!moving_block[i4][j6]) && (!listX[1+list_offset][0]->is_long_term)))
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]= 0;
                  dec_picture->mv  [LIST_0][i4][j4][1]= 0;
                  dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                }
                else
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]= pmvfw[0];
                  dec_picture->mv  [LIST_0][i4][j4][1]= pmvfw[1];
                  dec_picture->ref_idx[LIST_0][i4][j4] = fw_rFrame;
                }
              }
              else
              {
                dec_picture->ref_idx[LIST_0][i4][j4] = -1;                
                dec_picture->mv  [LIST_0][i4][j4][0]= 0;
                dec_picture->mv  [LIST_0][i4][j4][1]= 0;
              }
              
              if (bw_rFrame >=0)
              {
                if  (bw_rFrame==0 && ((!moving_block[i4][j6]) && (!listX[1+list_offset][0]->is_long_term)))
                {                  
                  
                  dec_picture->mv  [LIST_1][i4][j4][0]= 0;
                  dec_picture->mv  [LIST_1][i4][j4][1]= 0;
                  dec_picture->ref_idx[LIST_1][i4][j4] = bw_rFrame;
                  
                }
                else
                {
                  dec_picture->mv  [LIST_1][i4][j4][0]= pmvbw[0];
                  dec_picture->mv  [LIST_1][i4][j4][1]= pmvbw[1];
                  
                  dec_picture->ref_idx[LIST_1][i4][j4] = bw_rFrame;
                }               
              }
              else
              {                  
                dec_picture->mv  [LIST_1][i4][j4][0]=0;
                dec_picture->mv  [LIST_1][i4][j4][1]=0;
                dec_picture->ref_idx[LIST_1][i4][j4] = -1;
                
              }
              
              if (fw_rFrame < 0 && bw_rFrame < 0)
              {
                dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                dec_picture->ref_idx[LIST_1][i4][j4] = 0;
              }
              
              fw_refframe = (dec_picture->ref_idx[LIST_0][i4][j4]!=-1) ? dec_picture->ref_idx[LIST_0][i4][j4]:0;
              bw_refframe = (dec_picture->ref_idx[LIST_1][i4][j4]!=-1) ? dec_picture->ref_idx[LIST_1][i4][j4]:0;
              
              fw_ref_idx = fw_refframe;
              bw_ref_idx = bw_refframe;
              
              if      (dec_picture->ref_idx[LIST_1][i4][j4]==-1) direct_pdir = 0;
              else if (dec_picture->ref_idx[LIST_0][i4][j4]==-1) direct_pdir = 1;
              else                                               direct_pdir = 2;
              
            }
            else // Temporal Mode
            {
              
              int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
              int j6= imgblock_y + j;
              
              int refList = (co_located_ref_idx[LIST_0][i4][j6]== -1 ? LIST_1 : LIST_0);
              int ref_idx =  co_located_ref_idx[refList][i4][j6];


              if(ref_idx==-1) // co-located is intra mode
              {
                for(hv=0; hv<2; hv++)
                {
                  dec_picture->mv  [LIST_0][i4][j4][hv]=0;
                  dec_picture->mv  [LIST_1][i4][j4][hv]=0;                    
                }
                
                dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                dec_picture->ref_idx[LIST_1][i4][j4] = 0;
                
                fw_refframe = 0;
                fw_ref_idx = 0;
              }
              else // co-located skip or inter mode
              {
                int mapped_idx=0;
                int iref;          
                
                {
                  for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
                  {
                    if (dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][i4][j6])
                    {
                      mapped_idx=iref;
                      break;
                    }
                    else //! invalid index. Default to zero even though this case should not happen
                    {                        
                      mapped_idx=INVALIDINDEX;
                    }
                  }
                  if (INVALIDINDEX == mapped_idx)
                  {   
                    error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
                  }
                }
          
                fw_ref_idx = mapped_idx;
                mv_scale = img->mvscale[LIST_0 + list_offset][mapped_idx];

                //! In such case, an array is needed for each different reference.
//                if (mv_scale == 9999 || Co_located->is_long_term)
                if (mv_scale == 9999 || listX[LIST_0+list_offset][mapped_idx]->is_long_term)
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]=co_located_mv[refList][i4][j6][0];
                  dec_picture->mv  [LIST_0][i4][j4][1]=co_located_mv[refList][i4][j6][1];

                  dec_picture->mv  [LIST_1][i4][j4][0]=0;
                  dec_picture->mv  [LIST_1][i4][j4][1]=0;
                }
                else
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]=(mv_scale * co_located_mv[refList][i4][j6][0] + 128 ) >> 8;
                  dec_picture->mv  [LIST_0][i4][j4][1]=(mv_scale * co_located_mv[refList][i4][j6][1] + 128 ) >> 8;
                  
                  dec_picture->mv  [LIST_1][i4][j4][0]=dec_picture->mv  [LIST_0][i4][j4][0] - co_located_mv[refList][i4][j6][0] ;
                  dec_picture->mv  [LIST_1][i4][j4][1]=dec_picture->mv  [LIST_0][i4][j4][1] - co_located_mv[refList][i4][j6][1] ;
                }
                
                fw_refframe = dec_picture->ref_idx[LIST_0][i4][j4] = mapped_idx; //listX[1][0]->ref_idx[refList][i4][j4];
                bw_refframe = dec_picture->ref_idx[LIST_1][i4][j4] = 0;
                
                fw_ref_idx = fw_refframe;
                bw_ref_idx = bw_refframe;
              }
            }
            // store reference picture ID determined by direct mode
            dec_picture->ref_pic_id[LIST_0][i4][j4] = dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_0 + list_offset][dec_picture->ref_idx[LIST_0][i4][j4]];
            dec_picture->ref_pic_id[LIST_1][i4][j4] = dec_picture->ref_pic_num_dec[img->current_slice_nr][LIST_1 + list_offset][dec_picture->ref_idx[LIST_1][i4][j4]];  
          }
                 
          if (mv_mode==0 && img->direct_type )
          {
            if (dec_picture->ref_idx[LIST_0][i4][j4] >= 0)
            {
              
              vec1_x = i4*4*mv_mul + fw_mv_array[i4][j4][0];
              if (!curr_mb_field)
              {
                vec1_y = j4*4*mv_mul + fw_mv_array[i4][j4][1];
              }
              else
              {
                if (mb_nr%2 == 0)
                {
                  vec1_y = (img->block_y * 2 + joff) * mv_mul + fw_mv_array[i4][j4][1];
                }
                else
                {
                  vec1_y = ((img->block_y-4) * 2 + joff)* mv_mul + fw_mv_array[i4][j4][1];
                }
              }               
              get_block(fw_refframe, listX[0+list_offset], vec1_x, vec1_y, img, tmp_block);
            }
                  
            if (dec_picture->ref_idx[LIST_1][i4][j4] >= 0)
            {
              vec2_x = i4*4*mv_mul + bw_mv_array[i4][j4][0];
              if (!curr_mb_field)
              {
                vec2_y = j4*4*mv_mul + bw_mv_array[i4][j4][1];
              }
              else
              {
                if (mb_nr%2 == 0)
                {
                  vec2_y = (img->block_y * 2 + joff) * mv_mul + bw_mv_array[i4][j4][1];
                }
                else
                {
                  vec2_y = ((img->block_y-4) * 2 + joff)* mv_mul + bw_mv_array[i4][j4][1];
                }            
              }
              get_block(bw_refframe, listX[1+list_offset], vec2_x, vec2_y, img, tmp_blockbw);
            }              
          }
          else
          {
            vec1_x = i4*4*mv_mul + fw_mv_array[i4][j4][0];
            vec2_x = i4*4*mv_mul + bw_mv_array[i4][j4][0];
            
            if (!curr_mb_field)
            {
              vec1_y = j4*4*mv_mul + fw_mv_array[i4][j4][1];
              vec2_y = j4*4*mv_mul + bw_mv_array[i4][j4][1];
            }
            else
            {
              if (mb_nr%2 == 0)
              {
                vec1_y = (img->block_y * 2 + joff) * mv_mul + fw_mv_array[i4][j4][1];
                vec2_y = (img->block_y * 2 + joff) * mv_mul + bw_mv_array[i4][j4][1];
              }
              else
              {
                vec1_y = ((img->block_y-4) * 2 + joff)* mv_mul + fw_mv_array[i4][j4][1];
                vec2_y = ((img->block_y-4) * 2 + joff)* mv_mul + bw_mv_array[i4][j4][1];
              }
            }
            
            get_block(fw_refframe, listX[0+list_offset], vec1_x, vec1_y, img, tmp_block);
            get_block(bw_refframe, listX[1+list_offset], vec2_x, vec2_y, img, tmp_blockbw);
          }
          
          if (mv_mode==0 && img->direct_type && direct_pdir==0)
          {
            if (img->apply_weights)
            {
              if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
              {
                fw_ref_idx >>=1;
              }
              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = Clip1(((tmp_block[ii][jj] * img->wp_weight[0][fw_ref_idx][0]  
                                               + img->wp_round_luma)>>img->luma_log2_weight_denom) 
                                               + img->wp_offset[0][fw_refframe>>curr_mb_field][0]);
            }
            else
            {
              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
            }
          }
          else if (mv_mode==0 && img->direct_type && direct_pdir==1)
          {              
            if (img->apply_weights)
            {
              if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
              {
                fw_ref_idx >>=1;
                bw_ref_idx >>=1;
              }

              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = Clip1(((tmp_blockbw[ii][jj] * img->wp_weight[1][bw_ref_idx][0] 
                                               + img->wp_round_luma)>>img->luma_log2_weight_denom) 
                                               + img->wp_offset[1][bw_refframe>>curr_mb_field][0]);
            }
            else
            {
              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = tmp_blockbw[ii][jj];
            }
          }
          else if(img->apply_weights)
          {
            int alpha_fw, alpha_bw;
            int wt_list_offset = (active_pps->weighted_bipred_idc==2)?list_offset:0;

            if (mv_mode==0 && img->direct_type==0 )bw_ref_idx=0;    //temporal direct 

            if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
              (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
            {
              fw_ref_idx >>=1;
              bw_ref_idx >>=1;
            }

            alpha_fw = img->wbp_weight[0+wt_list_offset][fw_ref_idx][bw_ref_idx][0];
            alpha_bw = img->wbp_weight[1+wt_list_offset][fw_ref_idx][bw_ref_idx][0];

            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = (int)Clip1(((alpha_fw * tmp_block[ii][jj] + alpha_bw * tmp_blockbw[ii][jj]  
                                             + (1<<img->luma_log2_weight_denom)) >> (img->luma_log2_weight_denom+1)) 
                                             + ((img->wp_offset[wt_list_offset+0][fw_ref_idx][0] + img->wp_offset[wt_list_offset+1][bw_ref_idx][0] + 1) >>1));
          }
          else
          {
            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = (tmp_block[ii][jj]+tmp_blockbw[ii][jj]+1)/2;
          }
        }
      }
      if (smb && mv_mode!=IBLOCK)
      {
        itrans_sp(img,ioff,joff,i,j);
      }
      else
      {
        itrans   (img,ioff,joff,i,j);      // use DCT transform and make 4x4 block m7 from prediction block mpr
      }
        
      for(ii=0;ii<BLOCK_SIZE;ii++)
      {
        for(jj=0;jj<BLOCK_SIZE;jj++)
        {
          dec_picture->imgY[j4*BLOCK_SIZE+jj][i4*BLOCK_SIZE+ii]=img->m7[ii][jj]; // contruct picture from 4x4 blocks
        }
      }
    }
  }
  
  // chroma decoding *******************************************************
  for(uv=0;uv<2;uv++)
  {
    intra_prediction = IS_INTRA (currMB);
    
    if (intra_prediction)
    {
      intrapred_chroma_dec(img, uv);
    }
    
    for (j=4;j<6;j++)
    {
      joff=(j-4)*4;
      j4=img->pix_c_y+joff;
      for(i=0;i<2;i++)
      {
        ioff=i*4;
        i4=img->pix_c_x+ioff;
        
        mv_mode  = currMB->b8mode[2*(j-4)+i];
        pred_dir = currMB->b8pdir[2*(j-4)+i];
        assert (pred_dir<=2);

        if (!intra_prediction)
        {
          if (pred_dir != 2)
          {
            //--- FORWARD/BACKWARD PREDICTION ---
            mv_array = dec_picture->mv[LIST_0 + pred_dir];
            list = listX[0+list_offset+pred_dir];
            for(jj=0;jj<4;jj++)
            {
              jf=(j4+jj)/2;
              for(ii=0;ii<4;ii++)
              {
                if1=(i4+ii)/2;
                
                fw_refframe = ref_idx   = dec_picture->ref_idx[LIST_0+pred_dir][if1][jf];
                
                i1=(img->pix_c_x+ii+ioff)*f1+mv_array[if1][jf][0];

                if (!curr_mb_field)
                  j1=(img->pix_c_y+jj+joff)*f1+mv_array[if1][jf][1];
                else
                {
                  if (mb_nr%2 == 0) 
                    j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+mv_array[if1][jf][1];
                  else
                    j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +mv_array[if1][jf][1];
                }
                
                j1 += list[ref_idx]->chroma_vector_adjustment;
                
                ii0=max (0, min (i1>>3, img->width_cr-1));
                jj0=max (0, min (j1>>3, max_y_cr));
                ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                jj1=max (0, min ((j1>>3)+1, max_y_cr));
                
                if1=(i1 & f2);
                jf1=(j1 & f2);
                if0=f1-if1;
                jf0=f1-jf1;
                
                if (img->apply_weights)
                {
                  pred = (if0*jf0*list[ref_idx]->imgUV[uv][jj0][ii0]+
                          if1*jf0*list[ref_idx]->imgUV[uv][jj0][ii1]+
                          if0*jf1*list[ref_idx]->imgUV[uv][jj1][ii0]+
                          if1*jf1*list[ref_idx]->imgUV[uv][jj1][ii1]+f4)>>6;
                  if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                    (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
                  {
                    ref_idx >>=1;
                  }

                  img->mpr[ii+ioff][jj+joff] = Clip1(((img->wp_weight[pred_dir][ref_idx][uv+1] * pred  + img->wp_round_chroma)>>img->chroma_log2_weight_denom) + img->wp_offset[pred_dir][ref_idx][uv+1]);
                }
                else
                {
                  img->mpr[ii+ioff][jj+joff]=(if0*jf0*list[ref_idx]->imgUV[uv][jj0][ii0]+
                                              if1*jf0*list[ref_idx]->imgUV[uv][jj0][ii1]+
                                              if0*jf1*list[ref_idx]->imgUV[uv][jj1][ii0]+
                                              if1*jf1*list[ref_idx]->imgUV[uv][jj1][ii1]+f4)>>6;
                }
              }
            }
          }
          else
          {
            if (mv_mode != 0)
            {
              //===== BI-DIRECTIONAL PREDICTION =====
              fw_mv_array = dec_picture->mv[LIST_0];
              bw_mv_array = dec_picture->mv[LIST_1];
            }
            else
            {
              //===== DIRECT PREDICTION =====
              fw_mv_array = dec_picture->mv[LIST_0];
              bw_mv_array = dec_picture->mv[LIST_1];
            }
            
            for(jj=0;jj<4;jj++)
            {
              jf=(j4+jj)/2;
              for(ii=0;ii<4;ii++)
              {
                ifx=(i4+ii)/2;
                
                direct_pdir = 2;
                
                if (mv_mode != 0)
                {
                  fw_refframe = dec_picture->ref_idx[LIST_0][ifx][jf];
                  bw_refframe = dec_picture->ref_idx[LIST_1][ifx][jf];
                  
                  fw_ref_idx = fw_refframe;
                  bw_ref_idx = bw_refframe;
                }
                else
                {
                  if (!mv_mode && img->direct_type )
                  {
                    if (dec_picture->ref_idx[LIST_0][(ifx/2)*2][2*(jf/2)]!=-1)
                    {
                      fw_refframe = dec_picture->ref_idx[LIST_0][(ifx/2)*2][2*(jf/2)];
                      fw_ref_idx = fw_refframe;
                    }
                    if (dec_picture->ref_idx[LIST_1][(ifx/2)*2][2*(jf/2)]!=-1)
                    {
                      bw_refframe = dec_picture->ref_idx[LIST_1][(ifx/2)*2][2*(jf/2)];
                      bw_ref_idx = bw_refframe;
                    }
                    
                    if      (dec_picture->ref_idx[LIST_1][(ifx/2)*2][2*(jf/2)]==-1) direct_pdir = 0;
                    else if (dec_picture->ref_idx[LIST_0][(ifx/2)*2][2*(jf/2)]==-1) direct_pdir = 1;
                    else                                                            direct_pdir = 2;
                  }                
                  else
                  {
                    
                    fw_refframe = dec_picture->ref_idx[LIST_0][ifx][jf];
                    bw_refframe = dec_picture->ref_idx[LIST_1][ifx][jf];
                    
                    fw_ref_idx = fw_refframe;
                    bw_ref_idx = bw_refframe;
                  }
                }
                
                
                if (mv_mode==0 && img->direct_type )
                {
                  if (direct_pdir == 0 || direct_pdir == 2)
                  {
                    i1=(img->pix_c_x+ii+ioff)*f1+fw_mv_array[ifx][jf][0];                  
                    
                    if (!curr_mb_field)
                    {
                      j1=(img->pix_c_y+jj+joff)*f1+fw_mv_array[ifx][jf][1];
                    }
                    else
                    {
                      if (mb_nr%2 == 0) 
                        j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+fw_mv_array[ifx][jf][1];
                      else
                        j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +fw_mv_array[ifx][jf][1];
                    }
                
                    j1 += listX[0+list_offset][fw_refframe]->chroma_vector_adjustment;
                    
                    ii0=max (0, min (i1>>3, img->width_cr-1));
                    jj0=max (0, min (j1>>3, max_y_cr));
                    ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                    jj1=max (0, min ((j1>>3)+1, max_y_cr));
                    
                    if1=(i1 & f2);
                    jf1=(j1 & f2);
                    if0=f1-if1;
                    jf0=f1-jf1;
                    
                    fw_pred=(if0*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii0]+
                             if1*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii1]+
                             if0*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii0]+
                             if1*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  }
                  if (direct_pdir == 1 || direct_pdir == 2)
                  {
                    i1=(img->pix_c_x+ii+ioff)*f1+bw_mv_array[ifx][jf][0];
                    
                    if (!curr_mb_field)
                    {
                      j1=(img->pix_c_y+jj+joff)*f1+bw_mv_array[ifx][jf][1];
                    }
                    else
                    {
                      if (mb_nr%2 == 0) 
                        j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+bw_mv_array[ifx][jf][1];
                      else
                        j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +bw_mv_array[ifx][jf][1];
                    }
                    j1 += listX[1+list_offset][bw_refframe]->chroma_vector_adjustment;
                    
                    ii0=max (0, min (i1>>3, img->width_cr-1));
                    jj0=max (0, min (j1>>3, max_y_cr));
                    ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                    jj1=max (0, min ((j1>>3)+1, max_y_cr));
                    
                    if1=(i1 & f2);
                    jf1=(j1 & f2);
                    if0=f1-if1;
                    jf0=f1-jf1;
                
                    bw_pred=(if0*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii0]+
                             if1*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii1]+
                             if0*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii0]+
                             if1*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  }

                }
                else
                {
                  i1=(img->pix_c_x+ii+ioff)*f1+fw_mv_array[ifx][jf][0];

                  if (!curr_mb_field)
                  {
                    j1=(img->pix_c_y+jj+joff)*f1+fw_mv_array[ifx][jf][1];
                  }
                  else
                  {
                    if (mb_nr%2 == 0) 
                      j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+fw_mv_array[ifx][jf][1];
                    else
                      j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +fw_mv_array[ifx][jf][1];
                  }

                  j1 += listX[0+list_offset][fw_refframe]->chroma_vector_adjustment;
                  
                  ii0=max (0, min (i1>>3, img->width_cr-1));
                  jj0=max (0, min (j1>>3, max_y_cr));
                  ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                  jj1=max (0, min ((j1>>3)+1, max_y_cr));
                  
                  if1=(i1 & f2);
                  jf1=(j1 & f2);
                  if0=f1-if1;
                  jf0=f1-jf1;
                  
                  fw_pred=(if0*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii0]+
                           if1*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii1]+
                           if0*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii0]+
                           if1*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  
                  i1=(img->pix_c_x+ii+ioff)*f1+bw_mv_array[ifx][jf][0];
                  
                  if (!curr_mb_field)
                  {
                    j1=(img->pix_c_y+jj+joff)*f1+bw_mv_array[ifx][jf][1];
                  }
                  else
                  {
                    if (mb_nr%2 == 0) 
                      j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+bw_mv_array[ifx][jf][1];
                    else
                      j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +bw_mv_array[ifx][jf][1];
                  }

                  j1 += listX[1+list_offset][bw_refframe]->chroma_vector_adjustment;

                  ii0=max (0, min (i1>>3, img->width_cr-1));
                  jj0=max (0, min (j1>>3, max_y_cr));
                  ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                  jj1=max (0, min ((j1>>3)+1, max_y_cr));
                  
                  if1=(i1 & f2);
                  jf1=(j1 & f2);
                  if0=f1-if1;
                  jf0=f1-jf1;
                  
                  bw_pred=(if0*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii0]+
                           if1*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii1]+
                           if0*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii0]+
                           if1*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  
                }
                
                if (img->apply_weights)
                {
                  if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                    (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
                  {
                    fw_ref_idx >>=1;
                    bw_ref_idx >>=1;
                  }

                  if (img->direct_type && direct_pdir==1)
                  {
                    img->mpr[ii+ioff][jj+joff]= Clip1(((img->wp_weight[1][bw_ref_idx][uv+1] * bw_pred  + img->wp_round_chroma)>>img->chroma_log2_weight_denom) + img->wp_offset[1][bw_refframe>>curr_mb_field][uv+1]);   // Replaced with integer only operations
                  }
                  else if (img->direct_type && direct_pdir==0)
                  {
                    img->mpr[ii+ioff][jj+joff]=Clip1(((img->wp_weight[0][fw_ref_idx][uv+1] * fw_pred + img->wp_round_chroma)>>img->chroma_log2_weight_denom) + img->wp_offset[0][fw_refframe>>curr_mb_field][uv+1]);   // Replaced with integer only operations
                  }
                  else
                  {
                    
                    int wt_list_offset = (active_pps->weighted_bipred_idc==2)?list_offset:0;

                    int alpha_fw = img->wbp_weight[0+wt_list_offset][fw_ref_idx][bw_ref_idx][uv+1];
                    int alpha_bw = img->wbp_weight[1+wt_list_offset][fw_ref_idx][bw_ref_idx][uv+1];

                    img->mpr[ii+ioff][jj+joff]= Clip1(((alpha_fw * fw_pred + alpha_bw * bw_pred  + (1<<img->chroma_log2_weight_denom)) >> (img->chroma_log2_weight_denom + 1))+ ((img->wp_offset[wt_list_offset + 0][fw_ref_idx][uv+1] + img->wp_offset[wt_list_offset + 1][bw_ref_idx][uv+1] + 1)>>1) );
                  }
                }
                else
                {
                  if (img->direct_type && direct_pdir==1)
                  {
                    img->mpr[ii+ioff][jj+joff]=bw_pred;
                  }
                  else if (img->direct_type && direct_pdir==0)
                  {
                    img->mpr[ii+ioff][jj+joff]=fw_pred;
                  } 
                  else
                  {
                    img->mpr[ii+ioff][jj+joff]=(fw_pred + bw_pred + 1 )/2;
                  }
                }
              }
            }
          }
        }
        
        if (!smb)
        {
          itrans(img,ioff,joff,2*uv+i,j);
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              dec_picture->imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
    
    if(smb)
    {
      itrans_sp_chroma(img,2*uv);
      for (j=4;j<6;j++)
      {
        joff=(j-4)*4;
        j4=img->pix_c_y+joff;
        for(i=0;i<2;i++)
        {
          ioff=i*4;
          i4=img->pix_c_x+ioff;
          itrans(img,ioff,joff,2*uv+i,j);
          
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              dec_picture->imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
  }

  return 0;
}




/////////////////////////////////////<
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////