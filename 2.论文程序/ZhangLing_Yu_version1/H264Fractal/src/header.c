
/*!
 *************************************************************************************
 * \file header.c
 *
 * \brief
 *    H.264 Slice and Sequence headers
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *************************************************************************************
 */

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "elements.h"
// #include "header.h"
#include "rtp_.h"
#include "mbuffer.h"
#include "i_defines.h"
#include "vlc.h"
#include "parset.h"
#include "header.h"

// A little trick to avoid those horrible #if TRACE all over the source code
#if TRACE
#define SYMTRACESTRING(s) strncpy(sym->tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // do nothing
#endif

int * assignSE2partition[2] ;
int assignSE2partition_NoDP[SE_MAX_ELEMENTS] =
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
int assignSE2partition_DP[SE_MAX_ELEMENTS] =
  {  0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0 } ;



static int ref_pic_list_reordering();
static int dec_ref_pic_marking();
static int pred_weight_table();
/////////////////////////////////////>


extern StorablePicture *dec_picture;


extern unsigned CeilLog2(unsigned uiVal);
void dec_ref_pic_marking_dec(Bitstream *currStream);
/////////////////////////////////////<


/*!
 ********************************************************************************************
 * \brief 
 *    Write a slice header
 ********************************************************************************************
*/
int SliceHeader()
{
  int dP_nr = assignSE2partition_dec[input->partition_mode][SE_HEADER];
  DataPartition *partition = & ((img->currentSlice)->partArr[dP_nr]);
  Slice* currSlice = img->currentSlice;
  int len = 0;
  unsigned int field_pic_flag = 0, bottom_field_flag = 0;    // POC200301

  int num_bits_slice_group_change_cycle;
  float numtmp;	
	
  if (img->MbaffFrameFlag)
    len  = ue_v("SH: first_mb_in_slice", img->current_mb_nr >> 1,   partition);
  else
    len  = ue_v("SH: first_mb_in_slice", img->current_mb_nr,   partition);//return 1

  len += ue_v("SH: slice_type",        get_picture_type (),   partition);//return 7

  // Note: Encoder supports only one pic/seq parameter set, hence value is
  // hard coded to zero
//  len += ue_v("SH: pic_parameter_set_id" , 0 ,partition);
  len += ue_v("SH: pic_parameter_set_id" , active_pps->pic_parameter_set_id ,partition);//return 1

  // frame_num
//  if(input->no_frames_h264 >= 1<<(log2_max_frame_num_minus4+4))
//    error ("Too many frames.  Increase log2_max_frame_num_minus4",-999);  

  len += u_v (log2_max_frame_num_minus4 + 4,"SH: frame_num", img->frame_num, partition);//return 4

  if (!active_sps->frame_mbs_only_flag)//0
  {
    // field_pic_flag    u(1)
    field_pic_flag = (img->structure ==TOP_FIELD || img->structure ==BOTTOM_FIELD)?1:0;
    assert( field_pic_flag == img->fld_flag );
    len += u_1("SH: field_pic_flag", field_pic_flag, partition);

    if (field_pic_flag)
    {
      //bottom_field_flag     u(1)
      bottom_field_flag = (img->structure == BOTTOM_FIELD)?1:0;
      len += u_1("SH: bottom_field_flag" , bottom_field_flag ,partition);
    }
  }

  if (img->currentPicture->idr_flag)//1
  {
    // idr_pic_id
    // hard coded to zero because we don't have proper IDR handling at the moment
    len += ue_v ("SH: idr_pic_id", 0, partition);//return 1
  }

  // POC200301
  if (img->pic_order_cnt_type == 0)//1
  {
    if (active_sps->frame_mbs_only_flag)//1
    {
      img->pic_order_cnt_lsb = (img->toppoc & ~((((unsigned int)(-1)) << (log2_max_pic_order_cnt_lsb_minus4+4))) );//0
    }
    else
    {
      if (!field_pic_flag || img->structure == TOP_FIELD)
        img->pic_order_cnt_lsb = (img->toppoc & ~((((unsigned int)(-1)) << (log2_max_pic_order_cnt_lsb_minus4+4))) );
      else if ( img->structure == BOTTOM_FIELD )
        img->pic_order_cnt_lsb = (img->bottompoc & ~((((unsigned int)(-1)) << (log2_max_pic_order_cnt_lsb_minus4+4))) );
    }

    len += u_v (log2_max_pic_order_cnt_lsb_minus4+4, "SH: pic_order_cnt_lsb", img->pic_order_cnt_lsb, partition);//return 4

    if (img->pic_order_present_flag && !field_pic_flag)  // img->fld_flag//0
    {
      len += se_v ("SH: delta_pic_order_cnt_bottom", img->delta_pic_order_cnt_bottom, partition);
    }
  }
  if (img->pic_order_cnt_type == 1 && !img->delta_pic_order_always_zero_flag)//0
  {
    len += se_v ("SH: delta_pic_order_cnt[0]", img->delta_pic_order_cnt[0], partition);

    if (img->pic_order_present_flag && !field_pic_flag)  // img->fld_flag
    {
      len += se_v ("SH: delta_pic_order_cnt[1]", img->delta_pic_order_cnt[1], partition);
    }
  }

  // redundant slice info redundant_pic_cnt is missing here
  if (input->redundant_slice_flag)//0
  {
    len += ue_v ("SH: redundant_pic_cnt", img->redundant_pic_cnt, partition);
  }

  // Direct Mode Type selection for B pictures
  if (img->type==B_SLICE)//0
  {
    len +=  u_1 ("SH: direct_spatial_mv_pred_flag", img->direct_type, partition);  	
    //len +=  u_1 ("SH: direct_spatial_mv_pred_flag", input->direct_type, partition);
  }

  if ((img->type == P_SLICE) || (img->type == B_SLICE) || (img->type==SP_SLICE))//0
  {
    int override_flag;
    if ((img->type == P_SLICE) || (img->type==SP_SLICE))
    {
      override_flag = (img->num_ref_idx_l0_active != (active_pps->num_ref_idx_l0_active_minus1 +1)) ? 1 : 0;
    }
    else
    {
      override_flag = ((img->num_ref_idx_l0_active != (active_pps->num_ref_idx_l0_active_minus1 +1)) 
                      || (img->num_ref_idx_l1_active != (active_pps->num_ref_idx_l1_active_minus1 +1))) ? 1 : 0;
    }
    // num_ref_idx_active_override_flag here always 1
    len +=  u_1 ("SH: num_ref_idx_active_override_flag", override_flag, partition);
    
    if (override_flag) 
    {
      len += ue_v ("SH: num_ref_idx_l0_active_minus1", img->num_ref_idx_l0_active-1, partition);
      if (img->type==B_SLICE)
      {
        len += ue_v ("SH: num_ref_idx_l1_active_minus1", img->num_ref_idx_l1_active-1, partition);
      }
    }

  }
  len += ref_pic_list_reordering();//18+0

  //if (((img->type == P_SLICE || img->type == SP_SLICE) && input->WeightedPrediction) || 
  //   ((img->type == B_SLICE) && input->WeightedBiprediction == 1))
  if (((img->type == P_SLICE || img->type == SP_SLICE) && active_pps->weighted_pred_flag) || 
     ((img->type == B_SLICE) && active_pps->weighted_bipred_idc == 1))  //0
  {
    len += pred_weight_table();
  }

  if (img->nal_reference_idc)//0
    len += dec_ref_pic_marking();//18+2

  if(input->symbol_mode==CABAC && img->type!=I_SLICE /*&& img->type!=SI_IMG*/)//0
  {
    len += ue_v("SH: cabac_init_idc", img->model_number, partition);
  }

  // we transmit zero in the pps, so here the real QP
  //len += se_v("SH: slice_qp_delta", (img->qp - 26), partition);
  len += se_v("SH: slice_qp_delta", (currSlice->qp - 26 - active_pps->pic_init_qp_minus26), partition);  //20+5

  if (img->type==SP_SLICE /*|| img->type==SI_SLICE*/)//0
  {
    if (img->type==SP_SLICE) // Switch Flag only for SP pictures
    {
      len += u_1 ("SH: sp_for_switch_flag", 0, partition);   // 1 for switching SP, 0 for normal SP
    }
    len += se_v ("SH: slice_qs_delta", (img->qpsp - 26), partition );
  }
/*
  if (input->LFSendParameters)
  {
    len += ue_v("SH: disable_deblocking_filter_idc",input->LFDisableIdc, partition);  // Turn loop filter on/off on slice basis 

    if (input->LFDisableIdc!=1)
    {
      len += se_v ("SH: slice_alpha_c0_offset_div2", input->LFAlphaC0Offset / 2, partition);

      len += se_v ("SH: slice_beta_offset_div2", input->LFBetaOffset / 2, partition);
    }
  }
*/
  if (active_pps->deblocking_filter_control_present_flag)//0
  {
    len += ue_v("SH: disable_deblocking_filter_idc",img->LFDisableIdc, partition);  // Turn loop filter on/off on slice basis 

    if (img->LFDisableIdc!=1)
    {
      len += se_v ("SH: slice_alpha_c0_offset_div2", img->LFAlphaC0Offset / 2, partition);

      len += se_v ("SH: slice_beta_offset_div2", img->LFBetaOffset / 2, partition);
    }
  }

	
  if ( active_pps->num_slice_groups_minus1>0 &&
    active_pps->slice_group_map_type>=3 && active_pps->slice_group_map_type<=5)//0
  {
    numtmp=img->PicHeightInMapUnits*img->PicWidthInMbs/(float)(active_pps->slice_group_change_rate_minus1+1)+1;
    num_bits_slice_group_change_cycle = (int)ceil(log(numtmp)/log(2));
    
    //! img->slice_group_change_cycle can be changed before calling FmoInit()
    len += u_v (num_bits_slice_group_change_cycle, "SH: slice_group_change_cycle", img->slice_group_change_cycle, partition);
  }
  
  // NOTE: The following syntax element is actually part 
  //        Slice data partition A RBSP syntax

  if(input->partition_mode&&!img->currentPicture->idr_flag)
  {
    len += ue_v("DPA: slice_id", img->current_slice_nr, partition);
  }

  return len;//25
}

/*!
 ********************************************************************************************
 * \brief 
 *    writes the ref_pic_list_reordering syntax
 *    based on content of according fields in img structure
 *
 * \return
 *    number of bits used 
 ********************************************************************************************
*/
static int ref_pic_list_reordering()
{
  int dP_nr = assignSE2partition_dec[input->partition_mode][SE_HEADER];
  DataPartition *partition = &((img->currentSlice)->partArr[dP_nr]);
  Slice *currSlice = img->currentSlice;

  int i, len=0;

  if ((img->type!=I_SLICE) /*&&(img->type!=SI_IMG)*/ )
  {
    len += u_1 ("SH: ref_pic_list_reordering_flag_l0", currSlice->ref_pic_list_reordering_flag_l0, partition);
    if (currSlice->ref_pic_list_reordering_flag_l0)
    {
      i=-1;
      do
      {
        i++;
        len += ue_v ("SH: remapping_of_pic_num_idc", currSlice->remapping_of_pic_nums_idc_l0[i], partition);
        if (currSlice->remapping_of_pic_nums_idc_l0[i]==0 ||
            currSlice->remapping_of_pic_nums_idc_l0[i]==1)
        {
          len += ue_v ("SH: abs_diff_pic_num_minus1_l0", currSlice->abs_diff_pic_num_minus1_l0[i], partition);
        }
        else
        {
          if (currSlice->remapping_of_pic_nums_idc_l0[i]==2)
          {
            len += ue_v ("SH: long_term_pic_idx_l0", currSlice->long_term_pic_idx_l0[i], partition);
          }
        }

      } while (currSlice->remapping_of_pic_nums_idc_l0[i] != 3);
    }
  }

  if (img->type==B_SLICE)
  {
    len += u_1 ("SH: ref_pic_list_reordering_flag_l1", currSlice->ref_pic_list_reordering_flag_l1, partition);
    if (currSlice->ref_pic_list_reordering_flag_l1)
    {
      i=-1;
      do
      {
        i++;
        len += ue_v ("SH: remapping_of_pic_num_idc", currSlice->remapping_of_pic_nums_idc_l1[i], partition);
        if (currSlice->remapping_of_pic_nums_idc_l1[i]==0 ||
            currSlice->remapping_of_pic_nums_idc_l1[i]==1)
        {
          len += ue_v ("SH: abs_diff_pic_num_minus1_l1", currSlice->abs_diff_pic_num_minus1_l1[i], partition);
        }
        else
        {
          if (currSlice->remapping_of_pic_nums_idc_l1[i]==2)
          {
            len += ue_v ("SH: long_term_pic_idx_l1", currSlice->long_term_pic_idx_l1[i], partition);
          }
        }
      } while (currSlice->remapping_of_pic_nums_idc_l1[i] != 3);
    }
  }

  return len;
}


/*!
 ************************************************************************
 * \brief
 *    write the memory menagement control operations
 ************************************************************************
 */
static int dec_ref_pic_marking()
{
  int dP_nr = assignSE2partition_dec[input->partition_mode][SE_HEADER];
  DataPartition *partition = &((img->currentSlice)->partArr[dP_nr]);

  DecRefPicMarking_t *tmp_drpm;

  int val, len=0;

  if (img->currentPicture->idr_flag)
  {
    len += u_1("SH: no_output_of_prior_pics_flag", img->no_output_of_prior_pics_flag, partition);
    len += u_1("SH: long_term_reference_flag", img->long_term_reference_flag, partition);
  }
  else
  {
    img->adaptive_ref_pic_buffering_flag = (img->dec_ref_pic_marking_buffer!=NULL);

    len += u_1("SH: adaptive_ref_pic_buffering_flag", img->adaptive_ref_pic_buffering_flag, partition);

    if (img->adaptive_ref_pic_buffering_flag)
    {
      tmp_drpm = img->dec_ref_pic_marking_buffer;
      // write Memory Management Control Operation 
      do
      {
        if (tmp_drpm==NULL) error ("Error encoding MMCO commands", 500);
        
        val = tmp_drpm->memory_management_control_operation;
        len += ue_v("SH: memory_management_control_operation", val, partition);

        if ((val==1)||(val==3)) 
        {
          len += 1 + ue_v("SH: difference_of_pic_nums_minus1", tmp_drpm->difference_of_pic_nums_minus1, partition);
        }
        if (val==2)
        {
          len+= ue_v("SH: long_term_pic_num", tmp_drpm->long_term_pic_num, partition);
        }
        if ((val==3)||(val==6))
        {
          len+= ue_v("SH: long_term_frame_idx", tmp_drpm->long_term_frame_idx, partition);
        }
        if (val==4)
        {
          len += ue_v("SH: max_long_term_pic_idx_plus1", tmp_drpm->max_long_term_frame_idx_plus1, partition);
        }
        
        tmp_drpm=tmp_drpm->Next;
        
      } while (val != 0);
      
    }
  }
  return len;
}

/*!
 ************************************************************************
 * \brief
 *    write the memory menagement control operations
 ************************************************************************
 */
static int pred_weight_table()
{
  int dP_nr = assignSE2partition_dec[input->partition_mode][SE_HEADER];
  DataPartition *partition = &((img->currentSlice)->partArr[dP_nr]);

  int len = 0;
  int i,j;

  len += ue_v("SH: luma_log_weight_denom", luma_log_weight_denom, partition);
  
  len += ue_v("SH: chroma_log_weight_denom", chroma_log_weight_denom, partition);

  for (i=0; i< img->num_ref_idx_l0_active; i++)
  {
    if ( (wp_weight[0][i][0] != 1<<luma_log_weight_denom) || (wp_offset[0][i][0] != 0) )
    {
      len += u_1 ("SH: luma_weight_flag_l0", 1, partition);
      
      len += se_v ("SH: luma_weight_l0", wp_weight[0][i][0], partition);
        
      len += se_v ("SH: luma_offset_l0", wp_offset[0][i][0], partition);
    }
    else
    {
        len += u_1 ("SH: luma_weight_flag_l0", 0, partition);
    }

    if ( (wp_weight[0][i][1] != 1<<chroma_log_weight_denom) || (wp_offset[0][i][1] != 0) || 
     (wp_weight[0][i][2] != 1<<chroma_log_weight_denom) || (wp_offset[0][i][2] != 0)  )
    {
      len += u_1 ("chroma_weight_flag_l0", 1, partition);
      for (j=1; j<3; j++)
      {
        len += se_v ("chroma_weight_l0", wp_weight[0][i][j] ,partition);
      
        len += se_v ("chroma_offset_l0", wp_offset[0][i][j] ,partition);
      }
    }
    else
    {
      len += u_1 ("chroma_weight_flag_l0", 0, partition);
    }
  }

  if (img->type == B_SLICE)
  {
    for (i=0; i< img->num_ref_idx_l1_active; i++)
    {
      if ( (wp_weight[1][i][0] != 1<<luma_log_weight_denom) || (wp_offset[1][i][0] != 0) )
      {
        len += u_1 ("SH: luma_weight_flag_l1", 1, partition);
        
        len += se_v ("SH: luma_weight_l1", wp_weight[1][i][0], partition);
        
        len += se_v ("SH: luma_offset_l1", wp_offset[1][i][0], partition);
      }
      else
      {
        len += u_1 ("SH: luma_weight_flag_l1", 0, partition);
      }
      
      if ( (wp_weight[1][i][1] != 1<<chroma_log_weight_denom) || (wp_offset[1][i][1] != 0) || 
      (wp_weight[1][i][2] != 1<<chroma_log_weight_denom) || (wp_offset[1][i][2] != 0) )
      {
        len += u_1 ("chroma_weight_flag_l1", 1, partition);
        for (j=1; j<3; j++)
        {
          len += se_v ("chroma_weight_l1", wp_weight[1][i][j] ,partition);
          
          len += se_v ("chroma_offset_l1", wp_offset[1][i][j] ,partition);
        }
      }
      else
      {
        len += u_1 ("chroma_weight_flag_l1", 0, partition);
      }
    }
  }
  return len;
}
  
/********************************************************************************************
 ********************************************************************************************
 *
 * Local Support Functions
 *
 ********************************************************************************************
 ********************************************************************************************/



// StW Note: This function is a hack.  It would be cleaner if the encoder maintains
// the picture type in the given format.  Note further that I have yet to understand
// why the encoder needs to know whether a picture is predicted from one or more
// reference pictures.

/*!
 ************************************************************************
 * \brief
 *    Selects picture type and codes it to symbol
 ************************************************************************
 */
int get_picture_type()
{

  // set this value to zero for transmission without signaling 
  // that the whole picture has the same slice type
  int same_slicetype_for_whole_frame = 5;

  switch (img->type)
  {
  case I_SLICE:
    return 2 + same_slicetype_for_whole_frame;
    break;
  case P_SLICE:
    return 0 + same_slicetype_for_whole_frame;
    break;
  case B_SLICE:
    return 1 + same_slicetype_for_whole_frame;
    break;
  case SP_SLICE:
    return 3 + same_slicetype_for_whole_frame;
    break;
  default:
    error("Picture Type not supported!",1);
    break;
  }
   
  return 0;
}



/*!
 *****************************************************************************
 *
 * \brief 
 *    int Partition_BC_Header () write the Partition type B, C header
 *
 * \return
 *    Number of bits used by the partition header
 *
 * \par Parameters
 *    PartNo: Partition Number to which the header should be written
 *
 * \par Side effects
 *    Partition header as per VCEG-N72r2 is written into the appropriate 
 *    partition bit buffer
 *
 * \par Limitations/Shortcomings/Tweaks
 *    The current code does not support the change of picture parameters within
 *    one coded sequence, hence there is only one parameter set necessary.  This
 *    is hard coded to zero.
 *
 * \date
 *    October 24, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/


int Partition_BC_Header(int PartNo)
{
  DataPartition *partition = &((img->currentSlice)->partArr[PartNo]);
  SyntaxElement symbol, *sym = &symbol;

  int len = 0;

  assert (input->of_mode == PAR_OF_RTP);
  assert (PartNo > 0 && PartNo < img->currentSlice->max_part_nr);

  sym->type = SE_HEADER;         // This will be true for all symbols generated here
  sym->mapping = ue_linfo;       // Mapping rule: Simple code number to len/info
  sym->value2  = 0;

	//ZL 
	//changed according to the g050r1
	SYMTRACESTRING("RTP-PH: Slice ID");
  sym->value1 = img->current_slice_nr;
  len += writeSyntaxElement_UVLC (sym, partition);

	if(active_pps->redundant_pic_cnt_present_flag)
  {
  SYMTRACESTRING("RTP-PH: Picture ID");
  sym->value1 = img->currentSlice->picture_id;
  len += writeSyntaxElement_UVLC (sym, partition);
  }


  return len;
}
/////////////////////////////////////>

/////////////////////////////////////>


extern int UsedBits;
static void ref_pic_list_reordering_dec();

static void pred_weight_table_dec();

/////////////////////////////////////<

// /*!
//  ************************************************************************
//  * \brief
//  *    calculate Ceil(Log2(uiVal))
//  ************************************************************************
//  */
// unsigned CeilLog2( unsigned uiVal)
// {
//   unsigned uiTmp = uiVal-1;
//   unsigned uiRet = 0;
// 
//   while( uiTmp != 0 )
//   {
//     uiTmp >>= 1;
//     uiRet++;
//   }
//   return uiRet;
// }


/*!
 ************************************************************************
 * \brief
 *    read the first part of the header (only the pic_parameter_set_id)
 * \return
 *    Length of the first part of the slice header (in bits)
 ************************************************************************
 */
int FirstPartOfSliceHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition_dec[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;
  int tmp;

  UsedBits= partition->bitstream->frame_bitoffset; // was hardcoded to 31 for previous start-code. This is better.

  // Get first_mb_in_slice
  currSlice->start_mb_nr = ue_v_dec ("SH: first_mb_in_slice", currStream);

  tmp = ue_v_dec ("SH: slice_type", currStream);
  
  if (tmp>4) tmp -=5;

  img->type = currSlice->picture_type = (SliceType) tmp;

  currSlice->pic_parameter_set_id = ue_v_dec ("SH: pic_parameter_set_id", currStream);
  
  return UsedBits;
}

/*!
 ************************************************************************
 * \brief
 *    read the scond part of the header (without the pic_parameter_set_id 
 * \return
 *    Length of the second part of the Slice header in bits
 ************************************************************************
 */
int RestOfSliceHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition_dec[currSlice->dp_mode][SE_HEADER];//0
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;

  int val, len;

  img->frame_num = u_v_dec (active_sps->log2_max_frame_num_minus4 + 4, "SH: frame_num", currStream);//0,		frame_bitoffset	9 to 13

  /* Tian Dong: frame_num gap processing, if found */
  if (img->idr_flag)//1
  {
    img->pre_frame_num = img->frame_num;//0
    assert(img->frame_num == 0);
  }

  if (active_sps->frame_mbs_only_flag)//1
  {
    img->structure = FRAME;//0
    img->field_pic_flag=0;
  }
  else
  {
    // field_pic_flag   u(1)
    img->field_pic_flag = u_1_dec("SH: field_pic_flag", currStream);
    if (img->field_pic_flag)
    {
      // bottom_field_flag  u(1)
      img->bottom_field_flag = u_1_dec("SH: bottom_field_flag", currStream);

      img->structure = img->bottom_field_flag ? BOTTOM_FIELD : TOP_FIELD;
    }
    else
    {
      img->structure = FRAME;
      img->bottom_field_flag=0;
    }
  }

  currSlice->structure = img->structure;//0

  img->MbaffFrameFlag=(active_sps->mb_adaptive_frame_field_flag && (img->field_pic_flag==0));//0

  if (img->structure == FRAME       ) assert (img->field_pic_flag == 0);
  if (img->structure == TOP_FIELD   ) assert (img->field_pic_flag == 1 && img->bottom_field_flag == 0);
  if (img->structure == BOTTOM_FIELD) assert (img->field_pic_flag == 1 && img->bottom_field_flag == 1);

  if (img->idr_flag)//1
  {
	  img->idr_pic_id = ue_v_dec("SH: idr_pic_id", currStream);//0	frame_bitoffset	13 to 14

  }

  if (active_sps->pic_order_cnt_type == 0)//0==0
  {
//	frame_bitoffset	14 to 18 
    img->pic_order_cnt_lsb = u_v_dec(active_sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "SH: pic_order_cnt_lsb", currStream);//return 0;	frame_bitoffset	18 改变了frame_bitoffset
    if( active_pps->pic_order_present_flag  ==  1 &&  !img->field_pic_flag )
      img->delta_pic_order_cnt_bottom = se_v_dec("SH: delta_pic_order_cnt_bottom", currStream);
    else
      img->delta_pic_order_cnt_bottom = 0;  
  }
  if( active_sps->pic_order_cnt_type == 1 && !active_sps->delta_pic_order_always_zero_flag ) 
  {
    img->delta_pic_order_cnt[ 0 ] = se_v_dec("SH: delta_pic_order_cnt[0]", currStream);
    if( active_pps->pic_order_present_flag  ==  1  &&  !img->field_pic_flag )
      img->delta_pic_order_cnt[ 1 ] = se_v_dec("SH: delta_pic_order_cnt[1]", currStream);
  }else//1
  {
    if (active_sps->pic_order_cnt_type == 1)
    {
      img->delta_pic_order_cnt[ 0 ] = 0;
      img->delta_pic_order_cnt[ 1 ] = 0;
    }
  }
  
  //! redundant_pic_cnt is missing here
  if (active_pps->redundant_pic_cnt_present_flag)//0
  {
    img->redundant_pic_cnt = ue_v_dec ("SH: redundant_pic_cnt", currStream);
  }

  if(img->type==B_SLICE)
  {
    img->direct_type = u_1_dec ("SH: direct_spatial_mv_pred_flag", currStream);
  }

  img->num_ref_idx_l0_active = active_pps->num_ref_idx_l0_active_minus1 + 1;//5=4+1
  img->num_ref_idx_l1_active = active_pps->num_ref_idx_l1_active_minus1 + 1;//5=4+1

  if(img->type==P_SLICE || img->type == SP_SLICE || img->type==B_SLICE)
  {
    val = u_1_dec ("SH: num_ref_idx_override_flag", currStream);
    if (val)
    {
      img->num_ref_idx_l0_active = 1 + ue_v_dec ("SH: num_ref_idx_l0_active_minus1", currStream);
      
      if(img->type==B_SLICE)
      {
        img->num_ref_idx_l1_active = 1 + ue_v_dec ("SH: num_ref_idx_l1_active_minus1", currStream);
      }
    }
  }
  if (img->type!=B_SLICE)
  {
    img->num_ref_idx_l1_active = 0;
  }

  ref_pic_list_reordering_dec();

  img->apply_weights = ((active_pps->weighted_pred_flag && (currSlice->picture_type == P_SLICE || currSlice->picture_type == SP_SLICE) )
          || ((active_pps->weighted_bipred_idc > 0 ) && (currSlice->picture_type == B_SLICE)));//0

  if ((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
      (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE)))
  {
    pred_weight_table_dec();
  }


//  frame_bitoffset	18 to 20  暂时要着
  if (img->nal_reference_idc)//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZz
    dec_ref_pic_marking_dec(currStream);

  if (active_pps->entropy_coding_mode_flag && img->type!=I_SLICE && img->type!=SI_SLICE)
  {
    img->model_number = ue_v_dec("SH: cabac_init_idc", currStream);
  }
  else 
  {
    img->model_number = 0;
  }

  val = se_v_dec("SH: slice_qp_delta", currStream);
  currSlice->qp = img->qp = 26 + active_pps->pic_init_qp_minus26 + val;//改变了qp值！0+（-9）

  if(img->type==SP_SLICE || img->type == SI_SLICE) 
  {
    if(img->type==SP_SLICE)
    {
      img->sp_switch = u_1_dec ("SH: sp_for_switch_flag", currStream);
    }
    val = se_v_dec("SH: slice_qs_delta", currStream);
    img->qpsp = 26 + active_pps->pic_init_qs_minus26 + val;
  }

  if (active_pps->deblocking_filter_control_present_flag)
  {
    currSlice->LFDisableIdc = ue_v_dec ("SH: disable_deblocking_filter_idc", currStream);

    if (currSlice->LFDisableIdc!=1)
    {
      currSlice->LFAlphaC0Offset = 2 * se_v_dec("SH: slice_alpha_c0_offset_div2", currStream);
      currSlice->LFBetaOffset = 2 * se_v_dec("SH: slice_beta_offset_div2", currStream);
    }
    else
    {
      currSlice->LFAlphaC0Offset = currSlice->LFBetaOffset = 0;
    }
  }
  else 
  {
    currSlice->LFDisableIdc = currSlice->LFAlphaC0Offset = currSlice->LFBetaOffset = 0;
  }

  if (active_pps->num_slice_groups_minus1>0 && active_pps->slice_group_map_type>=3 &&
      active_pps->slice_group_map_type<=5)
  {
    len = (active_sps->pic_height_in_map_units_minus1+1)*(active_sps->pic_width_in_mbs_minus1+1)/ 
          (active_pps->slice_group_change_rate_minus1+1);
    if (((active_sps->pic_height_in_map_units_minus1+1)*(active_sps->pic_width_in_mbs_minus1+1))% 
          (active_pps->slice_group_change_rate_minus1+1))
          len +=1;

    len = CeilLog2(len+1);

    img->slice_group_change_cycle = u_v_dec (len, "SH: slice_group_change_cycle", currStream);
  }
  img->PicHeightInMbs = img->FrameHeightInMbs / ( 1 + img->field_pic_flag );
  img->PicSizeInMbs   = img->PicWidthInMbs * img->PicHeightInMbs;
  img->FrameSizeInMbs = img->PicWidthInMbs * img->FrameHeightInMbs;

  return UsedBits;
}


/*!
 ************************************************************************
 * \brief
 *    read the reference picture reordering information
 ************************************************************************
 */
static void ref_pic_list_reordering_dec()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition_dec[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;
  int i, val;

  alloc_ref_pic_list_reordering_buffer(currSlice);
  
  if (img->type!=I_SLICE && img->type!=SI_SLICE)
  {
    val = currSlice->ref_pic_list_reordering_flag_l0 = u_1_dec ("SH: ref_pic_list_reordering_flag_l0", currStream);
    
    if (val)
    {
      i=0;
      do
      {
        val = currSlice->remapping_of_pic_nums_idc_l0[i] = ue_v_dec("SH: remapping_of_pic_nums_idc_l0", currStream);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l0[i] = ue_v_dec("SH: abs_diff_pic_num_minus1_l0", currStream);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l0[i] = ue_v_dec("SH: long_term_pic_idx_l0", currStream);
          }
        }
        i++;
        // assert (i>img->num_ref_idx_l0_active);
      } while (val != 3);
    }
  }

  if (img->type==B_SLICE)
  {
    val = currSlice->ref_pic_list_reordering_flag_l1 = u_1_dec ("SH: ref_pic_list_reordering_flag_l1", currStream);
    
    if (val)
    {
      i=0;
      do
      {
        val = currSlice->remapping_of_pic_nums_idc_l1[i] = ue_v_dec("SH: remapping_of_pic_nums_idc_l1", currStream);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l1[i] = ue_v_dec("SH: abs_diff_pic_num_minus1_l1", currStream);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l1[i] = ue_v_dec("SH: long_term_pic_idx_l1", currStream);
          }
        }
        i++;
        // assert (i>img->num_ref_idx_l1_active);
      } while (val != 3);
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    read the weighted prediction tables
 ************************************************************************
 */
static void pred_weight_table_dec()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition_dec[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  Bitstream *currStream = partition->bitstream;
  int luma_weight_flag_l0, luma_weight_flag_l1, chroma_weight_flag_l0, chroma_weight_flag_l1;
  int i,j;

  img->luma_log2_weight_denom = ue_v_dec ("SH: luma_log2_weight_denom", currStream);
  img->wp_round_luma = img->luma_log2_weight_denom ? 1<<(img->luma_log2_weight_denom - 1): 0;
  
  img->chroma_log2_weight_denom = ue_v_dec ("SH: chroma_log2_weight_denom", currStream);
  img->wp_round_chroma = img->chroma_log2_weight_denom ? 1<<(img->chroma_log2_weight_denom - 1): 0;

  reset_wp_params(img);

  for (i=0; i<img->num_ref_idx_l0_active; i++)
  {
    luma_weight_flag_l0 = u_1_dec("SH: luma_weight_flag_l0", currStream);
    
    if (luma_weight_flag_l0)
    {
      img->wp_weight[0][i][0] = se_v_dec ("SH: luma_weight_l0", currStream);
      img->wp_offset[0][i][0] = se_v_dec ("SH: luma_offset_l0", currStream);
    }
    else
    {
      img->wp_weight[0][i][0] = 1<<img->luma_log2_weight_denom;
      img->wp_offset[0][i][0] = 0;
    }
    
    chroma_weight_flag_l0 = u_1_dec ("SH: chroma_weight_flag_l0", currStream);
    
    for (j=1; j<3; j++)
    {
      if (chroma_weight_flag_l0)
      {
        img->wp_weight[0][i][j] = se_v_dec("SH: chroma_weight_l0", currStream);
        img->wp_offset[0][i][j] = se_v_dec("SH: chroma_offset_l0", currStream);
      }
      else
      {
        img->wp_weight[0][i][j] = 1<<img->chroma_log2_weight_denom;
        img->wp_offset[0][i][j] = 0;
      }
    }
  }
  if ((img->type == B_SLICE) && active_pps->weighted_bipred_idc == 1)
  {
    for (i=0; i<img->num_ref_idx_l1_active; i++)
    {
      luma_weight_flag_l1 = u_1_dec("SH: luma_weight_flag_l1", currStream);
      
      if (luma_weight_flag_l1)
      {
        img->wp_weight[1][i][0] = se_v_dec ("SH: luma_weight_l1", currStream);
        img->wp_offset[1][i][0] = se_v_dec ("SH: luma_offset_l1", currStream);
      }
      else
      {
        img->wp_weight[1][i][0] = 1<<img->luma_log2_weight_denom;
        img->wp_offset[1][i][0] = 0;
      }
      
      chroma_weight_flag_l1 = u_1_dec ("SH: chroma_weight_flag_l1", currStream);
      
      for (j=1; j<3; j++)
      {
        if (chroma_weight_flag_l1)
        {
          img->wp_weight[1][i][j] = se_v_dec("SH: chroma_weight_l1", currStream);
          img->wp_offset[1][i][j] = se_v_dec("SH: chroma_offset_l1", currStream);
        }
        else
        {
          img->wp_weight[1][i][j] = 1<<img->chroma_log2_weight_denom;
          img->wp_offset[1][i][j] = 0;
        }
      }
    }
  }    
}


/*!
 ************************************************************************
 * \brief
 *    read the memory control operations
 ************************************************************************
 */
void dec_ref_pic_marking_dec(Bitstream *currStream)
{
  int val;

  DecRefPicMarking_t *tmp_drpm,*tmp_drpm2;

  // free old buffer content
  while (img->dec_ref_pic_marking_buffer)
  { 
    tmp_drpm=img->dec_ref_pic_marking_buffer;

    img->dec_ref_pic_marking_buffer=tmp_drpm->Next;
    free (tmp_drpm);
  } 

  if (img->idr_flag)
  {
    img->no_output_of_prior_pics_flag = u_1_dec("SH: no_output_of_prior_pics_flag", currStream);
    img->long_term_reference_flag = u_1_dec("SH: long_term_reference_flag", currStream);
  }
  else
  {
    img->adaptive_ref_pic_buffering_flag = u_1_dec("SH: adaptive_ref_pic_buffering_flag", currStream);
    if (img->adaptive_ref_pic_buffering_flag)
    {
      // read Memory Management Control Operation 
      do
      {
        tmp_drpm=(DecRefPicMarking_t*)calloc (1,sizeof (DecRefPicMarking_t));
        tmp_drpm->Next=NULL;
        
        val = tmp_drpm->memory_management_control_operation = ue_v_dec("SH: memory_management_control_operation", currStream);

        if ((val==1)||(val==3)) 
        {
          tmp_drpm->difference_of_pic_nums_minus1 = ue_v_dec("SH: difference_of_pic_nums_minus1", currStream);
        }
        if (val==2)
        {
          tmp_drpm->long_term_pic_num = ue_v_dec("SH: long_term_pic_num", currStream);
        }
          
        if ((val==3)||(val==6))
        {
          tmp_drpm->long_term_frame_idx = ue_v_dec("SH: long_term_frame_idx", currStream);
        }
        if (val==4)
        {
          tmp_drpm->max_long_term_frame_idx_plus1 = ue_v_dec("SH: max_long_term_pic_idx_plus1", currStream);
        }
        
        // add command
        if (img->dec_ref_pic_marking_buffer==NULL) 
        {
          img->dec_ref_pic_marking_buffer=tmp_drpm;
        }
        else
        {
          tmp_drpm2=img->dec_ref_pic_marking_buffer;
          while (tmp_drpm2->Next!=NULL) tmp_drpm2=tmp_drpm2->Next;
          tmp_drpm2->Next=tmp_drpm;
        }
        
      }while (val != 0);
      
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    To calculate the poc values
 *        based upon JVT-F100d2
 *  POC200301: Until Jan 2003, this function will calculate the correct POC
 *    values, but the management of POCs in buffered pictures may need more work.
 * \return
 *    none
 ************************************************************************
 */
void decode_poc(struct img_par *img)
{
  int i;
  // for POC mode 0:
  unsigned int        MaxPicOrderCntLsb = (1<<(active_sps->log2_max_pic_order_cnt_lsb_minus4+4));

  switch ( active_sps->pic_order_cnt_type )
  {
  case 0: // POC MODE 0
    // 1st
    if(img->idr_flag)
    {
      img->PrevPicOrderCntMsb = 0;
      img->PrevPicOrderCntLsb = 0;
    }
    else
    {
      if (img->last_has_mmco_5) 
      {
        if (img->last_pic_bottom_field)
        {
          img->PrevPicOrderCntMsb = 0;
          img->PrevPicOrderCntLsb = 0;
        }
        else
        {
          img->PrevPicOrderCntMsb = 0;
          img->PrevPicOrderCntLsb = img->toppoc;
        }
      }
    }
    // Calculate the MSBs of current picture
    if( img->pic_order_cnt_lsb  <  img->PrevPicOrderCntLsb  &&  
      ( img->PrevPicOrderCntLsb - img->pic_order_cnt_lsb )  >=  ( MaxPicOrderCntLsb / 2 ) )
      img->PicOrderCntMsb = img->PrevPicOrderCntMsb + MaxPicOrderCntLsb;
    else if ( img->pic_order_cnt_lsb  >  img->PrevPicOrderCntLsb  &&
      ( img->pic_order_cnt_lsb - img->PrevPicOrderCntLsb )  >  ( MaxPicOrderCntLsb / 2 ) )
      img->PicOrderCntMsb = img->PrevPicOrderCntMsb - MaxPicOrderCntLsb;
    else
      img->PicOrderCntMsb = img->PrevPicOrderCntMsb;
    
    // 2nd

    if ((!img->field_pic_flag)||(!img->bottom_field_flag ))
    {
      img->ThisPOC= img->toppoc = img->PicOrderCntMsb + img->pic_order_cnt_lsb;
    }

    if (!img->field_pic_flag)
    {
      img->bottompoc = img->toppoc + img->delta_pic_order_cnt_bottom;
    }
    else
    {
      if( img->bottom_field_flag ) 
      {
        img->ThisPOC= img->bottompoc = img->PicOrderCntMsb + img->pic_order_cnt_lsb;
      }
    }

    // last: some post-processing. 
      
//    if (img->newframe == 1)
    {
      if (!img->bottom_field_flag)
      {      
        img->framepoc = img->toppoc;
      }
      else
      {
        img->framepoc = img->bottompoc;
      }
    }

    if ( img->frame_num!=img->PreviousFrameNum)
      img->PreviousFrameNum=img->frame_num;

    if(!img->disposable_flag)
    {
      img->PrevPicOrderCntLsb = img->pic_order_cnt_lsb;
      img->PrevPicOrderCntMsb = img->PicOrderCntMsb;
    }

    break;

  case 1: // POC MODE 1
    // 1st
    if(img->idr_flag)
    {
      img->FrameNumOffset=0;     //  first pix of IDRGOP, 
      img->delta_pic_order_cnt[0]=0;                        //ignore first delta
      if(img->frame_num)  error("frame_num != 0 in idr pix", -1020);
    }
    else 
    {
      if (img->last_has_mmco_5)
      {
        img->PreviousFrameNumOffset = 0;
        img->PreviousFrameNum = 0;
      }
      if (img->frame_num<img->PreviousFrameNum)
      {             //not first pix of IDRGOP
        img->FrameNumOffset = img->PreviousFrameNumOffset + img->MaxFrameNum;
      }
      else 
      {
        img->FrameNumOffset = img->PreviousFrameNumOffset;
      }
    }

    // 2nd
    if(active_sps->num_ref_frames_in_pic_order_cnt_cycle) 
      img->AbsFrameNum = img->FrameNumOffset+img->frame_num;
    else 
      img->AbsFrameNum=0;
    if(img->disposable_flag && img->AbsFrameNum>0)
      img->AbsFrameNum--;

    // 3rd
    img->ExpectedDeltaPerPicOrderCntCycle=0;

    if(active_sps->num_ref_frames_in_pic_order_cnt_cycle)
      for(i=0;i<(int) active_sps->num_ref_frames_in_pic_order_cnt_cycle;i++)
        img->ExpectedDeltaPerPicOrderCntCycle += active_sps->offset_for_ref_frame[i];

    if(img->AbsFrameNum)
    {
      img->PicOrderCntCycleCnt = (img->AbsFrameNum-1)/active_sps->num_ref_frames_in_pic_order_cnt_cycle;
      img->FrameNumInPicOrderCntCycle = (img->AbsFrameNum-1)%active_sps->num_ref_frames_in_pic_order_cnt_cycle;
      img->ExpectedPicOrderCnt = img->PicOrderCntCycleCnt*img->ExpectedDeltaPerPicOrderCntCycle;
      for(i=0;i<=(int)img->FrameNumInPicOrderCntCycle;i++)
        img->ExpectedPicOrderCnt += active_sps->offset_for_ref_frame[i];
    }
    else 
      img->ExpectedPicOrderCnt=0;

    if(img->disposable_flag)
      img->ExpectedPicOrderCnt += active_sps->offset_for_non_ref_pic;

    if(img->field_pic_flag==0)
    {           //frame pix
      img->toppoc = img->ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
      img->bottompoc = img->toppoc + active_sps->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[1];
      img->ThisPOC = img->framepoc = (img->toppoc < img->bottompoc)? img->toppoc : img->bottompoc; // POC200301
    }
    else if (img->bottom_field_flag==0)
    {  //top field 
      img->ThisPOC = img->toppoc = img->ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
    } 
    else
    {  //bottom field
      img->ThisPOC = img->bottompoc = img->ExpectedPicOrderCnt + active_sps->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[0];
    }
    img->framepoc=img->ThisPOC;

    img->PreviousFrameNum=img->frame_num;
    img->PreviousFrameNumOffset=img->FrameNumOffset;
    
    break;


  case 2: // POC MODE 2
    if(img->idr_flag) // IDR picture
    {
      img->FrameNumOffset=0;     //  first pix of IDRGOP, 
      img->ThisPOC = img->framepoc = img->toppoc = img->bottompoc = 0;
      if(img->frame_num)  error("frame_num != 0 in idr pix", -1020);
    }
    else
    {
      if (img->last_has_mmco_5)
      {
        img->PreviousFrameNum = 0;
        img->PreviousFrameNumOffset = 0;
      }
      if (img->frame_num<img->PreviousFrameNum)
        img->FrameNumOffset = img->PreviousFrameNumOffset + img->MaxFrameNum;
      else 
        img->FrameNumOffset = img->PreviousFrameNumOffset;


      img->AbsFrameNum = img->FrameNumOffset+img->frame_num;
      if(img->disposable_flag)
        img->ThisPOC = (2*img->AbsFrameNum - 1);
      else
        img->ThisPOC = (2*img->AbsFrameNum);

      if (img->field_pic_flag==0)
        img->toppoc = img->bottompoc = img->framepoc = img->ThisPOC;
      else if (img->bottom_field_flag==0)
         img->toppoc = img->framepoc = img->ThisPOC;
      else img->bottompoc = img->framepoc = img->ThisPOC;
    }

    if (!img->disposable_flag)
      img->PreviousFrameNum=img->frame_num;
    img->PreviousFrameNumOffset=img->FrameNumOffset;
    break;


  default:
    //error must occurs
    assert( 1==0 );
    break;
  }
}

/*!
 ************************************************************************
 * \brief
 *    A little helper for the debugging of POC code
 * \return
 *    none
 ************************************************************************
 */
int dumppoc(struct img_par *img) {
    printf ("\nPOC locals...\n");
    printf ("toppoc                                %d\n", img->toppoc);
    printf ("bottompoc                             %d\n", img->bottompoc);
    printf ("frame_num                             %d\n", img->frame_num);
    printf ("field_pic_flag                        %d\n", img->field_pic_flag);
    printf ("bottom_field_flag                     %d\n", img->bottom_field_flag);
    printf ("POC SPS\n");
    printf ("log2_max_frame_num_minus4             %d\n", active_sps->log2_max_frame_num_minus4);         // POC200301
    printf ("log2_max_pic_order_cnt_lsb_minus4     %d\n", active_sps->log2_max_pic_order_cnt_lsb_minus4);
    printf ("pic_order_cnt_type                    %d\n", active_sps->pic_order_cnt_type);
    printf ("num_ref_frames_in_pic_order_cnt_cycle %d\n", active_sps->num_ref_frames_in_pic_order_cnt_cycle);
    printf ("delta_pic_order_always_zero_flag      %d\n", active_sps->delta_pic_order_always_zero_flag);
    printf ("offset_for_non_ref_pic                %d\n", active_sps->offset_for_non_ref_pic);
    printf ("offset_for_top_to_bottom_field        %d\n", active_sps->offset_for_top_to_bottom_field);
    printf ("offset_for_ref_frame[0]               %d\n", active_sps->offset_for_ref_frame[0]);
    printf ("offset_for_ref_frame[1]               %d\n", active_sps->offset_for_ref_frame[1]);
    printf ("POC in SLice Header\n");
    printf ("pic_order_present_flag                %d\n", active_pps->pic_order_present_flag);
    printf ("delta_pic_order_cnt[0]                %d\n", img->delta_pic_order_cnt[0]);
    printf ("delta_pic_order_cnt[1]                %d\n", img->delta_pic_order_cnt[1]);
    printf ("delta_pic_order_cnt[2]                %d\n", img->delta_pic_order_cnt[2]);
    printf ("idr_flag                              %d\n", img->idr_flag);
    printf ("MaxFrameNum                           %d\n", img->MaxFrameNum);

    return 0;
}

/*!
 ************************************************************************
 * \brief
 *    return the poc of img as per (8-1) JVT-F100d2
 *  POC200301
 ************************************************************************
 */
int picture_order(struct img_par *img)
{
  if (img->field_pic_flag==0) // is a frame
    return img->framepoc;
  else if (img->bottom_field_flag==0) // top field
    return img->toppoc;
  else // bottom field
    return img->bottompoc;
}









/////////////////////////////////////<














