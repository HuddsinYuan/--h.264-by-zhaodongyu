
/*!
 **************************************************************************************
 * \file
 *    parset.c
 * \brief
 *    Picture and Sequence Parameter set generation and handling
 *  \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 *
 **************************************************************************************
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <malloc.h>
#include <string.h>
 
#include "global.h"
// #include "contributors.h"
#include "parsetcommon.h"
#include "nalu.h"
#include "parset.h"
#include "fmo.h"
#include "vlc.h"
#include "mbuffer.h"
#include <cabac.h>
#include <erc_api.h>




// Local helpers
static int IdentifyProfile();
static int IdentifyLevel();
static int IdentifyNumRefFrames();
static int GenerateVUISequenceParameters();

extern ColocatedParams *Co_located;

/////////////////////////////////////>
extern int UsedBits;      // for internal statistics, is adjusted by se_v, ue_v, u_1
// extern ColocatedParams *Co_located;

seq_parameter_set_rbsp_t SeqParSet[MAXSPS];//全局变量++ MAXSPS=32，它是一个指向序列参数集的数组。
pic_parameter_set_rbsp_t PicParSet[MAXPPS];


extern StorablePicture* dec_picture;




/////////////////////////////////////<
/*! 
 *************************************************************************************
 * \brief
 *    generates a sequence and picture parameter set and stores these in global
 *    active_sps and active_pps
 *
 * \return
 *    A NALU containing the Sequence ParameterSet
 *
 *************************************************************************************
*/
void GenerateParameterSets ()
{
  seq_parameter_set_rbsp_t *sps = NULL; 
  pic_parameter_set_rbsp_t *pps = NULL;

  sps = AllocSPS();
  pps = AllocPPS();

  FillParameterSetStructures (sps, pps);
  
  active_sps = sps;
  active_pps = pps;
}

/*! 
*************************************************************************************
* \brief
*    frees global parameter sets active_sps and active_pps
*
* \return
*    A NALU containing the Sequence ParameterSet
*
*************************************************************************************
*/
void FreeParameterSets ()
{
  FreeSPS (active_sps);
  FreePPS (active_pps);
}

/*! 
*************************************************************************************
* \brief
*    int GenerateSeq_parameter_set_NALU ();
*
* \note
*    Uses the global variables through FillParameterSetStructures()
*
* \return
*    A NALU containing the Sequence ParameterSet
*
*************************************************************************************
*/

NALU_t *GenerateSeq_parameter_set_NALU ()
{
  NALU_t *n = AllocNALU(64000);
  int RBSPlen = 0;
  int NALUlen;
  byte rbsp[MAXRBSPSIZE];

  RBSPlen = GenerateSeq_parameter_set_rbsp (active_sps, rbsp);//8
  NALUlen = RBSPtoNALU (rbsp, n, RBSPlen, NALU_TYPE_SPS, NALU_PRIORITY_HIGHEST, 0, 1);//9
  n->startcodeprefix_len = 4;

  return n;
}


/*! 
*************************************************************************************
* \brief
*    NALU_t *GeneratePic_parameter_set_NALU ();
*
* \note
*    Uses the global variables through FillParameterSetStructures()
*
* \return
*    A NALU containing the Picture Parameter Set
*
*************************************************************************************
*/

NALU_t *GeneratePic_parameter_set_NALU()
{
  NALU_t *n = AllocNALU(64000);
  int RBSPlen = 0;
  int NALUlen;
  byte rbsp[MAXRBSPSIZE];

  RBSPlen = GeneratePic_parameter_set_rbsp (active_pps, rbsp);
  NALUlen = RBSPtoNALU (rbsp, n, RBSPlen, NALU_TYPE_PPS, NALU_PRIORITY_HIGHEST, 0, 1);
  n->startcodeprefix_len = 4;

  return n;
}

/*!
 ************************************************************************
 * \brief
 *    FillParameterSetStructures: extracts info from global variables and
 *    generates a picture and sequence parameter set structure
 *
 * \param sps
 *    Sequence parameter set to be filled
 * \param pps
 *    Picture parameter set to be filled
 * \par
 *    Function reads all kinds of values from several global variables,
 *    including input-> and image-> and fills in the sps and pps.  Many
 *    values are current hard-coded to defaults, especially most of the
 *    VUI stuff.  Currently, the sps and pps structures are fixed length
 *    with the exception of the fully flexible FMO map (mode 6).  This
 *    mode is not supported.  Hence, the function does not need to
 *    allocate memory for the FMOmap, the pointer slice_group_id is
 *    always NULL.  If one wants to implement FMO mode 6, one would need
 *    to malloc the memory for the map here, and the caller would need
 *    to free it after use.
 *
 * \par 
 *    Limitations
 *    Currently, the encoder does not support multiple parameter sets,
 *    primarily because the config file does not support it.  Hence the
 *    pic_parameter_set_id and the seq_parameter_set_id are always zero.
 *    If one day multiple parameter sets are implemented, it would
 *    make sense to break this function into two, one for the picture and
 *    one for the sequence.
 *    Currently, FMO is not supported
 *    The following pps and sps elements seem not to be used in the encoder 
 *    or decoder and, hence, a guessed default value is conveyed:
 *
 *    pps->num_ref_idx_l0_active_minus1 = img->num_ref_pic_active_fwd_minus1;
 *    pps->num_ref_idx_l1_active_minus1 = img->num_ref_pic_active_bwd_minus1;
 *    pps->chroma_qp_index_offset = 0;
 *    sps->required_frame_num_update_behaviour_flag = FALSE;
 *    sps->direct_temporal_constrained_flag = FALSE;
 *
 * \par
 *    Regarding the QP
 *    The previous software versions coded the absolute QP only in the 
 *    slice header.  This is kept, and the offset in the PPS is coded 
 *    even if we could save bits by intelligently using this field.
 *
 ************************************************************************
 */

void FillParameterSetStructures (seq_parameter_set_rbsp_t *sps, 
                                 pic_parameter_set_rbsp_t *pps)
{
  unsigned i;
  // *************************************************************************
  // Sequence Parameter Set
  // *************************************************************************
  assert (sps != NULL);
  assert (pps != NULL);
  // Profile and Level should be calculated using the info from the config
  // file.  Calculation is hidden in IndetifyProfile() and IdentifyLevel()
  sps->profile_idc=IdentifyProfile();//  77
  sps->level_idc=IdentifyLevel();//  30

  // needs to be set according to profile
  sps->constrained_set0_flag = 0;
  sps->constrained_set1_flag = 0;
  sps->constrained_set2_flag = 0;

  // Parameter Set ID hardcoded to zero
  sps->seq_parameter_set_id = 0;

  //! POC stuff:
  //! The following values are hard-coded in init_poc().  Apparently,
  //! the poc implementation covers only a subset of the poc functionality.
  //! Here, the same subset is implemented.  Changes in the POC stuff have
  //! also to be reflected here
  sps->log2_max_frame_num_minus4 = log2_max_frame_num_minus4;
  sps->log2_max_pic_order_cnt_lsb_minus4 = log2_max_pic_order_cnt_lsb_minus4;
  
  sps->pic_order_cnt_type = input->pic_order_cnt_type;
  sps->num_ref_frames_in_pic_order_cnt_cycle = img->num_ref_frames_in_pic_order_cnt_cycle; //1
  sps->delta_pic_order_always_zero_flag = img->delta_pic_order_always_zero_flag;
  sps->offset_for_non_ref_pic = img->offset_for_non_ref_pic;
  sps->offset_for_top_to_bottom_field = img->offset_for_top_to_bottom_field;

  for (i=0; i<img->num_ref_frames_in_pic_order_cnt_cycle; i++)
  {
    sps->offset_for_ref_frame[i] = img->offset_for_ref_frame[i];
  }
  // End of POC stuff

  // Number of Reference Frames
  sps->num_ref_frames = IdentifyNumRefFrames();//5

  //required_frame_num_update_behaviour_flag hardcoded to zero
  sps->gaps_in_frame_num_value_allowed_flag = FALSE1;    // double check

  sps->frame_mbs_only_flag = !(input->PicInterlace || input->MbInterlace);//1

  // Picture size, finally a simple one :-)
  sps->pic_width_in_mbs_minus1 = (input->imagewidth/16) -1;//39
  sps->pic_height_in_map_units_minus1 = ((input->imageheight/16)/ (2 - sps->frame_mbs_only_flag)) - 1;//29

  // a couple of flags, simple
  sps->mb_adaptive_frame_field_flag = (FRAME_CODING != input->MbInterlace);
  sps->direct_8x8_inference_flag = input->directInferenceFlag;

  // Sequence VUI not implemented, signalled as not present
  sps->vui_parameters_present_flag = FALSE1;
  
  {
    int PicWidthInMbs, PicHeightInMapUnits, FrameHeightInMbs;
    int width, height;
    PicWidthInMbs = (sps->pic_width_in_mbs_minus1 +1);
    PicHeightInMapUnits = (sps->pic_height_in_map_units_minus1 +1);
    FrameHeightInMbs = ( 2 - sps->frame_mbs_only_flag ) * PicHeightInMapUnits;
    
    width = PicWidthInMbs * MB_BLOCK_SIZE;
    height = FrameHeightInMbs * MB_BLOCK_SIZE;
    
    Co_located = alloc_colocated (width, height,sps->mb_adaptive_frame_field_flag);
    
  }
  // *************************************************************************
  // Picture Parameter Set 
  // *************************************************************************

  pps->seq_parameter_set_id = 0;
  pps->pic_parameter_set_id = 0;
  pps->entropy_coding_mode_flag = (input->symbol_mode==UVLC?0:1);

  // JVT-Fxxx (by Stephan Wenger, make this flag unconditional
  pps->pic_order_present_flag = img->pic_order_present_flag;


  // Begin FMO stuff
  pps->num_slice_groups_minus1 = input->num_slice_groups_minus1;

	
  //! Following set the parameter for different slice group types
  if (pps->num_slice_groups_minus1 > 0)
    switch (input->slice_group_map_type)
    {
    case 0:
			
      pps->slice_group_map_type = 0;
      for(i=0; i<=pps->num_slice_groups_minus1; i++)
      {
        pps->run_length_minus1[i]=input->run_length_minus1[i];
      }
			
      break;
    case 1:
      pps->slice_group_map_type = 1;
      break;
    case 2:
      // i loops from 0 to num_slice_groups_minus1-1, because no info for background needed
      pps->slice_group_map_type = 2;
      for(i=0; i<pps->num_slice_groups_minus1; i++)
      {
        pps->top_left[i] = input->top_left[i];
        pps->bottom_right[i] = input->bottom_right[i];      
      }
     break;
    case 3:
    case 4:
    case 5:
      pps->slice_group_map_type = input->slice_group_map_type;
			
      pps->slice_group_change_direction_flag = input->slice_group_change_direction_flag;
      pps->slice_group_change_rate_minus1 = input->slice_group_change_rate_minus1;
      break;
    case 6:
      pps->slice_group_map_type = 6;   
      pps->pic_size_in_map_units_minus1 = 
				((input->imageheight/MB_BLOCK_SIZE)/(2-sps->frame_mbs_only_flag))
				*(input->imagewidth/MB_BLOCK_SIZE) -1;
			
      for (i=0;i<=pps->pic_size_in_map_units_minus1; i++)
        pps->slice_group_id[i] = input->slice_group_id[i];
			
      break;
    default:
      printf ("Parset.c: slice_group_map_type invalid, default\n");
      assert (0==1);
    }
// End FMO stuff

  pps->num_ref_idx_l0_active_minus1 = sps->frame_mbs_only_flag ? (sps->num_ref_frames-1) : (2 * sps->num_ref_frames - 1) ;   // set defaults
  pps->num_ref_idx_l1_active_minus1 = sps->frame_mbs_only_flag ? (sps->num_ref_frames-1) : (2 * sps->num_ref_frames - 1) ;   // set defaults
  //pps->num_ref_idx_l1_active_minus1 = sps->frame_mbs_only_flag ? 0 : 1 ;   // set defaults

  
  pps->weighted_pred_flag = input->WeightedPrediction;
  pps->weighted_bipred_idc = input->WeightedBiprediction;

  pps->pic_init_qp_minus26 = 0;         // hard coded to zero, QP lives in the slice header
  pps->pic_init_qs_minus26 = 0;

  pps->chroma_qp_index_offset = input->chroma_qp_index_offset;      // double check: is this chroma fidelity thing already implemented???

  pps->deblocking_filter_control_present_flag = input->LFSendParameters;
  pps->constrained_intra_pred_flag = input->UseConstrainedIntraPred;
  
  pps->redundant_pic_cnt_present_flag = 0;

  // the picture vui consists currently of the cropping rectangle, which cannot
  // used by the current decoder and hence is never sent.
  sps->frame_cropping_flag = FALSE1;
};



/*! 
 *************************************************************************************
 * \brief
 *    int GenerateSeq_parameter_set_rbsp (seq_parameter_set_rbsp_t *sps, char *rbsp);
 *
 * \param sps
 *    sequence parameter structure
 * \param rbsp
 *    buffer to be filled with the rbsp, size should be at least MAXIMUMPARSETRBSPSIZE
 *
 * \return
 *    size of the RBSP in bytes
 *
 * \note
 *    Sequence Parameter VUI function is called, but the function implements
 *    an exit (-1)
 *************************************************************************************
 */
 
int GenerateSeq_parameter_set_rbsp (seq_parameter_set_rbsp_t *sps, char *rbsp)
{
  DataPartition *partition;
  int len = 0, LenInBytes;
  unsigned i;

  assert (rbsp != NULL);
  // In order to use the entropy coding functions from golomb.c we need 
  // to allocate a partition structure.  It will be freed later in this
  // function
  if ((partition=calloc(1,sizeof(DataPartition)))==NULL) no_mem_exit("SeqParameterSet:partition");
  if ((partition->bitstream=calloc(1, sizeof(Bitstream)))==NULL) no_mem_exit("SeqParameterSet:bitstream");
  // .. and use the rbsp provided (or allocated above) for the data
  partition->bitstream->streamBuffer = rbsp;
  partition->bitstream->bits_to_go = 8;

  len+=u_v  (8, "SPS: profile_idc",                             sps->profile_idc,              partition);//8

  len+=u_1  ("SPS: constrained_set0_flag",                      sps->constrained_set0_flag,    partition);//9
  len+=u_1  ("SPS: constrained_set1_flag",                      sps->constrained_set1_flag,    partition);//10
  len+=u_1  ("SPS: constrained_set2_flag",                      sps->constrained_set2_flag,    partition);//11
  len+=u_v  (5, "SPS: reserved_zero",                           0,                             partition);//16

  len+=u_v  (8, "SPS: level_idc",                               sps->level_idc,                partition);//24

  len+=ue_v ("SPS: seq_parameter_set_id",                    sps->seq_parameter_set_id,        partition);//25
  len+=ue_v ("SPS: log2_max_frame_num_minus4",               sps->log2_max_frame_num_minus4,   partition);//26
  len+=ue_v ("SPS: pic_order_cnt_type",                      sps->pic_order_cnt_type,          partition);//27
  // POC200301
  if (sps->pic_order_cnt_type == 0)
    len+=ue_v ("SPS: log2_max_pic_order_cnt_lsb_minus4",     sps->log2_max_pic_order_cnt_lsb_minus4,         partition);//28
  else if (sps->pic_order_cnt_type == 1)
  {
    len+=u_1  ("SPS: delta_pic_order_always_zero_flag",        sps->delta_pic_order_always_zero_flag,        partition);
    len+=se_v ("SPS: offset_for_non_ref_pic",                  sps->offset_for_non_ref_pic,                  partition);
    len+=se_v ("SPS: offset_for_top_to_bottom_field",          sps->offset_for_top_to_bottom_field,          partition);
    len+=ue_v ("SPS: num_ref_frames_in_pic_order_cnt_cycle",   sps->num_ref_frames_in_pic_order_cnt_cycle,   partition);
    for (i=0; i<sps->num_ref_frames_in_pic_order_cnt_cycle; i++)
      len+=se_v ("SPS: offset_for_ref_frame",                  sps->offset_for_ref_frame[i],                 partition);
  }
  len+=ue_v ("SPS: num_ref_frames",                          sps->num_ref_frames,                            partition);//33
  len+=u_1  ("SPS: gaps_in_frame_num_value_allowed_flag",    sps->gaps_in_frame_num_value_allowed_flag,      partition);//34
  len+=ue_v ("SPS: pic_width_in_mbs_minus1",                 sps->pic_width_in_mbs_minus1,                   partition);//43
  len+=ue_v ("SPS: pic_height_in_map_units_minus1",          sps->pic_height_in_map_units_minus1,            partition);//52
  len+=u_1  ("SPS: frame_mbs_only_flag",                     sps->frame_mbs_only_flag,                       partition);//53
  if (!sps->frame_mbs_only_flag)//0
  {
    len+=u_1  ("SPS: mb_adaptive_frame_field_flag",           sps->mb_adaptive_frame_field_flag,             partition);
  }
  len+=u_1  ("SPS: direct_8x8_inference_flag",                sps->direct_8x8_inference_flag,                partition);//54

  len+=u_1  ("SPS: frame_cropping_flag",                      sps->frame_cropping_flag,                      partition);//55
  if (sps->frame_cropping_flag)//0
  {
    len+=ue_v ("SPS: frame_cropping_rect_left_offset",          sps->frame_cropping_rect_left_offset,        partition);
    len+=ue_v ("SPS: frame_cropping_rect_right_offset",         sps->frame_cropping_rect_right_offset,       partition);
    len+=ue_v ("SPS: frame_cropping_rect_top_offset",           sps->frame_cropping_rect_top_offset,         partition);
    len+=ue_v ("SPS: frame_cropping_rect_bottom_offset",        sps->frame_cropping_rect_bottom_offset,      partition);
  }

  len+=u_1  ("SPS: vui_parameters_present_flag",                sps->vui_parameters_present_flag,            partition);//56
  if (sps->vui_parameters_present_flag)//0
    len+=GenerateVUISequenceParameters();    // currently a dummy, asserting

  SODBtoRBSP(partition->bitstream);     // copies the last couple of bits into the byte buffer
  
  LenInBytes=partition->bitstream->byte_pos;//8

  free (partition->bitstream);
  free (partition);
  
  return LenInBytes;
}


/*! 
 *************************************************************************************
 * \brief
 *    int GeneratePic_parameter_set_rbsp (pic_parameter_set_rbsp_t *sps, char *rbsp);
 *
 * \param pps
 *    picture parameter structure
 * \param rbsp
 *    buffer to be filled with the rbsp, size should be at least MAXIMUMPARSETRBSPSIZE
 *
 * \return
 *    size of the RBSP in bytes, negative in case of an error
 *
 * \note
 *    Picture Parameter VUI function is called, but the function implements
 *    an exit (-1)
 *************************************************************************************
 */
 
int GeneratePic_parameter_set_rbsp (pic_parameter_set_rbsp_t *pps, char *rbsp)
{
  DataPartition *partition;
  int len = 0, LenInBytes;
  unsigned i;
  unsigned NumberBitsPerSliceGroupId;

  assert (rbsp != NULL);

  // In order to use the entropy coding functions from golomb.c we need 
  // to allocate a partition structure.  It will be freed later in this
  // function
  if ((partition=calloc(1,sizeof(DataPartition)))==NULL) no_mem_exit("PicParameterSet:partition");
  if ((partition->bitstream=calloc(1, sizeof(Bitstream)))==NULL) no_mem_exit("PicParameterSet:bitstream");
  // .. and use the rbsp provided (or allocated above) for the data
  partition->bitstream->streamBuffer = rbsp;
  partition->bitstream->bits_to_go = 8;
  //sw paff
  pps->pic_order_present_flag = img->pic_order_present_flag;

  len+=ue_v ("PPS: pic_parameter_set_id",                    pps->pic_parameter_set_id,                      partition);
  len+=ue_v ("PPS: seq_parameter_set_id",                    pps->seq_parameter_set_id,                      partition);
  len+=u_1  ("PPS: entropy_coding_mode_flag",                pps->entropy_coding_mode_flag,                  partition);
  len+=u_1  ("PPS: pic_order_present_flag",                  pps->pic_order_present_flag,                    partition);
  len+=ue_v ("PPS: num_slice_groups_minus1",                 pps->num_slice_groups_minus1,                   partition);

  // FMO stuff
  if(pps->num_slice_groups_minus1 > 0 )
  {
    len+=ue_v ("PPS: slice_group_map_type",                 pps->slice_group_map_type,                   partition);
    if (pps->slice_group_map_type == 0)
      for (i=0; i<=pps->num_slice_groups_minus1; i++)
        len+=ue_v ("PPS: run_length_minus1[i]",                           pps->run_length_minus1[i],                             partition);
    else if (pps->slice_group_map_type==2)
      for (i=0; i<pps->num_slice_groups_minus1; i++)
      {

        len+=ue_v ("PPS: top_left[i]",                          pps->top_left[i],                           partition);
        len+=ue_v ("PPS: bottom_right[i]",                      pps->bottom_right[i],                       partition);
      }
    else if (pps->slice_group_map_type == 3 ||
             pps->slice_group_map_type == 4 ||
             pps->slice_group_map_type == 5) 
    {
      len+=u_1  ("PPS: slice_group_change_direction_flag",         pps->slice_group_change_direction_flag,         partition);
      len+=ue_v ("PPS: slice_group_change_rate_minus1",            pps->slice_group_change_rate_minus1,            partition);
    } 
    else if (pps->slice_group_map_type == 6)
    {
      if (pps->num_slice_groups_minus1>=4)
        NumberBitsPerSliceGroupId=3;
      else if (pps->num_slice_groups_minus1>=2)
        NumberBitsPerSliceGroupId=2;
      else if (pps->num_slice_groups_minus1>=1)
        NumberBitsPerSliceGroupId=1;
      else
        NumberBitsPerSliceGroupId=0;
        
      len+=ue_v ("PPS: pic_size_in_map_units_minus1",          pps->pic_size_in_map_units_minus1,             partition);
      for(i=0; i<=pps->pic_size_in_map_units_minus1; i++)
        len+= u_v  (NumberBitsPerSliceGroupId, "PPS: >slice_group_id[i]",   pps->slice_group_id[i],           partition);			
    }
  }
  // End of FMO stuff

  len+=ue_v ("PPS: num_ref_idx_l0_active_minus1",             pps->num_ref_idx_l0_active_minus1,              partition);
  len+=ue_v ("PPS: num_ref_idx_l1_active_minus1",             pps->num_ref_idx_l1_active_minus1,              partition);
  len+=u_1  ("PPS: weighted_pred_flag",                       pps->weighted_pred_flag,                        partition);
  len+=u_v  (2, "PPS: weighted_bipred_idc",                   pps->weighted_bipred_idc,                       partition);
  len+=se_v ("PPS: pic_init_qp_minus26",                      pps->pic_init_qp_minus26,                       partition);
  len+=se_v ("PPS: pic_init_qs_minus26",                      pps->pic_init_qs_minus26,                       partition);
  len+=se_v ("PPS: chroma_qp_index_offset",                   pps->chroma_qp_index_offset,                    partition);
  len+=u_1  ("PPS: deblocking_filter_control_present_flag",   pps->deblocking_filter_control_present_flag,    partition);
  len+=u_1  ("PPS: constrained_intra_pred_flag",              pps->constrained_intra_pred_flag,               partition);
  len+=u_1  ("PPS: redundant_pic_cnt_present_flag",           pps->redundant_pic_cnt_present_flag,            partition);

  SODBtoRBSP(partition->bitstream);     // copies the last couple of bits into the byte buffer
  
  LenInBytes=partition->bitstream->byte_pos;

  // Get rid of the helper structures
  free (partition->bitstream);
  free (partition);

  return LenInBytes;
}



/*! 
 *************************************************************************************
 * \brief
 *    Returns the Profile
 *
 * \return
 *    Profile according to Annex A
 *
 * \note
 *    Function is currently a dummy.  Should "calculate" the profile from those
 *    config file parameters.  E.g.
 *
 *    Profile = Baseline;
 *    if (CABAC Used || Interlace used) Profile=Main;
 *    if (!Cabac Used) && (Bframes | SPframes) Profile = Streaming;
 *
 *************************************************************************************
 */
int IdentifyProfile()
{
  return input->ProfileIDC;
};

/*! 
 *************************************************************************************
 * \brief
 *    Returns the Level
 *
 * \return
 *    Level according to Annex A
 *
 * \note
 *    This function is currently a dummy, but should calculate the level out of 
 *    the config file parameters (primarily the picture size)
 *************************************************************************************
 */
int IdentifyLevel()
{
  return input->LevelIDC;
};


/*! 
 *************************************************************************************
 * \brief
 *    Returns the number of reference frame buffers
 *
 * \return
 *    Number of reference frame buffers used
 *
 * \note
 *    This function currently maps to input->num_reference_frames.  With all this interlace
 *    stuff this may or may not be correct.  If you determine a problem with the
 *    memory management for Interlace, then this could be one possible problem.
 *    However, so far no problem have been determined by my limited testing of
 *    a stupid 1950's technology :-)  StW, 11/27/02
 *************************************************************************************
 */

int IdentifyNumRefFrames()
{
  if(input->num_reference_frames > 16)error("no ref frames too large",-100);
  
  return input->num_reference_frames;
}


/*! 
 *************************************************************************************
 * \brief
 *    Function body for VUI Parameter generation (to be done)
 *
 * \return
 *    exits with error message
 *************************************************************************************
 */
static int GenerateVUISequenceParameters()
{
  printf ("Sequence Parameter VUI not yet implemented, this should never happen, exit\n");
  exit (-1);
}
/////////////////////////////////////>

// fill sps with content of p

int InterpretSPS (DataPartition *p, seq_parameter_set_rbsp_t *sps)
{
	unsigned i;
	int reserved_zero;
	Bitstream *s = p->bitstream;
	
	assert (p != NULL);
	assert (p->bitstream != NULL);
	assert (p->bitstream->streamBuffer != 0);
	assert (sps != NULL);
	
	UsedBits = 0;
	
	sps->profile_idc                            = u_v_dec  (8, "SPS: profile_idc"                           , s);//77
	
	sps->constrained_set0_flag                  = u_1_dec  (   "SPS: constrained_set0_flag"                 , s);//0
	sps->constrained_set1_flag                  = u_1_dec  (   "SPS: constrained_set1_flag"                 , s);//0
	sps->constrained_set2_flag                  = u_1_dec  (   "SPS: constrained_set2_flag"                 , s);//0
	reserved_zero                               = u_v_dec  (5, "SPS: reserved_zero_5bits"                   , s);//0
	assert (reserved_zero==0);
	
	sps->level_idc                              = u_v_dec  (8, "SPS: level_idc"                             , s);//30
	
	
	sps->seq_parameter_set_id                   = ue_v_dec ("SPS: seq_parameter_set_id"                     , s);//0
	sps->log2_max_frame_num_minus4              = ue_v_dec ("SPS: log2_max_frame_num_minus4"                , s);//0
	sps->pic_order_cnt_type                     = ue_v_dec ("SPS: pic_order_cnt_type"                       , s);//0
	
	if (sps->pic_order_cnt_type == 0)
		sps->log2_max_pic_order_cnt_lsb_minus4 = ue_v_dec ("SPS: log2_max_pic_order_cnt_lsb_minus4"           , s);//0,s是bitstream
	else if (sps->pic_order_cnt_type == 1)
	{
		sps->delta_pic_order_always_zero_flag      = u_1_dec  ("SPS: delta_pic_order_always_zero_flag"       , s);
		sps->offset_for_non_ref_pic                = se_v_dec ("SPS: offset_for_non_ref_pic"                 , s);
		sps->offset_for_top_to_bottom_field        = se_v_dec ("SPS: offset_for_top_to_bottom_field"         , s);
		sps->num_ref_frames_in_pic_order_cnt_cycle = ue_v_dec ("SPS: num_ref_frames_in_pic_order_cnt_cycle"  , s);
		for(i=0; i<sps->num_ref_frames_in_pic_order_cnt_cycle; i++)
			sps->offset_for_ref_frame[i]               = se_v_dec ("SPS: offset_for_ref_frame[i]"              , s);
	}
	sps->num_ref_frames                        = ue_v_dec ("SPS: num_ref_frames"                         , s);//5
	sps->gaps_in_frame_num_value_allowed_flag  = u_1_dec  ("SPS: gaps_in_frame_num_value_allowed_flag"   , s);//0
	sps->pic_width_in_mbs_minus1               = ue_v_dec ("SPS: pic_width_in_mbs_minus1"                , s);//21
	sps->pic_height_in_map_units_minus1        = ue_v_dec ("SPS: pic_height_in_map_units_minus1"         , s);//17
	sps->frame_mbs_only_flag                   = u_1_dec  ("SPS: frame_mbs_only_flag"                    , s);//1
	if (!sps->frame_mbs_only_flag)
	{
		sps->mb_adaptive_frame_field_flag          = u_1_dec  ("SPS: mb_adaptive_frame_field_flag"           , s);
	}
	sps->direct_8x8_inference_flag             = u_1_dec  ("SPS: direct_8x8_inference_flag"              , s);//0
	sps->frame_cropping_flag                   = u_1_dec  ("SPS: frame_cropping_flag"                , s);//0
	
	if (sps->frame_cropping_flag)
	{
		sps->frame_cropping_rect_left_offset      = ue_v_dec ("SPS: frame_cropping_rect_left_offset"           , s);
		sps->frame_cropping_rect_right_offset     = ue_v_dec ("SPS: frame_cropping_rect_right_offset"          , s);
		sps->frame_cropping_rect_top_offset       = ue_v_dec ("SPS: frame_cropping_rect_top_offset"            , s);
		sps->frame_cropping_rect_bottom_offset    = ue_v_dec ("SPS: frame_cropping_rect_bottom_offset"         , s);
	}
	sps->vui_parameters_present_flag           = u_1_dec  ("SPS: vui_parameters_present_flag"            , s);//0
	if (sps->vui_parameters_present_flag)
	{
		printf ("VUI sequence parameters present but not supported, ignored, proceeding to next NALU\n");
	}
	sps->Valid = TRUE1;
	return UsedBits;//56
}

int InterpretPPS (DataPartition *p, pic_parameter_set_rbsp_t *pps)
{
	unsigned i;
	int NumberBitsPerSliceGroupId;
	Bitstream *s = p->bitstream;
	
	assert (p != NULL);
	assert (p->bitstream != NULL);
	assert (p->bitstream->streamBuffer != 0);
	assert (pps != NULL);
	
	UsedBits = 0;
	
	pps->pic_parameter_set_id                  = ue_v_dec ("PPS: pic_parameter_set_id"                   , s);
	pps->seq_parameter_set_id                  = ue_v_dec ("PPS: seq_parameter_set_id"                   , s);
	pps->entropy_coding_mode_flag              = u_1_dec  ("PPS: entropy_coding_mode_flag"               , s);
	
	//! Note: as per JVT-F078 the following bit is unconditional.  If F078 is not accepted, then
	//! one has to fetch the correct SPS to check whether the bit is present (hopefully there is
	//! no consistency problem :-(
	//! The current encoder code handles this in the same way.  When you change this, don't forget
	//! the encoder!  StW, 12/8/02
	pps->pic_order_present_flag                = u_1_dec  ("PPS: pic_order_present_flag"                 , s);
	
	pps->num_slice_groups_minus1               = ue_v_dec ("PPS: num_slice_groups_minus1"                , s);
	
	// FMO stuff begins here
	if (pps->num_slice_groups_minus1 > 0)
	{
		pps->slice_group_map_type               = ue_v_dec ("PPS: slice_group_map_type"                , s);
		if (pps->slice_group_map_type == 0)
		{
			for (i=0; i<=pps->num_slice_groups_minus1; i++)
				pps->run_length_minus1 [i]                  = ue_v_dec ("PPS: run_length_minus1 [i]"              , s);
		}
		else if (pps->slice_group_map_type == 2)
		{
			for (i=0; i<pps->num_slice_groups_minus1; i++)
			{
				//! JVT-F078: avoid reference of SPS by using ue(v) instead of u(v)
				pps->top_left [i]                          = ue_v_dec ("PPS: top_left [i]"                        , s);
				pps->bottom_right [i]                      = ue_v_dec ("PPS: bottom_right [i]"                    , s);
			}
		}
		else if (pps->slice_group_map_type == 3 ||
			pps->slice_group_map_type == 4 ||
			pps->slice_group_map_type == 5)
		{
			pps->slice_group_change_direction_flag     = u_1_dec  ("PPS: slice_group_change_direction_flag"      , s);
			pps->slice_group_change_rate_minus1        = ue_v_dec ("PPS: slice_group_change_rate_minus1"         , s);
		}
		else if (pps->slice_group_map_type == 6)
		{
			if (pps->num_slice_groups_minus1+1 >4)
				NumberBitsPerSliceGroupId = 3;
			else if (pps->num_slice_groups_minus1+1 > 2)
				NumberBitsPerSliceGroupId = 2;
			else
				NumberBitsPerSliceGroupId = 1;
			//! JVT-F078, exlicitly signal number of MBs in the map
			pps->num_slice_group_map_units_minus1      = ue_v_dec ("PPS: num_slice_group_map_units_minus1"               , s);
			for (i=0; i<=pps->num_slice_group_map_units_minus1; i++)
				pps->slice_group_id[i] = u_v_dec (NumberBitsPerSliceGroupId, "slice_group_id[i]", s);
		}
	}
	
// 	// End of FMO stuff
	
	pps->num_ref_idx_l0_active_minus1          = ue_v_dec ("PPS: num_ref_idx_l0_active_minus1"           , s);
	pps->num_ref_idx_l1_active_minus1          = ue_v_dec ("PPS: num_ref_idx_l1_active_minus1"           , s);
	pps->weighted_pred_flag                    = u_1_dec  ("PPS: weighted prediction flag"               , s);
	pps->weighted_bipred_idc                   = u_v_dec  ( 2, "PPS: weighted_bipred_idc"                , s);
	pps->pic_init_qp_minus26                   = se_v_dec ("PPS: pic_init_qp_minus26"                    , s);
	pps->pic_init_qs_minus26                   = se_v_dec ("PPS: pic_init_qs_minus26"                    , s);
	pps->chroma_qp_index_offset                = se_v_dec ("PPS: chroma_qp_index_offset"                 , s);
	pps->deblocking_filter_control_present_flag = u_1_dec ("PPS: deblocking_filter_control_present_flag" , s);
	pps->constrained_intra_pred_flag           = u_1_dec  ("PPS: constrained_intra_pred_flag"            , s);
	pps->redundant_pic_cnt_present_flag        = u_1_dec  ("PPS: redundant_pic_cnt_present_flag"         , s);
	
	pps->Valid = TRUE1;
	return UsedBits;
}


void PPSConsistencyCheck (pic_parameter_set_rbsp_t *pps)
{
	printf ("Consistency checking a picture parset, to be implemented\n");
	//  if (pps->seq_parameter_set_id invalid then do something)
}

void SPSConsistencyCheck (seq_parameter_set_rbsp_t *sps)
{
	printf ("Consistency checking a sequence parset, to be implemented\n");
}

void MakePPSavailable (int id, pic_parameter_set_rbsp_t *pps)
{
	assert (pps->Valid == TRUE1);
	
	if (PicParSet[id].Valid == TRUE1 && PicParSet[id].slice_group_id != NULL)
		free (PicParSet[id].slice_group_id);
	
	memcpy (&PicParSet[id], pps, sizeof (pic_parameter_set_rbsp_t));
	
	if ((PicParSet[id].slice_group_id = calloc (PicParSet[id].num_slice_group_map_units_minus1+1, sizeof(int))) == NULL)
		no_mem_exit ("MakePPSavailable: Cannot calloc slice_group_id");
	
	memcpy (PicParSet[id].slice_group_id, pps->slice_group_id, (pps->num_slice_group_map_units_minus1+1)*sizeof(int));
}

void MakeSPSavailable (int id, seq_parameter_set_rbsp_t *sps)
{
	assert (sps->Valid == TRUE1);
    //将读入的sps存入全局变量数组SeqParSet，以备后用
	memcpy (&SeqParSet[id], sps, sizeof (seq_parameter_set_rbsp_t));
}


void ProcessSPS (NALU_t *nalu)//解码器先将码流中的数据读入临时指针sps，
//之后存入全局变量数组SeqParSet，最后需要使用这些参数时，将SeqParSet中的数据读入active_sps，同理pps。
{
	DataPartition *dp = AllocPartition(1);
	seq_parameter_set_rbsp_t *sps = AllocSPS();// 返回值为指向序列参数集的指针
	int dummy;
	
	memcpy (dp->bitstream->streamBuffer, &nalu->buf[1], nalu->len-1);
	dp->bitstream->code_len = dp->bitstream->bitstream_length = RBSPtoSODB (dp->bitstream->streamBuffer, nalu->len-1);
	dp->bitstream->ei_flag = 0;
	dp->bitstream->read_len = dp->bitstream->frame_bitoffset = 0;
	//在此函数里将码流中的句法元素读入sps
	dummy = InterpretSPS (dp, sps);//log2_max_pic_order_cnt_lsb_minus4由0变为1   //58 
	
	if (active_sps)
	{
		if (sps->seq_parameter_set_id == active_sps->seq_parameter_set_id)
		{
			if (!sps_is_equal(sps, active_sps))
			{
				if (dec_picture)
				{
					// this may only happen on slice loss
					exit_picture();
				}
				active_sps=NULL;
			}
		}
	}
	// SPSConsistencyCheck (pps);
	MakeSPSavailable (sps->seq_parameter_set_id, sps);
	
	FreePartition (dp, 1);
	FreeSPS (sps);
}


void ProcessPPS (NALU_t *nalu)
{
	DataPartition *dp;
	pic_parameter_set_rbsp_t *pps;
	int dummy;
	
	dp = AllocPartition(1);
	pps = AllocPPS();
	memcpy (dp->bitstream->streamBuffer, &nalu->buf[1], nalu->len-1);
	dp->bitstream->code_len = dp->bitstream->bitstream_length = RBSPtoSODB (dp->bitstream->streamBuffer, nalu->len-1);
	dp->bitstream->ei_flag = 0;
	dp->bitstream->read_len = dp->bitstream->frame_bitoffset = 0;
	dummy = InterpretPPS (dp, pps);
	// PPSConsistencyCheck (pps);
	if (active_pps)
	{
		if (pps->pic_parameter_set_id == active_pps->pic_parameter_set_id)
		{
			if (!pps_is_equal(pps, active_pps))
			{
				if (dec_picture)
				{
					// this may only happen on slice loss
					exit_picture();
				}
				active_pps = NULL;
			}
		}
	}
	MakePPSavailable (pps->pic_parameter_set_id, pps);
	FreePartition (dp, 1);
	FreePPS (pps);
}

void activate_sps (seq_parameter_set_rbsp_t *sps)
{
	if (active_sps != sps)
	{
		if (dec_picture)
		{
			// this may only happen on slice loss
			exit_picture();
		}
		active_sps = sps;
		
		img -> MaxFrameNum = 1<<(sps->log2_max_frame_num_minus4+4);
		img -> PicWidthInMbs = (sps->pic_width_in_mbs_minus1 +1);
		img -> PicHeightInMapUnits = (sps->pic_height_in_map_units_minus1 +1);
		img -> FrameHeightInMbs = ( 2 - sps->frame_mbs_only_flag ) * img->PicHeightInMapUnits;
		img->FrameSizeInMbs = img->PicWidthInMbs * img->FrameHeightInMbs;
		
		img->width = img->PicWidthInMbs * MB_BLOCK_SIZE;
		img->width_cr = img->width /2;
		img->height = img->FrameHeightInMbs * MB_BLOCK_SIZE;
		img->height_cr = img->height / 2;
		
		init_global_buffers_dec();
		if (!img->no_output_of_prior_pics_flag)
		{
			flush_dpb();
		}
		init_dpb();
		
		if (NULL!=Co_located)
		{
// 			free_collocated(Co_located);
		}
		Co_located = alloc_colocated (img->width, img->height,sps->mb_adaptive_frame_field_flag);
		ercInit(img->width, img->height, 1);
	}
}

void activate_pps(pic_parameter_set_rbsp_t *pps)
{
	if (active_pps != pps)
	{
		if (dec_picture)
		{
			// this may only happen on slice loss
			exit_picture();
		}
		active_pps = pps;
	}
	
}
// 

void UseParameterSet (int PicParsetId)//此函数在解码IDR内有调用。
{
	seq_parameter_set_rbsp_t *sps = &SeqParSet[PicParSet[PicParsetId].seq_parameter_set_id];
	pic_parameter_set_rbsp_t *pps = &PicParSet[PicParsetId];
	int i;
	
	
	if (PicParSet[PicParsetId].Valid != TRUE1)
		printf ("Trying to use an invalid (uninitialized) Picture Parameter Set with ID %d, expect the unexpected...\n", PicParsetId);
	if (SeqParSet[PicParSet[PicParsetId].seq_parameter_set_id].Valid != TRUE1)
		printf ("PicParset %d references an invalid (uninitialized) Sequence Parameter Set with ID %d, expect the unexpected...\n", PicParsetId, PicParSet[PicParsetId].seq_parameter_set_id);
	
	sps =  &SeqParSet[PicParSet[PicParsetId].seq_parameter_set_id];//将SeqParSet中的数据读入active_sps
	
	
	// In theory, and with a well-designed software, the lines above
	// are everything necessary.  In practice, we need to patch many values
	// in img-> (but no more in inp-> -- these have been taken care of)
	
	// Sequence Parameter Set Stuff first
	
	//  printf ("Using Picture Parameter set %d and associated Sequence Parameter Set %d\n", PicParsetId, PicParSet[PicParsetId].seq_parameter_set_id);
	
	if (sps->pic_order_cnt_type < 0 || sps->pic_order_cnt_type > 2)  // != 1
	{
		printf ("invalid sps->pic_order_cnt_type = %d\n", sps->pic_order_cnt_type);
		error ("pic_order_cnt_type != 1", -1000);
	}
	
	if (sps->pic_order_cnt_type == 1)
	{
		if(sps->num_ref_frames_in_pic_order_cnt_cycle >= MAXnum_ref_frames_in_pic_order_cnt_cycle)
		{
			error("num_ref_frames_in_pic_order_cnt_cycle too large",-1011);
		}
	}
	
	activate_sps(sps);//
	activate_pps(pps);
	
	// currSlice->dp_mode is set by read_new_slice (NALU first byte available there)
	if (pps->entropy_coding_mode_flag == UVLC)
	{
		nal_startcode_follows = uvlc_startcode_follows;
		for (i=0; i<3; i++)
		{
			img->currentSlice->partArr[i].readSyntaxElement = readSyntaxElement_UVLC_dec;
		}
	}
	else
	{
		nal_startcode_follows = cabac_startcode_follows;
		for (i=0; i<3; i++)
		{
			img->currentSlice->partArr[i].readSyntaxElement = readSyntaxElement_CABAC;
		}
	}
}






/////////////////////////////////////<