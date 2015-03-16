#ifndef _CONFIGFILE_H_
#define _CONFIGFILE_H_
#include "global.h"
#include "i_global.h"

#define DEFAULTCONFIGFILENAME "encoder.cfg"

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
#define PROFILE_IDC     88
#define LEVEL_IDC       21

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////



typedef struct{
	char *TokenName;
	void *Place;
	int Type;
	double Default;
	int param_limits; //! 0: no limits, 1: both min and max, 2: only min (i.e. no negatives), 3: only two value
	double min_limit;
	double max_limit;
}Mapping;

InputParameters configinput;

#ifdef INCLUDED_BY_CONFIGFILE_C
Mapping Map[] = {
	{"Infile_c",               &configinput.infile_c,               1,     0.0,         0,       0.0,       0.0},
	{"Infile_r",               &configinput.infile_r,               1,     0.0,         0,       0.0,       0.0},
	{"Infile_l",               &configinput.infile_l,               1,     0.0,         0,       0.0,       0.0},
	{"Infile_c_plane",         &configinput.infile_c_plane,         1,     0.0,         0,       0.0,       0.0},
	{"Infile_r_plane",         &configinput.infile_r_plane,         1,     0.0,         0,       0.0,       0.0},
	{"Infile_l_plane",         &configinput.infile_l_plane,         1,     0.0,         0,       0.0,       0.0},

	{"Outfile_all_c",          &configinput.outfile_all_c,          1,     0.0,         0,       0.0,       0.0},
	{"Outfile_all_r",          &configinput.outfile_all_r,          1,     0.0,         0,       0.0,       0.0},
	{"Outfile_all_l",          &configinput.outfile_all_l,          1,     0.0,         0,       0.0,       0.0},
	{"Outfile_c_0",            &configinput.outfile_c[0],           1,     0.0,         0,       0.0,       0.0},
	{"Outfile_c_1",            &configinput.outfile_c[1],           1,     0.0,         0,       0.0,       0.0},
	{"Outfile_r_0",            &configinput.outfile_r[0],           1,     0.0,         0,       0.0,       0.0},
	{"Outfile_r_1",            &configinput.outfile_r[1],           1,     0.0,         0,       0.0,       0.0},
	{"Outfile_l_0",            &configinput.outfile_l[0],           1,     0.0,         0,       0.0,       0.0},
	{"Outfile_l_1",            &configinput.outfile_l[1],           1,     0.0,         0,       0.0,       0.0},

	{"Outfile_c_rec_0",        &configinput.outfile_c_rec[0],       1,     0.0,         0,       0.0,       0.0},
	{"Outfile_c_rec_1",        &configinput.outfile_c_rec[1],       1,     0.0,         0,       0.0,       0.0},
	{"Outfile_r_rec_0",        &configinput.outfile_r_rec[0],       1,     0.0,         0,       0.0,       0.0},
	{"Outfile_r_rec_1",        &configinput.outfile_r_rec[1],       1,     0.0,         0,       0.0,       0.0},
	{"Outfile_l_rec_0",        &configinput.outfile_l_rec[0],       1,     0.0,         0,       0.0,       0.0},
	{"Outfile_l_rec_1",        &configinput.outfile_l_rec[1],       1,     0.0,         0,       0.0,       0.0},
	{"Outall_c",               &configinput.outall_c,               1,     0.0,         0,       0.0,       0.0},
	{"Outall_264",             &configinput.outall_264,               1,     0.0,         0,       0.0,       0.0},
	{"Outall_r",               &configinput.outall_r,               1,     0.0,         0,       0.0,       0.0},
	{"Outall_l",               &configinput.outall_l,               1,     0.0,         0,       0.0,       0.0},
	{"IntraPeriod",            &configinput.IntraPeriod,            0,     0.0,         0,       0.0,       0.0},
    {"NumberBFrames",          &configinput.successive_Bframe,      0,     0.0,         0,       0.0,       0.0},
    {"DirectModeType",         &configinput.direct_type,            0,     0.0,         0,       0.0,       0.0},
    {"NumberReferenceFrames",  &configinput.num_reference_frames,   0,     0.0,         0,       0.0,       0.0},
    {"SymbolMode",             &configinput.symbol_mode,            0,     0.0,         0,       0.0,       0.0},
    {"OutFileMode",            &configinput.of_mode,                0,     0.0,         0,       0.0,       0.0},
    {"PartitionMode",          &configinput.partition_mode,         0,     0.0,         0,       0.0,       0.0},
	{"SliceMode",              &configinput.slice_mode,             0,     0.0,         0,       0.0,       0.0},
    {"SliceArgument",          &configinput.slice_argument,         0,     0.0,         0,       0.0,       0.0},


	{"QPFirstFrame",             &configinput.qp0,                     0,     0.0,         0,       0.0,       0.0},
    {"QPRemainingFrame",         &configinput.qpN,                     0,     0.0,         0,       0.0,       0.0},
    {"FrameSkip",                &configinput.jumpd,                   0,     0.0,         0,       0.0,       0.0},
    {"ChromaQPOffset",           &configinput.chroma_qp_index_offset,  0,     0.0,         0,       0.0,       0.0},    
    {"RandomIntraMBRefresh",     &configinput.RandomIntraMBRefresh,    0,     0.0,         0,       0.0,       0.0},
	{"InterSearch16x16",         &configinput.InterSearch16x16,        0,     0.0,         0,       0.0,       0.0},
    {"InterSearch16x8",          &configinput.InterSearch16x8 ,        0,     0.0,         0,       0.0,       0.0},
    {"InterSearch8x16",          &configinput.InterSearch8x16,         0,     0.0,         0,       0.0,       0.0},
    {"InterSearch8x8",           &configinput.InterSearch8x8 ,         0,     0.0,         0,       0.0,       0.0},
    {"InterSearch8x4",           &configinput.InterSearch8x4,          0,     0.0,         0,       0.0,       0.0},
    {"InterSearch4x8",           &configinput.InterSearch4x8,          0,     0.0,         0,       0.0,       0.0},
    {"InterSearch4x4",           &configinput.InterSearch4x4,          0,     0.0,         0,       0.0,       0.0},

    {"ProfileIDC",               &configinput.ProfileIDC,              0,     0.0,         0,       0.0,       0.0},
    {"LevelIDC",                 &configinput.LevelIDC,                0,     0.0,         0,       0.0,       0.0},

    {"OutputFile",               &configinput.outfile,                 1,     0.0,         0,       0.0,       0.0},
	{"ReconFile",                &configinput.ReconFile,               1,     0.0,         0,       0.0,       0.0},
    {"TraceFile",                &configinput.TraceFile,               1,     0.0,         0,       0.0,       0.0},
	{"RDOptimization",           &configinput.rdopt,                   0,     0.0,         0,       0.0,       0.0},
    {"FramesToBeEncoded",        &configinput.no_frames_h264,          0,     0.0,         0,       0.0,       0.0},


	{"I_Frame",                &configinput.I_frame,                0,     0,           2,       0,         0},
	{"I_Quality",              &configinput.i_quality,              0,     50,          1,       0,         100},
	{"ImageWidth",             &configinput.imagewidth,             0,     176.0,       1,       88.0,      2000.0},
	{"ImageHeight",            &configinput.imageheight,            0,     144.0,       1,       72.0,      2000.0},
	{"NoFrames",               &configinput.no_frames,              0,     100.0,       2,       1.0,       0.0},

	{"Right",                  &configinput.right,                  0,     0.0,         1,       0.0,       2.0},
	{"Left",                   &configinput.left,                   0,     0.0,         1,       0.0,       2.0},
	{"Num_Regions",            &configinput.num_regions,            0,     1.0,         1,       1.0,       3.0},
	{"Tol_16",                 &configinput.tol_16,                 2,     8.0,         1,       1.0,       20.0},
	{"Tol_8",                  &configinput.tol_8,                  2,     5.0,         1,       1.0,       10.0},
	{"Tol_4",                  &configinput.tol_4,                  2,     3.0,         1,       1.0,       10.0},
	{"Search_Range",           &configinput.search_range,           0,     8.0,         1,       3.0,       100.0},
    {"YUV_Format",             &configinput.yuv_format,             0,     0.0,         1,       0.0,       3.0},
	{"FrameRate",              &configinput.FrameRate,             0,     20.0,        1,       0.0,      100.0},
	{"DisplayEncoderParams",   &configinput.displayEncoderParams,   0,     0.0,         1,       0.0,       1.0},
	{"BitDepthLuma",           &configinput.bitDepthLuma,           0,     8.0,         1,       8.0,       12.0},



    {"ContextInitMethod",        &configinput.context_init_method,     0},
    {"FixedModelNumber",         &configinput.model_number,            0},

	{"num_slice_groups_minus1",           &configinput.num_slice_groups_minus1,           0,     0.0,         0,       0.0,       0.0},
    {"slice_group_map_type",              &configinput.slice_group_map_type,              0,     0.0,         0,       0.0,       0.0},
    {"slice_group_change_direction_flag", &configinput.slice_group_change_direction_flag, 0,     0.0,         0,       0.0,       0.0},
    {"slice_group_change_rate_minus1",    &configinput.slice_group_change_rate_minus1,    0,     0.0,         0,       0.0,       0.0},
    {"SliceGroupConfigFileName",          &configinput.SliceGroupConfigFileName,          1,     0.0,         0,       0.0,       0.0},
	{"UseRedundantSlice",        &configinput.redundant_slice_flag,    0},


#ifdef _FULL_SEARCH_RANGE_
    {"RestrictSearchRange",      &configinput.full_search,             0},
#endif
    {"RDOptimization",           &configinput.rdopt,                   0},
    {"LossRateA",                &configinput.LossRateA,               0},
    {"LossRateB",                &configinput.LossRateB,               0},
    {"LossRateC",                &configinput.LossRateC,               0},
    {"NumberOfDecoders",         &configinput.NoOfDecoders,            0},
    {"RestrictRefFrames",        &configinput.RestrictRef ,            0},
	
	{"UseConstrainedIntraPred",  &configinput.UseConstrainedIntraPred, 0,     0.0,         0,       0.0,       0.0},
#ifdef _ADAPT_LAST_GROUP_
    {"LastFrameNumber",          &configinput.last_frame,              0,     0.0,         0,       0.0,       0.0},
#endif
#ifdef _CHANGE_QP_
    {"ChangeQPI",                &configinput.qp02,                    0,     0.0,         0,       0.0,       0.0},
    {"ChangeQPP",                &configinput.qpN2,                    0,     0.0,         0,       0.0,       0.0},
    {"ChangeQPB",                &configinput.qpB2,                    0,     0.0,         0,       0.0,       0.0},
    {"ChangeQPStart",            &configinput.qp2start,                0,     0.0,         0,       0.0,       0.0},
#endif
#ifdef _LEAKYBUCKET_
    {"NumberofLeakyBuckets",     &configinput.NumberLeakyBuckets,      0,     0.0,         0,       0.0,       0.0},
    {"LeakyBucketRateFile",      &configinput.LeakyBucketRateFile,     1,     0.0,         0,       0.0,       0.0},
    {"LeakyBucketParamFile",     &configinput.LeakyBucketParamFile,    1,     0.0,         0,       0.0,       0.0},
#endif
    {"NumberFramesInEnhancementLayerSubSequence", &configinput.NumFramesInELSubSeq, 0,     0.0,         0,       0.0,       0.0},
    {"NumberOfFrameInSecondIGOP",&configinput.NumFrameIn2ndIGOP, 0},
    {"SparePictureOption",       &configinput.SparePictureOption,      0,     0.0,         0,       0.0,       0.0},
    {"SparePictureDetectionThr", &configinput.SPDetectionThreshold,    0,     0.0,         0,       0.0,       0.0},
    {"SparePicturePercentageThr",&configinput.SPPercentageThreshold,   0,     0.0,         0,       0.0,       0.0},
    {"PicOrderCntType",          &configinput.pic_order_cnt_type,      0,     0.0,         0,       0.0,       0.0},

    // Rate Control
    {"RateControlEnable",        &configinput.RCEnable,             0,     0.0,         0,       0.0,       0.0},
    {"Bitrate",                  &configinput.bit_rate,             0,     0.0,         0,       0.0,       0.0},
    {"InitialQP",                &configinput.SeinitialQP,          0,     0.0,         0,       0.0,       0.0},
    {"BasicUnit",                &configinput.basicunit,            0,     0.0,         0,       0.0,       0.0},
    {"ChannelType",              &configinput.channel_type,         0,     0.0,         0,       0.0,       0.0},


	{NULL,                     NULL,                               -1,     0.0,         0,       0.0,       0.0}
};
#endif


#ifndef INCLUDED_BY_CONFIGFILE_C
extern Mapping Map[];
#endif


void PatchInputNoFrames();
void Configure_h264 (int ac, char *av[]);

#endif