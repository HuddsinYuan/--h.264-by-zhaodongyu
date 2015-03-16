#ifndef _DEFINES_ENC_H_
#define _DEFINES_ENC_H_
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <string.h>
#include <math.h>

#include "i_global.h"

#define INIT_FRAME_RATE 30
#define FILE_NAME_SIZE 200
// #define byte unsigned char 
#define DEBUG 0
#define GREY_LEVELS 255
#define MAX_BITS 10

#define MIN_ALPHA -2.35 //alpha 最小值//不要更改，否则解码出错
#define MAX_ALPHA 4.0
#define MIN_BETA -60    //beta 最小值
#define MAX_BETA 255

#define LENGTH_QUAD_COUNT  16 //32 /* Number of bits used to store quad_count	                                                            */  
#define LENGTH_TRANS_COUNT 16 //32 /* Number of bits used to store trans_count	                                                            */ 
#define LENGTH_CODE_LENGTH 16 //32 /* Number of bits used to store code_length 	                                                            */ 
#define LENGTH_DATA_LENGTH 16 //32 /* Number of bits used to store data_length 	

#define BLOCK_SIZE_I 8 
#define POSITIVE_SCALING 0    /* A flag specifying use of positive scaling                                                                                                                   */
#define bound(a)   ((a) < 0.0 ? 0 : ((a)>255.0? 255 : a))
int imageWidth,imageHeight;
int search_range;
int search_mode;
int num_regions;
int **searchedpt;
byte i_quality;
int b_4_sub;
typedef struct
{
  short mv_x;
  short mv_y;
} MotionVector;

MotionVector UpCur,Lrefpre,Lrefcur,disparity,motion;    //////////////////视差估计中的最佳位置保存，用于下一宏块初始点的确定;
typedef struct transformation_node
{
	unsigned char block_type;    //
	unsigned char use;//OB时是否使用到当前块，当NOB时，该参数忽略
	unsigned char region;//OB时，如果使用到了当前块，该块数据哪个region，NOB时忽略,该值在码流中不保存用于压缩和解压缩过程
	unsigned char reference;          //0:参考块是C目；1：参考块是H目；2：参考块是M目；3：参考块是N目
	unsigned char partition;//块划分模式，0：没有划分；1：16x8；2：8x16；3：8x8；
	int x,y;
	double scale,offset;
    struct transformation_node *next;
}TRANS_NODE;

typedef struct 
{
  char infile_c[FILE_NAME_SIZE];
  char infile_r[FILE_NAME_SIZE];
  char infile_l[FILE_NAME_SIZE];
  char infile_c_plane[FILE_NAME_SIZE];
  char infile_r_plane[FILE_NAME_SIZE];
  char infile_l_plane[FILE_NAME_SIZE];

  char outfile_all_c[FILE_NAME_SIZE];
  char outfile_all_r[FILE_NAME_SIZE];
  char outfile_all_l[FILE_NAME_SIZE];
  char outfile_c[2][FILE_NAME_SIZE];
  char outfile_r[2][FILE_NAME_SIZE];
  char outfile_l[2][FILE_NAME_SIZE];


  char outfile_c_rec[2][FILE_NAME_SIZE];
  char outfile_r_rec[2][FILE_NAME_SIZE];
  char outfile_l_rec[2][FILE_NAME_SIZE];
  char outall_c[FILE_NAME_SIZE];
  char outall_264[FILE_NAME_SIZE];

  char outall_r[FILE_NAME_SIZE];
  char outall_l[FILE_NAME_SIZE];

  int I_frame;
  byte i_quality;
	
  int imagewidth;
  int imageheight;
  int no_frames;                              // number of frames to encode 需要编码的帧数
  int right;                                 
  int left;
  int num_regions;                            //* Sets Region-based (Set as 2) 基于对象编码的对象数目
  double tol_16;
  double tol_8;
  double tol_4;
  int search_range;
  int yuv_format;
  int frame_rate;
  int IntraPeriod;

  int bitDepthLuma;
  int displayEncoderParams;
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
  int ProfileIDC;               //!< profile idc
  int LevelIDC;                 //!< level idc

  int no_frames_h264;                //!< number of frames to be encoded
  int qp0;                      //!< QP of first frame         第一帧的QP
  int qpN;                      //!< QP of remaining frames    其余帧的QP
  int jumpd;                    //!< number of frames to skip in input sequence (e.g 2 takes frame 0,3,6,9...)
  int hadamard;                 /*!< 0: 'normal' SAD in 1/3 pixel search.  1: use 4x4 Haphazard transform and '
                                     Sum of absolute transform difference' in 1/3 pixel search                   */
//   int search_range;             /*!< search range - integer pel search and 16x16 blocks.  The search window is
//                                      generally around the predicted vector. Max vector is 2xmcrange.  For 8x8
//                                      and 4x4 block sizes the search range is 1/2 of that for 16x16 blocks.       */
  int num_reference_frames;     //!< number of reference frames to be used参考帧的数目
  int P_List0_refs;
  int B_List0_refs;
  int B_List1_refs;
//用分形的imagewidth,imageheight
//   int img_width;                //!< image width  (must be a multiple of 16 pels)
//   int img_height;               //!< image height (must be a multiple of 16 pels)
//   int yuv_format;               //!< GH: YUV format (0=4:0:0, 1=4:2:0, 2=4:2:2, 3=4:4:4,currently only 4:2:0 is supported)
  int color_depth;              //!< GH: YUV color depth per component in bit/pel (currently only 8 bit/pel is supported)
  int intra_upd;                /*!< For error robustness. 0: no special action. 1: One GOB/frame is intra coded
                                     as regular 'update'. 2: One GOB every 2 frames is intra coded etc.
                                     In connection with this intra update, restrictions is put on motion vectors
                                     to prevent errors to propagate from the past                                */
  int blc_size[8][2];           //!< array for different block sizes
  int slice_mode;               //!< Indicate what algorithm to use for setting slices
  int slice_argument;           //!< Argument to the specified slice algorithm
  int UseConstrainedIntraPred;  //!< 0: Inter MB pixels are allowed for intra prediction 1: Not allowed
  int  infile_header;           //!< If input file has a header set this to the length of the header
  char infile[100];             //!< YUV 4:2:0 input format
  char outfile[100];            //!< H.264 compressed output bitstream
  char ReconFile[100];          //!< Reconstructed Pictures
  char TraceFile[100];          //!< Trace Outputs
  int intra_period;             //!< Random Access period though intra

  int idr_enable;				//!< Encode intra slices as IDR
  int start_frame;				//!< Encode sequence starting from Frame start_frame

  // B pictures
  int successive_Bframe;        //!< number of B frames that will be used
  int qpB;                      //!< QP of B frames
  int direct_type;              //!< Direct Mode type to be used (0: Temporal, 1: Spatial)
  int directInferenceFlag;      //!< Direct Inference Flag

  // SP Pictures
  int sp_periodicity;           //!< The periodicity of SP-pictures
  int qpsp;                     //!< SP Picture QP for prediction error
  int qpsp_pred;                //!< SP Picture QP for predicted block

  int WeightedPrediction;        //!< Weighted prediciton for P frames (0: not used, 1: explicit)
  int WeightedBiprediction;      //!< Weighted prediciton for B frames (0: not used, 1: explicit, 2: implicit)
  int StoredBPictures;           //!< Stored (Reference) B pictures replace P pictures (0: not used, 1: used)

  int symbol_mode;              //!< Specifies the mode the symbols are mapped on bits
  int of_mode;                  //!< Specifies the mode of the output file
  int partition_mode;           //!< Specifies the mode of data partitioning

  int InterSearch16x16;
  int InterSearch16x8;
  int InterSearch8x16;
  int InterSearch8x8;
  int InterSearch8x4;
  int InterSearch4x8;
  int InterSearch4x4;

//   char PictureTypeSequence[MAXPICTURETYPESEQUENCELEN];
  int FrameRate;

  int chroma_qp_index_offset;
#ifdef _FULL_SEARCH_RANGE_
  int full_search;
#endif
#ifdef _ADAPT_LAST_GROUP_
  int last_frame;
#endif
#ifdef _CHANGE_QP_
  int qpN2, qpB2, qp2start;
  int qp02;
#endif
  int rdopt;
#ifdef _LEAKYBUCKET_
  int  NumberLeakyBuckets;
  char LeakyBucketRateFile[100];
  char LeakyBucketParamFile[100];
#endif

  int PicInterlace;           //!< picture adaptive frame/field
  int MbInterlace;            //!< macroblock adaptive frame/field

  int IntraBottom;            //!< Force Intra Bottom at GOP periods.

  int LossRateA;              //!< assumed loss probablility of partition A (or full slice), in per cent, used for loss-aware R/D optimization
  int LossRateB;              //!< assumed loss probablility of partition B, in per cent, used for loss-aware R/D 
  int LossRateC;              //!< assumed loss probablility of partition C, in per cent, used for loss-aware R/D 
  int NoOfDecoders;
  int RestrictRef;
  int NumFramesInELSubSeq;
  int NumFrameIn2ndIGOP;

  int RandomIntraMBRefresh;     //!< Number of pseudo-random intra-MBs per picture

  int LFSendParameters;
  int LFDisableIdc;
  int LFAlphaC0Offset;
  int LFBetaOffset;

  int SparePictureOption;
  int SPDetectionThreshold;
  int SPPercentageThreshold;

  // FMO
  char SliceGroupConfigFileName[100];    //!< Filename for config info fot type 0, 2, 6	
  int num_slice_groups_minus1;           //!< "FmoNumSliceGroups" in encoder.cfg, same as FmoNumSliceGroups, which should be erased later
  int slice_group_map_type; 

  int *top_left;                         //!< top_left and bottom_right store values indicating foregrounds
  int *bottom_right; 
  int *slice_group_id;                   //!< slice_group_id is for slice group type being 6  
  int *run_length_minus1;                //!< run_length_minus1 is for slice group type being 0

  int slice_group_change_direction_flag;
  int slice_group_change_rate_minus1;
  int slice_group_change_cycle;

  int redundant_slice_flag; //! whether redundant slices exist,  JVT-D101
  int pic_order_cnt_type;   // POC200301

  int context_init_method;
  int model_number;

  //! Rate Control on JVT standard 
  int RCEnable;    
  int bit_rate;
  int SeinitialQP;
  int basicunit;
  int channel_type;

  // FastME enable
  int FMEnable;

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
}InputParameters;

typedef struct  img_par
{
	int current_frame;
    int frmWidthInMbs;
    int frmHeightInMbs;
    int frmSizeInMbs;
/*	Macroblock  *mb_data;*/

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////

//   int number;                  //!< current image number to be encoded
  int pn;                      //!< picture number
  int nb_references;
  int current_mb_nr;

/////////////////////////////////////>

  int allrefzero;
//   int mpr[16][16];                            //!< predicted block
//   int mvscale[6][MAX_REFERENCE_PICTURES];
//   int m7[16][16];                             //!< final 4x4 block. Extended to 16x16 for ABT
  int cof[4][6][4][4];                        //!< correction coefficients from predicted
  int cofu[4];
//   int **ipredmode;                            //!< prediction type
//   int quad[256];
//   int ***nz_coeff;
  int **siblock;
//   int cod_counter;                            //!< Current count of number of skipped macroblocks in a row

  
  unsigned oldFrameSizeInMbs;

  int number;                                 //!< frame number

  unsigned num_dec_mb;
  int newframe;

//   int structure;                               //!< Identify picture structure type

  int idr_psnr_number;
  int psnr_number;

  time_t ltime_start;               // for time measurement
  time_t ltime_end;                 // for time measurement

#ifdef WIN32
  struct _timeb tstruct_start;
  struct _timeb tstruct_end;
#else
  struct timeb tstruct_start;
  struct timeb tstruct_end;
#endif
  
//   int **siblock;
//   int allrefzero;

/////////////////////////////////////<

  int total_number_mb;
  int current_slice_nr;
  int type;
  int structure;               //!< picture structure
  int num_reference_frames;    //!< number of reference frames to be used 使用的参考帧数目
  int max_num_references;      //!< maximum number of reference pictures that may occur
  int qp;                      //!< quant for the current frame           当前帧的QP
  int qpsp;                    //!< quant for the prediction frame of SP-frame
  int framerate;               //帧率
  int width;                   //!< Number of pels                       亮度宽度
  int width_cr;                //!< Number of pels chroma                色度宽度
  int height;                  //!< Number of lines                      亮度高度
  int height_cr;               //!< Number of lines  chroma              色度高度
  int subblock_x;              //!< current subblock horizontal          当前子块的水平坐标
  int subblock_y;              //!< current subblock vertical            当前子块的垂直坐标
  int is_intra_block;
  int is_v_block;
  int mb_y_upd;
  int mb_y_intra;              //!< which GOB to intra code
  int block_c_x;               //!< current block chroma vertical
  int **ipredmode;             //!< intra prediction mode              帧内预测模式
  int cod_counter;             //!< Current count of number of skipped macroblocks in a row
  int ***nz_coeff;             //!< number of coefficients per block (CAVLC) 每个块的系数

  int mb_x;                    //!< current MB horizontal       当前宏块水平坐标
  int mb_y;                    //!< current MB vertical         当前宏块垂直坐标
  int block_x;                 //!< current block horizontal    当前块水平坐标
  int block_y;                 //!< current block vertical      当前块垂直坐标
  int pix_x;                   //!< current pixel horizontal    当前像素水平坐标
  int pix_y;                   //!< current pixel vertical      当前像素垂直坐标
  int pix_c_x;                 //!< current pixel chroma horizontal    当前像素色度水平坐标
  int pix_c_y;                 //!< current pixel chroma vertical      当前像素色度垂直坐标

  int opix_x;                   //!< current original picture pixel horizontal   
  int opix_y;                   //!< current original picture pixel vertical
  int opix_c_x;                 //!< current original picture pixel chroma horizontal
  int opix_c_y;                 //!< current original picture pixel chroma vertical


  // some temporal buffers
  int mprr[9][16][16];         //!< all 9 prediction modes? // enlarged from 4 to 16 for ABT (is that neccessary?)

  int mprr_2[5][16][16];       //!< all 4 new intra prediction modes    
  int mprr_c[2][4][8][8];      //!< new chroma 8x8 intra prediction modes 色度8x8帧内预测模式
  int mpr[16][16];             //!< current best prediction mode          当前最佳的预测模式
  int m7[16][16];              //!< the diff pixel values between orginal image and prediction 原始图像与预测图像像素差值


//   int pre_16_fract[640][352];
  int ****cofAC;               //!< AC coefficients [8x8block][4x4block][level/run][scan_pos]  AC系数
  int ***cofDC;                //!< DC coefficients [yuv][level/run][scan_pos]                 DC系数

  Picture     *currentPicture; //!< The coded picture currently in the works (typically frame_pic, top_pic, or bottom_pic)
  Slice       *currentSlice;                                //!< pointer to current Slice data struct 指向当前条带数据结构的指针
  
  
/////////////////////////////////////>
  int sp_switch;                              //!< 1 for switching sp, 0 for normal sp

  Macroblock          *mb_data;                //!< array containing all MBs of a whole frame
  unsigned int field_pic_flag;
  unsigned int bottom_field_flag;
  int idr_flag;
  int idr_pic_id;

   // for signalling to the neighbour logic that this is a deblocker call
//    int DeblockCall;
   
  int **field_anchor;
   
//    DecRefPicMarking_t *dec_ref_pic_marking_buffer;                    //!< stores the memory management control operations
   
  int disposable_flag;                          //!< flag for disposable frame, 1:disposable
//    int num_ref_idx_l0_active;             //!< number of forward reference
//    int num_ref_idx_l1_active;             //!< number of backward reference
   
//    int slice_group_change_cycle;

   // for POC mode 0:
  signed int PrevPicOrderCntMsb;
  unsigned int PrevPicOrderCntLsb;
  signed int PicOrderCntMsb;

   // for POC mode 1:
  unsigned int AbsFrameNum;
  signed int ExpectedPicOrderCnt, PicOrderCntCycleCnt, FrameNumInPicOrderCntCycle;
  unsigned int PreviousFrameNum, FrameNumOffset;
  int ExpectedDeltaPerPicOrderCntCycle;
  int PreviousPOC, ThisPOC;
  int PreviousFrameNumOffset;


   //weighted prediction
  unsigned int luma_log2_weight_denom;
  unsigned int chroma_log2_weight_denom;
  int ***wp_weight;  // weight in [list][index][component] order
  int ***wp_offset;  // offset in [list][index][component] order
  int ****wbp_weight; //weight in [list][fw_index][bw_index][component] order
  int wp_round_luma;
  int wp_round_chroma;
  unsigned int apply_weights;

   int MaxFrameNum;
/////////////////////////////////////<
  
//   Macroblock    *mb_data;                                   //!< array containing all MBs of a whole frame
  SyntaxElement   MB_SyntaxElements[MAX_SYMBOLS_PER_MB];    //!< temporal storage for all chosen syntax elements of one MB

/////////////////////////////////////>
  
  
  int quad_dec[256];

  
  
/////////////////////////////////////<

  int *quad;               //!< Array containing square values,used for snr computation  */                                         /* Values are limited to 5000 for pixel differences over 70 (sqr(5000)).
  int *intra_block;

  int tr;
  int fld_type;                        //!< top or bottom field
  unsigned int fld_flag;                                
  int direct_intraP_ref[4][4];
  int pstruct_next_P;
  int imgtr_next_P_frm;
  int imgtr_last_P_frm;
  int imgtr_next_P_fld;
  int imgtr_last_P_fld;

  // B pictures
  int b_interval;
  int p_interval;
  int b_frame_to_code;
  int fw_mb_mode;
  int bw_mb_mode;

  int****** pred_mv;                 //!< motion vector predictors for all block types and all reference frames

  int****** all_mv;       //!< replaces local all_mv
  int LFDisableIdc;
  int LFAlphaC0Offset;
  int LFBetaOffset;

  int direct_type;              //!< Direct Mode type to be used (0: Temporal, 1: Spatial)

  int num_ref_idx_l0_active;
  int num_ref_idx_l1_active;

  int field_mode;     //!< For MB level field/frame -- field mode on flag
  int top_field;      //!< For MB level field/frame -- top field flag
  int mvscale[6][MAX_REFERENCE_PICTURES];
  int buf_cycle;
  int i16offset;

  int layer;             //!< which layer this picture belonged to
  int old_layer;         //!< old layer number
  int NoResidueDirect;

  int redundant_pic_cnt; // JVT-D101

  int MbaffFrameFlag;    //!< indicates frame with mb aff coding

  //the following should probably go in sequence parameters
  // unsigned int log2_max_frame_num_minus4;
  unsigned int pic_order_cnt_type;
  // for poc mode 0, POC200301
  // unsigned int log2_max_pic_order_cnt_lsb_minus4;  
  // for poc mode 1, POC200301
  unsigned int delta_pic_order_always_zero_flag;
           int offset_for_non_ref_pic;
           int offset_for_top_to_bottom_field;
  unsigned int num_ref_frames_in_pic_order_cnt_cycle;
           int offset_for_ref_frame[1];  // MAX_LENGTH_POC_CYCLE in decoder

  // POC200301
  //the following is for slice header syntax elements of poc
  // for poc mode 0.
  unsigned int pic_order_cnt_lsb;
           int delta_pic_order_cnt_bottom;
  // for poc mode 1.
           int delta_pic_order_cnt[2];


  // POC200301
  unsigned int field_picture;
    signed int toppoc;      //!< poc for this frame or field
    signed int bottompoc;   //!< for completeness - poc of bottom field of a frame (always = poc+1)
    signed int framepoc;    //!< min (toppoc, bottompoc)
//     signed int ThisPOC;     //!< current picture POC
  unsigned int frame_num;   //!< frame_num for this frame
  
  unsigned PicWidthInMbs;
  unsigned PicHeightInMapUnits;
  unsigned FrameHeightInMbs;
  unsigned PicHeightInMbs;
  unsigned PicSizeInMbs;
  unsigned FrameSizeInMbs;

  //the following should probably go in picture parameters
  unsigned int pic_order_present_flag; // ????????

  //the following are sent in the slice header
//  int delta_pic_order_cnt[2];
  int nal_reference_idc;

  int adaptive_ref_pic_buffering_flag;
  int no_output_of_prior_pics_flag;
  int long_term_reference_flag;

  DecRefPicMarking_t *dec_ref_pic_marking_buffer;

  int model_number;


  /*rate control*/
  int NumberofHeaderBits; 
  int NumberofTextureBits;
  int NumberofBasicUnitHeaderBits;
  int NumberofBasicUnitTextureBits;
  double TotalMADBasicUnit;
  int NumberofMBTextureBits;
  int NumberofMBHeaderBits;
  int NumberofCodedBFrame; 
  int NumberofCodedPFrame;
  int NumberofGOP;
  int TotalQpforPPicture;
  int NumberofPPicture;
  double MADofMB[6336];
  int BasicUnitQP;
  int TopFieldFlag;
  int FieldControl;
  int FieldFrame;
  int Frame_Total_Number_MB;
  int IFLAG;
  int NumberofCodedMacroBlocks;
  int BasicUnit;
  int write_macroblock;
  int bot_MB;
  int write_macroblock_frame;

  int DeblockCall;
        
  int last_pic_bottom_field;
  int last_has_mmco_5;
  int pre_frame_num;

  int slice_group_change_cycle;

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////

} ImageParameters;                       

extern InputParameters *input;
extern ImageParameters *img;

extern struct snr_par  *snr;

typedef struct  
{
	int frames;
	short Width;
	short Height;
    byte  stereo;    //几目压缩
	byte  regions;
	byte i_quality; //I帧质量
	byte i_interval;//I帧个数
	byte search_range;
	byte zero;
}compressionInfo;//8 bytes

typedef struct  
{
	byte *sourceframe_Y;
	byte *sourceframe_U;
	byte *sourceframe_V;
	byte *plane_Y;
	byte *plane_U;
	byte *plane_V;
	byte *plane_Y_domain;
	byte *plane_U_domain;
	byte *plane_V_domain;
}SourceFrame;


#define QUAN_A(x)\
{\
	int b,c;\
	b=(int)(x)%10;c=(int)(x)/10;\
	if(b>2&&b<8)\
	{b=5;}\
	else if(b>7)\
	{b=0;c+=1;}\
	else b=0;\
	x=(int)(c*10)+b;\
}//对x量化,个位数量化为0或5，十位数没变

#define QUAN_B(x)\
{\
	int b,c;\
	b=(int)(x)%10;c=(int)(x)/10;\
	if(b>=5)\
	{c++;}\
	x=(int)(c*10);\
}
static int imax(int a, int b)
{
  return ((a) > (b)) ? (a) : (b);
}
#endif




















