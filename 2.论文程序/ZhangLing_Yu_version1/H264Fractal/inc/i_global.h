
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////

/*!
 ************************************************************************
 *  \file
 *     global.h
 *  \brief
 *     global definitions for for H.264 decoder.
 ************************************************************************
 */

#ifndef I_GLOBAL_H_
#define I_GLOBAL_H_

#include "StdAfx.h"
#include <stdio.h>                              //!< for FILE
#include <time.h>
#include <sys/timeb.h>
#include "i_defines.h"
#include "parsetcommon.h"          

#include "context_ini.h"

#ifdef WIN32
   #define  snprintf _snprintf
#endif

typedef unsigned char   byte;       // "i_defines.h"已经定义了      !<  8 bit unsigned
typedef int             int32;
typedef unsigned int    u_int32;

#ifdef WIN32
  typedef __int64   int64;
# define INT64_MIN        (-9223372036854775807i64 - 1i64)
#else
  typedef long long int64;
# define INT64_MIN        (-9223372036854775807LL - 1LL)
#endif



/***********************************************************************
 * T y p e    d e f i n i t i o n s    f o r    T M L
 ***********************************************************************
 */

#define pel_t byte

//! Data Partitioning Modes
typedef enum
{
  PAR_DP_1,    //!< no data partitioning is supported
  PAR_DP_3,    //!< data partitioning with 3 partitions
} PAR_DP_TYPE;


//! Output File Types
typedef enum
{
  PAR_OF_ANNEXB,   //!< Current TML description
  PAR_OF_RTP,   //!< RTP Packet Output format
//  PAR_OF_IFF    //!< Interim File Format
} PAR_OF_TYPE;

//! Boolean Type原来就这样
/*typedef enum {
  FALSE,
  TRUE
} Boolean;
*/
typedef enum {
	FRAME_CODING,
	FIELD_CODING,
	ADAPTIVE_CODING
} CodingType;

//! definition of H.264 syntax elements//语法元素类型的定义
typedef enum {
  SE_HEADER,//头
  SE_PTYPE,//
  SE_MBTYPE,//宏块类型
  SE_REFFRAME,//参考帧
  SE_INTRAPREDMODE,//帧内预测模式
  SE_MVD,//MVD
  SE_CBP_INTRA,//帧内的CBP
  SE_LUM_DC_INTRA,
  SE_CHR_DC_INTRA,
  SE_LUM_AC_INTRA,
  SE_CHR_AC_INTRA,
  SE_CBP_INTER,//帧间的CBP
  SE_LUM_DC_INTER,
  SE_CHR_DC_INTER,
  SE_LUM_AC_INTER,
  SE_CHR_AC_INTER,
  SE_DELTA_QUANT_INTER,
  SE_DELTA_QUANT_INTRA,
  SE_BFRAME,//B帧
  SE_EOS,//结束符
  SE_MAX_ELEMENTS //!< number of maximum syntax elements, this MUST be the last one!
} SE_type;        // substituting the definitions in element.h


typedef enum {
  INTER_MB,
  INTRA_MB_4x4,
  INTRA_MB_16x16
} IntraInterDecision;

typedef enum {
  BITS_HEADER,//0
  BITS_TOTAL_MB,//1
  BITS_MB_MODE,//2
  BITS_INTER_MB,//3
  BITS_CBP_MB,//4
  BITS_COEFF_Y_MB,//5
  BITS_COEFF_UV_MB,//6
  BITS_DELTA_QUANT_MB,//7
  MAX_BITCOUNTER_MB//8
} BitCountType;

typedef enum {
  NO_SLICES,
  FIXED_MB,
  FIXED_RATE,
  CALLBACK1,
  FMO
} SliceMode;


typedef enum {
  UVLC,//0
  CABAC//1
} SymbolMode;//熵编码模式有两个，UVLC或CABAC

typedef enum {
 LIST_0=0,
 LIST_1=1
} Lists;


typedef enum {
  FRAME,
  TOP_FIELD,
  BOTTOM_FIELD
} PictureStructure;           //!< New enum for field processing


typedef enum {
  P_SLICE = 0,
  B_SLICE,
  I_SLICE,
  SP_SLICE,
  SI_SLICE
} SliceType;

/***********************************************************************
 * D a t a    t y p e s   f o r  C A B A C
 ***********************************************************************
 */

//! struct to characterize the state of the arithmetic coding engine
typedef struct
{
  unsigned int  Elow, Erange;
  unsigned int  Ebuffer;
  unsigned int  Ebits_to_go;
  unsigned int  Ebits_to_follow;
  byte          *Ecodestrm;
  int           *Ecodestrm_len;
//  int           *Ecodestrm_laststartcode;
  // storage in case of recode MB
  unsigned int  ElowS, ErangeS;
  unsigned int  EbufferS;
  unsigned int  Ebits_to_goS;
  unsigned int  Ebits_to_followS;
  byte          *EcodestrmS;
  int           *Ecodestrm_lenS;
  int           C, CS;
  int           E, ES;
  int           B, BS;
} EncodingEnvironment;

typedef EncodingEnvironment *EncodingEnvironmentPtr;

/////////////////////////////////////>


struct img_par;
// struct inp_par;
// struct stat_par;

//! struct to characterize the state of the arithmetic coding engine
typedef struct
{
	unsigned int    Dlow, Drange;
	unsigned int    Dvalue;
	unsigned int    Dbuffer;
	int             Dbits_to_go;
	byte            *Dcodestrm;
	int             *Dcodestrm_len;
} DecodingEnvironment;

typedef DecodingEnvironment *DecodingEnvironmentPtr;




// input parameters from configuration file
struct inp_par_dec
{
	char infile[100];                       //!< H.264 inputfile
	char outfile[100];                      //!< Decoded YUV 4:2:0 output
	char reffile[100];                      //!< Optional YUV 4:2:0 reference file for SNR measurement
	int FileFormat;                         //!< File format of the Input file, PAR_OF_ANNEXB or PAR_OF_RTP
	int dpb_size;                          //!< Frame buffer size
	int ref_offset;
	int poc_scale;
	
#ifdef _LEAKYBUCKET_
	unsigned long R_decoder;                //!< Decoder Rate in HRD Model
	unsigned long B_decoder;                //!< Decoder Buffer size in HRD model
	unsigned long F_decoder;                //!< Decoder Inital buffer fullness in HRD model
	char LeakyBucketParamFile[100];         //!< LeakyBucketParamFile
#endif
	
};

extern struct inp_par_dec *input_dec;



/////////////////////////////////////<

//! struct for context management
typedef struct
{
	unsigned short state;         // index into state-table CP  
	unsigned char  MPS;           // Least Probable Symbol 0/1 CP
	
	unsigned long  count;
	
} BiContextType;

typedef BiContextType *BiContextTypePtr;



/**********************************************************************
 * C O N T E X T S   F O R   T M L   S Y N T A X   E L E M E N T S
 **********************************************************************
 */

#define NUM_MB_TYPE_CTX  11
#define NUM_B8_TYPE_CTX  9
#define NUM_MV_RES_CTX   10
#define NUM_REF_NO_CTX   6
#define NUM_DELTA_QP_CTX 4
#define NUM_MB_AFF_CTX 4

typedef struct
{
	BiContextType mb_type_contexts [3][NUM_MB_TYPE_CTX];
	BiContextType b8_type_contexts [2][NUM_B8_TYPE_CTX];
	BiContextType mv_res_contexts  [2][NUM_MV_RES_CTX];
	BiContextType ref_no_contexts  [2][NUM_REF_NO_CTX];
	BiContextType delta_qp_contexts   [NUM_DELTA_QP_CTX];
	BiContextType mb_aff_contexts     [NUM_MB_AFF_CTX];
	
} MotionInfoContexts;

#define NUM_IPR_CTX    2
#define NUM_CIPR_CTX   4
#define NUM_CBP_CTX    4
#define NUM_BCBP_CTX   4
#define NUM_MAP_CTX   15
#define NUM_LAST_CTX  15
#define NUM_ONE_CTX    5
#define NUM_ABS_CTX    5

typedef struct
{
	BiContextType  ipr_contexts [NUM_IPR_CTX]; 
	BiContextType  cipr_contexts[NUM_CIPR_CTX]; 
	BiContextType  cbp_contexts [3][NUM_CBP_CTX];
	BiContextType  bcbp_contexts[NUM_BLOCK_TYPES][NUM_BCBP_CTX];
	BiContextType  map_contexts [NUM_BLOCK_TYPES][NUM_MAP_CTX];
	BiContextType  last_contexts[NUM_BLOCK_TYPES][NUM_LAST_CTX];
	BiContextType  one_contexts [NUM_BLOCK_TYPES][NUM_ONE_CTX];
	BiContextType  abs_contexts [NUM_BLOCK_TYPES][NUM_ABS_CTX];
	BiContextType  fld_map_contexts [NUM_BLOCK_TYPES][NUM_MAP_CTX];
	BiContextType  fld_last_contexts[NUM_BLOCK_TYPES][NUM_LAST_CTX];
} TextureInfoContexts;

//*********************** end of data type definition for CABAC *******************

typedef struct pix_pos
{
	int available;
	int mb_addr;
	int x;
	int y;
	int pos_x;
	int pos_y;
} PixelPos;

/*! Buffer structure for RMPNI commands */
typedef struct RMPNIbuffer_s
{
	int RMPNI;
	int Data;
	struct RMPNIbuffer_s *Next;
} RMPNIbuffer_t;

/*! Buffer structure for decoded referenc picture marking commands */
typedef struct DecRefPicMarking_s
{
	int memory_management_control_operation;
	int difference_of_pic_nums_minus1;
	int long_term_pic_num;
	int long_term_frame_idx;
	int max_long_term_frame_idx_plus1;
	struct DecRefPicMarking_s *Next;
} DecRefPicMarking_t;

//! Syntaxelement
typedef struct syntaxelement
{
	int                 type;           //!< type of syntax element for data part.
	int                 value1;         //!< numerical value of syntax element
	int                 value2;         //!< for blocked symbols, e.g. run/level
	int                 len;            //!< length of code
	int                 inf;            //!< info part of UVLC code
	unsigned int        bitpattern;     //!< UVLC bitpattern
	int                 context;        //!< CABAC context
	int                 k;              //!< CABAC context for coeff_count,uv
	
#if TRACE
#define             TRACESTRING_SIZE 100            //!< size of trace string
	char                tracestring[TRACESTRING_SIZE];  //!< trace string
#endif
	
	//!< for mapping of syntaxElement to UVLC
	void    (*mapping)(int value1, int value2, int* len_ptr, int* info_ptr);
	//!< used for CABAC: refers to actual coding method of each individual syntax element type
	void    (*writing)(struct syntaxelement *, EncodingEnvironmentPtr);
/////////////////////////////////////>
	
	//! used for CABAC: refers to actual coding method of each individual syntax element type
	void  (*reading)(struct syntaxelement *, struct inp_par_dec *, struct img_par *, DecodingEnvironmentPtr);
	
	//! for mapping of UVLC to syntaxElement
    void    (*mapping_dec)(int len, int info, int *value1, int *value2);
	
	
/////////////////////////////////////<
} SyntaxElement;

//! Macroblock
typedef struct macroblock
{
	int                 currSEnr;                   //!< number of current syntax element
	int                 slice_nr;
	int                 delta_qp;
	int                 qp ;
	int                 qpsp ;
	int                 bitcounter[MAX_BITCOUNTER_MB];
	
	struct macroblock   *mb_available_up;   //!< pointer to neighboring MB (CABAC)
	struct macroblock   *mb_available_left; //!< pointer to neighboring MB (CABAC)
	
	int                 mb_type;
	int                 mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];          //!< indices correspond to [forw,backw][block_y][block_x][x,y]
	int                 intra_pred_modes[BLOCK_MULTIPLE*BLOCK_MULTIPLE];
	int                 cbp ;
	int                 cbp_blk ;    //!< 1 bit set for every 4x4 block with coefs (not implemented for INTRA)
	int                 b8mode[4];
	int                 b8pdir[4];
	unsigned long       cbp_bits;//cbp_bits 则似乎是表示当前宏块的各个分割是否有残差。
/////////////////////////////////////>
	int                 ei_flag;
	int                 i16mode;
	int                 delta_quant;          //!< for rate control

/////////////////////////////////////<	
	int                 lf_disable;
	int                 lf_alpha_c0_offset;
	int                 lf_beta_offset;
	
	int                 c_ipred_mode;      //!< chroma intra prediction mode
	int                 IntraChromaPredModeFlag;
	
	int                 mb_field;
	
	int mbAddrA, mbAddrB, mbAddrC, mbAddrD;
	int mbAvailA, mbAvailB, mbAvailC, mbAvailD;
	
	// rate control
	double              actj;               // macroblock activity measure for macroblock j
	int                 prev_qp;
	int                 prev_delta_qp;
	int                 prev_cbp;
	int                 predict_qp;
	int                 predict_error;
	
	int                 LFDisableIdc;
	int                 LFAlphaC0Offset;
	int                 LFBetaOffset;
	
	int                 skip_flag;
//新添加的，用于存储分形参数四个，包括YUV三个分量；

	int    **x;
	int    **y;
	double **scale;
	double **offset;


} Macroblock;



//! Bitstream
typedef struct
{
	int             byte_pos;           //!< current position in bitstream;
	int             bits_to_go;         //!< current bitcounter
	byte            byte_buf;           //!< current buffer for last written byte
	int             stored_byte_pos;    //!< storage for position in bitstream;
	int             stored_bits_to_go;  //!< storage for bitcounter
	byte            stored_byte_buf;    //!< storage for buffer of last written byte
	
	byte            byte_buf_skip;      //!< current buffer for last written byte
	int             byte_pos_skip;      //!< storage for position in bitstream;
	int             bits_to_go_skip;    //!< storage for bitcounter
	
/////////////////////////////////////>
	int           ei_flag;            //!< error indication, 0: no error, else unspecified error

	// CABAC Decoding
	int           read_len;           //!< actual position in the codebuffer, CABAC only
	int           code_len;           //!< overall codebuffer length, CABAC only
	// UVLC Decoding
	int           frame_bitoffset;    //!< actual position in the codebuffer, bit-oriented, UVLC only
	int           bitstream_length;   //!< over codebuffer lnegth, byte oriented, UVLC only
	
/////////////////////////////////////<	
	byte            *streamBuffer;      //!< actual buffer for written bytes
	int             write_flag;         //!< Bitstream contains data and needs to be written
	
} Bitstream;


//! DataPartition
typedef struct datapartition
{

  Bitstream           *bitstream;
  EncodingEnvironment ee_cabac;
/////////////////////////////////////>
  
  DecodingEnvironment de_cabac;

  int     (*readSyntaxElement)(SyntaxElement *, struct img_par *, struct inp_par_dec *, struct datapartition *);
          /*!< virtual function;
               actual method depends on chosen data partition and
               entropy coding method  */

  
  
/////////////////////////////////////<
  int                 (*writeSyntaxElement)(SyntaxElement *, struct datapartition *);
                      /*!< virtual function;
                           actual method depends on chosen data partition and
                           entropy coding method  */
} DataPartition;


// ! Slice
typedef struct
{
	int                 picture_id;
	int                 qp;
	int                 picture_type; //!< picture type
	int                 start_mb_nr;
	int                 max_part_nr;  //!< number of different partitions
	
/////////////////////////////////////>
	int                 ei_flag;       //!< 0 if the partArr[0] contains valid information
	
	PictureStructure    structure;     //!< Identify picture structure type


// 	int                 max_part_nr;
	int                 dp_mode;       //!< data partioning mode
	int                 next_header;

	int                 pic_parameter_set_id;   //!<the ID of the picture parameter set the slice is reffering to
	
	int                 LFDisableIdc;     //!< Disable loop filter on slice
	int                 LFAlphaC0Offset;  //!< Alpha and C0 offset for filtering slice
	int                 LFBetaOffset;     //!< Beta offset for filtering slice

	
/////////////////////////////////////<
	int                 num_mb;       //!< number of MBs in the slice
	DataPartition       *partArr;     //!< array of partitions
	MotionInfoContexts  *mot_ctx;     //!< pointer to struct of context models for use in CABAC
	TextureInfoContexts *tex_ctx;     //!< pointer to struct of context models for use in CABAC
	
	// !KS: RMPNI buffer should be retired. just do some sore simple stuff
// 	RMPNIbuffer_t        *rmpni_buffer; //!< stores the slice temporary buffer remapping commands
	
	int                 ref_pic_list_reordering_flag_l0;
	int                 *remapping_of_pic_nums_idc_l0;
	int                 *abs_diff_pic_num_minus1_l0;
	int                 *long_term_pic_idx_l0;
	int                 ref_pic_list_reordering_flag_l1;
	int                 *remapping_of_pic_nums_idc_l1;
	int                 *abs_diff_pic_num_minus1_l1;
	int                 *long_term_pic_idx_l1;
	
	Boolean             (*slice_too_big)(int bits_slice); //!< for use of callback functions
	
	int                 field_ctx[3][2]; //GB
	
} Slice;


#define MAXSLICEPERPICTURE 100
typedef struct 
{
	int   no_slices;
	int   idr_flag;
	Slice *slices[MAXSLICEPERPICTURE];
	int bits_per_picture;
	float distortion_y;
	float distortion_u;
	float distortion_v;
}Picture;

Picture *top_pic;
Picture *bottom_pic;
Picture *frame_pic;

typedef struct
{
	// Size info
	int x_size, y_framesize, y_fieldsize;  
	char *yf, *uf, *vf;                    //!< frame representation
	char *yt, *ut, *vt;                    //!< top field
	char *yb, *ub, *vb;                    //!< bottom field
} Sourceframe;


// global picture format dependend buffers, mem allocation in image.c
// byte   **imgY_org_i;           //!< Reference luma image
// byte  ***imgUV_org_i;          //!< Reference croma image
//int    **refFrArr;           //!< Array for reference frames of each block
int    **img4Y_tmp;          //!< for quarter pel interpolation

int mb_type_16[900];//Mb的预测类型，zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

unsigned int log2_max_frame_num_minus4;
unsigned int log2_max_pic_order_cnt_lsb_minus4;


int UV_dct;
int  me_tot_time,me_time;
pic_parameter_set_rbsp_t *active_pps;
seq_parameter_set_rbsp_t *active_sps;//全局变量,主要为后续解码所用。

// B pictures
// motion vector : forward, backward, direct
int  mb_adaptive;     //!< For MB level field/frame coding tools
int  MBPairIsField;     //!< For MB level field/frame coding tools


//Weighted prediction
int ***wp_weight;  // weight in [list][index][component] order
int ***wp_offset;  // offset in [list][index][component] order
int ****wbp_weight;  // weight in [list][fwd_index][bwd_idx][component] order
int luma_log_weight_denom;
int chroma_log_weight_denom;
int wp_luma_round;
int wp_chroma_round;

// global picture format dependend buffers, mem allocation in image.c (field picture)
byte   **imgY_org_top;
byte   **imgY_org_bot;

byte  ***imgUV_org_top;
byte  ***imgUV_org_bot;

byte   **imgY_org_frm;
byte  ***imgUV_org_frm;

byte   **imgY_com;               //!< Encoded luma images
byte  ***imgUV_com;              //!< Encoded croma images

int   ***direct_ref_idx;         //!< direct mode reference index buffer
int    **direct_pdir;         //!< direct mode reference index buffer

// Buffers for rd optimization with packet losses, Dim. Kontopodis
byte **pixel_map;   //!< Shows the latest reference frame that is reliable for each pixel
byte **refresh_map; //!< Stores the new values for pixel_map  
int intras;         //!< Counts the intra updates in each frame.

int  Bframe_ctr, frame_no, nextP_tr_fld, nextP_tr_frm;
int  tot_time;

#define ET_SIZE 300      //!< size of error text buffer
char errortext[ET_SIZE]; //!< buffer for error message for exit with error()


//! Info for the "decoders-in-the-encoder" used for rdoptimization with packet losses
typedef struct
{
	int  **resY;             //!< Residue of Luminance
	byte ***decY;            //!< Decoded values at the simulated decoders
	byte ****decref;         //!< Reference frames of the simulated decoders
	byte ***decY_best;       //!< Decoded frames for the best mode for all decoders
	byte **RefBlock;
	byte **status_map;
	byte **dec_mb_mode;
} Decoders;
extern Decoders *decs;

//! SNRParameters
typedef struct snr_par
{
	float snr_y;               //!< current Y SNR
	float snr_u;               //!< current U SNR
	float snr_v;               //!< current V SNR
	float snr_y1;              //!< SNR Y(dB) first frame
	float snr_u1;              //!< SNR U(dB) first frame
	float snr_v1;              //!< SNR V(dB) first frame
	float snr_ya;              //!< Average SNR Y(dB) remaining frames
	float snr_ua;              //!< Average SNR U(dB) remaining frames
	float snr_va;              //!< Average SNR V(dB) remaining frames
} SNRParameters;

#define NUM_PIC_TYPE 5

//!< statistics 统计
typedef struct
{
	int   quant0;                 //!< quant for the first frame
	int   quant1;                 //!< average quant for the remaining frames
	float bitr;                   //!< bit rate for current frame, used only for output til terminal
	float bitr0;                  //!< stored bit rate for the first frame
	float bitrate;                //!< average bit rate for the sequence except first frame
	int   bit_ctr;                //!< counter for bit usage
	int   bit_ctr_0;              //!< stored bit use for the first frame
	int   bit_ctr_n;              //!< bit usage for the current frame
	int   bit_slice;              //!< number of bits in current slice
	int   bit_ctr_emulationprevention; //!< stored bits needed to prevent start code emulation
	
	// B pictures
	int   bit_ctr_P;
	int   bit_ctr_B;
	float bitrate_P;
	float bitrate_B;
	
	int   mode_use       [NUM_PIC_TYPE][MAXMODE]; //!< Macroblock mode usage for Intra frames
	int   bit_use_mode   [NUM_PIC_TYPE][MAXMODE]; //!< statistics of bit usage
	int   bit_use_stuffingBits[NUM_PIC_TYPE];
	int   bit_use_mb_type     [NUM_PIC_TYPE];
	int   bit_use_header      [NUM_PIC_TYPE];
	int   tmp_bit_use_cbp     [NUM_PIC_TYPE];
	int   bit_use_coeffY      [NUM_PIC_TYPE];
	int   bit_use_coeffC      [NUM_PIC_TYPE];
	int   bit_use_delta_quant [NUM_PIC_TYPE];
	
	int   em_prev_bits_frm;
	int   em_prev_bits_fld;
	int  *em_prev_bits;
	int   bit_ctr_parametersets;
	int   bit_ctr_parametersets_n;
} StatParameters;

//!< For MB level field/frame coding tools
//!< temporary structure to store MB data for field/frame coding
typedef struct
{
	double min_rdcost;
	
	int    rec_mbY[16][16];       // hold the Y component of reconstructed MB
	int    rec_mbU[8][8], rec_mbV[8][8]; 
	int    ****cofAC;
	int    ***cofDC;
	int    mb_type;
	int    b8mode[4], b8pdir[4];
	int    **ipredmode;
	int    intra_pred_modes[16];
	int    cbp, cbp_blk;    //cbp 表示 8x8 是否有残差,cbp_blk表示 4x4块是否有残差
	int    mode;
	int    ******pred_mv;        //!< predicted motion vectors
	int    ******all_mv;         //!< all modes motion vectors
	int    refar[2][4][4];       //!< reference frame array [list][x][y]
	int    i16offset;
	int    c_ipred_mode;
	int    qp;
	int    prev_qp;
	int    prev_delta_qp;
} RD_DATA;

RD_DATA *rdopt; 
RD_DATA rddata_top_frame_mb, rddata_bot_frame_mb; //!< For MB level field/frame coding tools
RD_DATA rddata_top_field_mb, rddata_bot_field_mb; //!< For MB level field/frame coding tools

// extern InputParameters *input;
// extern ImageParameters *img;
extern StatParameters *stat;

extern SNRParameters *snr;


// files
FILE *p_dec;                     //!< internal decoded image for debugging
FILE *p_stat;                    //!< status file for the last encoding session
FILE *p_log;                     //!< SNR file
// FILE *p_in;                      //!< YUV
FILE *p_trace;                   //!< Trace file


// B pictures
// motion vector : forward, backward, direct
int  mb_adaptive;     //!< For MB level field/frame coding tools
int  MBPairIsField;     //!< For MB level field/frame coding tools


/////////////////////////////////////>

// pic_parameter_set_rbsp_t *active_pps;
// seq_parameter_set_rbsp_t *active_sps;
// 
// global picture format dependend buffers, mem allocation in decod.c ******************
int  **refFrArr;                                //!< Array for reference frames of each block

byte **imgY_ref;                                //!< reference frame find snr
byte ***imgUV_ref;

int  ReMapRef[20];

// B pictures
int  Bframe_ctr;
int  frame_no;

int  g_nFrame;

// For MB level frame/field coding
int  TopFieldForSkip_Y[16][16];
int  TopFieldForSkip_UV[2][16][16];


#define ET_SIZE 300      //!< size of error text buffer
char errortext[ET_SIZE]; //!< buffer for error message for exit with error()



/////////////////////////////////////<

/***********************************************************************
 * P r o t o t y p e s   f o r    T M L
 ***********************************************************************
 */

void intrapred_luma(int CurrPixX,int CurrPixY, int *left_available, int *up_available, int *all_available);
// void init();
// int  find_sad(int hadamard, int m7[16][16]);
int  dct_luma(int pos_mb1,int pos_mb2,int *cnt_nonz,int);
int  dct_luma_sp(int pos_mb1,int pos_mb2,int *cnt_nonz);
void copyblock_sp(int pos_mb1,int pos_mb2);
int  dct_chroma(int uv,int i11);
int  dct_chroma_sp(int uv,int i11);
// int  motion_search(int isi);
int  sign(int a,int b);
void intrapred_chroma(int,int,int uv);
void intrapred_luma_16x16();
int  find_sad_16x16(int *intra_mode);

int dct_luma_16x16(int);

void init_poc();

void init_img();
void report();
void information_init();
int  get_picture_type();
int  clip1a(int a);

void  LumaPrediction4x4 (int, int, int, int, int, int, int);
int   SATD (int*, int);


pel_t* FastLineX (int, pel_t*, int, int, int, int);
pel_t* UMVLineX  (int, pel_t*, int, int, int, int);

void LumaResidualCoding ();
void LumaResidualCoding_fract (int,int);
void ChromaResidualCoding (int*);
void ChromaResidualCoding_fract (int*);

void IntraChromaPrediction8x8 (int*, int*, int*);
int  writeMBHeader   (int rdopt); 
// void DeblockFrame(ImageParameters *img, byte **, byte ***);

extern int*   refbits;
extern int**** motion_cost;

void  Get_Direct_Motion_Vectors ();
void  PartitionMotionSearch     (int, int, double);
int   BIDPartitionCost          (int, int, int, int, int);
int   LumaResidualCoding8x8     (int*, int*, int, int, int, int, int, int);
int   LumaResidualCoding8x8_fract(int*, int*, int);
int   writeLumaCoeff8x8         (int, int);
int   writeMotionVector8x8      (int  i0, int  j0, int  i1, int  j1, int  refframe,int  list_idx, int  mv_mode) ;
// int   writeMotionVector8x8_fract(int  i0, int  j0, int  i1, int  j1, int  refframe,int  list_idx, int  mv_mode)
                               ;
int   writeReferenceFrame       (int, int, int, int, int);
int   writeScaleAndOffset       (int con);
int   writeScaleAndOffset_new   (int con);
int   writeXAndY                (int con);


// int   writeAbpCoeffIndex        (int, int, int, int);
int   writeIntra4x4Modes        (int);
int   writeChromaIntraPredMode  ();
// 
void estimate_weighting_factor_B_slice();
void estimate_weighting_factor_P_slice();

int  Get_Direct_Cost8x8 (int, double);
int  Get_Direct_CostMB  (double);
int  B8Mode2Value (int b8mode, int b8pdir);

int GetSkipCostMB (double lambda);
void FindSkipModeMotionVector ();
// 
// 
// // dynamic mem allocation
int  init_global_buffers();
void free_global_buffers();
void no_mem_exit  (char *where);

int  get_mem_mv  (int*******);
void free_mem_mv (int******);
void free_img    ();

int  get_mem_ACcoeff  (int*****);
int  get_mem_DCcoeff  (int****);
void free_mem_ACcoeff (int****);
void free_mem_DCcoeff (int***);

int  decide_fld_frame(float snr_frame_Y, float snr_field_Y, int bit_field, int bit_frame, double lambda_picture);
void combine_field();
// 
Picture *malloc_picture();
void     free_picture (Picture *pic);
// 
int   encode_one_slice(int SLiceGroupId, Picture *pic);   //! returns the number of MBs in the slice

void  start_macroblock(int mb_addr, int mb_field);
void  set_MB_parameters (int mb_addr);           //! sets up img-> according to input-> and currSlice->

int   writeMotionInfo2NAL ();
void   writeXYAndSOInfo2NAL_fract (int con);
void  terminate_macroblock(Boolean *end_of_slice, Boolean *recode_macroblock);
int   slice_too_big(int rlc_bits);
void  write_one_macroblock(int eos_bit);
void write_one_macroblock_fract (int eos_bit);
void  proceed2nextMacroblock();
void  proceed2nextFrame();

void free_slice_list(Picture *currPic);

void report_stats_on_error();

#if TRACE
void  trace2out(SyntaxElement *se);
#endif


void error(char *text, int code);
int  start_sequence();
int  terminate_sequence();
int  start_slice();
int  terminate_slice();

// // B pictures
// int  get_fwMV(int *min_fw_sad, int tot_intra_sad);
// void get_bwMV(int *min_bw_sad);
// void get_bid(int *bid_sad, int fw_predframe_no);
// void get_dir(int *dir_sad);
// void compare_sad(int tot_intra_sad, int fw_sad, int bw_sad, int bid_sad, int dir_sad, int);
// int  BlkSize2CodeNumber(int blc_size_h, int blc_size_v);
// 
// void InitMotionVectorSearchModule();
// 
int  field_flag_inference();
// 
void set_mbaff_parameters();  // For MB AFF
void writeVlcByteAlign(Bitstream* currStream);


int   writeLumaCoeff4x4_CABAC(int, int, int);
int   writeCBPandLumaCoeff();
int   writeChromaCoeff();
// int   writeMB_bits_for_4x4_luma   (int, int, int);
// int   writeMB_bits_for_16x16_luma ();
// int   writeMB_bits_for_luma       (int);
// int   writeMB_bits_for_DC_chroma  (int);
// int   writeMB_bits_for_AC_chroma  (int);
// int   writeMB_bits_for_CBP        ();
// 
// int   SingleUnifiedMotionSearch   (int, int, int**, int***, int*****, int, int*****, double);

// //============= rate-distortion optimization ===================
void  clear_rdopt      ();
void  init_rdopt       ();
void  RD_Mode_Decision ();
// //============= rate-distortion opt with packet losses ===========
void decode_one_macroblock();
void decode_one_mb (int, Macroblock*);
void decode_one_b8block (int, int, int, int, int);
void Get_Reference_Block(byte **imY, int block_y, int block_x, int mvhor, int mvver, byte **out);
byte Get_Reference_Pixel(byte **imY, int y, int x);
// int Half_Upsample(byte **imY, int j, int i);
void DecOneForthPix(byte **dY, byte ***dref);
// void compute_residue(int mode);
void compute_residue_b8block (int, int);
void compute_residue_mb (int);
// void UpdateDecoders();
void Build_Status_Map(byte **s_map);
void Error_Concealment(byte **inY, byte **s_map, byte ***refY);
void Conceal_Error(byte **inY, int mb_y, int mb_x, byte ***refY, byte **s_map);
//============= restriction of reference frames based on the latest intra-refreshes==========
void UpdatePixelMap();

//============= fast full integer search =======================
#ifdef _FAST_FULL_ME_
void  ClearFastFullIntegerSearch    ();
void  ResetFastFullIntegerSearch    ();
#endif

void process_2nd_IGOP();
void SetImgType();

// Tian Dong: for IGOPs
extern Boolean In2ndIGOP;
extern int start_frame_no_in_this_IGOP;
extern int start_tr_in_this_IGOP;
extern int FirstFrameIn2ndIGOP;
#define IMG_NUMBER ( img->current_frame - start_frame_no_in_this_IGOP )
#define PAYLOAD_TYPE_IDERP 8

void AllocNalPayloadBuffer();
void FreeNalPayloadBuffer();
void SODBtoRBSP(Bitstream *currStream);
int RBSPtoEBSP(byte *streamBuffer, int begin_bytepos, int end_bytepos, int min_num_bytes);
int Bytes_After_Header;

// JVT-D101: the bit for redundant_pic_cnt in slice header may be changed, 
// therefore the bit position in the bitstream must be stored.
int rpc_bytes_to_go;
int rpc_bits_to_go;
void modify_redundant_pic_cnt(unsigned char *streamBuffer);
// End JVT-D101

// Fast ME enable
int BlockMotionSearch (int,int,int,int,int,int,double);
void i_encode_one_macroblock (void);

int i_coding_main();




/////////////////////////////////////>

typedef struct old_slice_par
{
	unsigned field_pic_flag;
	unsigned bottom_field_flag;
	unsigned frame_num;
	int nal_ref_idc;
	unsigned pic_oder_cnt_lsb;
	int delta_pic_oder_cnt_bottom;
	int delta_pic_order_cnt[2];
	int idr_flag;
	int idr_pic_id;
	int pps_id;
} OldSliceParams;

extern OldSliceParams old_slice;

// files
FILE *p_out;                    //!< pointer to output YUV file
//FILE *p_out2;                    //!< pointer to debug output YUV file
FILE *p_ref;                    //!< pointer to input original reference YUV file file
FILE *p_log;                    //!< SNR file

#if TRACE
FILE *p_trace;
#endif

// prototypes
void init_conf_dec(struct inp_par_dec *inp, char *config_filename);
void report_dec(struct inp_par_dec *inp, struct img_par *img, struct snr_par *snr);
void init_dec(struct img_par *img);

void malloc_slice_dec(struct inp_par_dec *inp, struct img_par *img);
void free_slice_dec(struct inp_par_dec *inp, struct img_par *img);

int  decode_one_frame(struct img_par *img,struct inp_par_dec *inp, struct snr_par *snr);
void init_picture(struct img_par *img, struct inp_par_dec *inp);
void exit_picture();

int  read_new_slice();
void decode_one_slice(struct img_par *img,struct inp_par_dec *inp);

void start_macroblock_dec(struct img_par *img,struct inp_par_dec *inp, int CurrentMBInScanOrder);
int  read_one_macroblock(struct img_par *img,struct inp_par_dec *inp);
void read_ipred_modes(struct img_par *img,struct inp_par_dec *inp);
int  decode_one_macroblock_dec(struct img_par *img,struct inp_par_dec *inp);
int  exit_macroblock(struct img_par *img,struct inp_par_dec *inp, int eos_bit);
void decode_ipcm_mb(struct img_par *img);


void readMotionInfoFromNAL (struct img_par *img,struct inp_par_dec *inp);
void readCBPandCoeffsFromNAL(struct img_par *img,struct inp_par_dec *inp);
void readIPCMcoeffsFromNAL(struct img_par *img, struct inp_par_dec *inp, struct datapartition *dP);


void copyblock_sp_dec(struct img_par *img,int block_x,int block_y);
void itrans_sp_chroma(struct img_par *img,int ll);
void itrans(struct img_par *img,int ioff,int joff,int i0,int j0);
void itrans_sp(struct img_par *img,int ioff,int joff,int i0,int j0);
int  intrapred(struct img_par *img,int ioff,int joff,int i4,int j4);
void itrans_2(struct img_par *img);
int  intrapred_luma_16x16_dec(struct img_par *img,int predmode);
void intrapred_chroma_dec(struct img_par *img, int uv);
// int  sign(int a , int b);

// SLICE function pointers
int  (*nal_startcode_follows) ();

// NAL functions TML/CABAC bitstream
int  uvlc_startcode_follows();
int  cabac_startcode_follows();
void free_Partition(Bitstream *currStream);

// ErrorConcealment
// void reset_ec_flags();

// void error(char *text, int code);
int  is_new_picture();
void init_old_slice();

// dynamic mem allocation
int  init_global_buffers_dec();
void free_global_buffers_dec();

void frame_postprocessing(struct img_par *img, struct inp_par_dec *inp);
void field_postprocessing(struct img_par *img, struct inp_par_dec *inp);
// int  bottom_field_picture(struct img_par *img,struct inp_par_dec *inp);
void decode_slice(struct img_par *img,struct inp_par_dec *inp, int current_header);

#define PAYLOAD_TYPE_IDERP 8
int RBSPtoSODB(byte *streamBuffer, int last_byte_pos);
int EBSPtoRBSP(byte *streamBuffer, int end_bytepos, int begin_bytepos);

// For MB level frame/field coding
void init_super_macroblock(struct img_par *img,struct inp_par_dec *inp);
void exit_super_macroblock(struct img_par *img,struct inp_par_dec *inp);
int  decode_super_macroblock(struct img_par *img,struct inp_par_dec *inp);
void decode_one_Copy_topMB(struct img_par *img,struct inp_par_dec *inp);

void SetOneRefMV(struct img_par* img);
int peekSyntaxElement_UVLC(SyntaxElement *sym, struct img_par *img, struct inp_par_dec *inp, struct datapartition *dP);

void fill_wp_params(struct img_par *img);

void reset_wp_params(struct img_par *img);

void FreePartition (DataPartition *dp, int n);
DataPartition *AllocPartition();

void tracebits2(const char *trace_str, int len, int info);

void init_decoding_engine_IPCM(struct img_par *img);
void readIPCMBytes_CABAC(SyntaxElement *sym, Bitstream *currStream);


/////////////////////////////////////<


Slice *malloc_slice();
void  free_slice(Slice *slice);
void  init_slice(int start_mb_addr);
void set_ref_pic_num();
void PSNR();



#endif


//////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////