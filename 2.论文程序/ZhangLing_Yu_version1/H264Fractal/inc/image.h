#ifndef _IMAGE_H_
#define _IMAGE_H_
#include "global.h"
#include "compute.h"
#include "huffman.h"
#include "mbuffer.h"
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
// #include "mbuffer.h"

extern StorablePicture *enc_picture;
extern StorablePicture *enc_frame_picture;
extern StorablePicture *enc_top_picture;
extern StorablePicture *enc_bottom_picture;
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////

FILE *fp_in,*fp_plane,*fp_out_all_rec,**fp_out_region_rec,h264_out;
int num_regions;
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// long **xy_p,**s_p,**o_p,**redual_ple,**redual_pru,**redual_ple_u,**redual_ple_v,**redual_pru_u,**redual_pru_v;;// ¡Ÿ ±±‰¡ø
HCODE ***huff_ptr;
// void start_oneframe(char);
// int encode_oneframe (char);

void ReadOneFrame(int curFrame,SourceFrame sourceframe);
void CopyFrameToImgOrg(SourceFrame sourceframe,byte **,byte ***);

void encode_Oneframe();

/*void decode_Oneframe();*/
/*void encode_one_macroblock(int CurMb,TRANS_NODE **trans,int con);
void decode_one_macroblock(int,TRANS_NODE**,int);
double full_search(int,int,int,int,int);
int bound_chk(int m, int n, int cx, int cy,int,int,int);
void encode_block_8(int b_x,int b_y,int,int,TRANS_NODE**,int);
void encode_block_4(int b_x,int b_y,int,int);

int encode_block_rect(int b_x,int b_y,int CurMb,int Curblock_8,TRANS_NODE **trans,int con,int mode,int depth);
void decode_block_rect(int b_x,int b_y,int CurMb,int Curblock_8,TRANS_NODE **trans,int con,int mode,int,int depth);

void decode_block_8(int b_x,int b_y,int,int,TRANS_NODE **trans,int con,int i);
void decode_block_4(int b_x,int b_y,int,int,int);
*/
void write_trans();

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
int encode_one_frame ();
Boolean dummy_slice_too_big(int bits_slice);
void copy_rdopt_data (int field_type);    //!< For MB level field/frame coding tools

void UnifiedOneForthPix (StorablePicture *s);


/////////////////////////////////////>

extern StorablePicture *dec_picture;



void find_snr_dec(struct snr_par *snr, StorablePicture *p, FILE *p_ref);
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE]);
// int  picture_order(struct img_par *img);



/////////////////////////////////////<



/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////


#endif