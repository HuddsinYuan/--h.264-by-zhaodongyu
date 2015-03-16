#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "defines_enc.h"
//#include "image.h"
#include "windows.h"//专门为了计算时间而定的

InputParameters *input;
ImageParameters *img;
compressionInfo info;

SourceFrame sourceframe;
/*SourceFrame sourceframe_ref; */   
       
int imageWidth,imageHeight;
int search_range;
int num_regions;
DWORD total_bit;
DWORD tol_time;

int N_ple,N_pru,N_ple_u,N_ple_v,N_pru_u,N_pru_v;//计算数字的种类的变量。
double psnr_avg1,psnr_avg2,psnr_avg3;
double snr_avg1,snr_avg2,snr_avg3;

// int frame_type;

byte i_quality;

byte   **imgY_org;           //!< Reference luma image	//++ 待编码原始图像的亮度像素值
byte  ***imgUV_org;          //!< Reference croma image	//++ 待编码原始图像的色度像素值

byte   **plane_Y;           //!< Reference luma image	
byte   ***plane_UV;           //!< Reference luma image	

byte   **plane_Y_domain;           //!< Reference luma image
byte   ***plane_UV_domain;           //!< Reference luma image

byte   **plane_Y_domain_r;           //!< Reference luma image
byte   ***plane_UV_domain_r;           //!< Reference luma image

byte   **plane_Y_domain_l;           //!< Reference luma image
byte   ***plane_UV_domain_l;           //!< Reference luma image

byte   **plane_Y_domain_temp;           //!< Reference luma image  //临时变量 不分配内存
byte   ***plane_UV_domain_temp;           //!< Reference luma image

byte **imgY_ref;             //参考图像
byte ***imgUV_ref;

byte **imgY_reff;             //参考图像

byte **imgY_ref_c;             //中间目参考图像
byte ***imgUV_ref_c;

byte **imgY_ref_r;             //右目参考图像
byte ***imgUV_ref_r;

byte **imgY_ref_h;             //水平
byte ***imgUV_ref_h;

byte **imgY_ref_m;             //垂直
byte ***imgUV_ref_m;

byte **imgY_ref_n;             //四分之一
byte ***imgUV_ref_n;

byte **imgY_ref_l;             //左目参考图像
byte ***imgUV_ref_l;

byte **imgY_ref_temp;             //临时参考图像变量，不分配内存
byte ***imgUV_ref_temp;

byte **imgY_rec;             //重建图像像素值
byte ***imgUV_rec;

byte ***imgY_rec_region;             //重建图像像素值，用于对象
byte ****imgUV_rec_region;

// byte **imgY_rec2;             //重建图像像素值
// byte ***imgUV_rec2;

TRANS_NODE **trans_Y,**trans_U,**trans_V;
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// int **x_trans,**y_trans,**s_trans,**o_trans;//huffman code作统计用
// int **redual_transle,**redual_transru,**redual_transle_u,**redual_transle_v,**redual_transru_u,**redual_transru_v;//huffman code作统计用
// int huff_search_range;
// int trans_count[2],trans_count_x[2],trans_count_y[2],trans_count_s[2],trans_count_o[2],trans_count_re[2],trans_count_re_u[2],trans_count_re_v[2];
// TRANS_NODE **trans_Y_R,**trans_U_R,**trans_V_R;
// TRANS_NODE **trans_Y_L,**trans_U_L,**trans_V_L;

int obj,region;
int no;

double dsum1,dsum2,rsum1,rsum2,rdsum;
int current_macroblock;
////////////////////////////////////////////////////////////////////
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// unsigned char *code[10][2],*code_ptr[10][2];//码流，2考虑到基于对象
// long code_length[10][2];
byte partition[2][6000];    //块划分的模式,对于图像尺寸较大的，6000也许较小，需要设置较大的数
int partition_length[2];  //块划分的个数
//byte trans_use[2][6000];   //当OB时，记录每一个trans是否被使用到
//byte trans_region[2][6000]; //当OB时，记录每个trans的块属性，即是否是边界块
byte trans_reference[2][6000];//当stereo时，记录每一个trans的参考帧属性
////////////////////////////////////////////////////////////////////
FILE *fp_in_c,*fp_in_r,*fp_in_l,//输入的文件
     *fp_plane_c,*fp_plane_r,*fp_plane_l,//与输入文件对应的ALPHA平面
     *fp_out_c_rec[2],*fp_out_r_rec[2],*fp_out_l_rec[2], //输出重建图像，OB
	 *fp_out_all_c_rec,*fp_out_all_r_rec,*fp_out_all_l_rec,*fp_out_264;//输出的整体重建图像
///bit文件
FILE *fp_out_c[2],*fp_out_r[2],*fp_out_l[2],  //输出的压缩文件，OB
     *fp_out_all_c,*fp_out_all_r,*fp_out_all_l;//输出的压缩文件NOB

FILE *fp_out[2];//temp

//16x16块的range像素和  中目的变量
double **sum_16_org,**sum_8_org,**sum_4_org;
double **sum_16_U_org,**sum_8_U_org,**sum_4_U_org;
double **sum_16_V_org,**sum_8_V_org,**sum_4_V_org;

//16x8块range像素和
double **sum_16_8_org,**sum_8_16_org,**sum_8_4_org,**sum_4_8_org;
double **sum_16_8_U_org,**sum_8_16_U_org,**sum_8_4_U_org,**sum_4_8_U_org;
double **sum_16_8_V_org,**sum_8_16_V_org,**sum_8_4_V_org,**sum_4_8_V_org;

//16x16块的range像素平方和
double **sum2_16_org,**sum2_8_org,**sum2_4_org;
double **sum2_16_U_org,**sum2_8_U_org,**sum2_4_U_org;
double **sum2_16_V_org,**sum2_8_V_org,**sum2_4_V_org;
//16x8块的range像素平方和
double **sum2_16_8_org,**sum2_8_16_org,**sum2_8_4_org,**sum2_4_8_org;
double **sum2_16_8_U_org,**sum2_8_16_U_org,**sum2_8_4_U_org,**sum2_4_8_U_org;
double **sum2_16_8_V_org,**sum2_8_16_V_org,**sum2_8_4_V_org,**sum2_4_8_V_org;

//16x16 8x8 4x4的domain像素和
double **sum_16_ref,**sum_8_ref,**sum_4_ref;
double **sum_16_U_ref,**sum_8_U_ref,**sum_4_U_ref;
double **sum_16_V_ref,**sum_8_V_ref,**sum_4_V_ref;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref,**sum_8_16_ref,**sum_8_4_ref,**sum_4_8_ref;
double **sum_16_8_U_ref,**sum_8_16_U_ref,**sum_8_4_U_ref,**sum_4_8_U_ref;
double **sum_16_8_V_ref,**sum_8_16_V_ref,**sum_8_4_V_ref,**sum_4_8_V_ref;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref,**sum2_8_ref,**sum2_4_ref;
double **sum2_16_U_ref,**sum2_8_U_ref,**sum2_4_U_ref;
double **sum2_16_V_ref,**sum2_8_V_ref,**sum2_4_V_ref;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref,**sum2_8_16_ref,**sum2_8_4_ref,**sum2_4_8_ref;
double **sum2_16_8_U_ref,**sum2_8_16_U_ref,**sum2_8_4_U_ref,**sum2_4_8_U_ref;
double **sum2_16_8_V_ref,**sum2_8_16_V_ref,**sum2_8_4_V_ref,**sum2_4_8_V_ref;


/////////////////////////////////////////////////////////////////////////////用于多目压缩
// //右目

//16x16 8x8 4x4的domain像素和
double **sum_16_ref_r,**sum_8_ref_r,**sum_4_ref_r;
double **sum_16_U_ref_r,**sum_8_U_ref_r,**sum_4_U_ref_r;
double **sum_16_V_ref_r,**sum_8_V_ref_r,**sum_4_V_ref_r;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref_r,**sum_8_16_ref_r,**sum_8_4_ref_r,**sum_4_8_ref_r;
double **sum_16_8_U_ref_r,**sum_8_16_U_ref_r,**sum_8_4_U_ref_r,**sum_4_8_U_ref_r;
double **sum_16_8_V_ref_r,**sum_8_16_V_ref_r,**sum_8_4_V_ref_r,**sum_4_8_V_ref_r;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref_r,**sum2_8_ref_r,**sum2_4_ref_r;
double **sum2_16_U_ref_r,**sum2_8_U_ref_r,**sum2_4_U_ref_r;
double **sum2_16_V_ref_r,**sum2_8_V_ref_r,**sum2_4_V_ref_r;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref_r,**sum2_8_16_ref_r,**sum2_8_4_ref_r,**sum2_4_8_ref_r;
double **sum2_16_8_U_ref_r,**sum2_8_16_U_ref_r,**sum2_8_4_U_ref_r,**sum2_4_8_U_ref_r;
double **sum2_16_8_V_ref_r,**sum2_8_16_V_ref_r,**sum2_8_4_V_ref_r,**sum2_4_8_V_ref_r;
/////////////////////////////////////////////////////////////////////////////

//                                   左目

//16x16 8x8 4x4的domain像素和
double **sum_16_ref_l,**sum_8_ref_l,**sum_4_ref_l;
double **sum_16_U_ref_l,**sum_8_U_ref_l,**sum_4_U_ref_l;
double **sum_16_V_ref_l,**sum_8_V_ref_l,**sum_4_V_ref_l;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref_l,**sum_8_16_ref_l,**sum_8_4_ref_l,**sum_4_8_ref_l;
double **sum_16_8_U_ref_l,**sum_8_16_U_ref_l,**sum_8_4_U_ref_l,**sum_4_8_U_ref_l;
double **sum_16_8_V_ref_l,**sum_8_16_V_ref_l,**sum_8_4_V_ref_l,**sum_4_8_V_ref_l;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref_l,**sum2_8_ref_l,**sum2_4_ref_l;
double **sum2_16_U_ref_l,**sum2_8_U_ref_l,**sum2_4_U_ref_l;
double **sum2_16_V_ref_l,**sum2_8_V_ref_l,**sum2_4_V_ref_l;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref_l,**sum2_8_16_ref_l,**sum2_8_4_ref_l,**sum2_4_8_ref_l;
double **sum2_16_8_U_ref_l,**sum2_8_16_U_ref_l,**sum2_8_4_U_ref_l,**sum2_4_8_U_ref_l;
double **sum2_16_8_V_ref_l,**sum2_8_16_V_ref_l,**sum2_8_4_V_ref_l,**sum2_4_8_V_ref_l;
/////////////////////////////////////////////////////////////////////////////

//分数像素目
//16x16 8x8 4x4的domain像素和
double **sum_16_ref_H,**sum_8_ref_H,**sum_4_ref_H;
double **sum_16_U_ref_H,**sum_8_U_ref_H,**sum_4_U_ref_H;
double **sum_16_V_ref_H,**sum_8_V_ref_H,**sum_4_V_ref_H;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref_H,**sum_8_16_ref_H,**sum_8_4_ref_H,**sum_4_8_ref_H;
double **sum_16_8_U_ref_H,**sum_8_16_U_ref_H,**sum_8_4_U_ref_H,**sum_4_8_U_ref_H;
double **sum_16_8_V_ref_H,**sum_8_16_V_ref_H,**sum_8_4_V_ref_H,**sum_4_8_V_ref_H;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref_H,**sum2_8_ref_H,**sum2_4_ref_H;
double **sum2_16_U_ref_H,**sum2_8_U_ref_H,**sum2_4_U_ref_H;
double **sum2_16_V_ref_H,**sum2_8_V_ref_H,**sum2_4_V_ref_H;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref_H,**sum2_8_16_ref_H,**sum2_8_4_ref_H,**sum2_4_8_ref_H;
double **sum2_16_8_U_ref_H,**sum2_8_16_U_ref_H,**sum2_8_4_U_ref_H,**sum2_4_8_U_ref_H;
double **sum2_16_8_V_ref_H,**sum2_8_16_V_ref_H,**sum2_8_4_V_ref_H,**sum2_4_8_V_ref_H;
////////////////////////////////////////////////////////////////

//16x16 8x8 4x4的domain像素和
double **sum_16_ref_M,**sum_8_ref_M,**sum_4_ref_M;
double **sum_16_U_ref_M,**sum_8_U_ref_M,**sum_4_U_ref_M;
double **sum_16_V_ref_M,**sum_8_V_ref_M,**sum_4_V_ref_M;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref_M,**sum_8_16_ref_M,**sum_8_4_ref_M,**sum_4_8_ref_M;
double **sum_16_8_U_ref_M,**sum_8_16_U_ref_M,**sum_8_4_U_ref_M,**sum_4_8_U_ref_M;
double **sum_16_8_V_ref_M,**sum_8_16_V_ref_M,**sum_8_4_V_ref_M,**sum_4_8_V_ref_M;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref_M,**sum2_8_ref_M,**sum2_4_ref_M;
double **sum2_16_U_ref_M,**sum2_8_U_ref_M,**sum2_4_U_ref_M;
double **sum2_16_V_ref_M,**sum2_8_V_ref_M,**sum2_4_V_ref_M;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref_M,**sum2_8_16_ref_M,**sum2_8_4_ref_M,**sum2_4_8_ref_M;
double **sum2_16_8_U_ref_M,**sum2_8_16_U_ref_M,**sum2_8_4_U_ref_M,**sum2_4_8_U_ref_M;
double **sum2_16_8_V_ref_M,**sum2_8_16_V_ref_M,**sum2_8_4_V_ref_M,**sum2_4_8_V_ref_M;
/////////////////////////////////////////////
//16x16 8x8 4x4的domain像素和
double **sum_16_ref_N,**sum_8_ref_N,**sum_4_ref_N;
double **sum_16_U_ref_N,**sum_8_U_ref_N,**sum_4_U_ref_N;
double **sum_16_V_ref_N,**sum_8_V_ref_N,**sum_4_V_ref_N;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref_N,**sum_8_16_ref_N,**sum_8_4_ref_N,**sum_4_8_ref_N;
double **sum_16_8_U_ref_N,**sum_8_16_U_ref_N,**sum_8_4_U_ref_N,**sum_4_8_U_ref_N;
double **sum_16_8_V_ref_N,**sum_8_16_V_ref_N,**sum_8_4_V_ref_N,**sum_4_8_V_ref_N;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref_N,**sum2_8_ref_N,**sum2_4_ref_N;
double **sum2_16_U_ref_N,**sum2_8_U_ref_N,**sum2_4_U_ref_N;
double **sum2_16_V_ref_N,**sum2_8_V_ref_N,**sum2_4_V_ref_N;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref_N,**sum2_8_16_ref_N,**sum2_8_4_ref_N,**sum2_4_8_ref_N;
double **sum2_16_8_U_ref_N,**sum2_8_16_U_ref_N,**sum2_8_4_U_ref_N,**sum2_4_8_U_ref_N;
double **sum2_16_8_V_ref_N,**sum2_8_16_V_ref_N,**sum2_8_4_V_ref_N,**sum2_4_8_V_ref_N;

//////////////////////////////////////////////////////////////临时变量，不分配内存

//16x16 8x8 4x4的domain像素和
double **sum_16_ref_temp,**sum_8_ref_temp,**sum_4_ref_temp;
double **sum_16_U_ref_temp,**sum_8_U_ref_temp,**sum_4_U_ref_temp;
double **sum_16_V_ref_temp,**sum_8_V_ref_temp,**sum_4_V_ref_temp;
//16x8 8x16 8x4 4x8的domain像素和
double **sum_16_8_ref_temp,**sum_8_16_ref_temp,**sum_8_4_ref_temp,**sum_4_8_ref_temp;
double **sum_16_8_U_ref_temp,**sum_8_16_U_ref_temp,**sum_8_4_U_ref_temp,**sum_4_8_U_ref_temp;
double **sum_16_8_V_ref_temp,**sum_8_16_V_ref_temp,**sum_8_4_V_ref_temp,**sum_4_8_V_ref_temp;
//16x16 8x8 4x4的domain像素平方和
double **sum2_16_ref_temp,**sum2_8_ref_temp,**sum2_4_ref_temp;
double **sum2_16_U_ref_temp,**sum2_8_U_ref_temp,**sum2_4_U_ref_temp;
double **sum2_16_V_ref_temp,**sum2_8_V_ref_temp,**sum2_4_V_ref_temp;
//16x8 8x16 8x4 4x8的domain像素平方和
double **sum2_16_8_ref_temp,**sum2_8_16_ref_temp,**sum2_8_4_ref_temp,**sum2_4_8_ref_temp;
double **sum2_16_8_U_ref_temp,**sum2_8_16_U_ref_temp,**sum2_8_4_U_ref_temp,**sum2_4_8_U_ref_temp;
double **sum2_16_8_V_ref_temp,**sum2_8_16_V_ref_temp,**sum2_8_4_V_ref_temp,**sum2_4_8_V_ref_temp;

char currentVideo;
int num_16,num_8,num_4,nomatch,skip,num_16_trans,num_rect1,num_rect2,num_all,num_rl;

/*void writeFile(FILE *,);*/
#endif
