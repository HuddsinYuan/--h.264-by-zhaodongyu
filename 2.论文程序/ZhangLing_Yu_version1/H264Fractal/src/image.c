#include <stdio.h>
#include <stdlib.h>


#include "image.h"
#include "block_enc.h"
#include "block_dec.h"
#include "i_global.h"

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////
#include <assert.h>
#include "global.h"
#include "refbuf.h"
#include "mbuffer.h"
#include "header.h"
#include "intrarefresh.h"
#include "fmo.h"
#include "sei.h"
#include "memalloc.h"
#include "nalu.h"
#include "ratectl.h"
#include "mb_access.h"
#include "rtp_.h"
#include <annexb.h>
#include <parset.h>

#include "code.h"

#include <limits.h>
/////////////////////////////////////>


#include "erc_globals.h"
#include "biaridecod.h"
// #include "errorconcealment.h"
#include "context_ini.h"
#include "cabac.h"
// #include "loopfilter.h"
#include "vlc.h"
#include "erc_api.h"

/////////////////////////////////////<

FILE *trans_show_Y,*trans_show_UV;

/////////////////////////////////////>
extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern frame erc_recfr;
extern int erc_mvperMB;
extern struct img_par *erc_img;
extern tol_time;

extern StorablePicture **listX[6];
extern ColocatedParams *Co_located;

StorablePicture *dec_picture;

OldSliceParams old_slice;


/////////////////////////////////////<
static void init_frame ();
static void init_field ();

static void put_buffer_frame();
static void put_buffer_top();
static void put_buffer_bot();
static void find_distortion();
static void find_snr();

extern void DeblockFrame(ImageParameters *img, byte **, byte ***);
void org_enc_show();

static void ReportFirstframe(int tmp_time, int me_time);
static void ReportIntra(int tmp_time, int me_time);
static void ReportSP(int tmp_time, int me_time);
static void ReportBS(int tmp_time, int me_time);
static void ReportP(int tmp_time, int me_time);
static void ReportB(int tmp_time, int me_time);
static void ReportNALNonVLCBits(int tmp_time, int me_time);

void code_a_picture(Picture *pic);
void frame_picture (Picture *frame);
void field_picture (Picture *top, Picture *bottom);


/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////

void ReadOneFrame(int curFrame,SourceFrame sourceframe)
{
    int hsize=input->imagewidth;
    int vsize=input->imageheight;

	fseek(fp_in, curFrame*vsize*hsize*3/2, 0);
	fread(sourceframe.sourceframe_Y, sizeof(byte), vsize*hsize, fp_in);	
	fread(sourceframe.sourceframe_U, sizeof(byte), vsize*hsize/4, fp_in);
    fread(sourceframe.sourceframe_V, sizeof(byte), vsize*hsize/4, fp_in);
	if(num_regions==2)
	{
		fseek(fp_plane, curFrame*vsize*hsize*3/2, 0);
		fread(sourceframe.plane_Y, sizeof(byte), vsize*hsize, fp_plane);	
	    fread(sourceframe.plane_U, sizeof(byte), vsize*hsize/4, fp_plane);
        fread(sourceframe.plane_V, sizeof(byte), vsize*hsize/4, fp_plane);
	}
}

void CopyFrameToImgOrg(SourceFrame sourceframe,byte **img_y,byte ***img_uv)
{
    int i,j;
	int width,height;

	for (height=input->imageheight,i=0;i<height;i++)
		for (width=input->imagewidth,j=0;j<width;j++)
			img_y[i][j] = sourceframe.sourceframe_Y[i*width+j];//读Y分量

	for (height=input->imageheight/2,i=0;i<height;i++)
		for (width=input->imagewidth/2,j=0;j<width;j++)
		{
			img_uv[0][i][j] = sourceframe.sourceframe_U[i*width+j];//读U分量
			img_uv[1][i][j] = sourceframe.sourceframe_V[i*width+j];//读V分量
		}

	if(num_regions==2)//基于对象//no
	{
		for (height=input->imageheight,i=0;i<height;i++)
			for (width=input->imagewidth,j=0;j<width;j++)
			{
//				plane_Y_domain_temp[i][j]=sourceframe.plane_Y_domain[i*width+j];
				plane_Y[i][j] = sourceframe.plane_Y[i*width+j];
			}
				
		for (height=input->imageheight/2,i=0;i<height;i++)
			for (width=input->imagewidth/2,j=0;j<width;j++)
			{
//				plane_UV_domain_temp[0][i][j] = sourceframe.plane_U_domain[i*width+j];
//				plane_UV_domain_temp[1][i][j] = sourceframe.plane_V_domain[i*width+j];
				plane_UV[0][i][j] = sourceframe.plane_U[i*width+j];
				plane_UV[1][i][j] = sourceframe.plane_V[i*width+j];
			}
	}
}


void start_oneframe(char video)
{
	int i,j,k;
	int scaleNum=(int)((MAX_ALPHA-MIN_ALPHA)*100)/5+1;//128，应该是比例因子的范围吧？？
	int offsetNum=(MAX_BETA-MIN_BETA)/5+1;//64，应该是偏移因子的范围吧？？？
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 	trans_count[0]=trans_count[1]=trans_count_re[0]=trans_count_re_u[0]=trans_count_re_v[0]=0;
// 	for (i=0;i<num_regions;i++)
// 	{
// 		memset(x_trans[i],0,imageWidth*imageHeight/16);
// 		memset(y_trans[i],0,imageWidth*imageHeight/16);
// 		memset(s_trans[i],0,imageWidth*imageHeight/16);
// 		memset(o_trans[i],0,imageWidth*imageHeight/16);
// 		memset(redual_transle[i],0,50);
// 		memset(redual_transru[i],0,10);
// 
// 		memset(redual_transle_u[i],0,16);
// 		memset(redual_transru_u[i],0,6);
// 		memset(redual_transle_v[i],0,16);
// 		memset(redual_transru_v[i],0,6);
// 
// 	}
// 
// 	for(j=0;j<10;j++)///////////ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 		for (i=0;i<num_regions;i++)
// 		{
// 			code[j][i]=code_ptr[j][i];//后面这个在哪里赋值的
// 			code_length[j][i] = 0;
// 			for(k=0; k<(imageWidth*imageHeight*MAX_BITS/(8*4)); k++)
// 				code[j][i][k]=0;
// 		}
// 
	partition_length[0]=partition_length[1]=0;

	num_regions=input->num_regions;//1
    


	switch(video)
	{
	case 'C':
		fp_in=fp_in_c;fp_plane=fp_plane_c;
		if (num_regions>1)//no
		{
			fp_out[0]=fp_out_c[0];
			fp_out[1]=fp_out_c[1];
		}
		else//yes
			fp_out[0]=fp_out_all_c;
		
		fp_out_all_rec=fp_out_all_c_rec;//？
		fp_out_region_rec=fp_out_c_rec;//？
// 		imgY_ref_temp = imgY_ref;//
// 		imgUV_ref_temp = imgUV_ref;//
		imgY_ref_temp = imgY_ref;//
		imgUV_ref_temp = imgUV_ref;//
		plane_Y_domain_temp = plane_Y_domain;
		plane_UV_domain_temp = plane_UV_domain;

		sum_16_ref_temp=sum_16_ref;sum_8_ref_temp=sum_8_ref;sum_4_ref_temp=sum_4_ref;
		sum_16_U_ref_temp=sum_16_U_ref;sum_8_U_ref_temp=sum_8_U_ref;sum_4_U_ref_temp=sum_4_U_ref;
		sum_16_V_ref_temp=sum_16_V_ref;sum_8_V_ref_temp=sum_8_V_ref;sum_4_V_ref_temp=sum_4_V_ref;
		
		sum_16_8_ref_temp=sum_16_8_ref;sum_8_16_ref_temp=sum_8_16_ref;sum_8_4_ref_temp=sum_8_4_ref;sum_4_8_ref_temp=sum_4_8_ref;
		sum_16_8_U_ref_temp=sum_16_8_U_ref;sum_8_16_U_ref_temp=sum_8_16_U_ref;sum_8_4_U_ref_temp=sum_8_4_U_ref;sum_4_8_U_ref_temp=sum_4_8_U_ref;
		sum_16_8_V_ref_temp=sum_16_8_V_ref;sum_8_16_V_ref_temp=sum_8_16_V_ref;sum_8_4_V_ref_temp=sum_8_4_V_ref;sum_4_8_V_ref_temp=sum_4_8_V_ref;

        sum2_16_ref_temp=sum2_16_ref;sum2_8_ref_temp=sum2_8_ref;sum2_4_ref_temp=sum2_4_ref;
        sum2_16_U_ref_temp=sum2_16_U_ref;sum2_8_U_ref_temp=sum2_8_U_ref;sum2_4_U_ref_temp=sum2_4_U_ref;
		sum2_16_V_ref_temp=sum2_16_V_ref;sum2_8_V_ref_temp=sum2_8_V_ref;sum2_4_V_ref_temp=sum2_4_V_ref;

        sum2_16_8_ref_temp=sum2_16_8_ref;sum2_8_16_ref_temp=sum2_8_16_ref;sum2_8_4_ref_temp=sum2_8_4_ref;sum2_4_8_ref_temp=sum2_4_8_ref;
        sum2_16_8_U_ref_temp=sum2_16_8_U_ref;sum2_8_16_U_ref_temp=sum2_8_16_U_ref;sum2_8_4_U_ref_temp=sum2_8_4_U_ref;sum2_4_8_U_ref_temp=sum2_4_8_U_ref;
        sum2_16_8_V_ref_temp=sum2_16_8_V_ref;sum2_8_16_V_ref_temp=sum2_8_16_V_ref;sum2_8_4_V_ref_temp=sum2_8_4_V_ref;sum2_4_8_V_ref_temp=sum2_4_8_V_ref;
		break;
	case 'R':
		fp_in=fp_in_r;fp_plane=fp_plane_r;
		if (num_regions>1)
		{
			fp_out[0]=fp_out_r[0];
			fp_out[1]=fp_out_r[1];
		}
		else
			fp_out[0]=fp_out_all_r;
		
		fp_out_all_rec=fp_out_all_r_rec;
		fp_out_region_rec=fp_out_r_rec;
		imgY_ref_temp = imgY_ref_r;
		imgUV_ref_temp = imgUV_ref_r;
		plane_Y_domain_temp = plane_Y_domain_r;
		plane_UV_domain_temp = plane_UV_domain_r;

		sum_16_ref_temp=sum_16_ref_r;sum_8_ref_temp=sum_8_ref_r;sum_4_ref_temp=sum_4_ref_r;
		sum_16_U_ref_temp=sum_16_U_ref_r;sum_8_U_ref_temp=sum_8_U_ref_r;sum_4_U_ref_temp=sum_4_U_ref_r;
		sum_16_V_ref_temp=sum_16_V_ref_r;sum_8_V_ref_temp=sum_8_V_ref_r;sum_4_V_ref_temp=sum_4_V_ref_r;
		
		sum_16_8_ref_temp=sum_16_8_ref_r;sum_8_16_ref_temp=sum_8_16_ref_r;sum_8_4_ref_temp=sum_8_4_ref_r;sum_4_8_ref_temp=sum_4_8_ref_r;
		sum_16_8_U_ref_temp=sum_16_8_U_ref_r;sum_8_16_U_ref_temp=sum_8_16_U_ref_r;sum_8_4_U_ref_temp=sum_8_4_U_ref_r;sum_4_8_U_ref_temp=sum_4_8_U_ref_r;
		sum_16_8_V_ref_temp=sum_16_8_V_ref_r;sum_8_16_V_ref_temp=sum_8_16_V_ref_r;sum_8_4_V_ref_temp=sum_8_4_V_ref_r;sum_4_8_V_ref_temp=sum_4_8_V_ref_r;

        sum2_16_ref_temp=sum2_16_ref_r;sum2_8_ref_temp=sum2_8_ref_r;sum2_4_ref_temp=sum2_4_ref_r;
        sum2_16_U_ref_temp=sum2_16_U_ref_r;sum2_8_U_ref_temp=sum2_8_U_ref_r;sum2_4_U_ref_temp=sum2_4_U_ref_r;
		sum2_16_V_ref_temp=sum2_16_V_ref_r;sum2_8_V_ref_temp=sum2_8_V_ref_r;sum2_4_V_ref_temp=sum2_4_V_ref_r;

        sum2_16_8_ref_temp=sum2_16_8_ref_r;sum2_8_16_ref_temp=sum2_8_16_ref_r;sum2_8_4_ref_temp=sum2_8_4_ref_r;sum2_4_8_ref_temp=sum2_4_8_ref_r;
        sum2_16_8_U_ref_temp=sum2_16_8_U_ref_r;sum2_8_16_U_ref_temp=sum2_8_16_U_ref_r;sum2_8_4_U_ref_temp=sum2_8_4_U_ref_r;sum2_4_8_U_ref_temp=sum2_4_8_U_ref_r;
        sum2_16_8_V_ref_temp=sum2_16_8_V_ref_r;sum2_8_16_V_ref_temp=sum2_8_16_V_ref_r;sum2_8_4_V_ref_temp=sum2_8_4_V_ref_r;sum2_4_8_V_ref_temp=sum2_4_8_V_ref_r;
		break;
	case 'L':
		fp_in=fp_in_l;fp_plane=fp_plane_l;
		if (num_regions>1)
		{
			fp_out[0]=fp_out_l[0];
			fp_out[1]=fp_out_l[1];
		}
		else
			fp_out[0]=fp_out_all_l;
		
		fp_out_all_rec=fp_out_all_l_rec;
		fp_out_region_rec=fp_out_l_rec;
		imgY_ref_temp = imgY_ref_l;
		imgUV_ref_temp = imgUV_ref_l;
		plane_Y_domain_temp = plane_Y_domain_l;
		plane_UV_domain_temp = plane_UV_domain_l;
		
		sum_16_ref_temp=sum_16_ref_l;sum_8_ref_temp=sum_8_ref_l;sum_4_ref_temp=sum_4_ref_l;
		sum_16_U_ref_temp=sum_16_U_ref_l;sum_8_U_ref_temp=sum_8_U_ref_l;sum_4_U_ref_temp=sum_4_U_ref_l;
		sum_16_V_ref_temp=sum_16_V_ref_l;sum_8_V_ref_temp=sum_8_V_ref_l;sum_4_V_ref_temp=sum_4_V_ref_l;
		
		sum_16_8_ref_temp=sum_16_8_ref_l;sum_8_16_ref_temp=sum_8_16_ref_l;sum_8_4_ref_temp=sum_8_4_ref_l;sum_4_8_ref_temp=sum_4_8_ref_l;
		sum_16_8_U_ref_temp=sum_16_8_U_ref_l;sum_8_16_U_ref_temp=sum_8_16_U_ref_l;sum_8_4_U_ref_temp=sum_8_4_U_ref_l;sum_4_8_U_ref_temp=sum_4_8_U_ref_l;
		sum_16_8_V_ref_temp=sum_16_8_V_ref_l;sum_8_16_V_ref_temp=sum_8_16_V_ref_l;sum_8_4_V_ref_temp=sum_8_4_V_ref_l;sum_4_8_V_ref_temp=sum_4_8_V_ref_l;
		
        sum2_16_ref_temp=sum2_16_ref_l;sum2_8_ref_temp=sum2_8_ref_l;sum2_4_ref_temp=sum2_4_ref_l;
        sum2_16_U_ref_temp=sum2_16_U_ref_l;sum2_8_U_ref_temp=sum2_8_U_ref_l;sum2_4_U_ref_temp=sum2_4_U_ref_l;
		sum2_16_V_ref_temp=sum2_16_V_ref_l;sum2_8_V_ref_temp=sum2_8_V_ref_l;sum2_4_V_ref_temp=sum2_4_V_ref_l;
		
        sum2_16_8_ref_temp=sum2_16_8_ref_l;sum2_8_16_ref_temp=sum2_8_16_ref_l;sum2_8_4_ref_temp=sum2_8_4_ref_l;sum2_4_8_ref_temp=sum2_4_8_ref_l;
        sum2_16_8_U_ref_temp=sum2_16_8_U_ref_l;sum2_8_16_U_ref_temp=sum2_8_16_U_ref_l;sum2_8_4_U_ref_temp=sum2_8_4_U_ref_l;sum2_4_8_U_ref_temp=sum2_4_8_U_ref_l;
        sum2_16_8_V_ref_temp=sum2_16_8_V_ref_l;sum2_8_16_V_ref_temp=sum2_8_16_V_ref_l;sum2_8_4_V_ref_temp=sum2_8_4_V_ref_l;sum2_4_8_V_ref_temp=sum2_4_8_V_ref_l;
		break;
	case 'H':
// 		fp_in=fp_in_l;fp_plane=fp_plane_l;
// 		if (num_regions>1)
// 		{
// 			fp_out[0]=fp_out_l[0];
// 			fp_out[1]=fp_out_l[1];
// 		}
// 		else
// 			fp_out[0]=fp_out_all_l;
		
// 		fp_out_all_rec=fp_out_all_l_rec;
// 		fp_out_region_rec=fp_out_l_rec;
		imgY_ref_temp = imgY_ref_h;
		imgUV_ref_temp = imgUV_ref_h;
// 		plane_Y_domain_temp = plane_Y_domain_l;
// 		plane_UV_domain_temp = plane_UV_domain_l;
// 		
		sum_16_ref_temp=sum_16_ref_H;sum_8_ref_temp=sum_8_ref_H;sum_4_ref_temp=sum_4_ref_H;
		sum_16_U_ref_temp=sum_16_U_ref_H;sum_8_U_ref_temp=sum_8_U_ref_H;sum_4_U_ref_temp=sum_4_U_ref_H;
		sum_16_V_ref_temp=sum_16_V_ref_H;sum_8_V_ref_temp=sum_8_V_ref_H;sum_4_V_ref_temp=sum_4_V_ref_H;
		
		sum_16_8_ref_temp=sum_16_8_ref_H;sum_8_16_ref_temp=sum_8_16_ref_H;sum_8_4_ref_temp=sum_8_4_ref_H;sum_4_8_ref_temp=sum_4_8_ref_H;
		sum_16_8_U_ref_temp=sum_16_8_U_ref_H;sum_8_16_U_ref_temp=sum_8_16_U_ref_H;sum_8_4_U_ref_temp=sum_8_4_U_ref_H;sum_4_8_U_ref_temp=sum_4_8_U_ref_H;
		sum_16_8_V_ref_temp=sum_16_8_V_ref_H;sum_8_16_V_ref_temp=sum_8_16_V_ref_H;sum_8_4_V_ref_temp=sum_8_4_V_ref_H;sum_4_8_V_ref_temp=sum_4_8_V_ref_H;
		
        sum2_16_ref_temp=sum2_16_ref_H;sum2_8_ref_temp=sum2_8_ref_H;sum2_4_ref_temp=sum2_4_ref_H;
        sum2_16_U_ref_temp=sum2_16_U_ref_H;sum2_8_U_ref_temp=sum2_8_U_ref_H;sum2_4_U_ref_temp=sum2_4_U_ref_H;
		sum2_16_V_ref_temp=sum2_16_V_ref_H;sum2_8_V_ref_temp=sum2_8_V_ref_H;sum2_4_V_ref_temp=sum2_4_V_ref_H;
		
        sum2_16_8_ref_temp=sum2_16_8_ref_H;sum2_8_16_ref_temp=sum2_8_16_ref_H;sum2_8_4_ref_temp=sum2_8_4_ref_H;sum2_4_8_ref_temp=sum2_4_8_ref_H;
        sum2_16_8_U_ref_temp=sum2_16_8_U_ref_H;sum2_8_16_U_ref_temp=sum2_8_16_U_ref_H;sum2_8_4_U_ref_temp=sum2_8_4_U_ref_H;sum2_4_8_U_ref_temp=sum2_4_8_U_ref_H;
        sum2_16_8_V_ref_temp=sum2_16_8_V_ref_H;sum2_8_16_V_ref_temp=sum2_8_16_V_ref_H;sum2_8_4_V_ref_temp=sum2_8_4_V_ref_H;sum2_4_8_V_ref_temp=sum2_4_8_V_ref_H;
		break;
	case 'M':
		// 		fp_in=fp_in_l;fp_plane=fp_plane_l;
		// 		if (num_regions>1)
		// 		{
		// 			fp_out[0]=fp_out_l[0];
		// 			fp_out[1]=fp_out_l[1];
		// 		}
		// 		else
		// 			fp_out[0]=fp_out_all_l;
		
		// 		fp_out_all_rec=fp_out_all_l_rec;
		// 		fp_out_region_rec=fp_out_l_rec;
		imgY_ref_temp = imgY_ref_m;
		imgUV_ref_temp = imgUV_ref_m;
		// 		plane_Y_domain_temp = plane_Y_domain_l;
		// 		plane_UV_domain_temp = plane_UV_domain_l;
		// 		
		sum_16_ref_temp=sum_16_ref_M;sum_8_ref_temp=sum_8_ref_M;sum_4_ref_temp=sum_4_ref_M;
		sum_16_U_ref_temp=sum_16_U_ref_M;sum_8_U_ref_temp=sum_8_U_ref_M;sum_4_U_ref_temp=sum_4_U_ref_M;
		sum_16_V_ref_temp=sum_16_V_ref_M;sum_8_V_ref_temp=sum_8_V_ref_M;sum_4_V_ref_temp=sum_4_V_ref_M;
		
		sum_16_8_ref_temp=sum_16_8_ref_M;sum_8_16_ref_temp=sum_8_16_ref_M;sum_8_4_ref_temp=sum_8_4_ref_M;sum_4_8_ref_temp=sum_4_8_ref_M;
		sum_16_8_U_ref_temp=sum_16_8_U_ref_M;sum_8_16_U_ref_temp=sum_8_16_U_ref_M;sum_8_4_U_ref_temp=sum_8_4_U_ref_M;sum_4_8_U_ref_temp=sum_4_8_U_ref_M;
		sum_16_8_V_ref_temp=sum_16_8_V_ref_M;sum_8_16_V_ref_temp=sum_8_16_V_ref_M;sum_8_4_V_ref_temp=sum_8_4_V_ref_M;sum_4_8_V_ref_temp=sum_4_8_V_ref_M;
		
        sum2_16_ref_temp=sum2_16_ref_M;sum2_8_ref_temp=sum2_8_ref_M;sum2_4_ref_temp=sum2_4_ref_M;
        sum2_16_U_ref_temp=sum2_16_U_ref_M;sum2_8_U_ref_temp=sum2_8_U_ref_M;sum2_4_U_ref_temp=sum2_4_U_ref_M;
		sum2_16_V_ref_temp=sum2_16_V_ref_M;sum2_8_V_ref_temp=sum2_8_V_ref_M;sum2_4_V_ref_temp=sum2_4_V_ref_M;
		
        sum2_16_8_ref_temp=sum2_16_8_ref_M;sum2_8_16_ref_temp=sum2_8_16_ref_M;sum2_8_4_ref_temp=sum2_8_4_ref_M;sum2_4_8_ref_temp=sum2_4_8_ref_M;
        sum2_16_8_U_ref_temp=sum2_16_8_U_ref_M;sum2_8_16_U_ref_temp=sum2_8_16_U_ref_M;sum2_8_4_U_ref_temp=sum2_8_4_U_ref_M;sum2_4_8_U_ref_temp=sum2_4_8_U_ref_M;
        sum2_16_8_V_ref_temp=sum2_16_8_V_ref_M;sum2_8_16_V_ref_temp=sum2_8_16_V_ref_M;sum2_8_4_V_ref_temp=sum2_8_4_V_ref_M;sum2_4_8_V_ref_temp=sum2_4_8_V_ref_M;
		break;
	case 'N':
		// 		fp_in=fp_in_l;fp_plane=fp_plane_l;
		// 		if (num_regions>1)
		// 		{
		// 			fp_out[0]=fp_out_l[0];
		// 			fp_out[1]=fp_out_l[1];
		// 		}
		// 		else
		// 			fp_out[0]=fp_out_all_l;
		
		// 		fp_out_all_rec=fp_out_all_l_rec;
		// 		fp_out_region_rec=fp_out_l_rec;
		imgY_ref_temp = imgY_ref_n;
		imgUV_ref_temp = imgUV_ref_n;
		// 		plane_Y_domain_temp = plane_Y_domain_l;
		// 		plane_UV_domain_temp = plane_UV_domain_l;
		// 		
		sum_16_ref_temp=sum_16_ref_N;sum_8_ref_temp=sum_8_ref_N;sum_4_ref_temp=sum_4_ref_N;
		sum_16_U_ref_temp=sum_16_U_ref_N;sum_8_U_ref_temp=sum_8_U_ref_N;sum_4_U_ref_temp=sum_4_U_ref_N;
		sum_16_V_ref_temp=sum_16_V_ref_N;sum_8_V_ref_temp=sum_8_V_ref_N;sum_4_V_ref_temp=sum_4_V_ref_N;
		
		sum_16_8_ref_temp=sum_16_8_ref_N;sum_8_16_ref_temp=sum_8_16_ref_N;sum_8_4_ref_temp=sum_8_4_ref_N;sum_4_8_ref_temp=sum_4_8_ref_N;
		sum_16_8_U_ref_temp=sum_16_8_U_ref_N;sum_8_16_U_ref_temp=sum_8_16_U_ref_N;sum_8_4_U_ref_temp=sum_8_4_U_ref_N;sum_4_8_U_ref_temp=sum_4_8_U_ref_N;
		sum_16_8_V_ref_temp=sum_16_8_V_ref_N;sum_8_16_V_ref_temp=sum_8_16_V_ref_N;sum_8_4_V_ref_temp=sum_8_4_V_ref_N;sum_4_8_V_ref_temp=sum_4_8_V_ref_N;
		
        sum2_16_ref_temp=sum2_16_ref_N;sum2_8_ref_temp=sum2_8_ref_N;sum2_4_ref_temp=sum2_4_ref_N;
        sum2_16_U_ref_temp=sum2_16_U_ref_N;sum2_8_U_ref_temp=sum2_8_U_ref_N;sum2_4_U_ref_temp=sum2_4_U_ref_N;
		sum2_16_V_ref_temp=sum2_16_V_ref_N;sum2_8_V_ref_temp=sum2_8_V_ref_N;sum2_4_V_ref_temp=sum2_4_V_ref_N;
		
        sum2_16_8_ref_temp=sum2_16_8_ref_N;sum2_8_16_ref_temp=sum2_8_16_ref_N;sum2_8_4_ref_temp=sum2_8_4_ref_N;sum2_4_8_ref_temp=sum2_4_8_ref_N;
        sum2_16_8_U_ref_temp=sum2_16_8_U_ref_N;sum2_8_16_U_ref_temp=sum2_8_16_U_ref_N;sum2_8_4_U_ref_temp=sum2_8_4_U_ref_N;sum2_4_8_U_ref_temp=sum2_4_8_U_ref_N;
        sum2_16_8_V_ref_temp=sum2_16_8_V_ref_N;sum2_8_16_V_ref_temp=sum2_8_16_V_ref_N;sum2_8_4_V_ref_temp=sum2_8_4_V_ref_N;sum2_4_8_V_ref_temp=sum2_4_8_V_ref_N;
		break;
	}
//huffman//////ZZZZZZZZZZZZZZZZZZZ	
// 	for(i=0;i<num_regions;i++) 
// 	{
// 		for (j=0;j<(1<<huff_search_range);j++)
// 		    xy_p[i][j] = 0;
// 		for (j=0;j<scaleNum;j++)
// 			s_p[i][j] = 0;
// 		for (j=0;j<offsetNum;j++)
// 			o_p[i][j] = 0;
// 		for (j=0;j<50;j++)
// 			redual_ple[i][j] = 0;
// 		for (j=0;j<10;j++)
// 			redual_pru[i][j] = 0;
// 		for (j=0;j<10;j++)
// 			redual_ple_u[i][j] = 0;
// 		for (j=0;j<10;j++)
// 			redual_pru_u[i][j] = 0;
// 		for (j=0;j<10;j++)
// 			redual_ple_v[i][j] = 0;
// 		for (j=0;j<10;j++)
// 			redual_pru_v[i][j] = 0;
// 
// 	}
}

int encode_oneframe(int frameNum)
{    
	int x,y;
	int len;

//>encode_one_frame

	time_t ltime1;
	time_t ltime2;

#ifdef WIN32
	struct _timeb tstruct1;
	struct _timeb tstruct2;
#else
	struct timeb tstruct1;
	struct timeb tstruct2;
#endif


	int tmp_time;
	int bits_frm = 0, bits_fld = 0;
	float dis_frm = 0, dis_frm_y = 0, dis_frm_u = 0, dis_frm_v = 0;
	float dis_fld = 0, dis_fld_y = 0, dis_fld_u = 0, dis_fld_v = 0;
	
	//Rate control
	int pic_type, bits = 0; 
	
	me_time=0;

#ifdef WIN32
	_ftime (&tstruct1);           // start time ms
#else
	ftime (&tstruct1);
#endif
	time (&ltime1);               // start time s

	//Rate control 
	img->write_macroblock = 0;

	put_buffer_frame ();      // sets the pointers to the frame structures 
	                          // (and not to one of the field structures)
	init_frame ();//初始化一帧
// 	FrameNumberInFile = 
	CalculateFrameNumber();

	ReadOneFrame(frameNum,sourceframe);  //读取一帧图像
	CopyFrameToImgOrg(sourceframe,imgY_org,imgUV_org);//一维变二维,将原始图像的Y,U,V分量赋给数组img_y[i][j]，img_uv[0][i][j]
    compute_range_Sum();//事先计算range的像素和与像素平方和，节省计算时间

	// Set parameters for directmode and Deblocking filter 为直接模式和滤波设置参数
	img->direct_type     = input->direct_type;//1
	img->LFDisableIdc    = input->LFDisableIdc;//0
	img->LFAlphaC0Offset = input->LFAlphaC0Offset;//0
	img->LFBetaOffset    = input->LFBetaOffset;//0
	
    if (input->PicInterlace == FIELD_CODING)//场编码
    {
      //Rate control
      img->FieldControl=1;

      img->field_picture = 1;  // we encode fields
      field_picture (top_pic, bottom_pic);
      img->fld_flag = 1;
    }
    else//非场编码
    {
      //Rate control
      img->FieldControl=0;

      // For frame coding, turn MB level field/frame coding flag on
      if (input->MbInterlace)  //0                                    
        mb_adaptive = 1;

      img->field_picture = 0; // we encode a frame//下面的码率控制删除了

    }
    if( active_sps->frame_mbs_only_flag)//yes,仅允许宏块帧编码
		img->TopFieldFlag=0;

//>frame_picture
	img->structure = FRAME;
	img->PicSizeInMbs = img->FrameSizeInMbs;//一帧包含的宏块总数

// 	enc_picture  = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);//添加后会大幅度的降低UV的PSNR，暂时先留着
	//下面的添加用处不大，几乎就没有吧。
// 	img->ThisPOC=enc_picture->poc=img->framepoc;//2
// 	enc_picture->top_poc    = img->toppoc;//2
// 	enc_picture->bottom_poc = img->bottompoc;//2
// 	
// 	enc_picture->frame_poc = img->framepoc;//2
// 	
// 	enc_picture->pic_num = img->frame_num;//1
// 	enc_picture->coded_frame = 1;
// 	
// 	enc_picture->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);//0

// 	enc_frame_picture  = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);//添加后会大幅度的降低UV的PSNR，暂时先留着
// 	img->ThisPOC=enc_frame_picture->poc=img->framepoc;//2
// 	enc_frame_picture->top_poc    = img->toppoc;//2
// 	enc_frame_picture->bottom_poc = img->bottompoc;//2
// 	
// 	enc_frame_picture->frame_poc = img->framepoc;//2
// 	
// 	enc_frame_picture->pic_num = img->frame_num;//1
// 	enc_frame_picture->coded_frame = 1;
// 	
// 	enc_frame_picture->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);//0
// 	
// 	enc_picture=enc_frame_picture;
	 
	stat->em_prev_bits_frm = 0;//stat是什么作用？？？
	stat->em_prev_bits = &stat->em_prev_bits_frm;//0
	img->fld_flag = 0;

//>code_a_picture
    img->currentPicture = frame_pic;
	img->currentPicture->idr_flag = ((!IMG_NUMBER) && (!(img->structure==BOTTOM_FIELD))) || (input->idr_enable && (img->type == I_SLICE || img->type==SP_SLICE || img->type==SI_SLICE)&& (!(img->structure==BOTTOM_FIELD)));//等于0
	frame_pic->no_slices = 0;
	frame_pic->distortion_u = frame_pic->distortion_v = frame_pic->distortion_y = 0.0;
    frame_pic->idr_flag =0;
	// The slice_group_change_cycle can be changed here.
	// FmoInit() is called before coding each picture, frame or field
	img->slice_group_change_cycle=1;
	FmoInit(img, active_pps, active_sps);///////////////////////////////////与片组相关，确定宏块到片组的映射关系。
	
	FmoStartPicture ();           //! picture level initialization of FMO///得到每一个片中的第一个宏块地址。

//>encode_one_slice
// 	Boolean end_of_slice = FALSE;
// 	Boolean recode_macroblock;

	img->cod_counter = 0;
	init_slice (0);////////////////////////////为当前要编码的slice申请一个Slice类型结构体，并进行初始化
	Bytes_After_Header = img->currentSlice->partArr[0].bitstream->byte_pos;//等于0
	len = start_slice ();/////////////////写slice_header//24
	// Rate control
	img->NumberofHeaderBits +=len;//2096+24=2120
	// Update statistics
	stat->bit_slice += len;//0+24
	stat->bit_use_header[img->type] += len;//0+24
	encode_Oneframe();  //分形编码一个P帧
 	decode_Oneframe();//分形解码一个P帧
	process_oneframe();//对分形的残差进行变换量化反量化反变换，并且熵编码分形语法元素，写入码流
	terminate_slice();
	stat->bit_slice = 0;
//<encode_one_slice
	DeblockFrame(img, enc_picture->imgY, enc_picture->imgUV);//---zl---去方块滤波的加入---
    // 	org_enc_show(1);
//<code_a_picture

//<frame_picture
    writeout_picture (frame_pic);//
	if (frame_pic)
		free_slice_list(frame_pic);
// 	FreeSourceframe (srcframe);//P帧的时候没有申请这个空间了
	if (img->fld_flag)
		stat->bit_ctr_emulationprevention += stat->em_prev_bits_fld;
	else
		stat->bit_ctr_emulationprevention += stat->em_prev_bits_frm;
	
	if (img->type != B_SLICE)
	{
		img->pstruct_next_P = img->fld_flag;
	}

// 	PSNR();
	find_snr();//fract
	time (&ltime2);               // end time sec
#ifdef WIN32
	_ftime (&tstruct2);           // end time ms
#else
	ftime (&tstruct2);            // end time ms
#endif
	tmp_time = (ltime2 * 1000 + tstruct2.millitm) - (ltime1 * 1000 + tstruct1.millitm);
	tot_time = tot_time + tmp_time;
		switch (img->type)
		{
		case I_SLICE:
			stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
			ReportIntra(tmp_time,me_time);
			//ReportIntra(tmp_time);
			break;
		case SP_SLICE:
			stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
			ReportSP(tmp_time,me_time);
			//ReportSP(tmp_time);
			break;
		case B_SLICE:
			stat->bit_ctr_B += stat->bit_ctr - stat->bit_ctr_n;
			if (img->nal_reference_idc>0)
				ReportBS(tmp_time,me_time);
			//ReportBS(tmp_time);
			else
				ReportB(tmp_time,me_time);
			//ReportB(tmp_time);
			
			break;
		default:      // P
// 			stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
			ReportP(tmp_time,me_time);
			//ReportP(tmp_time);
		}
// 	}
	stat->bit_ctr_n = stat->bit_ctr;
	stat->bit_ctr_parametersets_n=0;
// 	FreeSourceframe (srcframe);
//<encode_one_frame

	//为下一帧做准备，存储上一帧的yuv分量的四个分形参数到pre_

	proceed2nextFrame();

}
void proceed2nextFrame()//
{
	int i,j;
	for (j=0;j<img->height/BLOCK_SIZE;j++)
	{
		for (i=0;i<(img->width/BLOCK_SIZE)*3/2;i++)
		{
			enc_picture->pre_ref_pic_id[0][i][j]=enc_picture->ref_pic_id[0][i][j];//。
			enc_picture->pre_ref_pic_id[1][i][j]=enc_picture->ref_pic_id[1][i][j];//。
		}
	}
//要不要把mv也存储到一个变量比如说pre_mv中？这样就可以根据前后帧两个块的类型不同来确定如何编码。暂时不要了吧，以后再说。
}

//---------------------------------------------------------分形解码一个P帧,生成重建帧（滤波后）和参考帧（不滤波）---------------------------------------------------------------------
void decode_Oneframe()
{
	int i,i1,j1,k;

	if(num_regions>1)//no  //num_regions=1
	{
		for(i=0;i<num_regions;i++)
		{
			memset(imgY_rec_region[i][0],0,input->imageheight*input->imagewidth);
			memset(imgUV_rec_region[0][i][0],0,input->imageheight*input->imagewidth/4);
			memset(imgUV_rec_region[1][i][0],0,input->imageheight*input->imagewidth/4);
		}	
	}
	for (current_macroblock=0;current_macroblock<img->frmSizeInMbs;current_macroblock++)//遍历每一个宏块
	{
		decode_one_macroblock(current_macroblock,trans_Y,1);//解码一个宏块的Y分量
	}	
	
	if(num_regions>1)//no   //num_regions=1
	for(i=0;i<num_regions;i++)
	{
		for(j1=0;j1<input->imageheight;j1++)
			for(i1=0;i1<input->imagewidth;i1++)
				if(imgY_rec_region[i][j1][i1]!=0)
					imgY_rec[j1][i1]=imgY_rec_region[i][j1][i1];
	}

	for (current_macroblock=0;current_macroblock<img->frmSizeInMbs/4;current_macroblock++)
	{
		decode_one_macroblock(current_macroblock,trans_U,2);//解码一个宏块的U分量
	}

	if(num_regions>1)//no
	for(i=0;i<num_regions;i++)
	{
		for(j1=0;j1<input->imageheight/2;j1++)
			for(i1=0;i1<input->imagewidth/2;i1++)
				if(imgUV_rec_region[0][i][j1][i1]!=0)
					imgUV_rec[0][j1][i1]=imgUV_rec_region[0][i][j1][i1];
	}

	for (current_macroblock=0;current_macroblock<img->frmSizeInMbs/4;current_macroblock++)
	{
		decode_one_macroblock(current_macroblock,trans_V,3);//解码一个宏块的V分量
	}

	if(num_regions>1)//no 
	for(i=0;i<num_regions;i++)
	{
		for(j1=0;j1<input->imageheight/2;j1++)
			for(i1=0;i1<input->imagewidth/2;i1++)
				if(imgUV_rec_region[1][i][j1][i1]!=0)
					imgUV_rec[1][j1][i1]=imgUV_rec_region[1][i][j1][i1];
	}

	if(num_regions>1)//no
		for(i=0;i<num_regions;i++)
		{
			fwrite(imgY_rec_region[i][0],1,input->imageheight*input->imagewidth,fp_out_region_rec[i]);
			fwrite(imgUV_rec_region[0][i][0],1,input->imageheight*input->imagewidth/4,fp_out_region_rec[i]);
			fwrite(imgUV_rec_region[1][i][0],1,input->imageheight*input->imagewidth/4,fp_out_region_rec[i]);
		}

///////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////
	if (num_regions==1)//yes   作用：对位于16*16块边界像素进行滤波
	{
		for (k=0;k<input->imageheight;k++)
			for(i=0;i<input->imagewidth;i++)
			{ 
				if (((i+1)%16==0)&&(i>2)&&(i<input->imagewidth-2))//将原来的具体数据改为input->imagewidth
				{
					imgY_reff[0][k*input->imagewidth+i]=(imgY_rec[0][k*input->imagewidth+i-1]+3*imgY_rec[0][k*input->imagewidth+i]
						+imgY_rec[0][k*input->imagewidth+i+1])/5;//16*16块边界像素进行滤波
				}
				else
					if (((i)%16==0)&&(i>2)&&(i<input->imagewidth-2))
					{
						imgY_reff[0][k*input->imagewidth+i]=(imgY_rec[0][k*input->imagewidth+i-1]+3*imgY_rec[0][k*input->imagewidth+i]
							+imgY_rec[0][k*input->imagewidth+i+1])/5;//16*16块边界像素进行滤波
					}else
						if (((k+1)%16==0)&&(k>2)&&(k<input->imageheight-2))
						{
							imgY_reff[0][k*input->imagewidth+i]=(imgY_rec[0][(k-1)*input->imagewidth+i]+3*imgY_rec[0][k*input->imagewidth+i]
								+imgY_rec[0][(k+1)*input->imagewidth+i])/5;//16*16块边界像素进行滤波
						}else
							if (((k)%16==0)&&(k>2)&&(k<input->imageheight-2))
							{
								imgY_reff[0][k*input->imagewidth+i]=(imgY_rec[0][(k-1)*input->imagewidth+i]+3*imgY_rec[0][k*input->imagewidth+i]
									+imgY_rec[0][(k+1)*input->imagewidth+i])/5;//16*16块边界像素进行滤波
							}else
				imgY_reff[k][i]=imgY_rec[k][i];//不在块边界的不用滤波
			}
	}
// 	if (1==num_regions)
// 	{
// 		memcpy(imgY_rec[0],imgY_reff[0],input->imageheight*input->imagewidth);
// 	}

//写YUV图像
// 	fwrite(imgY_rec[0],1,input->imageheight*input->imagewidth,fp_out_all_rec);
// 	fwrite(imgUV_rec[0][0],1,input->imageheight*input->imagewidth/4,fp_out_all_rec);
// 	fwrite(imgUV_rec[1][0],1,input->imageheight*input->imagewidth/4,fp_out_all_rec);

	if (num_regions>1)//no
	{
		memcpy(plane_Y_domain_temp[0],plane_Y[0],input->imageheight*input->imagewidth);
		memcpy(plane_UV_domain_temp[0][0],plane_UV[0][0],input->imageheight*input->imagewidth/4);
		memcpy(plane_UV_domain_temp[1][0],plane_UV[1][0],input->imageheight*input->imagewidth/4);
	}
	
//制作参考帧，为分形解码的图像
	memcpy(imgY_ref_temp[0],imgY_rec[0],input->imageheight*input->imagewidth);//参考帧Y分量
	memcpy(imgUV_ref_temp[0][0],imgUV_rec[0][0],input->imageheight*input->imagewidth/4);//参考帧U分量
	memcpy(imgUV_ref_temp[1][0],imgUV_rec[1][0],input->imageheight*input->imagewidth/4);//参考帧V分量

}


/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////


void dct_luma_fractal_4(int block4_x,int block4_y,int block16_x,int block16_y)//对分形的残差进行变换量化反量化反变换
{
	int diff[16],i,j,k,dummy,nonzero;
    int block16_x_n,block16_y_n,block_x,block_y;
	int pix_x,pix_y;

	block16_x_n=img->width/16;
	block16_y_n=img->height/16;

	pix_x=block16_x*16+block4_x*4;
	pix_y=block16_y*16+block4_y*4;
	
	img->pix_x=block16_x*16;
	img->pix_y=block16_y*16;

	block_x=block4_x*4;
	block_y=block4_y*4;
	UV_dct = 0;
	for ( j=0; j<4; j++)
		for ( i=0; i<4; i++)
		{
			img->mpr[i+block_x][j+block_y] = imgY_rec[j+pix_y][i+pix_x];
			img->m7[i][j] = imgY_org[pix_y+j][pix_x+i] - imgY_rec[j+pix_y][i+pix_x];
		}
	dummy = 0;
// 	   printf("\n");

    nonzero = dct_luma(block_x, block_y, &dummy, 1);
// 	   printf("\n");

}

void dct_luma_fractal_16(int block16_x,int block16_y)//对分形的残差进行变换量化反量化反变换
{
	int block4_x,block4_y,block8x8,dummy,nonzero,block16x16;
	int*  ACLevel;
	int*  ACRun;
	int  mb_x, mb_y, i, j, k;
	int  b8, b4,level,run;
	int block16_x_n,block16_y_n,block_x,block_y;
	int pix_x,pix_y,sum_cnt_nonz;
  
	block16_x_n=img->width/16;
	block16_y_n=img->height/16;

	img->opix_x=img->pix_x=block16_x*16;
	img->opix_y=img->pix_y=block16_y*16;

// 	img->current_mb_nr=block16_y*block16_x_n+block16_x;
// 	Macroblock *currMB = &img->mb_data[img->current_mb_nr];
// 	currMB->cbp     = 0 ;
// 	currMB->cbp_blk = 0 ;
// 	sum_cnt_nonz    = 0 ;
// 	for (block16x16=0;block16x16<4;block16x16++)
// 	{
// 
// 		for (block8x8=0;block8x8<4;block8x8++)
// 		{
// 			sum_cnt_nonz += LumaResidualCoding8x8_fract (&(currMB->cbp),&(currMB->cbp_blk),block8x8);
// 		}
//     	if (sum_cnt_nonz <= 5 )
// 		{
// 			currMB->cbp     &= 0xfffff0 ;
// 			currMB->cbp_blk &= 0xff0000 ;
// 			for (i=0; i < MB_BLOCK_SIZE; i++)
// 			{
// 				for (j=0; j < MB_BLOCK_SIZE; j++)
// 				{
// 					enc_picture->imgY[img->pix_y+j][img->pix_x+i]=img->mpr[i][j];
// 				}
// 			}
// 		}
// 
// 	}

// 	int a4=0;
// 	printf("\n");
	
	for(block4_y=0;block4_y<4;block4_y++)
		for (block4_x=0;block4_x<4;block4_x++)
		{
// 			a4++;
// 			printf("a4=%d",a4);
			dct_luma_fractal_4(block4_x,block4_y, block16_x,block16_y);
				
		}
// 	printf("\n");

}

void process_oneframe()//对分形的残差进行变换量化反量化反变换，并且熵编码分形语法元素，写入码流
//在这里面完成h264 encode_one_slice()的功能，比如添加片头（此处的片就是一帧的理解吧）的处理，和片尾的处理等
{
	int block16_x_n,block16_y_n,block16_x,block16_y,j,i;
	int sum_cnt_nonz;
	int CurrentMbAddr;
	Boolean end_of_slice = FALSE1;
	Boolean recode_macroblock;
	int NumberOfCodedMBs = 0,dummy;

	int MVbit;
	MVbit=0;
// 	int a16=0;
	block16_x_n=img->width/16;//一个图像行包含多少的宏块
	block16_y_n=img->height/16;//一个图像列包含包少个宏块
	N_ple=0,N_pru=0,N_ple_u=0,N_ple_v=0,N_pru_u=0,N_pru_v=0;//计算数字的种类的变量。
	UV_dct = 0;
	
	CurrentMbAddr=0;
	while (end_of_slice == FALSE1) // loop over macroblocks 对每个slice进行编码，用帧或场方式
	{

		recode_macroblock = FALSE1;
		rdopt = &rddata_top_frame_mb;   // store data in top frame MB 
		block16_y=CurrentMbAddr/block16_x_n;//将当前宏块的光栅扫描地址转换为二维x坐标(以宏块为单位)
		block16_x=CurrentMbAddr%block16_x_n;//将当前宏块的光栅扫描地址转换为二维y坐标(以宏块为单位)
		
//		start_macroblock (CurrentMbAddr, FALSE1);
		LumaResidualCoding_fract(block16_x,block16_y);//fract-DCT变换量化反量化反变换---
//实现264中set_stored_macroblock_parameters ()的功能中，
//+++    (1)保存重建图像（残差编码的时候保存了）
//+++    (2)变换系数和cbp（残差保存了）  
//+++    (3)宏块模式 （在宏块编码的时候保存了） 
//+++    currMB->mb_type = mode(best_mode);  
//+++    (4)参考帧保存  （不用保存）
//+++    rdopt.c 1511行,对前向参考,后向参考帧分别进行了保存  
//+++    (5)帧内预测模式  （这个没有保存，因为是帧间，不知道有影响不）
//+++    currMB->c_ipred_mode = best_c_imode;  
//+++    (6)为当前块保存运动向量 （在宏块编码的时候保存了）  
//+++    SetMotionVectorsMB (currMB, bframe);  

		dummy = 0;
		set_MB_parameters (CurrentMbAddr);//只计算了当前宏块的坐标，没有干什么
		ChromaResidualCoding_fract (&dummy);//---zl---DCT变换量化重排序---

		img->current_mb_nr=CurrentMbAddr;//当前宏块地知
// MVbit+=
		write_one_macroblock_fract (1);//熵编码一个宏块的分形语法元素，写入码流
		terminate_macroblock (&end_of_slice, &recode_macroblock);
		
		// printf ("encode_one_slice: mb %d,  slice %d,   bitbuf bytepos %d EOS %d\n", 
		//       img->current_mb_nr, img->current_slice_nr, 
		//       img->currentSlice->partArr[0].bitstream->byte_pos, end_of_slice);
		
		if (recode_macroblock == FALSE1)       // The final processing of the macroblock has been done
		{
			CurrentMbAddr = FmoGetNextMBNr (CurrentMbAddr);//下一个宏块地址
			if (CurrentMbAddr == -1)   // end of slice
			{
				// printf ("FMO End of Slice Group detected, current MBs %d, force end of slice\n", NumberOfCodedMBs+1);
				end_of_slice = TRUE1;
			}
			NumberOfCodedMBs++;       // only here we are sure that the coded MB is actually included in the slice//在当前帧中已编码的宏块个数
			proceed2nextMacroblock (CurrentMbAddr);
		}
		else//no
		{
			//!Go back to the previous MB to recode it
			img->current_mb_nr = FmoGetPreviousMBNr(img->current_mb_nr);
			if(img->current_mb_nr == -1 )   // The first MB of the slice group  is too big,
// which means it's impossible to encode picture using current slice bits restriction
			{
				printf (errortext, ET_SIZE, "Error encoding first MB with spcified parameter, bits of current MB may be too big");
				error (errortext, 300);
			}
		}
	}
}

void dct_luma_fractal_4_uv(int block4_x,int block4_y,int block16_x,int block16_y)//对分形的残差进行变换量化反量化反变换
{
	int diff[16],i,j,k,dummy,nonzero;
    int block16_x_n,block16_y_n,block_x,block_y;
	int pix_x,pix_y;
	
	block16_x_n=img->width/32;
	block16_y_n=img->height/32;
	
	pix_x=block16_x*16+block4_x*4;
	pix_y=block16_y*16+block4_y*4;
	
	img->pix_x=block16_x*16;
	img->pix_y=block16_y*16;
	
	block_x=block4_x*4;
	block_y=block4_y*4;
	UV_dct = 1;
	for (j=0; j<4; j++)
		for (  i=0; i<4; i++, k++)
		{
			img->mpr[i+block_x][j+block_y] = imgUV_rec[0][j+pix_y][i+pix_x];
			img->m7[i][j] = imgUV_org[0][pix_y+j][pix_x+i] - imgUV_rec[0][j+pix_y][i+pix_x];
		}
	dummy = 0;
	nonzero = dct_luma(block_x, block_y, &dummy, 1);

	UV_dct = 2;
	for (j=0; j<4; j++)
		for (  i=0; i<4; i++, k++)
		{
			img->mpr[i+block_x][j+block_y] = imgUV_rec[1][j+pix_y][i+pix_x];
			img->m7[i][j] = imgUV_org[1][pix_y+j][pix_x+i] - imgUV_rec[1][j+pix_y][i+pix_x];
		}
	dummy = 0;
	nonzero = dct_luma(block_x, block_y, &dummy, 1);
		
}

void dct_luma_fractal_16_uv(int block16_x,int block16_y)//对分形的残差进行变换量化反量化反变换
{
	int block4_x,block4_y,dummy,nonzero;
	int  mb_x, mb_y, i, j, k;
	int  b8, b4,level,run;
	int block16_x_n,block16_y_n,block_x,block_y;
	int pix_x,pix_y;
	
	for(block4_y=0;block4_y<4;block4_y++)
		for (block4_x=0;block4_x<4;block4_x++)
		{
			dct_luma_fractal_4_uv(block4_x,block4_y, block16_x,block16_y);
		}
}

void dct_luma_fractal_uv()
{
	int block16_x_n,block16_y_n,block16_x,block16_y;
	block16_x_n=img->width/32;//22
	block16_y_n=img->height/32;//18
   	for (block16_y=0;block16_y<block16_y_n;block16_y++)
		for (block16_x=0;block16_x<block16_x_n;block16_x++)
		{
			dct_luma_fractal_16_uv(block16_x,block16_y);
		}
}
void tran_show(int CurrentMbAddr,int mode)
{
	int i,j,k;
	if (1==mode)//Y分量
	{
		if (0==trans_Y[0][CurrentMbAddr].partition)//16*16划分
		{
			fprintf(trans_show_Y,"\nluma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_Y[0][CurrentMbAddr].partition,
				trans_Y[0][CurrentMbAddr].block_type,trans_Y[0][CurrentMbAddr].x,
			trans_Y[0][CurrentMbAddr].y,trans_Y[0][CurrentMbAddr].reference,trans_Y[0][CurrentMbAddr].offset/5,
			trans_Y[0][CurrentMbAddr].scale*20);
		}else if(1==trans_Y[0][CurrentMbAddr].partition || 2==trans_Y[0][CurrentMbAddr].partition)//16*8或8*16划分
		{
			for (i=0;i<2;i++)//一个16*16宏块包含两个矩形块
			{
				fprintf(trans_show_Y,"\nluma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_Y[0][CurrentMbAddr].next[i].partition,
					trans_Y[0][CurrentMbAddr].next[i].block_type,trans_Y[0][CurrentMbAddr].next[i].x,
					trans_Y[0][CurrentMbAddr].next[i].y,trans_Y[0][CurrentMbAddr].next[i].reference,trans_Y[0][CurrentMbAddr].next[i].offset/5,
					trans_Y[0][CurrentMbAddr].next[i].scale*20);
			}
		}else//8*8划分
		{
			for (i=0;i<4;i++)//一个宏块包含4个8*8方块
			{
			
				if (0==trans_Y[0][CurrentMbAddr].next[i].partition )//8*8划分
				{
					fprintf(trans_show_Y,"\nluma,mode=3.%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_Y[0][CurrentMbAddr].next[i].partition,
						trans_Y[0][CurrentMbAddr].next[i].block_type,trans_Y[0][CurrentMbAddr].next[i].x,
						trans_Y[0][CurrentMbAddr].next[i].y,trans_Y[0][CurrentMbAddr].next[i].reference,trans_Y[0][CurrentMbAddr].next[i].offset/5,
						trans_Y[0][CurrentMbAddr].next[i].scale*20);
				}else if (1==trans_Y[0][CurrentMbAddr].next[i].partition || 2==trans_Y[0][CurrentMbAddr].next[i].partition)//8*4或4*8划分
				{
					for (j=0;j<2;j++)//一个8*8包含两个矩形小块
					{
						fprintf(trans_show_Y,"\nluma,mode=3.%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_Y[0][CurrentMbAddr].next[i].partition,
							trans_Y[0][CurrentMbAddr].next[i].block_type,trans_Y[0][CurrentMbAddr].next[i].next[j].x,
							trans_Y[0][CurrentMbAddr].next[i].next[j].y,trans_Y[0][CurrentMbAddr].next[i].next[j].reference,trans_Y[0][CurrentMbAddr].next[i].next[j].offset/5,
							trans_Y[0][CurrentMbAddr].next[i].next[j].scale*20);
					}
				}else//4*4划分
				{
					for (k=0;k<4;k++)//一个8*8块包含4个4*4小方块
					{
						fprintf(trans_show_Y,"\nluma,mode=3.%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_Y[0][CurrentMbAddr].next[i].partition,
							trans_Y[0][CurrentMbAddr].next[i].block_type,trans_Y[0][CurrentMbAddr].next[i].next[k].x,
							trans_Y[0][CurrentMbAddr].next[i].next[k].y,trans_Y[0][CurrentMbAddr].next[i].next[k].reference,trans_Y[0][CurrentMbAddr].next[i].next[k].offset/5,
							trans_Y[0][CurrentMbAddr].next[i].next[k].scale*20);
					}

				}
			}
		}

	}else if (2==mode)//U分量
	{
		if (0==trans_U[0][CurrentMbAddr].partition)
		{
			fprintf(trans_show_UV,"\n2、chroma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_U[0][CurrentMbAddr].partition,
				trans_U[0][CurrentMbAddr].block_type,trans_U[0][CurrentMbAddr].x,
				trans_U[0][CurrentMbAddr].y,trans_U[0][CurrentMbAddr].reference,trans_U[0][CurrentMbAddr].offset/5,
				trans_U[0][CurrentMbAddr].scale*20);

		}else if(1==trans_U[0][CurrentMbAddr].partition || 2==trans_U[0][CurrentMbAddr].partition)
		{
			for (i=0;i<2;i++)
			{
				fprintf(trans_show_UV,"\n2、chroma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_U[0][CurrentMbAddr].partition,
					trans_U[0][CurrentMbAddr].block_type,trans_U[0][CurrentMbAddr].next[i].x,
					trans_U[0][CurrentMbAddr].next[i].y,trans_U[0][CurrentMbAddr].next[i].reference,trans_U[0][CurrentMbAddr].next[i].offset/5,
					trans_U[0][CurrentMbAddr].next[i].scale*20);
			}
		}else
		{
			for (i=0;i<4;i++)
			{
				fprintf(trans_show_UV,"\n2、chroma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_U[0][CurrentMbAddr].partition,
					trans_U[0][CurrentMbAddr].block_type,trans_U[0][CurrentMbAddr].next[i].x,
					trans_U[0][CurrentMbAddr].next[i].y,trans_U[0][CurrentMbAddr].next[i].reference,trans_U[0][CurrentMbAddr].next[i].offset/5,
					trans_U[0][CurrentMbAddr].next[i].scale*20);
			}
		}
	}else//V分量
	{
		if (0==trans_V[0][CurrentMbAddr].partition)
		{
			fprintf(trans_show_UV,"\n3、chroma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_V[0][CurrentMbAddr].partition,
				trans_V[0][CurrentMbAddr].block_type,trans_V[0][CurrentMbAddr].x,
				trans_V[0][CurrentMbAddr].y,trans_V[0][CurrentMbAddr].reference,trans_V[0][CurrentMbAddr].offset/5,
				trans_V[0][CurrentMbAddr].scale*20);
		}else if(1==trans_V[0][CurrentMbAddr].partition || 2==trans_V[0][CurrentMbAddr].partition)
		{
			for (i=0;i<2;i++)
			{
				fprintf(trans_show_UV,"\n3、chroma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_V[0][CurrentMbAddr].partition,
					trans_V[0][CurrentMbAddr].next[i].block_type,trans_V[0][CurrentMbAddr].next[i].x,
					trans_V[0][CurrentMbAddr].next[i].y,trans_V[0][CurrentMbAddr].next[i].reference,trans_V[0][CurrentMbAddr].next[i].offset/5,
					trans_V[0][CurrentMbAddr].next[i].scale*20);
			}
		}else
		{
			for (i=0;i<4;i++)
			{
				fprintf(trans_show_UV,"\n3、chroma,mode=%d\n%d %2d %2d  %d  %3.2f  %3.3f",trans_V[0][CurrentMbAddr].partition,
					trans_V[0][CurrentMbAddr].next[i].block_type,trans_V[0][CurrentMbAddr].next[i].x,
					trans_V[0][CurrentMbAddr].next[i].y,trans_V[0][CurrentMbAddr].next[i].reference,trans_V[0][CurrentMbAddr].next[i].offset/5,
					trans_V[0][CurrentMbAddr].next[i].scale*20);
			}
		}
	}

}
void encode_Oneframe()//分形编码一个P帧，Y,U,V分量单独编码
{
	int i,i1,j1,k,CurrentMbAddr;
	int block16_x_n,block16_y_n,block16_x,block16_y,j;
	int len;

	for (current_macroblock=0;current_macroblock<img->frmSizeInMbs;current_macroblock++)//循环每个宏块分形编码Y分量
	{ 

		CurrentMbAddr=img->current_mb_nr=current_macroblock;//当前宏块地址
// 		start_macroblock (current_macroblock, FALSE1);
// 		img->mb_data[img->current_mb_nr].mb_type=0;//P8*8
// 		p_encode_one_macroblock(current_macroblock,trans_Y,1);

		start_macroblock (CurrentMbAddr, FALSE1);//初始化当前宏块的MV,MVD,参考索引，滤波，bitcounter等
		encode_one_macroblock(current_macroblock,trans_Y,1);//分形编码一个宏块的Y分量，并将分形参数保存入trans
	}
	// 	printf("\nS_P=%d\n",S_P);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	for (current_macroblock=0;current_macroblock<img->frmSizeInMbs/4;current_macroblock++)//循环每个宏块分形编码U分量
	{
		encode_one_macroblock(current_macroblock,trans_U,2);
	}
	for (current_macroblock=0;current_macroblock<img->frmSizeInMbs/4;current_macroblock++)//循环每个宏块分形编码V分量
	{
		encode_one_macroblock(current_macroblock,trans_V,3);
	}
///张玲/////保存数据///////	
	if (img->current_frame==1)//是分形编码的第一帧
	{
		if ((trans_show_Y=fopen("trans_show_Y.txt","w"))!=NULL)//写Y分量的分形参数
		{
	    	fprintf(trans_show_Y,"\n////////////////++++++++++++++current_frame = %2d+++++++++++++++/////////////////\n",img->current_frame);
			for (CurrentMbAddr=0;CurrentMbAddr<img->FrameSizeInMbs;CurrentMbAddr++)//写每个宏块的参数
			{
				fprintf(trans_show_Y,"\nCurrentMb=%3d",CurrentMbAddr);//写当前宏块地址
				tran_show(CurrentMbAddr,1);//写Y分量的分形参数
			}
		}
		fclose(trans_show_Y);
		if ((trans_show_UV=fopen("trans_show_UV.txt","w"))!=NULL)//写UV分量的分形参数
		{
			fprintf(trans_show_UV,"\n//////////////+++++++++++++++current_frame = %2d+++++++++++++++/////////////////\n",img->current_frame);
			for (CurrentMbAddr=0;CurrentMbAddr<img->FrameSizeInMbs/4;CurrentMbAddr++)
			{
				fprintf(trans_show_UV,"\nCurrentMb=%3d",CurrentMbAddr);//写当前宏块地址
				tran_show(CurrentMbAddr,2);//写U分量的分形参数
				tran_show(CurrentMbAddr,3);//写V分量的分形参数
				
			}
		}
		fclose(trans_show_UV);
	}else//不是分形编码的第一帧
	{
		if ((trans_show_Y=fopen("trans_show_Y.txt","a"))!=NULL)
		{
			fprintf(trans_show_Y,"\n////////////////++++++++++++++current_frame = %2d+++++++++++++++/////////////////\n",img->current_frame);
			for (CurrentMbAddr=0;CurrentMbAddr<img->FrameSizeInMbs;CurrentMbAddr++)
			{
				fprintf(trans_show_Y,"\nCurrentMb=%3d",CurrentMbAddr);//写当前宏块地址
				tran_show(CurrentMbAddr,1);//写Y分量的分形参数
			}
		}
		fclose(trans_show_Y);
		if ((trans_show_UV=fopen("trans_show_UV.txt","a"))!=NULL)
		{
	    	fprintf(trans_show_UV,"\n//////////////+++++++++++++++current_frame = %2d+++++++++++++++/////////////////\n",img->current_frame);
			for (CurrentMbAddr=0;CurrentMbAddr<img->FrameSizeInMbs/4;CurrentMbAddr++)
			{
				fprintf(trans_show_UV,"\nCurrentMb=%3d",CurrentMbAddr);//写当前宏块地址
				tran_show(CurrentMbAddr,2);//写U分量的分形参数
				tran_show(CurrentMbAddr,3);//写V分量的分形参数
			}
		}
		fclose(trans_show_UV);

	}
		
		
//将解码分解到每块需要添加的东东。
// 	memcpy(imgY_ref_temp[0],imgY_rec[0],input->imageheight*input->imagewidth);
// 	memcpy(imgUV_ref_temp[0][0],imgUV_rec[0][0],input->imageheight*input->imagewidth/4);
// 	memcpy(imgUV_ref_temp[1][0],imgUV_rec[1][0],input->imageheight*input->imagewidth/4);



	num_all=num_16+num_8+num_4+num_rect1;
	//printf("  %d    %d    %d    %d/%d   %d",num_16,num_8,num_4,num_rect1,num_rect2,num_all);
	// 	printf("      (%d)   ",num_rl);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
}

void org_enc_show(int n)
{

	int i,j,mb_x,mb_y;
	mb_x=(n%img->width)*16;
	mb_y=(n/img->width)*16;
	printf("\nog(en) 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 ");
	for (j=mb_y;j<mb_y+16;j++)
	{
		printf("\n %2d ",j);
		for (i=mb_x;i<mb_x+16;i++)
		{

			printf(" %d",imgY_org[j][i]);
		}
		printf("\n(%2d)",j);
		for (i=mb_x;i<mb_x+16;i++)
		{
			
			printf(" %d",enc_picture->imgY[j][i]);
		}
	}
}

static int  writeout_picture(Picture *pic);

static int  picture_structure_decision(Picture *frame, Picture *top, Picture *bot);
static void distortion_fld (float *dis_fld_y, float *dis_fld_u, float *dis_fld_v);

static void field_mode_buffer(int bit_field, float snr_field_y, float snr_field_u, float snr_field_v);
static void frame_mode_buffer (int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v);



static void copy_motion_vectors_MB();

static void CopyFrameToOldImgOrgVariables (Sourceframe *sf);
static void CopyTopFieldToOldImgOrgVariables (Sourceframe *sf);
static void CopyBottomFieldToOldImgOrgVariables (Sourceframe *sf);
static Sourceframe *AllocSourceframe (int xs, int ys);
static void FreeSourceframe (Sourceframe *sf);
static void i_ReadOneFrame (int FrameNoInFile, int HeaderSize, int xs, int ys, Sourceframe *sf);
static void writeUnit(Bitstream* currStream ,int partition);


#ifdef _ADAPT_LAST_GROUP_
int *last_P_no;
int *last_P_no_frm;
int *last_P_no_fld;
#endif





static int CalculateFrameNumber();  // Calculates the next frame number
static int FrameNumberInFile;       // The current frame number in the input file
static Sourceframe *srcframe;


StorablePicture *enc_picture;
StorablePicture *enc_frame_picture;
StorablePicture *enc_top_picture;
StorablePicture *enc_bottom_picture;
//Rate control
int    QP;

const int ONE_FOURTH_TAP[3][2] =
{
	{20,20},
	{-5,-4},
	{ 1, 0},
};


void MbAffPostProc()
{
	byte temp[16][32];
	
	byte ** imgY  = enc_picture->imgY;
	byte ***imgUV = enc_picture->imgUV;
	
	int i, x, y, x0, y0, uv;
	for (i=0; i<(int)img->PicSizeInMbs; i+=2)
	{
		if (enc_picture->mb_field[i])
		{
			get_mb_pos(i, &x0, &y0);
			for (y=0; y<(2*MB_BLOCK_SIZE);y++)
				for (x=0; x<MB_BLOCK_SIZE; x++)
					temp[x][y] = imgY[y0+y][x0+x];
				
				for (y=0; y<MB_BLOCK_SIZE;y++)
					for (x=0; x<MB_BLOCK_SIZE; x++)
					{
						imgY[y0+(2*y)][x0+x]   = temp[x][y];
						imgY[y0+(2*y+1)][x0+x] = temp[x][y+MB_BLOCK_SIZE];
					}
					
					x0 = x0/2;
					y0 = y0/2;
					
					for (uv=0; uv<2; uv++)
					{
						for (y=0; y<(2*MB_BLOCK_SIZE/2);y++)
							for (x=0; x<MB_BLOCK_SIZE/2; x++)
								temp[x][y] = imgUV[uv][y0+y][x0+x];
							
							for (y=0; y<MB_BLOCK_SIZE/2;y++)
								for (x=0; x<MB_BLOCK_SIZE/2; x++)
								{
									imgUV[uv][y0+(2*y)][x0+x]   = temp[x][y];
									imgUV[uv][y0+(2*y+1)][x0+x] = temp[x][y+MB_BLOCK_SIZE/2];
								}
					}
		}
	}
}

/*!
************************************************************************
* \brief
*    Encodes a picture
*
*    This is the main picture coding loop.. It is called by all this
*    frame and field coding stuff after the img-> elements have been
*    set up.  Not sure whether it is useful for MB-adaptive frame/field
*    coding
************************************************************************
*/
void code_a_picture(Picture *pic)
{
	int NumberOfCodedMBs = 0;
	int SliceGroup = 0;
	int j;
	
	img->currentPicture = pic;
	img->currentPicture->idr_flag = ((!IMG_NUMBER) && (!(img->structure==BOTTOM_FIELD))) || (input->idr_enable && (img->type == I_SLICE || img->type==SP_SLICE || img->type==SI_SLICE)&& (!(img->structure==BOTTOM_FIELD)));//图像IDR标志位
	//1
	pic->no_slices = 0;
	pic->distortion_u = pic->distortion_v = pic->distortion_y = 0.0;
	
	// restrict list 1 size
	img->num_ref_idx_l0_active = max(1, (img->type==B_SLICE ? active_pps->num_ref_idx_l0_active_minus1 + 1: active_pps->num_ref_idx_l0_active_minus1 +1 )); //3
	img->num_ref_idx_l1_active = (img->type==B_SLICE ? active_pps->num_ref_idx_l1_active_minus1 + 1 : 0);//0
	
	
	// generate reference picture lists
	init_lists(img->type, img->structure);//(2,0)
	
	// assign list 0 size from list size
	img->num_ref_idx_l0_active = listXsize[0];//0
	img->num_ref_idx_l1_active = listXsize[1];//0
	
	//if (!img->MbaffFrameFlag)
	{
		if ((img->type == P_SLICE || img->type == SP_SLICE) && input->P_List0_refs)//no
		{
			img->num_ref_idx_l0_active = min(img->num_ref_idx_l0_active, input->P_List0_refs);
			listXsize[0] = min(listXsize[0], input->P_List0_refs);  
		}
		if (img->type == B_SLICE )//no
		{
			
			if (input->B_List0_refs)
			{
				img->num_ref_idx_l0_active = min(img->num_ref_idx_l0_active, input->B_List0_refs);
				listXsize[0] = min(listXsize[0], input->B_List0_refs);  
			}
			if (input->B_List1_refs)
			{
				
				img->num_ref_idx_l1_active = min(img->num_ref_idx_l1_active, input->B_List1_refs);
				listXsize[1] = min(listXsize[1], input->B_List1_refs);  
			}
		}
	} 
	
	
	//if (img->MbaffFrameFlag)
	if (img->structure==FRAME)//yes
		init_mbaff_lists();

	if (img->type != I_SLICE && (input->WeightedPrediction == 1 || (input->WeightedBiprediction > 0 && (img->type == B_SLICE))))//no
	{
		if (img->type==P_SLICE || img->type==SP_SLICE)
			estimate_weighting_factor_P_slice ();
		else
			estimate_weighting_factor_B_slice ();
	}
	
	
	RandomIntraNewPicture ();     //! Allocates forced INTRA MBs (even for fields!)
	
	// The slice_group_change_cycle can be changed here.
	// FmoInit() is called before coding each picture, frame or field
	img->slice_group_change_cycle=1;
	FmoInit(img, active_pps, active_sps);///////////////////////////////////与片组相关，确定宏块到片组的映射关系。
	
	FmoStartPicture ();           //! picture level initialization of FMO///得到每一个片中的第一个宏块地址。
	
	while (NumberOfCodedMBs < img->total_number_mb)   //0<396    // loop over slices以片为单位循环编码所有片
	{
		// Encode one SLice Group
		while (!FmoSliceGroupCompletelyCoded (SliceGroup)) //以片组为单位循环处理每个片组中的所有片数据
		{
			// Encode the current slice
			NumberOfCodedMBs += encode_one_slice (SliceGroup, pic);//重要的编码一个条带的函数!!!!!!!!!!!!!!!!!!!
			FmoSetLastMacroblockInSlice (img->current_mb_nr);
			// Proceed to next slice
			img->current_slice_nr++;
			stat->bit_slice = 0;
		}
		// Proceed to next SliceGroup
		SliceGroup++;
	}
	FmoEndPicture ();

////注意看经过不
	if (input->rdopt == 2 && (img->type != B_SLICE))
		for (j = 0; j < input->NoOfDecoders; j++)
			DeblockFrame (img, decs->decY_best[j], NULL);
		
	DeblockFrame (img, enc_picture->imgY, enc_picture->imgUV);
		
	if (img->MbaffFrameFlag)
		MbAffPostProc();
		
}



/*!
************************************************************************
* \brief
*    Encodes one frame
************************************************************************
*/
int encode_one_frame ()
{
	static int prev_frame_no = 0; // POC200301
	static int consecutive_non_reference_pictures = 0; // POC200301
	
#ifdef _LEAKYBUCKET_  //                             ZZZZZZZZZZZZZZZZZZZZZZZZZ
	extern long Bit_Buffer[10000];
	extern unsigned long total_frame_buffer;
#endif
	
	time_t ltime1;
	time_t ltime2;
	
#ifdef WIN32
	struct _timeb tstruct1;
	struct _timeb tstruct2;
#else
	struct timeb tstruct1;
	struct timeb tstruct2;
#endif
	
	int tmp_time;
	int bits_frm = 0, bits_fld = 0;
	float dis_frm = 0, dis_frm_y = 0, dis_frm_u = 0, dis_frm_v = 0;
	float dis_fld = 0, dis_fld_y = 0, dis_fld_u = 0, dis_fld_v = 0;
	
	//Rate control
	int pic_type, bits = 0; 
	
	me_time=0;
	
#ifdef WIN32
	_ftime (&tstruct1);           // start time ms
#else
	ftime (&tstruct1);
#endif
	time (&ltime1);               // start time s
	
	//Rate control 
	img->write_macroblock = 0;
	/*
	//Shankar Regunathan (Oct 2002)
	//Prepare Panscanrect SEI payload
	UpdatePanScanRectInfo ();
	//Prepare Arbitrarydata SEI Payload
	UpdateUser_data_unregistered ();
	//Prepare Registered data SEI Payload
	UpdateUser_data_registered_itu_t_t35 ();
	//Prepare RandomAccess SEI Payload
	UpdateRandomAccess ();
	*/
	
	put_buffer_frame ();      // sets the pointers to the frame structures 
                            	// (and not to one of the field structures)
	init_frame ();
	FrameNumberInFile = CalculateFrameNumber();//0
	
	srcframe = AllocSourceframe (img->width, img->height);
	i_ReadOneFrame (FrameNumberInFile, input->infile_header, img->width, img->height, srcframe);
	CopyFrameToOldImgOrgVariables (srcframe);
	
	// Set parameters for directmode and Deblocking filter
	img->direct_type     = input->direct_type;//1
	img->LFDisableIdc    = input->LFDisableIdc;//0
	img->LFAlphaC0Offset = input->LFAlphaC0Offset;//0
	img->LFBetaOffset    = input->LFBetaOffset;//0
	
	if (img->type == B_SLICE)//no
		Bframe_ctr++;         // Bframe_ctr only used for statistics, should go to stat->
	
	if (input->PicInterlace == FIELD_CODING)//no
	{
		//Rate control
		img->FieldControl=1;
		
		img->field_picture = 1;  // we encode fields
		field_picture (top_pic, bottom_pic);
		img->fld_flag = 1;
	}
	else//yes
	{
		//Rate control
		img->FieldControl=0;
		
		// For frame coding, turn MB level field/frame coding flag on
		if (input->MbInterlace)//no
			mb_adaptive = 1;
		
		img->field_picture = 0; // we encode a frame
		
		//Rate control
		if(input->RCEnable)//no
		{ 
		/*update the number of MBs in the basic unit for MB adaptive 
			f/f coding*/
			if((input->MbInterlace)&&(input->basicunit<img->Frame_Total_Number_MB)\
				&&(img->type==P_SLICE)&&(img->IFLAG==0))
				img->BasicUnit=input->basicunit*2;
			else
				img->BasicUnit=input->basicunit;
			
			rc_init_pict(1,0,1); 
			img->qp  = updateQuantizationParameter(0); 
			
			
			pic_type = img->type;
			QP =0;
		}
		
		if( active_sps->frame_mbs_only_flag)//yes
			img->TopFieldFlag=0;
		
		frame_picture (frame_pic);///////////////////////////////重要的编码一个I帧的函数！！！！！
		
		// For field coding, turn MB level field/frame coding flag off
		if (input->MbInterlace)
			mb_adaptive = 0;
		
		if (input->PicInterlace == ADAPTIVE_CODING)
		{
			//Rate control
			img->FieldControl=1;
			img->write_macroblock = 0;
			img->bot_MB = 0;
			
			img->field_picture = 1;  // we encode fields
			field_picture (top_pic, bottom_pic);
			
			//! Note: the distortion for a field coded picture is stored in the top field
			//! the distortion values in the bottom field are dummies
			dis_fld = top_pic->distortion_y + top_pic->distortion_u + top_pic->distortion_v;
			dis_frm = frame_pic->distortion_y + frame_pic->distortion_u + frame_pic->distortion_v;
			
			img->fld_flag = picture_structure_decision (frame_pic, top_pic, bottom_pic);
			update_field_frame_contexts (img->fld_flag);
			
			//Rate control
			if(img->fld_flag==0)
				img->FieldFrame=1;
			/*the current choice is field coding*/
			else
				img->FieldFrame=0;
		}
		else
			
			img->fld_flag = 0;
	}
	
	if (img->fld_flag)
		stat->bit_ctr_emulationprevention += stat->em_prev_bits_fld;
	else
		stat->bit_ctr_emulationprevention += stat->em_prev_bits_frm;
	
	if (img->type != B_SLICE)
	{
		img->pstruct_next_P = img->fld_flag;
	}
	
	// Here, img->structure may be either FRAME or BOTTOM FIELD depending on whether AFF coding is used
	// The picture structure decision changes really only the fld_flag
	
	if (img->fld_flag)            // field mode (use field when fld_flag=1 only)
	{
		field_mode_buffer (bits_fld, dis_fld_y, dis_fld_u, dis_fld_v);
		writeout_picture (top_pic);
		writeout_picture (bottom_pic);
	}
	else                          //frame mode
	{
		frame_mode_buffer (bits_frm, dis_frm_y, dis_frm_u, dis_frm_v);
		writeout_picture (frame_pic);
	}
	
	if (frame_pic)
		free_slice_list(frame_pic);
	if (top_pic)
		free_slice_list(top_pic);
	if (bottom_pic)
		free_slice_list(bottom_pic);
	
		/*
		// Tian Dong (Sept 2002)
		// in frame mode, the newly reconstructed frame has been inserted to the mem buffer
		// and it is time to prepare the spare picture SEI payload.
		if (input->InterlaceCodingOption == FRAME_CODING
		&& input->SparePictureOption && img->type != B_SLICE)
		CalculateSparePicture ();
	*/
	
	//Rate control
	if(input->RCEnable)
	{
		bits = stat->bit_ctr-stat->bit_ctr_n;
		rc_update_pict_frame(bits);
	}
	
	/*
    
	  if (input->InterlaceCodingOption == FRAME_CODING)
	  {
	  if (input->rdopt == 2 && img->type != B_SLICE)
      UpdateDecoders ();      // simulate packet losses and move decoded image to reference buffers
	  
		if (input->RestrictRef)
		UpdatePixelMap ();
		}
	*/
	
	find_snr ();//里面有重建图像
	
	time (&ltime2);               // end time sec
#ifdef WIN32
	_ftime (&tstruct2);           // end time ms
#else
	ftime (&tstruct2);            // end time ms
#endif
	
	tmp_time = (ltime2 * 1000 + tstruct2.millitm) - (ltime1 * 1000 + tstruct1.millitm);
	tot_time = tot_time + tmp_time;
	
	if (input->PicInterlace == ADAPTIVE_CODING)    //  =3   
	{
		if (img->fld_flag)
		{
			// store bottom field
			store_picture_in_dpb(enc_bottom_picture);
			free_storable_picture(enc_frame_picture);
		}
		else
		{
			// replace top with frame
			replace_top_pic_with_frame(enc_frame_picture);
			free_storable_picture(enc_bottom_picture);
		}
	}
	else
	{
		if (img->fld_flag)
		{
			store_picture_in_dpb(enc_bottom_picture);
		}
		else
		{
			store_picture_in_dpb(enc_frame_picture);
		}
	}
	
	
 #ifdef _LEAKYBUCKET_   //                  ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	// Store bits used for this frame and increment counter of no. of coded frames
	Bit_Buffer[total_frame_buffer] = stat->bit_ctr - stat->bit_ctr_n;
	total_frame_buffer++;
#endif
	
	// POC200301: Verify that POC coding type 2 is not used if more than one consecutive 
	// non-reference frame is requested or if decoding order is different from output order
	if (img->pic_order_cnt_type == 2)
	{
		if (!img->nal_reference_idc) consecutive_non_reference_pictures++;
		else consecutive_non_reference_pictures = 0;
		
		if (frame_no < prev_frame_no || consecutive_non_reference_pictures>1)
			error("POC type 2 cannot be applied for the coding pattern where the encoding /decoding order of pictures are different from the output order.\n", -1);
		prev_frame_no = frame_no;
	}
	
	if (stat->bit_ctr_parametersets_n!=0)
		ReportNALNonVLCBits(tmp_time, me_time);
	
	if (IMG_NUMBER == 0)
		ReportFirstframe(tmp_time,me_time);
    //ReportFirstframe(tmp_time);
	else
	{
		//Rate control
		if(input->RCEnable)
		{
			if((!input->PicInterlace)&&(!input->MbInterlace))
				bits=stat->bit_ctr-stat->bit_ctr_n;
			else
			{
				bits = stat->bit_ctr -Pprev_bits; // used for rate control update */
				Pprev_bits = stat->bit_ctr;
			}
		}
		
		switch (img->type)
		{
		case I_SLICE:
			stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
			ReportIntra(tmp_time,me_time);
			//ReportIntra(tmp_time);
			break;
		case SP_SLICE:
			stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
			ReportSP(tmp_time,me_time);
			//ReportSP(tmp_time);
			break;
		case B_SLICE:
			stat->bit_ctr_B += stat->bit_ctr - stat->bit_ctr_n;
			if (img->nal_reference_idc>0)
				ReportBS(tmp_time,me_time);
			//ReportBS(tmp_time);
			else
				ReportB(tmp_time,me_time);
			//ReportB(tmp_time);
			
			break;
		default:      // P
			stat->bit_ctr_P += stat->bit_ctr - stat->bit_ctr_n;
			ReportP(tmp_time,me_time);
			//ReportP(tmp_time);
		}
	}
	stat->bit_ctr_n = stat->bit_ctr;
	
	//Rate control
	if(input->RCEnable) 
	{
		rc_update_pict(bits);
		/*update the parameters of quadratic R-D model*/
		if((img->type==P_SLICE)&&(active_sps->frame_mbs_only_flag))
			updateRCModel();
		else if((img->type==P_SLICE)&&(!active_sps->frame_mbs_only_flag)\
			&&(img->IFLAG==0))
			updateRCModel();
	}
	
	stat->bit_ctr_parametersets_n=0;
	
	FreeSourceframe (srcframe);
	
	if (IMG_NUMBER == 0)
		return 0;
	else
		return 1;
}


/*!
************************************************************************
* \brief
*    Encodes a frame picture
************************************************************************
*/
void frame_picture (Picture *frame)
{
	
	img->structure = FRAME;//0
	img->PicSizeInMbs = img->FrameSizeInMbs;//396
	
	enc_frame_picture  = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);//分配内存
	img->ThisPOC=enc_frame_picture->poc=img->framepoc;//0
	enc_frame_picture->top_poc    = img->toppoc;//0
	enc_frame_picture->bottom_poc = img->bottompoc;//0
	
	enc_frame_picture->frame_poc = img->framepoc;//0
	
	enc_frame_picture->pic_num = img->frame_num;//0
	enc_frame_picture->coded_frame = 1;
	
	enc_frame_picture->MbaffFrameFlag = img->MbaffFrameFlag = (input->MbInterlace != FRAME_CODING);//0
	
	enc_picture=enc_frame_picture;
	
	stat->em_prev_bits_frm = 0;
	stat->em_prev_bits = &stat->em_prev_bits_frm;
	
	if (img->MbaffFrameFlag)//no
	{
		CopyTopFieldToOldImgOrgVariables (srcframe);
		CopyBottomFieldToOldImgOrgVariables (srcframe);
	}
	
	img->fld_flag = 0;
	code_a_picture(frame);/////////////////////////重要的编码一个I帧的函数
	
	frame->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);
	
	if (img->structure==FRAME)
	{
		find_distortion (snr, img);      
		frame->distortion_y = snr->snr_y;
		frame->distortion_u = snr->snr_u;
		frame->distortion_v = snr->snr_v;
	}
}


/*!
************************************************************************
* \brief
*    Encodes a field picture, consisting of top and bottom field
************************************************************************
*/
void field_picture (Picture *top, Picture *bottom)
{
	//Rate control
	int old_pic_type;              // picture type of top field used for rate control    
	int TopFieldBits;
	
	//Rate control
	old_pic_type = img->type;
	
	stat->em_prev_bits_fld = 0;
	stat->em_prev_bits = &stat->em_prev_bits_fld;
	img->current_frame *= 2;
	img->buf_cycle *= 2;
	img->height = input->imageheight / 2;
	img->height_cr = input->imageheight / 4;
	img->fld_flag = 1;
	img->PicSizeInMbs = img->FrameSizeInMbs/2;
	// Top field
	
	//  img->bottom_field_flag = 0;
	enc_top_picture  = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);
	enc_top_picture->poc=img->toppoc;
	enc_top_picture->frame_poc = img->toppoc;
	enc_top_picture->pic_num = img->frame_num;
	enc_top_picture->coded_frame = 0;
	enc_top_picture->MbaffFrameFlag = img->MbaffFrameFlag = FALSE1;
	img->ThisPOC = img->toppoc;
	
	img->structure = TOP_FIELD;
	enc_picture = enc_top_picture;
	put_buffer_top ();
	init_field ();
	if (img->type == B_SLICE)       //all I- and P-frames
		nextP_tr_fld--;
	
	CopyTopFieldToOldImgOrgVariables (srcframe);
	
	img->fld_flag = 1;
	//  img->bottom_field_flag = 0;
	
	//Rate control
	if(input->RCEnable)
	{
		img->BasicUnit=input->basicunit;
		
		if(input->PicInterlace==FIELD_CODING)
			rc_init_pict(0,1,1); 
		else
			rc_init_pict(0,1,0);
		
		img->qp  = updateQuantizationParameter(1); 
	}
	img->TopFieldFlag=1;
	
	code_a_picture(top_pic);
	enc_picture->structure = 1;
    
	store_picture_in_dpb(enc_top_picture);
	
	top->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);
	
	//Rate control
	TopFieldBits=top->bits_per_picture;
	
	//  Bottom field
	//  img->bottom_field_flag = 0;
	enc_bottom_picture  = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);
	enc_bottom_picture->poc=img->bottompoc;
	enc_bottom_picture->frame_poc = img->bottompoc;
	enc_bottom_picture->pic_num = img->frame_num;
	enc_bottom_picture->coded_frame = 0;
	enc_bottom_picture->MbaffFrameFlag = img->MbaffFrameFlag = FALSE1;
	img->ThisPOC = img->bottompoc;
	img->structure = BOTTOM_FIELD;
	enc_picture = enc_bottom_picture;
	put_buffer_bot ();
	img->current_frame++;
	
	init_field ();
	
	if (img->type == B_SLICE)       //all I- and P-frames
		nextP_tr_fld++;             //check once coding B field
	
	if (img->type == I_SLICE && input->IntraBottom!=1)
		img->type = P_SLICE;
	
	CopyBottomFieldToOldImgOrgVariables (srcframe);
	img->fld_flag = 1;
	//  img->bottom_field_flag = 1;
	
	//Rate control
	if(input->RCEnable)  setbitscount(TopFieldBits);
	if(input->RCEnable)
	{
		rc_init_pict(0,0,0); 
		img->qp  = updateQuantizationParameter(0); 
	}
	img->TopFieldFlag=0;
	
	enc_picture->structure = 2;
	code_a_picture(bottom_pic);
	
	bottom->bits_per_picture = 8 * ((((img->currentSlice)->partArr[0]).bitstream)->byte_pos);
	
	// the distortion for a field coded frame (consisting of top and bottom field)
	// lives in the top->distortion varaibles, thye bottom-> are dummies
	distortion_fld (&top->distortion_y, &top->distortion_u, &top->distortion_v);
	
}


/*!
************************************************************************
* \brief
*    This function write out a picture
* \return
*    0 if OK,                                                         \n
*    1 in case of error
*
************************************************************************
*/
static int writeout_picture(Picture *pic)
{
	Bitstream *currStream;
	int partition, slice;
	Slice *currSlice;

	img->currentPicture=pic;
	
	for (slice=0; slice<pic->no_slices; slice++)//1
	{
		currSlice = pic->slices[slice];
// 		printf("\nwriteout_picture:slice=%d\n",slice);
		for (partition=0; partition<currSlice->max_part_nr; partition++)//1
		{
			currStream = (currSlice->partArr[partition]).bitstream;
			assert (currStream->bits_to_go == 8);    //! should always be the case, the 
			//! byte alignment is done in terminate_slice
// 			printf("\nwriteout_picture:partition=%d, currStream->byte_pos=%d",partition,currSlice->partArr[partition].bitstream->byte_pos);
			writeUnit (currSlice->partArr[partition].bitstream,partition);
			
		}           // partition loop
	}           // slice loop
	return 0;   
}


/*!
************************************************************************
* \brief
*    Distortion Field
************************************************************************
*/
static void distortion_fld (float *dis_fld_y, float *dis_fld_u, float *dis_fld_v)
{
	
	img->current_frame /= 2;
	img->buf_cycle /= 2;
	img->height = input->imageheight;
	img->height_cr = input->imageheight / 2;
	img->total_number_mb =
		(img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
	
	combine_field ();
	
	imgY_org= imgY_org_frm;
	imgUV_org= imgUV_org_frm;
	
	find_distortion (snr, img);   // find snr from original frame picture
	
	*dis_fld_y = snr->snr_y;
	*dis_fld_u = snr->snr_u;
	*dis_fld_v = snr->snr_v;
}


/*!
************************************************************************
* \brief
*    Picture Structure Decision
************************************************************************
*/
static int picture_structure_decision (Picture *frame, Picture *top, Picture *bot)
{
	double lambda_picture;
	int spframe = (img->type == SP_SLICE);
	int bframe = (img->type == B_SLICE);
	float snr_frame, snr_field;
	int bit_frame, bit_field;
	
	lambda_picture = 0.85 * pow (2, (img->qp - SHIFT_QP) / 3.0) * (bframe || spframe ? 4 : 1);
	
	snr_frame = frame->distortion_y + frame->distortion_u + frame->distortion_v;
	//! all distrortions of a field picture are accumulated in the top field
	snr_field = top->distortion_y + top->distortion_u + top->distortion_v;
	bit_field = top->bits_per_picture + bot->bits_per_picture;
	bit_frame = frame->bits_per_picture;
	
	return decide_fld_frame (snr_frame, snr_field, bit_field, bit_frame, lambda_picture);
}


/*!
************************************************************************
* \brief
*    Field Mode Buffer
************************************************************************
*/
static void field_mode_buffer (int bit_field, float snr_field_y, float snr_field_u, float snr_field_v)
{
	put_buffer_frame ();
	
	snr->snr_y = snr_field_y;
	snr->snr_u = snr_field_u;
	snr->snr_v = snr_field_v;
}


/*!
************************************************************************
* \brief
*    Frame Mode Buffer
************************************************************************
*/
static void frame_mode_buffer (int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v)
{
	put_buffer_frame ();
	
	if ((input->PicInterlace != FRAME_CODING)||(input->MbInterlace != FRAME_CODING))
	{
		img->height = img->height / 2;
		img->height_cr = img->height_cr / 2;
		img->current_frame *= 2;
		
		put_buffer_top ();
		
		img->current_frame++;
		put_buffer_bot ();
		
		img->current_frame /= 2;         // reset the img->current_frame to field
		img->height = input->imageheight;
		img->height_cr = input->imageheight / 2;
		img->total_number_mb =
			(img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
		
		snr->snr_y = snr_frame_y;
		snr->snr_u = snr_frame_u;
		snr->snr_v = snr_frame_v;
		put_buffer_frame ();
		
	}
}


/*!
************************************************************************
* \brief
*    mmco initializations should go here
************************************************************************
*/
static void init_dec_ref_pic_marking_buffer()
{
	img->dec_ref_pic_marking_buffer=NULL;
}


/*!
************************************************************************
* \brief
*    Initializes the parameters for a new frame
************************************************************************
*/
static void init_frame ()
{
	int i;
	int prevP_no, nextP_no;
	
	last_P_no = last_P_no_frm;
	
	img->current_mb_nr = 0;
	img->current_slice_nr = 0;
	stat->bit_slice = 0;
	
	img->mb_y = img->mb_x = 0;
	img->block_y = img->pix_y = img->pix_c_y = 0; 
	img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;
	
	// The 'slice_nr' of each macroblock is set to -1 here, to guarantee the correct encoding 
	// with FMO (if no FMO, encoding is correct without following assignment), 
	// for which MBs may not be encoded with scan order
	for(i=0;i< ((img->width/MB_BLOCK_SIZE)*(img->height/MB_BLOCK_SIZE));i++)//将每一个宏块的slice_nr设置为i-1，有什么用呢？？？
		img->mb_data[i].slice_nr=-1;
	
	if (img->type != B_SLICE)//非B条带
	{
		img->tr = start_tr_in_this_IGOP + IMG_NUMBER * (input->jumpd + 1);//等于当前帧号
		
		img->imgtr_last_P_frm = img->imgtr_next_P_frm;//0
		img->imgtr_next_P_frm = img->tr;//将当前帧号赋给下一个P帧？？？
		
#ifdef _ADAPT_LAST_GROUP_
		if (input->last_frame && img->current_frame + 1 == input->no_frames_h264)//no
			img->tr = input->last_frame;
#endif
		
		if (IMG_NUMBER != 0 && input->successive_Bframe != 0)     // B pictures to encode //no
			nextP_tr_frm = img->tr;
		
		
		if(!input->RCEnable)                  // 非码率控制
		{
			if (img->type == I_SLICE)//I条带
#ifdef _CHANGE_QP_
				if (input->qp2start > 0 && img->tr >= input->qp2start)
					img->qp = input->qp02;
				else
#endif    
					img->qp = input->qp0;   // set quant. parameter for I-frame
			else//非I条带
			{
#ifdef _CHANGE_QP_
					if (input->qp2start > 0 && img->tr >= input->qp2start)//no
						img->qp = input->qpN2;
					else
#endif
						img->qp = input->qpN;//输入的QP
					
					if (img->type == SP_SLICE)//SP条带
					{
						img->qp = input->qpsp;
						img->qpsp = input->qpsp_pred;
					}   
			}
		}
		
		img->mb_y_intra = img->mb_y_upd;  //  img->mb_y_intra indicates which GOB to intra code for this frame
		
		if (input->intra_upd > 0) // if error robustness, find next GOB to update //no
		{
			img->mb_y_upd = (IMG_NUMBER / input->intra_upd) % (img->height / MB_BLOCK_SIZE);
		}
	}//end of 非B条带
	else
	{
		img->p_interval = input->jumpd + 1;
		prevP_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * img->p_interval;
		nextP_no = start_tr_in_this_IGOP + (IMG_NUMBER) * img->p_interval;
		
#ifdef _ADAPT_LAST_GROUP_
		last_P_no[0] = prevP_no;
		for (i = 1; i < img->buf_cycle; i++)
			last_P_no[i] = last_P_no[i - 1] - img->p_interval;
		
		if (input->last_frame && img->current_frame + 1 == input->no_frames_h264)
		{
			nextP_no = input->last_frame;
			img->p_interval = nextP_no - prevP_no;
		}
#endif
		
		img->b_interval =
			(int) ((float) (input->jumpd + 1) / (input->successive_Bframe + 1.0) +
			0.49999);
		
		img->tr = prevP_no + img->b_interval * img->b_frame_to_code;      // from prev_P
		
		if (img->tr >= nextP_no)
			img->tr = nextP_no - 1;
		//Rate control
		if(!input->RCEnable)                  // without using rate control
		{    
#ifdef _CHANGE_QP_
			if (input->qp2start > 0 && img->tr >= input->qp2start)
				img->qp = input->qpB2;
			else
#endif
				img->qp = input->qpB;
		}
		
	}
	
// 	UpdateSubseqInfo (img->layer);        // Tian Dong (Sept 2002)                注释掉了
// 	UpdateSceneInformation (0, 0, 0, -1); // JVT-D099, scene information SEI, nothing included by default   注释掉了
	
	//! Commented out by StW, needs fixing in SEI.h to keep the trace file clean
	//  PrepareAggregationSEIMessage ();
	
	img->total_number_mb = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);//宏块总数
	
	img->no_output_of_prior_pics_flag = 0;
	img->long_term_reference_flag = 0;
	
	init_dec_ref_pic_marking_buffer();
}

/*!
************************************************************************
* \brief
*    Initializes the parameters for a new field
************************************************************************
*/
static void init_field ()
{
	int i;
	int prevP_no, nextP_no;
	
	last_P_no = last_P_no_fld;
	
	img->current_mb_nr = 0;
	img->current_slice_nr = 0;
	stat->bit_slice = 0;
	
	input->jumpd *= 2;
	input->successive_Bframe *= 2;
	img->current_frame /= 2;
	img->buf_cycle /= 2;
	
	img->mb_y = img->mb_x = 0;
	img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
	img->block_x = img->pix_x = img->block_c_x = img->pix_c_x = 0;        // define horizontal positions
	
	if (img->type != B_SLICE)
    {
		img->tr = img->current_frame * (input->jumpd + 2) + img->fld_type;
		
		if (!img->fld_type)
        {
			img->imgtr_last_P_fld = img->imgtr_next_P_fld;
			img->imgtr_next_P_fld = img->tr;
        }
		
#ifdef _ADAPT_LAST_GROUP_
		if (input->last_frame && img->current_frame + 1 == input->no_frames_h264)
			img->tr = input->last_frame;
#endif
		if (img->current_frame != 0 && input->successive_Bframe != 0)    // B pictures to encode
			nextP_tr_fld = img->tr;
		
		//Rate control
		if(!input->RCEnable)                  // without using rate control
		{
			if (img->type == I_SLICE)
#ifdef _CHANGE_QP_
				if (input->qp2start > 0 && img->tr >= input->qp2start)
					img->qp = input->qp02;
				else
#endif    
					img->qp = input->qp0;   // set quant. parameter for I-frame
				else
				{
#ifdef _CHANGE_QP_
					if (input->qp2start > 0 && img->tr >= input->qp2start)
						img->qp = input->qpN2;
					else
#endif
						img->qp = input->qpN;
					if (img->type == SP_SLICE)
					{
						img->qp = input->qpsp;
						img->qpsp = input->qpsp_pred;
					}
				}
		}
		
		img->mb_y_intra = img->mb_y_upd;  //  img->mb_y_intra indicates which GOB to intra code for this frame
		
		if (input->intra_upd > 0) // if error robustness, find next GOB to update
        {
			img->mb_y_upd =
				(img->current_frame / input->intra_upd) % (img->width / MB_BLOCK_SIZE);
        }
    }
	else
    {
		img->p_interval = input->jumpd + 2;
		prevP_no = (img->current_frame - 1) * img->p_interval + img->fld_type;
		nextP_no = img->current_frame * img->p_interval + img->fld_type;
#ifdef _ADAPT_LAST_GROUP_
		if (!img->fld_type)       // top field
        {
			last_P_no[0] = prevP_no + 1;
			last_P_no[1] = prevP_no;
			for (i = 1; i <= img->buf_cycle; i++)
            {
				last_P_no[2 * i] = last_P_no[2 * i - 2] - img->p_interval;
				last_P_no[2 * i + 1] = last_P_no[2 * i - 1] - img->p_interval;
            }
        }
		else                      // bottom field
        {
			last_P_no[0] = nextP_no - 1;
			last_P_no[1] = prevP_no;
			for (i = 1; i <= img->buf_cycle; i++)
            {
				last_P_no[2 * i] = last_P_no[2 * i - 2] - img->p_interval;
				last_P_no[2 * i + 1] = last_P_no[2 * i - 1] - img->p_interval;
            }
        }
		
		if (input->last_frame && img->current_frame + 1 == input->no_frames_h264)
        {
			nextP_no = input->last_frame;
			img->p_interval = nextP_no - prevP_no;
        }
#endif
		
		img->b_interval =
			(int) ((float) (input->jumpd + 1) / (input->successive_Bframe + 1.0) + 0.49999);
		
		img->tr = prevP_no + (img->b_interval + 1) * img->b_frame_to_code;        // from prev_P
		if (img->tr >= nextP_no)
			img->tr = nextP_no - 1; // ?????
		//Rate control
		if(!input->RCEnable)                  // without using rate control
		{
#ifdef _CHANGE_QP_
			if (input->qp2start > 0 && img->tr >= input->qp2start)
				img->qp = input->qpB2;
			else
#endif
				img->qp = input->qpB;
		}
    }
	input->jumpd /= 2;
	input->successive_Bframe /= 2;
	img->buf_cycle *= 2;
	img->current_frame = 2 * img->current_frame + img->fld_type;
	img->total_number_mb = (img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);
}


#define Clip(min,max,val) (((val)<(min))?(min):(((val)>(max))?(max):(val)))
/*!
************************************************************************
* \brief
*    Generate Full Pel Representation
************************************************************************
*/
static void GenerateFullPelRepresentation (pel_t ** Fourthpel,
                                           pel_t * Fullpel, int xsize,
                                           int ysize)
{
	int x, y;
	
	for (y = 0; y < ysize; y++)
		for (x = 0; x < xsize; x++)
			PutPel_11 (Fullpel, y, x, FastPelY_14 (Fourthpel, y * 4, x * 4, ysize, xsize), xsize);
}


/*!
************************************************************************
* \brief
*    Upsample 4 times, store them in out4x.  Color is simply copied
*
* \par Input:
*    srcy, srcu, srcv, out4y, out4u, out4v
*
* \par Side Effects_
*    Uses (writes) img4Y_tmp.  This should be moved to a static variable
*    in this module
************************************************************************/
void UnifiedOneForthPix (StorablePicture *s)//产生１／４像素
{
	int is;
	int i, j, j4;
	int ie2, je2, jj, maxy;
	
	byte **out4Y;
	byte  *ref11;
	byte  **imgY = s->imgY;
	
	int img_width =s->size_x;
	int img_height=s->size_y;
	
	// don't upsample twice
	if (s->imgY_ups || s->imgY_11)
		return;
	
	s->imgY_11 = malloc ((s->size_x * s->size_y) * sizeof (byte));
	if (NULL == s->imgY_11)
		no_mem_exit("alloc_storable_picture: s->imgY_11");
	
	get_mem2D (&(s->imgY_ups), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);//图象边缘扩展并水平垂直各上采样4倍
	
	if (input->WeightedPrediction || input->WeightedBiprediction)//加权预测
	{
		s->imgY_11_w = malloc ((s->size_x * s->size_y) * sizeof (byte));
		get_mem2D (&(s->imgY_ups_w), (2*IMG_PAD_SIZE + s->size_y)*4, (2*IMG_PAD_SIZE + s->size_x)*4);
	}
	out4Y = s->imgY_ups;
	ref11 = s->imgY_11;
	//水平位置内插１／２像素
	for (j = -IMG_PAD_SIZE; j < s->size_y + IMG_PAD_SIZE; j++)
	{
		for (i = -IMG_PAD_SIZE; i < s->size_x + IMG_PAD_SIZE; i++)
		{
			jj = max (0, min (s->size_y - 1, j));
			// 	  ONE_FOURTH_TAP 2维数组ONE_FOURTH_TAP[0][1]ONE_FOURTH_TAP[1][2]ONE_FOURTH_TAP[2][2]无用
			// 		  ONE_FOURTH_TAP[0][0]=20 ONE_FOURTH_TAP[1][0]=-5 ONE_FOURTH_TAP[2][0]=1对应标准求半像素公式(E-5F+20G+20H-5I+J)/32的系数
			// 		  E F GbH I J其中b是半像素位置，其他大写字母是整数像素位置，依据距离远近加权，得到预测系数，当然我认为这些tap系数可以修改，有兴趣自己研究
			is =
				(ONE_FOURTH_TAP[0][0] *               //20G+20H G=(jj,i) H=(jj,i-1)
				(imgY[jj][max (0, min (s->size_x - 1, i))] +
                imgY[jj][max (0, min (s->size_x - 1, i + 1))]) +
				ONE_FOURTH_TAP[1][0] *                //-5F-5I
				(imgY[jj][max (0, min (s->size_x - 1, i - 1))] +
                imgY[jj][max (0, min (s->size_x - 1, i + 2))]) +
				ONE_FOURTH_TAP[2][0] *                 //E+J
				(imgY[jj][max (0, min (s->size_x - 1, i - 2))] +
                imgY[jj][max (0, min (s->size_x - 1, i + 3))]));
            img4Y_tmp[j + IMG_PAD_SIZE][(i + IMG_PAD_SIZE) * 2] = imgY[jj][max (0, min (s->size_x - 1, i))] * 1024;    // 整数像素1/1 pix pos
			img4Y_tmp[j + IMG_PAD_SIZE][(i + IMG_PAD_SIZE) * 2 + 1] = is * 32;  // 1/2 pix pos这里本来该除32的，却乘了32,定标法而已　看后面除以1024就明白
		}
	}
	//垂直位置内插１／２像素　参照以上，联系下DA的超采样问题，这里依然可能有内插导致频率混叠问题，所以加权也是低通滤波问题，有兴趣继续研究，点到为止
	for (i = 0; i < (s->size_x + 2 * IMG_PAD_SIZE) * 2; i++)
	{
		for (j = 0; j < s->size_y + 2 * IMG_PAD_SIZE; j++)
		{
			j4 = j * 4;
			maxy = s->size_y + 2 * IMG_PAD_SIZE - 1;
			// change for TML4, use 6 TAP vertical filter
			is =
				(ONE_FOURTH_TAP[0][0] *
				(img4Y_tmp[j][i] + img4Y_tmp[min (maxy, j + 1)][i]) +
				ONE_FOURTH_TAP[1][0] * (img4Y_tmp[max (0, j - 1)][i] +
				img4Y_tmp[min (maxy, j + 2)][i]) +
				ONE_FOURTH_TAP[2][0] * (img4Y_tmp[max (0, j - 2)][i] +
				img4Y_tmp[min (maxy, j + 3)][i])) / 32;
			
			PutPel_14 (out4Y, (j - IMG_PAD_SIZE) * 4, (i - IMG_PAD_SIZE * 2) * 2, (pel_t) max (0, min (255, (int) ((img4Y_tmp[j][i] + 512) / 1024))));  // 1/2 pix
			PutPel_14 (out4Y, (j - IMG_PAD_SIZE) * 4 + 2, (i - IMG_PAD_SIZE * 2) * 2, (pel_t) max (0, min (255, (int) ((is + 512) / 1024))));   // 1/2 pix
		}
	}
	
	/* 1/4 pix *///1/4像素内插 a=(G+b)/2 G整数像素点　b半像素点
	/* luma */
	ie2 = (s->size_x + 2 * IMG_PAD_SIZE - 1) * 4;
	je2 = (s->size_y + 2 * IMG_PAD_SIZE - 1) * 4;
	
	for (j = 0; j < je2 + 4; j += 2)
		for (i = 0; i < ie2 + 3; i += 2)
		{
			/*  '-'  *///水平内插
			PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4, i - IMG_PAD_SIZE * 4 + 1,//- IMG_PAD_SIZE * 4基于减去填充考虑得到在图象内的点
				(pel_t) (max (0, min (255,(int) (FastPelY_14 (out4Y, j - IMG_PAD_SIZE * 4,//G=(i - IMG_PAD_SIZE * 4,j - IMG_PAD_SIZE * 4)
				i - IMG_PAD_SIZE * 4, img_height, img_width) + FastPelY_14 (out4Y,
				j - IMG_PAD_SIZE * 4, min (ie2 + 2, i + 2) - IMG_PAD_SIZE * 4, img_height, img_width)+1) / 2))));//b=(i + 2- IMG_PAD_SIZE * 4,j - IMG_PAD_SIZE * 4)
		}
		for (i = 0; i < ie2 + 4; i++)
		{
			for (j = 0; j < je2 + 3; j += 2)
			{
				if (i % 2 == 0)
				{
					/*  '|'  *///垂直内插
					PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4 + 1, i - IMG_PAD_SIZE * 4,
						(pel_t) (max (0, min (255, (int) (FastPelY_14 (out4Y, j - IMG_PAD_SIZE * 4,//G=(i - IMG_PAD_SIZE * 4,j - IMG_PAD_SIZE * 4)
						i - IMG_PAD_SIZE * 4, img_height, img_width) + FastPelY_14 (out4Y,
						min (je2 + 2, j + 2) - IMG_PAD_SIZE * 4, i - IMG_PAD_SIZE * 4, img_height, img_width)+1) / 2))));//b=(i - IMG_PAD_SIZE * 4,j + 2 - IMG_PAD_SIZE * 4)
				}
				else if ((j % 4 == 0 && i % 4 == 1) || (j % 4 == 2 && i % 4 == 3))
				{
					/*  '/'  *///４５度内插
					PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4 + 1, i - IMG_PAD_SIZE * 4,
						(pel_t) (max (0, min (255, (int) (FastPelY_14 (out4Y, j - IMG_PAD_SIZE * 4,
						min (ie2 + 2, i + 1) - IMG_PAD_SIZE * 4, img_height, img_width) + FastPelY_14 (out4Y,
						min (je2 + 2, j + 2) - IMG_PAD_SIZE * 4, i - IMG_PAD_SIZE * 4 - 1, img_height, img_width) + 1) / 2))));
				}
				else
				{
					/*  '\'  *///反４５度内插
					PutPel_14 (out4Y, j - IMG_PAD_SIZE * 4 + 1, i - IMG_PAD_SIZE * 4,
						(pel_t) (max (0, min (255, (int) (FastPelY_14 (out4Y, j - IMG_PAD_SIZE * 4,
						i - IMG_PAD_SIZE * 4 - 1, img_height, img_width) + FastPelY_14 (out4Y,
						min (je2 + 2, j + 2) - IMG_PAD_SIZE * 4, 
						min (ie2 + 2, i + 1) - IMG_PAD_SIZE * 4, img_height, img_width) + 1) / 2))));
				}
			}
		}
		
		/*  Chroma: *///色度只是复制整数点，无内插 
					  /*    for (j = 0; j < img->height_cr; j++)
					  {
					  memcpy (outU[j], imgU[j], img->width_cr);       // just copy 1/1 pix, interpolate "online" 
					  memcpy (outV[j], imgV[j], img->width_cr);
					  }
		*/
		// Generate 1/1th pel representation (used for integer pel MV search)
		GenerateFullPelRepresentation (out4Y, ref11, s->size_x, s->size_y);//吧out4Y复制给ref11
		
}


/*!
************************************************************************
* \brief
*    Find SNR for all three components
************************************************************************
*/
static void find_snr ()
{
	int i, j;
	int diff_y, diff_u, diff_v;
	int impix;
	
	//  Calculate  PSNR for Y, U and V.
	
	//     Luma.
	impix = img->height * img->width;
	
	if (img->fld_flag != 0)//no
	{
		
		diff_y = 0;
		for (i = 0; i < img->width; ++i)
		{
			for (j = 0; j < img->height; ++j)
			{
				diff_y += img->quad[imgY_org[j][i] - imgY_com[j][i]];
			}
		}
		
		//     Chroma.
		diff_u = 0;
		diff_v = 0;
		
		for (i = 0; i < img->width_cr; i++)
		{
			for (j = 0; j < img->height_cr; j++)
			{
				diff_u += img->quad[imgUV_org[0][j][i] - imgUV_com[0][j][i]];
				diff_v += img->quad[imgUV_org[1][j][i] - imgUV_com[1][j][i]];
			}
		}
	}
	else//yes
	{ 
		imgY_org = imgY_org_frm;
		imgUV_org = imgUV_org_frm;
		
		if(input->PicInterlace==ADAPTIVE_CODING)
		{
			enc_picture = enc_frame_picture;
		}  
		
		diff_y = 0;
		for (i = 0; i < img->width; ++i)
		{
			for (j = 0; j < img->height; ++j)
			{
				diff_y += img->quad[imgY_org[j][i] - enc_picture->imgY[j][i]];//enc_picture->imgY[j][i]是重建图像， img->quad是平方的意思
			}
		}
		
		//     Chroma.
		diff_u = 0;
		diff_v = 0;
		
		for (i = 0; i < img->width_cr; i++)
		{
			for (j = 0; j < img->height_cr; j++)
			{
				diff_u += img->quad[imgUV_org[0][j][i] - enc_picture->imgUV[0][j][i]];
				diff_v += img->quad[imgUV_org[1][j][i] - enc_picture->imgUV[1][j][i]];
			}
		}
	}
#if ZEROSNR
	if (diff_y == 0)
		diff_y = 1;
	if (diff_u == 0)
		diff_u = 1;
	if (diff_v == 0)
		diff_v = 1; 
#endif
	
	//  Collecting SNR statistics
	if (diff_y != 0)
	{
		snr->snr_y = (float) (10 * log10 (65025 * (float) impix / (float) diff_y));         // luma snr for current frame
		snr->snr_u = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_u)));   // u croma snr for current frame, 1/4 of luma samples
		snr->snr_v = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_v)));   // v croma snr for current frame, 1/4 of luma samples
	}
	
	if (img->current_frame == 0)
	{
		snr->snr_y1 = (float) (10 * log10 (65025 * (float) impix / (float) diff_y));        // keep luma snr for first frame
		snr->snr_u1 = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_u)));  // keep croma u snr for first frame
		snr->snr_v1 = (float) (10 * log10 (65025 * (float) impix / (float) (4 * diff_v)));  // keep croma v snr for first frame
		snr->snr_ya = snr->snr_y1;
		snr->snr_ua = snr->snr_u1;
		snr->snr_va = snr->snr_v1;
	}
	// B pictures
	else
	{
		snr->snr_ya = (float) (snr->snr_ya * (img->current_frame + Bframe_ctr) + snr->snr_y) / (img->current_frame + Bframe_ctr + 1); // average snr lume for all frames inc. first
		snr->snr_ua = (float) (snr->snr_ua * (img->current_frame + Bframe_ctr) + snr->snr_u) / (img->current_frame + Bframe_ctr + 1); // average snr u croma for all frames inc. first
		snr->snr_va = (float) (snr->snr_va * (img->current_frame + Bframe_ctr) + snr->snr_v) / (img->current_frame + Bframe_ctr + 1); // average snr v croma for all frames inc. first
		
	}
}

/*!
************************************************************************
* \brief
*    Find distortion for all three components
************************************************************************
*/
static void find_distortion ()
{
	int i, j;
	int diff_y, diff_u, diff_v;
	int impix;
	
	//  Calculate  PSNR for Y, U and V.
	
	//     Luma.
	impix = img->height * img->width;
	
	if (img->structure!=FRAME)
	{
		
		diff_y = 0;
		for (i = 0; i < img->width; ++i)
		{
			for (j = 0; j < img->height; ++j)
			{
				diff_y += img->quad[abs (imgY_org[j][i] - imgY_com[j][i])];
			}
		}
		
		//     Chroma.
		
		diff_u = 0;
		diff_v = 0;
		
		for (i = 0; i < img->width_cr; i++)
		{
			for (j = 0; j < img->height_cr; j++)
			{
				diff_u += img->quad[abs (imgUV_org[0][j][i] - imgUV_com[0][j][i])];
				diff_v += img->quad[abs (imgUV_org[1][j][i] - imgUV_com[1][j][i])];
			}
		}
	}else
	{
		imgY_org   = imgY_org_frm;
		imgUV_org = imgUV_org_frm;
		
		diff_y = 0;
		for (i = 0; i < img->width; ++i)
		{
			for (j = 0; j < img->height; ++j)
			{
				diff_y += img->quad[abs (imgY_org[j][i] - enc_picture->imgY[j][i])];
			}
		}
		
		//     Chroma.
		
		diff_u = 0;
		diff_v = 0;
		
		for (i = 0; i < img->width_cr; i++)
		{
			for (j = 0; j < img->height_cr; j++)
			{
				diff_u += img->quad[abs (imgUV_org[0][j][i] - enc_picture->imgUV[0][j][i])];
				diff_v += img->quad[abs (imgUV_org[1][j][i] - enc_picture->imgUV[1][j][i])];
			}
		}
	}
	// Calculate real PSNR at find_snr()
	snr->snr_y= (float)diff_y;
	snr->snr_u= (float)diff_u;
	snr->snr_v= (float)diff_v;
}


/*!
 ************************************************************************
 * \brief
 *    Just a placebo
 ************************************************************************
 */
Boolean dummy_slice_too_big (int bits_slice)
{
  return FALSE1;
}





/*! 
***************************************************************************
// For MB level field/frame coding
***************************************************************************
*/
void copy_rdopt_data (int bot_block)
{
	int mb_nr = img->current_mb_nr;
	Macroblock *currMB = &img->mb_data[mb_nr];
	int i, j, k, l;
	
	int bframe = (img->type == B_SLICE);
	int mode;
	int b8mode, b8pdir;
	
	int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
	
	mode             = rdopt->mode;
	currMB->mb_type  = rdopt->mb_type;   // copy mb_type 
	currMB->cbp      = rdopt->cbp;   // copy cbp
	currMB->cbp_blk  = rdopt->cbp_blk;   // copy cbp_blk
	img->i16offset   = rdopt->i16offset;
	
	currMB->prev_qp=rdopt->prev_qp;
	currMB->prev_delta_qp=rdopt->prev_delta_qp;
	currMB->qp=rdopt->qp;
	
	currMB->c_ipred_mode = rdopt->c_ipred_mode;
	
	for (i = 0; i < 6; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 2; k++)
				for (l = 0; l < 18; l++)
					img->cofAC[i][j][k][l] = rdopt->cofAC[i][j][k][l];
				
				for (i = 0; i < 3; i++)
					for (k = 0; k < 2; k++)
						for (l = 0; l < 18; l++)
							img->cofDC[i][k][l] = rdopt->cofDC[i][k][l];
						
						for (j = 0; j < 4; j++)
							for (i = 0; i < 4; i++)
							{
								enc_picture->ref_idx[LIST_0][img->block_x + i][img->block_y + j] = rdopt->refar[LIST_0][j][i];
								enc_picture->ref_pic_id [LIST_0][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]];
								if (bframe)
								{
									enc_picture->ref_idx[LIST_1][img->block_x + i][img->block_y + j] = rdopt->refar[LIST_1][j][i];
									enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j]];
								}
							}
							
							//===== reconstruction values =====
							for (j = 0; j < 16; j++)
								for (i = 0; i < 16; i++)
								{
									enc_picture->imgY[img->pix_y + j][img->pix_x + i] = rdopt->rec_mbY[j][i];
								}
								
								
								for (j = 0; j < 8; j++)
									for (i = 0; i < 8; i++)
									{
										enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x + i] = rdopt->rec_mbU[j][i];
										enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x + i] = rdopt->rec_mbV[j][i];
									}
									
									
									for (i = 0; i < 4; i++)
									{
										currMB->b8mode[i] = rdopt->b8mode[i];
										currMB->b8pdir[i] = rdopt->b8pdir[i];
									}
									
									//==== intra prediction modes ====
									if (mode == P8x8)
									{
										for (k = 0, j = img->block_y; j < img->block_y + 4; j++)
											for (i = img->block_x; i < img->block_x + 4; i++, k++)
											{
												img->ipredmode[i][j]        = rdopt->ipredmode[i][j];
												currMB->intra_pred_modes[k] = rdopt->intra_pred_modes[k];
											}
									}
									else if (mode != I4MB)
									{
										for (k = 0, j = img->block_y; j < img->block_y + 4; j++)
											for (     i = img->block_x; i < img->block_x + 4; i++, k++)
											{
												img->ipredmode[i][j]        = DC_PRED;
												currMB->intra_pred_modes[k] = DC_PRED;
											}
									}
									else if (mode == I4MB)
									{
										for (k = 0, j = img->block_y; j < img->block_y + 4; j++)
											for (     i = img->block_x; i < img->block_x + 4; i++, k++)
											{
												img->ipredmode[i][j]        = rdopt->ipredmode[i][j];
												currMB->intra_pred_modes[k] = rdopt->intra_pred_modes[k];
												
											}
											
									}
									
									if (img->MbaffFrameFlag)
									{
										// motion vectors
										copy_motion_vectors_MB ();
										
										
										if (!IS_INTRA(currMB))
										{
											for (j = 0; j < 4; j++)
												for (i = 0; i < 4; i++)
												{
													b8mode = currMB->b8mode[i/2+2*(j/2)];
													b8pdir = currMB->b8pdir[i/2+2*(j/2)];
													
													if (b8pdir!=1)
													{
														enc_picture->mv[LIST_0][i+img->block_x][j+img->block_y][0] = rdopt->all_mv[i][j][LIST_0][rdopt->refar[LIST_0][j][i]][b8mode][0];
														enc_picture->mv[LIST_0][i+img->block_x][j+img->block_y][1] = rdopt->all_mv[i][j][LIST_0][rdopt->refar[LIST_0][j][i]][b8mode][1];
													}
													else
													{
														enc_picture->mv[LIST_0][i+img->block_x][j+img->block_y][0] = 0;
														enc_picture->mv[LIST_0][i+img->block_x][j+img->block_y][1] = 0;
													}
													if (bframe)
													{
														if (b8pdir!=0)
														{
															enc_picture->mv[LIST_1][i+img->block_x][j+img->block_y][0] = rdopt->all_mv[i][j][LIST_1][rdopt->refar[LIST_1][j][i]][b8mode][0];
															enc_picture->mv[LIST_1][i+img->block_x][j+img->block_y][1] = rdopt->all_mv[i][j][LIST_1][rdopt->refar[LIST_1][j][i]][b8mode][1];
														}
														else
														{
															enc_picture->mv[LIST_1][i+img->block_x][j+img->block_y][0] = 0;
															enc_picture->mv[LIST_1][i+img->block_x][j+img->block_y][1] = 0;
														}
													}
												}
										}
										else
										{
											for (j = 0; j < 4; j++)
												for (i = 0; i < 4; i++)
												{
													enc_picture->mv[LIST_0][i+img->block_x][j+img->block_y][0] = 0;
													enc_picture->mv[LIST_0][i+img->block_x][j+img->block_y][1] = 0;
													
													if (bframe)
													{
														enc_picture->mv[LIST_1][i+img->block_x][j+img->block_y][0] = 0;
														enc_picture->mv[LIST_1][i+img->block_x][j+img->block_y][1] = 0;
													}
												}
										}
									}
									
}                             // end of copy_rdopt_data

static void copy_motion_vectors_MB ()
{
	int i,j,k,l;
	
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			for (k = 0; k < img->max_num_references; k++)
			{
				for (l = 0; l < 9; l++)
				{
					img->all_mv[i][j][LIST_0][k][l][0] = rdopt->all_mv[i][j][LIST_0][k][l][0];
					img->all_mv[i][j][LIST_0][k][l][1] = rdopt->all_mv[i][j][LIST_0][k][l][1];
					
					img->all_mv[i][j][LIST_1][k][l][0] = rdopt->all_mv[i][j][LIST_1][k][l][0];
					img->all_mv[i][j][LIST_1][k][l][1] = rdopt->all_mv[i][j][LIST_1][k][l][1];
					
					img->pred_mv[i][j][LIST_0][k][l][0] = rdopt->pred_mv[i][j][LIST_0][k][l][0];
					img->pred_mv[i][j][LIST_0][k][l][1] = rdopt->pred_mv[i][j][LIST_0][k][l][1];
					
					img->pred_mv[i][j][LIST_1][k][l][0] = rdopt->pred_mv[i][j][LIST_1][k][l][0];
					img->pred_mv[i][j][LIST_1][k][l][1] = rdopt->pred_mv[i][j][LIST_1][k][l][1];
					
				}
			}
		}
	}
}


static void ReportNALNonVLCBits(int tmp_time, int me_time)
{
	//! Need to add type (i.e. SPS, PPS, SEI etc).
    printf ("%04d(NVB)%8d \n", frame_no, stat->bit_ctr_parametersets_n);
	
}
static void ReportFirstframe(int tmp_time,int me_time)
{
	//Rate control
	int bits;
	printf ("%04d(IDR)%8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n,0,
		img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
		img->fld_flag ? "FLD" : "FRM", intras);
	
	//Rate control
	if(input->RCEnable)
	{
		if((!input->PicInterlace)&&(!input->MbInterlace))
			bits = stat->bit_ctr-stat->bit_ctr_n; // used for rate control update 
		else
		{
			bits = stat->bit_ctr - Iprev_bits; // used for rate control update 
			Iprev_bits = stat->bit_ctr;
		}
	}
	total_bit+=(stat->bit_ctr - stat->bit_ctr_n);
	tol_time+=tmp_time;
	stat->bitr0 = stat->bitr;
	stat->bit_ctr_0 = stat->bit_ctr;
	stat->bit_ctr = 0;
	psnr_avg1+=snr->snr_y;
	psnr_avg2+=snr->snr_u;
	psnr_avg3+=snr->snr_v;

}


static void ReportIntra(int tmp_time, int me_time)
{
	frame_no=img->current_frame;

	if (img->currentPicture->idr_flag == 1)
		printf ("%04d(IDR)%8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n, 0,
		img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
		img->fld_flag ? "FLD" : "FRM", intras); 
	else
		printf ("%04d(I)  %8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n, 0,
		img->qp, snr->snr_y, snr->snr_u, snr->snr_v, tmp_time, me_time,
		img->fld_flag ? "FLD" : "FRM", intras);
	psnr_avg1+=snr->snr_y;
	psnr_avg2+=snr->snr_u;
	psnr_avg3+=snr->snr_v;
	total_bit+=(stat->bit_ctr - stat->bit_ctr_n);
	tol_time+=tmp_time;

}

static void ReportSP(int tmp_time, int me_time)
{
	printf ("%04d(SP) %8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n, active_pps->weighted_pred_flag, img->qp, snr->snr_y,
		snr->snr_u, snr->snr_v, tmp_time, me_time,
		img->fld_flag ? "FLD" : "FRM", intras);
}

static void ReportBS(int tmp_time, int me_time)
{
	printf ("%04d(BS) %8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d %1d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n, active_pps->weighted_bipred_idc, img->qp, snr->snr_y,
		snr->snr_u, snr->snr_v, tmp_time, me_time,
		img->fld_flag ? "FLD" : "FRM", intras,img->direct_type);
}

static void ReportB(int tmp_time, int me_time)
{
    printf ("%04d(B)  %8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d %1d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n, active_pps->weighted_bipred_idc,img->qp,
		snr->snr_y, snr->snr_u, snr->snr_v, tmp_time,me_time,
		img->fld_flag ? "FLD" : "FRM",intras,img->direct_type);
}


static void ReportP(int tmp_time, int me_time)
{            
    printf ("%04d(P)  %8d %1d %2d %7.3f %7.3f %7.3f  %7d   %5d     %3s   %3d\n",
		frame_no, stat->bit_ctr - stat->bit_ctr_n, active_pps->weighted_pred_flag, img->qp, snr->snr_y,
		snr->snr_u, snr->snr_v, tmp_time, me_time,
		img->fld_flag ? "FLD" : "FRM", intras);
	total_bit+=(stat->bit_ctr - stat->bit_ctr_n);
	tol_time+=tmp_time;
	
}

/*!
************************************************************************
* \brief
*    Copies contents of a Sourceframe structure into the old-style
*    variables imgY_org_frm and imgUV_org_frm.  No other side effects
* \param sf
*    the source frame the frame is to be taken from
************************************************************************
*/

static void CopyFrameToOldImgOrgVariables (Sourceframe *sf)
{
	int x, y;
	for (y=0; y<sf->y_framesize; y++)
		for (x=0; x<sf->x_size; x++)
			imgY_org_frm[y][x] = sf->yf[y*sf->x_size+x];
	for (y=0; y<sf->y_framesize/2; y++)
		for (x=0; x<sf->x_size/2; x++)
		{
			char ab=sf->uf[y*sf->x_size/2+x];
			imgUV_org_frm[0][y][x] = sf->uf[y*sf->x_size/2+x];
			imgUV_org_frm[1][y][x] = sf->vf[y*sf->x_size/2+x];
		}
// 	for (y=0; y<1; y++)
// 	{
// 		for (x=0; x<sf->x_size/4; x++)
// 		{
// 			if (x%2==0 && x%8!=0)
// 			{
// 				printf("%d(%d,%d)   ",imgY_org_frm[y][x],imgUV_org_frm[0][y/2][x/2],imgUV_org_frm[1][y/2][x/2]);
// 			} 
// 			else if (x%2!=0 && x%8!=0)
// 			{
// 				printf("%d  ",imgY_org_frm[y][x]);
// 			}
// 			else
// 			{
// 				printf("%d(%d,%d)\n",imgY_org_frm[y][x],imgUV_org_frm[0][y/2][x/2],imgUV_org_frm[1][y/2][x/2]);
// 			}
// 		}
// 	}
// 	printf("\n");
// 
}


/*!
************************************************************************
* \brief
*    Copies contents of a Sourceframe structure into the old-style
*    variables imgY_org_top and imgUV_org_top.  No other side effects
* \param sf
*    the source frame the field is to be taken from
************************************************************************
*/

static void CopyTopFieldToOldImgOrgVariables (Sourceframe *sf)
{
	int x, y;
	
	for (y=0; y<sf->y_fieldsize; y++)
		for (x=0; x<sf->x_size; x++)
			imgY_org_top [y][x] = sf->yt[y*sf->x_size+x];
		for (y=0; y<sf->y_fieldsize/2; y++)
			for (x=0;x<sf->x_size/2; x++)
			{
				imgUV_org_top[0][y][x] = sf->ut[y*sf->x_size/2+x];
				imgUV_org_top[1][y][x] = sf->vt[y*sf->x_size/2+x];
			}
}
/*!
************************************************************************
* \brief
*    Copies contents of a Sourceframe structure into the old-style
*    variables imgY_org_bot and imgUV_org_bot.  No other side effects
* \param sf
*    the source frame the field is to be taken from
************************************************************************
*/

static void CopyBottomFieldToOldImgOrgVariables (Sourceframe *sf)
{
	int x, y;
	
	for (y=0; y<sf->y_fieldsize; y++)
		for (x=0; x<sf->x_size; x++)
			imgY_org_bot [y][x] = sf->yb[y*sf->x_size+x];
		for (y=0; y<sf->y_fieldsize/2; y++)
			for (x=0;x<sf->x_size/2; x++)
			{
				imgUV_org_bot[0][y][x] = sf->ub[y*sf->x_size/2+x];
				imgUV_org_bot[1][y][x] = sf->vb[y*sf->x_size/2+x];
			}
}


/*!
************************************************************************
* \brief
*    Allocates Sourceframe structure
* \param xs
*    horizontal size of frame in pixels
* \param ys
*    vertical size of frame in pixels, must be divisible by 2
* \return
*    pointer to initialized source frame structure
************************************************************************
*/

static Sourceframe *AllocSourceframe (int xs, int ys)
{
	Sourceframe *sf = NULL;
	const unsigned int bytes_y = xs*ys;
	const unsigned int bytes_uv = (xs*ys)/4;
	
	if ((sf = calloc (1, sizeof (Sourceframe))) == NULL) no_mem_exit ("ReadOneFrame: sf");
	if (sf->yf == NULL) if ((sf->yf = calloc (1, bytes_y)) == NULL) no_mem_exit ("ReadOneFrame: sf->yf");
	if (sf->yt == NULL) if ((sf->yt = calloc (1, bytes_y/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->yt");
	if (sf->yb == NULL) if ((sf->yb = calloc (1, bytes_y/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->yb");
	if (sf->uf == NULL) if ((sf->uf = calloc (1, bytes_uv)) == NULL) no_mem_exit ("ReadOneFrame: sf->uf");
	if (sf->ut == NULL) if ((sf->ut = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->ut");
	if (sf->ub == NULL) if ((sf->ub = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->ub");
	if (sf->vf == NULL) if ((sf->vf = calloc (1, bytes_uv)) == NULL) no_mem_exit ("ReadOneFrame: sf->vf");
	if (sf->vt == NULL) if ((sf->vt = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->vt");
	if (sf->vb == NULL) if ((sf->vb = calloc (1, bytes_uv/2)) == NULL) no_mem_exit ("ReadOneFrame: sf->vb");
	sf->x_size = xs;
	sf->y_framesize = ys;
	sf->y_fieldsize = ys/2;
	
	return sf;
}
 

/*!
************************************************************************
* \brief
*    Frees Sourceframe structure
* \param sf
*    pointer to Sourceframe previoously allocated with ALlocSourceframe()
* \return
*    none
************************************************************************
*/

static void FreeSourceframe (Sourceframe *sf)
{
	if (sf!=NULL) 
	{
		if (sf->yf != NULL) free (sf->yf);
		if (sf->yt != NULL) free (sf->yt);
		if (sf->yb != NULL) free (sf->yb);
		if (sf->uf != NULL) free (sf->uf);
		if (sf->ut != NULL) free (sf->ut);
		if (sf->ub != NULL) free (sf->ub);
		if (sf->vf != NULL) free (sf->vf);
		if (sf->vt != NULL) free (sf->vt);
		if (sf->vb != NULL) free (sf->vb);
		free (sf);
	}
}

/*!
************************************************************************
* \brief
*    Calculates the absolute frame number in the source file out
*    of various variables in img-> and input->
* \return
*    frame number in the file to be read
* \par side effects
*    global variable frame_no updated -- dunno, for what this one is necessary
************************************************************************
*/
static int CalculateFrameNumber()
{
	if (img->type == B_SLICE)
		frame_no = start_tr_in_this_IGOP + (IMG_NUMBER - 1) * (input->jumpd + 1) + img->b_interval * img->b_frame_to_code;
	else
    {
		frame_no = start_tr_in_this_IGOP + IMG_NUMBER * (input->jumpd + 1);
#ifdef _ADAPT_LAST_GROUP_
		if (input->last_frame && img->current_frame + 1 == input->no_frames_h264)
			frame_no = input->last_frame;
#endif
    }
	return frame_no;
}


/*!
************************************************************************
* \brief
*    Generate Field Component from Frame Components by copying
* \param src
*    source frame component
* \param top
*    destination top field component
* \param bot
*    destination bottom field component
* \param xs
*    horizontal size of frame in pixels
* \param ys
*    vertical size of frame in pixels, must be divisible by 2
************************************************************************
*/
static void GenerateFieldComponent (char *src, char *top, char *bot, int xs, int ys)
{
	int fieldline;
	assert (ys % 2 == 0);
	
	for (fieldline = 0; fieldline < ys/2; fieldline++)
	{
		memcpy (&top[xs * fieldline], &src[xs * (fieldline * 2 + 0)], xs);
		memcpy (&bot[xs * fieldline], &src[xs * (fieldline * 2 + 1)], xs);
	}
}


/*!
************************************************************************
* \brief
*    Reads one new frame from file
* \param FrameNoInFile
*    Frame number in the source file
* \param HeaderSize
*    Number of bytes in the source file to be skipped
* \param xs
*    horizontal size of frame in pixels, must be divisible by 16
* \param ys
*    vertical size of frame in pixels, must be divisible by 16 or
*    32 in case of MB-adaptive frame/field coding
* \param sf
*    Sourceframe structure to which the frame is written
************************************************************************
*/
static void i_ReadOneFrame (int FrameNoInFile, int HeaderSize, int xs, int ys, Sourceframe *sf)
{
	  int i;
								  
	  const unsigned int bytes_y = xs*ys;
	  const unsigned int bytes_uv = (xs*ys)/4;
	  const int framesize_in_bytes = bytes_y + 2*bytes_uv;
								  
//	  assert (xs % MB_BLOCK_SIZE == 0);
//	  assert (ys % MB_BLOCK_SIZE == 0);
	  assert (/*fp_in_c*/fp_in != NULL);
	  assert (sf != NULL);
	  assert (sf->yf != NULL);
								  
	  assert (FrameNumberInFile == FrameNoInFile);
	  // printf ("ReadOneFrame: frame_no %d xs %d ys %d\n", FrameNoInFile, xs, ys);
								  
	  if (fseek (/*fp_in_c*/fp_in, HeaderSize, SEEK_SET) != 0)
		  error ("i_ReadOneFrame: cannot fseek to (Header size) in fp_in_c", -1);
								  
	  // the reason for the following loop is to support source files bigger than
	  // MAXINT.  In most operating systems, including Windows, it is possible to
	  // fseek to file positions bigger than MAXINT by using this relative seeking
      // technique.  StW, 12/30/02
	  // Skip starting frames
	  for (i=0; i<input->start_frame; i++)
	    if (fseek (/*fp_in_c*/fp_in, framesize_in_bytes, SEEK_CUR) != 0) 
		{
		  printf ("i_ReadOneFrame: cannot advance file pointer in fp_in_c beyond frame %d, looping to picture zero\n", i);
		  if (fseek (/*fp_in_c*/fp_in, HeaderSize, SEEK_SET) != 0)
    	    report_stats_on_error();
		  exit (-1);
		 } 
							  
		  for (i=0; i<FrameNoInFile; i++)
			 if (fseek (/*fp_in_c*/fp_in, framesize_in_bytes, SEEK_CUR) != 0) 
			{
        	  printf ("i_ReadOneFrame: cannot advance file pointer in fp_in_c beyond frame %d, looping to picture zero\n", i);
     		  if (fseek (/*fp_in_c*/fp_in, HeaderSize, SEEK_SET) != 0)
			  error ("i_ReadOneFrame: cannot fseek to (Header size) in fp_in_c", -1);
			 }
										  
		  // Here we are at the correct position for the source frame in the file.  Now
		  // read it.
		  if (fread (sf->yf, 1, bytes_y, /*fp_in_c*/fp_in) != bytes_y)
		  {
		     printf ("i_ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_y);
		     report_stats_on_error();
			 exit (-1);
		  }
		  if (fread (sf->uf, 1, bytes_uv, /*fp_in_c*/fp_in) != bytes_uv)
		  {
	    	  printf ("i_ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_uv);
			  report_stats_on_error();
			  exit (-1);
		  }
		  if (fread (sf->vf, 1, bytes_uv, /*fp_in_c*/fp_in) != bytes_uv)
		  {
			  printf ("i_ReadOneFrame: cannot read %d bytes from input file, unexpected EOF?, exiting", bytes_uv);
			  report_stats_on_error();
			  exit (-1);
		  }
										  
		  // Complete frame is read into sf->?f, now setup 
		  // top and bottom field (sf->?t and sf->?b)
										  
	    	  GenerateFieldComponent (sf->yf, sf->yt, sf->yb, xs, ys);
			  GenerateFieldComponent (sf->uf, sf->ut, sf->ub, xs/2, ys/2);
			  GenerateFieldComponent (sf->vf, sf->vt, sf->vb, xs/2, ys/2);
		  }
							  
							  
		  /*!
		  ************************************************************************
		  * \brief
		  *    point to frame coding variables 
		  ************************************************************************
*/
static void put_buffer_frame()
{
	imgY_org  = imgY_org_frm;
	imgUV_org = imgUV_org_frm;  
}

/*!
************************************************************************
* \brief
*    point to top field coding variables 
************************************************************************
*/
static void put_buffer_top()
{
	img->fld_type = 0;
	
	imgY_org = imgY_org_top;
	imgUV_org = imgUV_org_top;
}

		  /*!
		  ************************************************************************
		  * \brief
		  *    point to bottom field coding variables 
		  ************************************************************************
*/
static void put_buffer_bot()
{
	img->fld_type = 1;
	
	imgY_org = imgY_org_bot;
	imgUV_org = imgUV_org_bot;
}

/*!
************************************************************************
* \brief
*    Writes a NAL unit of a partition or slice
************************************************************************
*/

static void writeUnit(Bitstream* currStream,int partition)
{
	NALU_t *nalu;
	assert (currStream->bits_to_go == 8);
	nalu = AllocNALU(img->width*img->height*4);
	nalu->startcodeprefix_len = 2+(img->current_mb_nr == 0?ZEROBYTES_SHORTSTARTCODE+1:ZEROBYTES_SHORTSTARTCODE);
// 	printf ("\nwriteUnit:nalu->startcodeprefix_len %d\n", nalu->startcodeprefix_len);
	nalu->len = currStream->byte_pos +1;            // add one for the first byte of the NALU
// 	printf ("\nwriteUnit:currStream->byte_pos %d,nalu->len %d",currStream->byte_pos, nalu->len);
	memcpy (&nalu->buf[1], currStream->streamBuffer, nalu->len-1);
	if (img->currentPicture->idr_flag)//
	{
		nalu->nal_unit_type = NALU_TYPE_IDR;
		nalu->nal_reference_idc = NALU_PRIORITY_HIGHEST;
	}
	else if (img->type == B_SLICE)
	{
		//different nal header for different partitions
		if(input->partition_mode == 0)
		{
			nalu->nal_unit_type = NALU_TYPE_SLICE;
		}
		else
		{
			nalu->nal_unit_type = NALU_TYPE_DPA +  partition;
		}
		
		if (img->nal_reference_idc !=0)
		{
			nalu->nal_reference_idc = NALU_PRIORITY_HIGH;
		}
		else
		{
			nalu->nal_reference_idc = NALU_PRIORITY_DISPOSABLE;
		}
	}
	else   // non-b frame, non IDR slice
	{
		//different nal header for different partitions
		if(input->partition_mode == 0)//PartitionMode=0
		{
			nalu->nal_unit_type = NALU_TYPE_SLICE;
		}
		else
		{
			nalu->nal_unit_type = NALU_TYPE_DPA +  partition;
		}
		if (img->nal_reference_idc !=0)
		{
			nalu->nal_reference_idc = NALU_PRIORITY_HIGH;
		}
		else
		{
			nalu->nal_reference_idc = NALU_PRIORITY_DISPOSABLE;
		}
	}
	nalu->forbidden_bit = 0;
	stat->bit_ctr += WriteNALU (nalu);
	
	FreeNALU(nalu);
}

/////////////////////////////////////>

/*!
***********************************************************************
* \brief
*    decodes one I- or P-frame
*
***********************************************************************
*/

int decode_one_frame(struct img_par *img,struct inp_par_dec *inp, struct snr_par *snr)
{
	int current_header;
	Slice *currSlice = img->currentSlice;
	
	img->current_slice_nr = 0;
	img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
	currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
	img->num_dec_mb = 0;
	img->newframe = 1;
	
	while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
	{
		current_header = read_new_slice();//改变了qp的值
		
		if (current_header == EOS)//2=1
		{
			exit_picture();
			return EOS;
		}
		
		decode_slice(img, inp, current_header);//查img来源路线
		
		img->newframe = 0;
		img->current_slice_nr++;
	}
	
// 	exit_picture();
	
	return (SOP);
}


/*!
************************************************************************
* \brief
*    initialize ref_pic_num array
************************************************************************
*/
// void set_ref_pic_num_dec()
// {
// 	int i,j;
// 	
// 	int slice_id=img->current_slice_nr;
// 	
// 	for (i=0;i<listXsize[LIST_0];i++)
// 	{
// 		dec_picture->ref_pic_num        [slice_id][LIST_0][i]=listX[LIST_0][i]->poc * 2 + ((listX[LIST_0][i]->structure==BOTTOM_FIELD)?1:0) ; 
// 		dec_picture->frm_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->frame_poc * 2; 
// 		dec_picture->top_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->top_poc * 2; 
// 		dec_picture->bottom_ref_pic_num [slice_id][LIST_0][i]=listX[LIST_0][i]->bottom_poc * 2 + 1; 
// 		//printf("POCS %d %d %d %d ",listX[LIST_0][i]->frame_poc,listX[LIST_0][i]->bottom_poc,listX[LIST_0][i]->top_poc,listX[LIST_0][i]->poc);
// 		//printf("refid %d %d %d %d\n",(int) dec_picture->frm_ref_pic_num[LIST_0][i],(int) dec_picture->top_ref_pic_num[LIST_0][i],(int) dec_picture->bottom_ref_pic_num[LIST_0][i],(int) dec_picture->ref_pic_num[LIST_0][i]);
// 	}
// 	
// 	for (i=0;i<listXsize[LIST_1];i++)
// 	{
// 		dec_picture->ref_pic_num        [slice_id][LIST_1][i]=listX[LIST_1][i]->poc  *2 + ((listX[LIST_1][i]->structure==BOTTOM_FIELD)?1:0);
// 		dec_picture->frm_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->frame_poc * 2; 
// 		dec_picture->top_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->top_poc * 2; 
// 		dec_picture->bottom_ref_pic_num [slice_id][LIST_1][i]=listX[LIST_1][i]->bottom_poc * 2 + 1; 
// 	}
// 	
// 	if (img->structure==FRAME)
// 		for (j=2;j<6;j++)
// 			for (i=0;i<listXsize[j];i++)
// 			{    
// 				dec_picture->ref_pic_num        [slice_id][j][i] = listX[j][i]->poc * 2 + ((listX[j][i]->structure==BOTTOM_FIELD)?1:0);
// 				dec_picture->frm_ref_pic_num    [slice_id][j][i] = listX[j][i]->frame_poc * 2 ;
// 				dec_picture->top_ref_pic_num    [slice_id][j][i] = listX[j][i]->top_poc * 2 ;
// 				dec_picture->bottom_ref_pic_num [slice_id][j][i] = listX[j][i]->bottom_poc * 2 + 1;
// 			}
// 			
// }
// 

/*!
************************************************************************
* \brief
*    Find PSNR for all three components.Compare decoded frame with
*    the original sequence. Read inp->jumpd frames to reflect frame skipping.
************************************************************************
*/
void find_snr_dec(
			  struct snr_par  *snr,   //!< pointer to snr parameters
			  StorablePicture *p,     //!< picture to be compared
			  FILE *p_ref)            //!< open reference YUV file
{
	int i,j;
	int diff_y,diff_u,diff_v;
	int uv;
	int  status;
	
	// calculate frame number
	int  psnrPOC = active_sps->mb_adaptive_frame_field_flag ? p->poc /(input_dec->poc_scale) : p->poc/(3-input_dec->poc_scale);
	//  int  psnrPOC = p->MbaffFrameFlag ? p->poc /(input->poc_scale) : p->poc/(3-input->poc_scale);
	
	// KS: Code below might work better if you have fields and a large (>1) poc offset between them
	//  int  poc_diff=max(1,(p->bottom_poc - p->top_poc));
	//  int  psnrPOC = active_sps->mb_adaptive_frame_field_flag ? p->poc /(input->poc_scale*poc_diff) : p->poc/((3-input->poc_scale)*poc_diff);
	
	if (psnrPOC==0 && img->psnr_number)
		img->idr_psnr_number=img->psnr_number + 1;
	img->psnr_number=max(img->psnr_number,img->idr_psnr_number+psnrPOC);
	
	frame_no = img->idr_psnr_number+psnrPOC;
	
	
	rewind(p_ref);
	
	for (i=0; i<frame_no; i++)
	{
		status = fseek (p_ref, (long) p->size_y* (long) (p->size_x*3/2), SEEK_CUR);
		if (status != 0)
		{
			printf(errortext, ET_SIZE, "Error in seeking frame number: %d", frame_no);
			fprintf(stderr, "%s", errortext);
			return;
			//      snprintf(errortext, ET_SIZE, "Error in seeking frame number: %d", frame_no);
			//      error(errortext, 500);
		}
	}
	
	for (j=0; j < p->size_y; j++)
		for (i=0; i < p->size_x; i++)
			imgY_ref[j][i]=fgetc(p_ref);
		
		for (uv=0; uv < 2; uv++)
			for (j=0; j < p->size_y_cr ; j++)
				for (i=0; i < p->size_x_cr; i++)
					imgUV_ref[uv][j][i]=fgetc(p_ref);
				
				img->quad_dec[0]=0;
				diff_y=0;
				for (j=0; j < p->size_y; ++j)
				{
					for (i=0; i < p->size_x; ++i)
					{
						diff_y += img->quad_dec[abs(p->imgY[j][i]-imgY_ref[j][i])];
					}
				}
				
				// Chroma
				diff_u=0;
				diff_v=0;
				
				for (j=0; j < p->size_y_cr; ++j)
				{
					for (i=0; i < p->size_x_cr; ++i)
					{
						diff_u += img->quad_dec[abs(imgUV_ref[0][j][i]-p->imgUV[0][j][i])];
						diff_v += img->quad_dec[abs(imgUV_ref[1][j][i]-p->imgUV[1][j][i])];
					}
				}
				
				/*  if (diff_y == 0)
				diff_y = 1;
				if (diff_u == 0)
				diff_u = 1;
				if (diff_v == 0)
				diff_v = 1; */
				
				// Collecting SNR statistics
				if (diff_y != 0)
					snr->snr_y=(float)(10*log10(65025*(float)(p->size_x)*(p->size_y)/(float)diff_y));        // luma snr for current frame
				else
					snr->snr_y=0;
				if (diff_u != 0)
					snr->snr_u=(float)(10*log10(65025*(float)(p->size_x)*(p->size_y)/(float)(4*diff_u)));    //  chroma snr for current frame
				else
					snr->snr_u=0;
				if (diff_v != 0)
					snr->snr_v=(float)(10*log10(65025*(float)(p->size_x)*(p->size_y)/(float)(4*diff_v)));    //  chroma snr for current frame
				else
					snr->snr_v=0;
				
				if (img->number == 0) // first
				{
					snr->snr_ya=snr->snr_y1=snr->snr_y;                                                        // keep luma snr for first frame
					snr->snr_ua=snr->snr_u1=snr->snr_u;                                                        // keep chroma snr for first frame
					snr->snr_va=snr->snr_v1=snr->snr_v;                                                        // keep chroma snr for first frame
					
				}
				else
				{
					snr->snr_ya=(float)(snr->snr_ya*(img->number+Bframe_ctr)+snr->snr_y)/(img->number+Bframe_ctr+1); // average snr chroma for all frames
					snr->snr_ua=(float)(snr->snr_ua*(img->number+Bframe_ctr)+snr->snr_u)/(img->number+Bframe_ctr+1); // average snr luma for all frames
					snr->snr_va=(float)(snr->snr_va*(img->number+Bframe_ctr)+snr->snr_v)/(img->number+Bframe_ctr+1); // average snr luma for all frames
				} 
}


/*!
************************************************************************
* \brief
*    Interpolation of 1/4 subpixel
************************************************************************
*/
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
{
	
	int dx, dy;
	int x, y;
	int i, j;
	int maxold_x,maxold_y;
	int result;
	int pres_x;
	int pres_y; 
	int tmp_res[4][9];
	static const int COEF[6] = {    1, -5, 20, 20, -5, 1  };
	
	dx = x_pos&3;
	dy = y_pos&3;
	x_pos = (x_pos-dx)/4;
	y_pos = (y_pos-dy)/4;
	
	maxold_x = dec_picture->size_x-1;
	maxold_y = dec_picture->size_y-1;
	
	if (dec_picture->mb_field[img->current_mb_nr])
		maxold_y = dec_picture->size_y/2 - 1;
	
	if (dx == 0 && dy == 0) {  /* fullpel position */
		for (j = 0; j < BLOCK_SIZE; j++)
			for (i = 0; i < BLOCK_SIZE; i++)
				block[i][j] = list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
	}
	else { /* other positions */
		
		if (dy == 0) { /* No vertical interpolation */
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (result = 0, x = -2; x < 4; x++)
						result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
					block[i][j] = max(0, min(255, (result+16)/32));
				}
			}
			
			if ((dx&1) == 1) {
				for (j = 0; j < BLOCK_SIZE; j++)
					for (i = 0; i < BLOCK_SIZE; i++)
						block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+dx/2))] +1 )/2;
			}
		}
		else if (dx == 0) {  /* No horizontal interpolation */
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (result = 0, y = -2; y < 4; y++)
						result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
					block[i][j] = max(0, min(255, (result+16)/32));
				}
			}
			
			if ((dy&1) == 1) {
				for (j = 0; j < BLOCK_SIZE; j++)
					for (i = 0; i < BLOCK_SIZE; i++)
						block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+dy/2))][max(0,min(maxold_x,x_pos+i))] +1 )/2;
			}
		}
		else if (dx == 2) {  /* Vertical & horizontal interpolation */
			
			for (j = -2; j < BLOCK_SIZE+3; j++) {
				for (i = 0; i < BLOCK_SIZE; i++)
					for (tmp_res[i][j+2] = 0, x = -2; x < 4; x++)
						tmp_res[i][j+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
			}
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (result = 0, y = -2; y < 4; y++)
						result += tmp_res[i][j+y+2]*COEF[y+2];
					block[i][j] = max(0, min(255, (result+512)/1024));
				} 
			}
			
			if ((dy&1) == 1) {
				for (j = 0; j < BLOCK_SIZE; j++)
					for (i = 0; i < BLOCK_SIZE; i++)
						block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[i][j+2+dy/2]+16)/32)) +1 )/2;
			}
		}
		else if (dy == 2) {  /* Horizontal & vertical interpolation */
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = -2; i < BLOCK_SIZE+3; i++)
					for (tmp_res[j][i+2] = 0, y = -2; y < 4; y++)
						tmp_res[j][i+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
			}
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (result = 0, x = -2; x < 4; x++)
						result += tmp_res[j][i+x+2]*COEF[x+2];
					block[i][j] = max(0, min(255, (result+512)/1024));
				}
			}
			
			if ((dx&1) == 1) {
				for (j = 0; j < BLOCK_SIZE; j++)
					for (i = 0; i < BLOCK_SIZE; i++)
						block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[j][i+2+dx/2]+16)/32))+1)/2;
			}
		}
		else {  /* Diagonal interpolation */
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = 0; i < BLOCK_SIZE; i++) {
					pres_y = dy == 1 ? y_pos+j : y_pos+j+1;
					pres_y = max(0,min(maxold_y,pres_y));
					for (result = 0, x = -2; x < 4; x++)
						result += list[ref_frame]->imgY[pres_y][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
					block[i][j] = max(0, min(255, (result+16)/32));
				}
			}
			
			for (j = 0; j < BLOCK_SIZE; j++) {
				for (i = 0; i < BLOCK_SIZE; i++) {
					pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
					pres_x = max(0,min(maxold_x,pres_x));
					for (result = 0, y = -2; y < 4; y++)
						result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][pres_x]*COEF[y+2];
					block[i][j] = (block[i][j] + max(0, min(255, (result+16)/32)) +1 ) / 2;
				}
			}
			
		}
  }
  
}


// void reorder_lists(int currSliceType, Slice * currSlice)
// {
// 	
// 	if ((currSliceType != I_SLICE)&&(currSliceType != SI_SLICE))
// 	{
// 		if (currSlice->ref_pic_list_reordering_flag_l0)
// 		{
// 			reorder_ref_pic_list(listX[0], &listXsize[0], 
// 				img->num_ref_idx_l0_active - 1, 
// 				currSlice->remapping_of_pic_nums_idc_l0, 
// 				currSlice->abs_diff_pic_num_minus1_l0, 
// 				currSlice->long_term_pic_idx_l0);
// 		}
// 		if (NULL == listX[0][img->num_ref_idx_l0_active-1])
// 		{
// 			error("number of entries in list 0 smaller than num_ref_idx_l0_active_minus1",500);
// 		}
// 		// that's a definition
// 		listXsize[0] = img->num_ref_idx_l0_active;
// 	}
// 	if (currSliceType == B_SLICE)
// 	{
// 		if (currSlice->ref_pic_list_reordering_flag_l1)
// 		{
// 			reorder_ref_pic_list(listX[1], &listXsize[1], 
// 				img->num_ref_idx_l1_active - 1, 
// 				currSlice->remapping_of_pic_nums_idc_l1, 
// 				currSlice->abs_diff_pic_num_minus1_l1, 
// 				currSlice->long_term_pic_idx_l1);
// 		}
// 		if (NULL == listX[1][img->num_ref_idx_l1_active-1])
// 		{
// 			error("number of entries in list 1 smaller than num_ref_idx_l1_active_minus1",500);
// 		}
// 		// that's a definition
// 		listXsize[1] = img->num_ref_idx_l1_active;
// 	}
// 	
// 	free_ref_pic_list_reordering_buffer(currSlice);
// }
// 
// 
// /*!
// ************************************************************************
// * \brief
// *    initialize ref_pic_num array
// ************************************************************************
// */
// void set_ref_pic_num()
// {
// 	int i,j;
// 	
// 	int slice_id=img->current_slice_nr;
// 	
// 	for (i=0;i<listXsize[LIST_0];i++)
// 	{
// 		dec_picture->ref_pic_num        [slice_id][LIST_0][i]=listX[LIST_0][i]->poc * 2 + ((listX[LIST_0][i]->structure==BOTTOM_FIELD)?1:0) ; 
// 		dec_picture->frm_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->frame_poc * 2; 
// 		dec_picture->top_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->top_poc * 2; 
// 		dec_picture->bottom_ref_pic_num [slice_id][LIST_0][i]=listX[LIST_0][i]->bottom_poc * 2 + 1; 
// 		//printf("POCS %d %d %d %d ",listX[LIST_0][i]->frame_poc,listX[LIST_0][i]->bottom_poc,listX[LIST_0][i]->top_poc,listX[LIST_0][i]->poc);
// 		//printf("refid %d %d %d %d\n",(int) dec_picture->frm_ref_pic_num[LIST_0][i],(int) dec_picture->top_ref_pic_num[LIST_0][i],(int) dec_picture->bottom_ref_pic_num[LIST_0][i],(int) dec_picture->ref_pic_num[LIST_0][i]);
// 	}
// 	
// 	for (i=0;i<listXsize[LIST_1];i++)
// 	{
// 		dec_picture->ref_pic_num        [slice_id][LIST_1][i]=listX[LIST_1][i]->poc  *2 + ((listX[LIST_1][i]->structure==BOTTOM_FIELD)?1:0);
// 		dec_picture->frm_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->frame_poc * 2; 
// 		dec_picture->top_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->top_poc * 2; 
// 		dec_picture->bottom_ref_pic_num [slice_id][LIST_1][i]=listX[LIST_1][i]->bottom_poc * 2 + 1; 
// 	}
// 	
// 	if (img->structure==FRAME)
// 		for (j=2;j<6;j++)
// 			for (i=0;i<listXsize[j];i++)
// 			{    
// 				dec_picture->ref_pic_num        [slice_id][j][i] = listX[j][i]->poc * 2 + ((listX[j][i]->structure==BOTTOM_FIELD)?1:0);
// 				dec_picture->frm_ref_pic_num    [slice_id][j][i] = listX[j][i]->frame_poc * 2 ;
// 				dec_picture->top_ref_pic_num    [slice_id][j][i] = listX[j][i]->top_poc * 2 ;
// 				dec_picture->bottom_ref_pic_num [slice_id][j][i] = listX[j][i]->bottom_poc * 2 + 1;
// 			}
// 			
// }
// 
// 
/*!
************************************************************************
* \brief
*    Reads new slice from bit_stream
************************************************************************
*/
int read_new_slice()
{
	NALU_t *nalu = AllocNALU(MAX_CODED_FRAME_SIZE);
	int current_header;
	int ret;
	int BitsUsedByHeader;
	Slice *currSlice = img->currentSlice;
	Bitstream *currStream;
	
	int slice_id_a, slice_id_b, slice_id_c;
	int redundant_pic_cnt_b, redundant_pic_cnt_c;
	long ftell_position, expected_slice_type;
	
	//  int i;
	expected_slice_type = NALU_TYPE_DPA;//2
	
	while (1)
	{
		ftell_position = ftell(bits);//0//13//22
		
		if (input_dec->FileFormat == PAR_OF_ANNEXB)//0=0
			ret=GetAnnexbNALU (nalu);//13//9//7329//到这里都一样了。
		else
			ret=GetRTPNALU (nalu);
		
		//In some cases, zero_byte shall be present. If current NALU is a VCL NALU, we can't tell
		//whether it is the first VCL NALU at this point, so only non-VCL NAL unit is checked here.
		CheckZeroByteNonVCL(nalu, &ret);
		
		NALUtoRBSP(nalu);
		//    printf ("nalu->len %d\n", nalu->len);
		
		if (ret < 0)//13//9
			printf ("Error while getting the NALU in file format %s, exit\n", input_dec->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
		if (ret == 0)//13//9
		{
			//      printf ("read_new_slice: returning %s\n", "EOS");
			if(expected_slice_type != NALU_TYPE_DPA)
			{
				/* oops... we found the next slice, go back! */
				fseek(bits, ftell_position, SEEK_SET);
				FreeNALU(nalu);
				return current_header;
			}
			else
				return EOS;
		}
		
		// Got a NALU
		if (nalu->forbidden_bit)
		{
			printf ("Found NALU w/ forbidden_bit set, bit error?  Let's try...\n");
		}
		
		switch (nalu->nal_unit_type)
		{
		case NALU_TYPE_SLICE:
// 			printf("NALU_TYPE_SLICE\n");
		case NALU_TYPE_IDR:
// 			printf("NALU_TYPE_IDR\n");
			img->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
			img->nal_reference_idc = nalu->nal_reference_idc;
			img->disposable_flag = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
			currSlice->dp_mode = PAR_DP_1;
			currSlice->max_part_nr = 1;
			currSlice->ei_flag = 0;
			currStream = currSlice->partArr[0].bitstream;
			currStream->ei_flag = 0;
			currStream->frame_bitoffset = currStream->read_len = 0;
			memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
			currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
			
			// Some syntax of the Slice Header depends on the parameter set, which depends on
			// the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
			// of the slice header first, then setup the active parameter sets, and then read
			// the rest of the slice header
			BitsUsedByHeader = FirstPartOfSliceHeader();//9
			UseParameterSet (currSlice->pic_parameter_set_id);
			BitsUsedByHeader+= RestOfSliceHeader ();//改变qp的值   37=9+28返回值不一样
			
			FmoInit_dec (active_pps, active_sps);
			
			if(is_new_picture())
			{
				init_picture(img, input);
				
				current_header = SOP;
				//check zero_byte if it is also the first NAL unit in the access unit
				CheckZeroByteVCL(nalu, &ret);
			}
			else
				current_header = SOS;
			
			init_lists(img->type, img->currentSlice->structure);
// 			reorder_lists (img->type, img->currentSlice);
			
			if (img->structure==FRAME)
			{
				init_mbaff_lists();
			}
			
			/*        if (img->frame_num==1) // write a reference list
			{
			count ++;
			if (count==1)
            for (i=0; i<listXsize[0]; i++)
			write_picture(listX[0][i], p_out2);
			}
			*/
			
			// From here on, active_sps, active_pps and the slice header are valid
			if (img->MbaffFrameFlag)
				img->current_mb_nr = currSlice->start_mb_nr << 1;
			else
				img->current_mb_nr = currSlice->start_mb_nr;
			
			if (active_pps->entropy_coding_mode_flag)
			{
				int ByteStartPosition = currStream->frame_bitoffset/8;
				if (currStream->frame_bitoffset%8 != 0) 
				{
					ByteStartPosition++;
				}
				arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
			}
			// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
			FreeNALU(nalu);
			return current_header;
			break;
		case NALU_TYPE_DPA:
// 			printf("NALU_TYPE_DPA\n");
			//! The state machine here should follow the same ideas as the old readSliceRTP()
			//! basically:
			//! work on DPA (as above)
			//! read and process all following SEI/SPS/PPS/PD/Filler NALUs
			//! if next video NALU is dpB, 
			//!   then read and check whether it belongs to DPA, if yes, use it
			//! else
			//!   ;   // nothing
			//! read and process all following SEI/SPS/PPS/PD/Filler NALUs
			//! if next video NALU is dpC
			//!   then read and check whether it belongs to DPA (and DPB, if present), if yes, use it, done
			//! else
			//!   use the DPA (and the DPB if present)
			
			/* 
            LC: inserting the code related to DP processing, mainly copying some of the parts
            related to NALU_TYPE_SLICE, NALU_TYPE_IDR.
			*/
			
			if(expected_slice_type != NALU_TYPE_DPA)
			{
				/* oops... we found the next slice, go back! */
				fseek(bits, ftell_position, SEEK_SET);
				FreeNALU(nalu);
				return current_header;
			}
			
			img->idr_flag          = (nalu->nal_unit_type == NALU_TYPE_IDR);
			img->nal_reference_idc = nalu->nal_reference_idc;
			img->disposable_flag   = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
			currSlice->dp_mode     = PAR_DP_3;
			currSlice->max_part_nr = 3;
			currSlice->ei_flag     = 0;
			currStream             = currSlice->partArr[0].bitstream;
			currStream->ei_flag    = 0;
			currStream->frame_bitoffset = currStream->read_len = 0;
			memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
			currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
			
			BitsUsedByHeader     = FirstPartOfSliceHeader();
			UseParameterSet (currSlice->pic_parameter_set_id);
			BitsUsedByHeader    += RestOfSliceHeader ();
			
			FmoInit_dec (active_pps, active_sps);
			
			if(is_new_picture())
			{
				init_picture(img, input);
				current_header = SOP;
				CheckZeroByteVCL(nalu, &ret);
			}
			else
				current_header = SOS;
			
			
			init_lists(img->type, img->currentSlice->structure);
// 			reorder_lists (img->type, img->currentSlice);
			
			if (img->structure==FRAME)
			{
				init_mbaff_lists();
			}
			
			// From here on, active_sps, active_pps and the slice header are valid
			if (img->MbaffFrameFlag)
				img->current_mb_nr = currSlice->start_mb_nr << 1;
			else
				img->current_mb_nr = currSlice->start_mb_nr;
			
			
				/* 
				LC:
				Now I need to read the slice ID, which depends on the value of 
				redundant_pic_cnt_present_flag (pag.49). 
			*/
			
			slice_id_a  = ue_v_dec("NALU:SLICE_A slice_idr", currStream);
			if (active_pps->entropy_coding_mode_flag)
			{
				int ByteStartPosition = currStream->frame_bitoffset/8;
				if (currStream->frame_bitoffset%8 != 0) 
				{
					ByteStartPosition++;
				}
				arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
			}
			// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
			break;
		case NALU_TYPE_DPB:
// 			printf("NALU_TYPE_DPB\n");
			/* LC: inserting the code related to DP processing */
			
			currStream             = currSlice->partArr[1].bitstream;
			currStream->ei_flag    = 0;
			currStream->frame_bitoffset = currStream->read_len = 0;
			memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
			currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
			
			slice_id_b  = ue_v_dec("NALU:SLICE_B slice_idr", currStream);
			if (active_pps->redundant_pic_cnt_present_flag)
				redundant_pic_cnt_b = ue_v_dec("NALU:SLICE_B redudand_pic_cnt", currStream);
			else
				redundant_pic_cnt_b = 0;
			
			/*  LC: Initializing CABAC for the current data stream. */
			
			if (active_pps->entropy_coding_mode_flag)
			{
				int ByteStartPosition = currStream->frame_bitoffset/8;
				if (currStream->frame_bitoffset % 8 != 0) 
					ByteStartPosition++;
				
				arideco_start_decoding (&currSlice->partArr[1].de_cabac, currStream->streamBuffer, 
					ByteStartPosition, &currStream->read_len, img->type);
				
			}
			
			/* LC: resilience code to be inserted */
			/*         FreeNALU(nalu); */
			/*         return current_header; */
			
			break;
		case NALU_TYPE_DPC:
// 			printf("NALU_TYPE_DPC\n");
			/* LC: inserting the code related to DP processing */
			currStream             = currSlice->partArr[2].bitstream;
			currStream->ei_flag    = 0;
			currStream->frame_bitoffset = currStream->read_len = 0;
			memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
			currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
			
			slice_id_c  = ue_v_dec("NALU:SLICE_C slice_idr", currStream);
			if (active_pps->redundant_pic_cnt_present_flag)
				redundant_pic_cnt_c = ue_v_dec("NALU:SLICE_C redudand_pic_cnt", currStream);
			else
				redundant_pic_cnt_c = 0;
			
			/* LC: Initializing CABAC for the current data stream. */
			
			if (active_pps->entropy_coding_mode_flag)
			{
				int ByteStartPosition = currStream->frame_bitoffset/8;
				if (currStream->frame_bitoffset % 8 != 0) 
					ByteStartPosition++;
				
				arideco_start_decoding (&currSlice->partArr[2].de_cabac, currStream->streamBuffer, 
					ByteStartPosition, &currStream->read_len, img->type);
			}
			
			/* LC: resilience code to be inserted */
			
			FreeNALU(nalu);
			return current_header;
			
			break;
		case NALU_TYPE_SEI:
// 			printf("NALU_TYPE_SEI\n");
// 			printf ("read_new_slice: Found NALU_TYPE_SEI, len %d\n", nalu->len);
//  			InterpretSEIMessage(nalu->buf,nalu->len,img);//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZz
			break;
		case NALU_TYPE_PPS:
// 			printf("NALU_TYPE_PPS\n");
			ProcessPPS(nalu);
			break;
			
		case NALU_TYPE_SPS:
// 			printf("NALU_TYPE_SPS\n");
			ProcessSPS(nalu);
			break;
		case NALU_TYPE_AUD:
// 			printf("NALU_TYPE_AUD\n");
			//        printf ("read_new_slice: Found 'Access Unit Delimiter' NAL unit, len %d, ignored\n", nalu->len);
			break;
		case NALU_TYPE_EOSEQ:
// 			printf("NALU_TYPE_EOSEQ\n");
			//        printf ("read_new_slice: Found 'End of Sequence' NAL unit, len %d, ignored\n", nalu->len);
			break;
		case NALU_TYPE_EOSTREAM:
// 			printf("NALU_TYPE_EOSTREAM\n");
			//        printf ("read_new_slice: Found 'End of Stream' NAL unit, len %d, ignored\n", nalu->len);
			break;
		case NALU_TYPE_FILL:
// 			printf("NALU_TYPE_FILL\n");
			printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", nalu->len);
			printf ("Skipping these filling bits, proceeding w/ next NALU\n");
			break;
		default:
			printf("default");
			printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", nalu->nal_unit_type, nalu->len);
    }
  }
  FreeNALU(nalu);
  
  return  current_header;
}


/*!
************************************************************************
* \brief
*    Initializes the parameters for a new picture
************************************************************************
*/
void init_picture(struct img_par *img, struct inp_par_dec *inp)
{
	int i,k,l;
	if (dec_picture)
	{
		// this may only happen on slice loss
		exit_picture();
	}
	
	if (img->frame_num != img->pre_frame_num && img->frame_num != (img->pre_frame_num + 1) % img->MaxFrameNum) 
	{
		if (active_sps->gaps_in_frame_num_value_allowed_flag == 0)
		{
			/* Advanced Error Concealment would be called here to combat unintentional loss of pictures. */
			error("An unintentional loss of pictures occurs! Exit\n", 100);
		}
		fill_frame_num_gap_dec(img);//ZZZZZZZZZZZZZZZZZZ看看经过不，不经过就不添加了
	}
	img->pre_frame_num = img->frame_num;
	img->num_dec_mb = 0;
	
	//calculate POCZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	decode_poc(img);
	//  dumppoc (img);
	
	if (img->structure==FRAME ||img->structure==TOP_FIELD)
	{
#ifdef WIN32
		_ftime (&(img->tstruct_start));             // start time ms
#else
		ftime (&(img->tstruct_start));              // start time ms
#endif
		time( &(img->ltime_start));                // start time s
	}
	
	dec_picture = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);
	dec_picture->top_poc=img->toppoc;
	dec_picture->bottom_poc=img->bottompoc;
	dec_picture->frame_poc=img->framepoc;
	
	// reset all variables of the error concealment instance before decoding of every frame.
	// here the third parameter should, if perfectly, be equal to the number of slices per frame.
	// using little value is ok, the code will allocate more memory if the slice number is larger
	ercReset(erc_errorVar, img->PicSizeInMbs, img->PicSizeInMbs, dec_picture->size_x);
	erc_mvperMB = 0;
	
	switch (img->structure )
	{
	case TOP_FIELD:
		{
			dec_picture->poc=img->toppoc;
			img->number *= 2;
			break;
		}
	case BOTTOM_FIELD:
		{
			dec_picture->poc=img->bottompoc;
			img->number++;
			break;
		}
	case FRAME:
		{
			dec_picture->poc=img->framepoc;
			break;
		}
	default:
		error("img->structure not initialized", 235);
	}
	
	img->current_slice_nr=0;
	
//  	if (img->type > SI_SLICE)////ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 	{
// 		set_ec_flag(SE_PTYPE);
// 		img->type = P_SLICE;  // concealed element
// 	}
	
	// CAVLC init
	for (i=0;i < (int)img->PicSizeInMbs; i++)
		for (k=0;k<4;k++)
			for (l=0;l<6;l++)
				img->nz_coeff[i][k][l]=-1;  // CAVLC
			
			if(active_pps->constrained_intra_pred_flag)
			{
				for (i=0; i<(int)img->PicSizeInMbs; i++)
				{
					img->intra_block[i] = 1;
				}
			}
			
			// Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
			// TO set Macroblock Map (mark all MBs as 'have to be concealed')
			for(i=0; i<(int)img->PicSizeInMbs; i++)
			{
				img->mb_data[i].slice_nr = -1; 
				img->mb_data[i].ei_flag = 1;
			}
			
			img->mb_y = img->mb_x = 0;
			img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
			img->block_x = img->pix_x = img->pix_c_x = 0; // define horizontal positions
			
			dec_picture->slice_type = img->type;
			dec_picture->used_for_reference = (img->nal_reference_idc != 0);
			dec_picture->idr_flag = img->idr_flag;
			dec_picture->no_output_of_prior_pics_flag = img->no_output_of_prior_pics_flag;
			dec_picture->long_term_reference_flag = img->long_term_reference_flag;
			dec_picture->adaptive_ref_pic_buffering_flag = img->adaptive_ref_pic_buffering_flag;
			
			dec_picture->dec_ref_pic_marking_buffer = img->dec_ref_pic_marking_buffer;
			img->dec_ref_pic_marking_buffer = NULL;
			
			dec_picture->MbaffFrameFlag = img->MbaffFrameFlag;
			dec_picture->PicWidthInMbs = img->PicWidthInMbs;
			dec_picture->pic_num = img->frame_num;
			dec_picture->frame_num = img->frame_num;
			dec_picture->coded_frame = (img->structure==FRAME);
			
			dec_picture->frame_mbs_only_flag = active_sps->frame_mbs_only_flag;
			dec_picture->frame_cropping_flag = active_sps->frame_cropping_flag;
			
			if (dec_picture->frame_cropping_flag)
			{
				dec_picture->frame_cropping_rect_left_offset   = active_sps->frame_cropping_rect_left_offset;
				dec_picture->frame_cropping_rect_right_offset  = active_sps->frame_cropping_rect_right_offset;
				dec_picture->frame_cropping_rect_top_offset    = active_sps->frame_cropping_rect_top_offset;
				dec_picture->frame_cropping_rect_bottom_offset = active_sps->frame_cropping_rect_bottom_offset;
			}
			
}

/*!
************************************************************************
* \brief
*    finish decoding of a picture, conceal errors and store it 
*    into the DPB
************************************************************************
*/
void exit_picture()//dec
{
	int ercStartMB;
	int ercSegment;
	frame recfr;
// 	unsigned int i;
	int structure, frame_poc, slice_type, refpic;
	
	int tmp_time;                   // time used by decoding the last frame
	
	// return if the last picture has already been finished
	if (dec_picture==NULL)
	{
		return;
	}
	
	//deblocking for frame or field
// 	DeblockPicture( img, dec_picture );
	
	if (dec_picture->MbaffFrameFlag)
		MbAffPostProc();
	
	recfr.yptr = &dec_picture->imgY[0][0];
	recfr.uptr = &dec_picture->imgUV[0][0][0];
	recfr.vptr = &dec_picture->imgUV[1][0][0];
	
	//! this is always true at the beginning of a picture
	ercStartMB = 0;
	ercSegment = 0;
	
	//! mark the start of the first segment
// 	if (!dec_picture->MbaffFrameFlag)
// 	{
// 		ercStartSegment(0, ercSegment, 0 , erc_errorVar);
// 		//! generate the segments according to the macroblock map
// 		for(i = 1; i<dec_picture->PicSizeInMbs; i++)
// 		{
// 			if(img->mb_data[i].ei_flag != img->mb_data[i-1].ei_flag)
// 			{
// 				ercStopSegment(i-1, ercSegment, 0, erc_errorVar); //! stop current segment
// 				
// 				//! mark current segment as lost or OK
// 				if(img->mb_data[i-1].ei_flag)
// 					ercMarkCurrSegmentLost(dec_picture->size_x, erc_errorVar);
// 				else
// 					ercMarkCurrSegmentOK(dec_picture->size_x, erc_errorVar);
// 				
// 				ercSegment++;  //! next segment
// 				ercStartSegment(i, ercSegment, 0 , erc_errorVar); //! start new segment
// 				ercStartMB = i;//! save start MB for this segment 
// 			}
// 		}
// 		//! mark end of the last segment
// 		ercStopSegment(dec_picture->PicSizeInMbs-1, ercSegment, 0, erc_errorVar);
// 		if(img->mb_data[i-1].ei_flag)
// 			ercMarkCurrSegmentLost(dec_picture->size_x, erc_errorVar);
// 		else
// 			ercMarkCurrSegmentOK(dec_picture->size_x, erc_errorVar);
// 		
// 		//! call the right error concealment function depending on the frame type.
// 		erc_mvperMB /= dec_picture->PicSizeInMbs;
// 		
// 		erc_img = img;
// 		if(dec_picture->slice_type == I_SLICE || dec_picture->slice_type == SI_SLICE) // I-frame
// 			ercConcealIntraFrame(&recfr, dec_picture->size_x, dec_picture->size_y, erc_errorVar);
// 		else
// 			ercConcealInterFrame(&recfr, erc_object_list, dec_picture->size_x, dec_picture->size_y, erc_errorVar);
// 	}
	
	if (img->structure == FRAME)         // buffer mgt. for frame mode
		frame_postprocessing(img, input);
	else
		field_postprocessing(img, input);   // reset all interlaced variables
	
	structure  = dec_picture->structure;
	slice_type = dec_picture->slice_type;
	frame_poc  = dec_picture->frame_poc;
	refpic     = dec_picture->used_for_reference;
	
	store_picture_in_dpb_dec(dec_picture);
// 	dec_picture=NULL;ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	
	if (img->last_has_mmco_5)
	{
		img->pre_frame_num = 0;
	}
	
	if ((structure==FRAME)||structure==BOTTOM_FIELD)
	{
		
#ifdef WIN32
		_ftime (&(img->tstruct_end));             // start time ms
#else
		ftime (&(img->tstruct_end));              // start time ms
#endif
		
		time( &(img->ltime_end));                // start time s
		
		tmp_time=(img->ltime_end*1000+img->tstruct_end.millitm) - (img->ltime_start*1000+img->tstruct_start.millitm);
		tot_time=tot_time + tmp_time;
		
		if(slice_type == I_SLICE) // I picture
			fprintf(stdout,"%3d(I)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
			frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
		else if(slice_type == P_SLICE) // P pictures
			fprintf(stdout,"%3d(P)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
			frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
		else if(slice_type == SP_SLICE) // SP pictures
			fprintf(stdout,"%3d(SP) %3d %5d %7.4f %7.4f %7.4f %5d\n",
			frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
		else if (slice_type == SI_SLICE)
			fprintf(stdout,"%3d(SI) %3d %5d %7.4f %7.4f %7.4f %5d\n",
			frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
		else if(refpic) // stored B pictures
			fprintf(stdout,"%3d(BS) %3d %5d %7.4f %7.4f %7.4f %5d\n",
			frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
		else // B pictures
			fprintf(stdout,"%3d(B)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
			frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
		
		fflush(stdout);
		
		if(slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic)   // I or P pictures
			img->number++;
		else
			Bframe_ctr++;    // B pictures
		
		g_nFrame++;
	}
	
	img->current_mb_nr = -4712;   // impossible value for debugging, StW
	img->current_slice_nr = 0;
	
}

/*!
************************************************************************
* \brief
*    write the encoding mode and motion vectors of current 
*    MB to the buffer of the error concealment module.
************************************************************************
*/

void ercWriteMBMODEandMV(struct img_par *img,struct inp_par_dec *inp)
{
	extern objectBuffer_t *erc_object_list;
	int i, ii, jj, currMBNum = img->current_mb_nr;
	int mbx = xPosMB(currMBNum,dec_picture->size_x), mby = yPosMB(currMBNum,dec_picture->size_x);
	objectBuffer_t *currRegion, *pRegion;
	Macroblock *currMB = &img->mb_data[currMBNum];
	int***  mv;
	
	currRegion = erc_object_list + (currMBNum<<2);
	
	if(img->type != B_SLICE) //non-B frame
	{
		for (i=0; i<4; i++)
		{
			pRegion             = currRegion + i;
			pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
			currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  :
			currMB->b8mode[i]==0      ? REGMODE_INTER_COPY :
			currMB->b8mode[i]==1      ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);
			if (currMB->b8mode[i]==0 || currMB->b8mode[i]==IBLOCK)  // INTRA OR COPY
			{
				pRegion->mv[0]    = 0;
				pRegion->mv[1]    = 0;
				pRegion->mv[2]    = 0;
			}
			else
			{
				ii              = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
				jj              = 4*mby + (i/2)*2;
				if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS
				{
					pRegion->mv[0]  = (dec_picture->mv[LIST_0][ii][jj][0] + dec_picture->mv[LIST_0][ii+1][jj][0] + dec_picture->mv[LIST_0][ii][jj+1][0] + dec_picture->mv[LIST_0][ii+1][jj+1][0] + 2)/4;
					pRegion->mv[1]  = (dec_picture->mv[LIST_0][ii][jj][1] + dec_picture->mv[LIST_0][ii+1][jj][1] + dec_picture->mv[LIST_0][ii][jj+1][1] + dec_picture->mv[LIST_0][ii+1][jj+1][1] + 2)/4;
				}
				else // 16x16, 16x8, 8x16, 8x8
				{
					pRegion->mv[0]  = dec_picture->mv[LIST_0][ii][jj][0];
					pRegion->mv[1]  = dec_picture->mv[LIST_0][ii][jj][1];
					//          pRegion->mv[0]  = dec_picture->mv[LIST_0][4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][0];
					//          pRegion->mv[1]  = dec_picture->mv[LIST_0][4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][1];
				}
				erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
				pRegion->mv[2]    = dec_picture->ref_idx[LIST_0][ii][jj];
			}
		}
	}
	else  //B-frame
	{
		for (i=0; i<4; i++)
		{
			ii                  = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
			jj                  = 4*mby + (i/2)*2;
			pRegion             = currRegion + i;
			pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
			currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
			if (currMB->mb_type==I16MB || currMB->b8mode[i]==IBLOCK)  // INTRA
			{
				pRegion->mv[0]    = 0;
				pRegion->mv[1]    = 0;
				pRegion->mv[2]    = 0;
			}
			else
			{
				int idx = (dec_picture->ref_idx[0][ii][jj]<0)?1:0;
				//        int idx = (currMB->b8mode[i]==0 && currMB->b8pdir[i]==2 ? LIST_0 : currMB->b8pdir[i]==1 ? LIST_1 : LIST_0);
				//        int idx = currMB->b8pdir[i]==0 ? LIST_0 : LIST_1;
				mv                = dec_picture->mv[idx];
				pRegion->mv[0]    = (mv[ii][jj][0] + mv[ii+1][jj][0] + mv[ii][jj+1][0] + mv[ii+1][jj+1][0] + 2)/4;
				pRegion->mv[1]    = (mv[ii][jj][1] + mv[ii+1][jj][1] + mv[ii][jj+1][1] + mv[ii+1][jj+1][1] + 2)/4;
				erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
				
				pRegion->mv[2]  = (dec_picture->ref_idx[idx][ii][jj]);
				/*        
				if (currMB->b8pdir[i]==0 || (currMB->b8pdir[i]==2 && currMB->b8mode[i]!=0)) // forward or bidirect
				{
				pRegion->mv[2]  = (dec_picture->ref_idx[LIST_0][ii][jj]);
				///???? is it right, not only "img->fw_refFrArr[jj][ii-4]"
				}
				else
				{
				pRegion->mv[2]  = (dec_picture->ref_idx[LIST_1][ii][jj]);
				//          pRegion->mv[2]  = 0;
				}
				*/
			}
		}
	}
}

/*!
************************************************************************
* \brief
*    set defaults for old_slice
*    NAL unit of a picture"
************************************************************************
*/
void init_old_slice()
{
	old_slice.field_pic_flag = 0;
	
	old_slice.pps_id = INT_MAX;
	
	old_slice.frame_num = INT_MAX;
	
	old_slice.nal_ref_idc = INT_MAX;
	
	old_slice.idr_flag = 0;
	
	old_slice.pic_oder_cnt_lsb          = UINT_MAX;
	old_slice.delta_pic_oder_cnt_bottom = INT_MAX;
	
	old_slice.delta_pic_order_cnt[0] = INT_MAX;
	old_slice.delta_pic_order_cnt[1] = INT_MAX;
	
}

/*!
************************************************************************
* \brief
*    save slice parameters that are needed for checking of "first VCL
*    NAL unit of a picture"
************************************************************************
/*!
************************************************************************
* \brief
*    Initializes the parameters for a new slice and
*     allocates the memory for the coded slice in the Picture structure
*  \par Side effects:
*      Adds slice/partition header symbols to the symbol buffer
*      increments Picture->no_slices, allocates memory for the
*      slice, sets img->currSlice
************************************************************************
*/
static void init_slice (int start_mb_addr)
{
	int i;
	Picture *currPic = img->currentPicture;
	Slice *curr_slice;
	DataPartition *dataPart;
	Bitstream *currStream;
	
	img->current_mb_nr = start_mb_addr;
	
	// Allocate new Slice in the current Picture, and set img->currentSlice
	assert (currPic != NULL);
	currPic->no_slices++;
	if (currPic->no_slices >= MAXSLICEPERPICTURE)
		error ("Too many slices per picture, increase MAXLSICESPERPICTURE in global.h.", -1);
	currPic->slices[currPic->no_slices-1] = malloc_slice();//frame_pic.slices初始化（大多参数）
	curr_slice = currPic->slices[currPic->no_slices-1];
	img->currentSlice = curr_slice;
	
	curr_slice->picture_id = img->tr % 256;//0//2????
	curr_slice->qp = img->qp;//28
	curr_slice->start_mb_nr = start_mb_addr;//0
	curr_slice->slice_too_big = dummy_slice_too_big;
	
	for (i = 0; i < curr_slice->max_part_nr; i++)
	{
		dataPart = &(curr_slice->partArr[i]);
		if (input->symbol_mode == UVLC)
			dataPart->writeSyntaxElement = writeSyntaxElement_UVLC;
		else
			dataPart->writeSyntaxElement = writeSyntaxElement_CABAC;
		
		currStream = dataPart->bitstream;
		currStream->bits_to_go = 8;
		currStream->byte_pos = 0;
		currStream->byte_buf = 0;
	}
}



void exit_slice()
{
	
	old_slice.pps_id = img->currentSlice->pic_parameter_set_id;
	
	old_slice.frame_num = img->frame_num;
	
	old_slice.field_pic_flag = img->field_pic_flag;
	
	if(img->field_pic_flag)
	{
		old_slice.bottom_field_flag = img->bottom_field_flag;
	}
	
	old_slice.nal_ref_idc   = img->nal_reference_idc;
	
	old_slice.idr_flag = img->idr_flag;
	if (img->idr_flag)
	{
		old_slice.idr_pic_id = img->idr_pic_id;
	}
	
	if (active_sps->pic_order_cnt_type == 0)
	{
		old_slice.pic_oder_cnt_lsb          = img->pic_order_cnt_lsb;
		old_slice.delta_pic_oder_cnt_bottom = img->delta_pic_order_cnt_bottom;
	}
	
	if (active_sps->pic_order_cnt_type == 1)
	{
		old_slice.delta_pic_order_cnt[0] = img->delta_pic_order_cnt[0];
		old_slice.delta_pic_order_cnt[1] = img->delta_pic_order_cnt[1];
	}
}

/*!
************************************************************************
* \brief
*    detect if current slice is "first VCL NAL unit of a picture"
************************************************************************
*/
int is_new_picture()
{
	int result=0;
	
	result |= (old_slice.pps_id != img->currentSlice->pic_parameter_set_id);
	
	result |= (old_slice.frame_num != img->frame_num);
	
	result |= (old_slice.field_pic_flag != img->field_pic_flag);
	
	if(img->field_pic_flag && old_slice.field_pic_flag)
	{
		result |= (old_slice.bottom_field_flag != img->bottom_field_flag);
	}
	
	result |= (old_slice.nal_ref_idc   != img->nal_reference_idc);
	
	result |= ( old_slice.idr_flag != img->idr_flag);
	
	if (img->idr_flag && old_slice.idr_flag)
	{
		result |= (old_slice.idr_pic_id != img->idr_pic_id);
	}
	
	if (active_sps->pic_order_cnt_type == 0)
	{
		result |=  (old_slice.pic_oder_cnt_lsb          != img->pic_order_cnt_lsb);
		result |=  (old_slice.delta_pic_oder_cnt_bottom != img->delta_pic_order_cnt_bottom);
	}
	
	if (active_sps->pic_order_cnt_type == 1)
	{
		result |= (old_slice.delta_pic_order_cnt[0] != img->delta_pic_order_cnt[0]);
		result |= (old_slice.delta_pic_order_cnt[1] != img->delta_pic_order_cnt[1]);
	}
	
	return result;
}


// /*!
// ************************************************************************
// * \brief
// *    write the encoding mode and motion vectors of current 
// *    MB to the buffer of the error concealment module.
// ************************************************************************
// */
// 
// void ercWriteMBMODEandMV(struct img_par *img,struct inp_par_dec *inp)
// {
// 	extern objectBuffer_t *erc_object_list;
// 	int i, ii, jj, currMBNum = img->current_mb_nr;
// 	int mbx = xPosMB(currMBNum,dec_picture->size_x), mby = yPosMB(currMBNum,dec_picture->size_x);
// 	objectBuffer_t *currRegion, *pRegion;
// 	Macroblock *currMB = &img->mb_data[currMBNum];
// 	int***  mv;
// 	
// 	currRegion = erc_object_list + (currMBNum<<2);
// 	
// 	if(img->type != B_SLICE) //non-B frame
// 	{
// 		for (i=0; i<4; i++)
// 		{
// 			pRegion             = currRegion + i;
// 			pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
// 			currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  :
// 			currMB->b8mode[i]==0      ? REGMODE_INTER_COPY :
// 			currMB->b8mode[i]==1      ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);
// 			if (currMB->b8mode[i]==0 || currMB->b8mode[i]==IBLOCK)  // INTRA OR COPY
// 			{
// 				pRegion->mv[0]    = 0;
// 				pRegion->mv[1]    = 0;
// 				pRegion->mv[2]    = 0;
// 			}
// 			else
// 			{
// 				ii              = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
// 				jj              = 4*mby + (i/2)*2;
// 				if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS
// 				{
// 					pRegion->mv[0]  = (dec_picture->mv[LIST_0][ii][jj][0] + dec_picture->mv[LIST_0][ii+1][jj][0] + dec_picture->mv[LIST_0][ii][jj+1][0] + dec_picture->mv[LIST_0][ii+1][jj+1][0] + 2)/4;
// 					pRegion->mv[1]  = (dec_picture->mv[LIST_0][ii][jj][1] + dec_picture->mv[LIST_0][ii+1][jj][1] + dec_picture->mv[LIST_0][ii][jj+1][1] + dec_picture->mv[LIST_0][ii+1][jj+1][1] + 2)/4;
// 				}
// 				else // 16x16, 16x8, 8x16, 8x8
// 				{
// 					pRegion->mv[0]  = dec_picture->mv[LIST_0][ii][jj][0];
// 					pRegion->mv[1]  = dec_picture->mv[LIST_0][ii][jj][1];
// 					//          pRegion->mv[0]  = dec_picture->mv[LIST_0][4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][0];
// 					//          pRegion->mv[1]  = dec_picture->mv[LIST_0][4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][1];
// 				}
// 				erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
// 				pRegion->mv[2]    = dec_picture->ref_idx[LIST_0][ii][jj];
// 			}
// 		}
// 	}
// 	else  //B-frame
// 	{
// 		for (i=0; i<4; i++)
// 		{
// 			ii                  = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
// 			jj                  = 4*mby + (i/2)*2;
// 			pRegion             = currRegion + i;
// 			pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
// 			currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
// 			if (currMB->mb_type==I16MB || currMB->b8mode[i]==IBLOCK)  // INTRA
// 			{
// 				pRegion->mv[0]    = 0;
// 				pRegion->mv[1]    = 0;
// 				pRegion->mv[2]    = 0;
// 			}
// 			else
// 			{
// 				int idx = (dec_picture->ref_idx[0][ii][jj]<0)?1:0;
// 				//        int idx = (currMB->b8mode[i]==0 && currMB->b8pdir[i]==2 ? LIST_0 : currMB->b8pdir[i]==1 ? LIST_1 : LIST_0);
// 				//        int idx = currMB->b8pdir[i]==0 ? LIST_0 : LIST_1;
// 				mv                = dec_picture->mv[idx];
// 				pRegion->mv[0]    = (mv[ii][jj][0] + mv[ii+1][jj][0] + mv[ii][jj+1][0] + mv[ii+1][jj+1][0] + 2)/4;
// 				pRegion->mv[1]    = (mv[ii][jj][1] + mv[ii+1][jj][1] + mv[ii][jj+1][1] + mv[ii+1][jj+1][1] + 2)/4;
// 				erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
// 				
// 				pRegion->mv[2]  = (dec_picture->ref_idx[idx][ii][jj]);
// 				/*        
// 				if (currMB->b8pdir[i]==0 || (currMB->b8pdir[i]==2 && currMB->b8mode[i]!=0)) // forward or bidirect
// 				{
// 				pRegion->mv[2]  = (dec_picture->ref_idx[LIST_0][ii][jj]);
// 				///???? is it right, not only "img->fw_refFrArr[jj][ii-4]"
// 				}
// 				else
// 				{
// 				pRegion->mv[2]  = (dec_picture->ref_idx[LIST_1][ii][jj]);
// 				//          pRegion->mv[2]  = 0;
// 				}
// 				*/
// 			}
// 		}
// 	}
// }
// 
/*!
************************************************************************
* \brief
*    decodes one slice
************************************************************************
*/
void decode_one_slice(struct img_par *img,struct inp_par_dec *inp)
{
	
	Boolean end_of_slice = FALSE1;
	int read_flag;
	img->cod_counter=-1;
	
// 	set_ref_pic_num_dec();
	
	if (img->type == B_SLICE)
		compute_collocated(Co_located, listX);
	
	//reset_ec_flags();
	
	while (end_of_slice == FALSE1) // loop over macroblocks
	{
		
#if TRACE
		fprintf(p_trace,"\n*********** POC: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->ThisPOC, img->current_mb_nr, img->current_slice_nr, img->type);
#endif
		
		// Initializes the current macroblock
		start_macroblock_dec(img,inp, img->current_mb_nr);
		// Get the syntax elements from the NAL
		read_flag = read_one_macroblock(img,inp); //1   看img的赋值。
		decode_one_macroblock_dec(img,inp);
		
		if(img->MbaffFrameFlag && dec_picture->mb_field[img->current_mb_nr])
		{
			img->num_ref_idx_l0_active >>= 1;
			img->num_ref_idx_l1_active >>= 1;
		}
		
		ercWriteMBMODEandMV(img,inp);
		
		end_of_slice=exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2));//MBAmap的信息
	}
	
	exit_slice();
	//reset_ec_flags();
	
}


void decode_slice(struct img_par *img,struct inp_par_dec *inp, int current_header)
{
	Slice *currSlice = img->currentSlice;
	
	if (active_pps->entropy_coding_mode_flag)
	{
		init_contexts_dec (img);
		cabac_new_slice();
	}
	
	if ( (active_pps->weighted_bipred_idc > 0  && (img->type == B_SLICE)) || (active_pps->weighted_pred_flag && img->type !=I_SLICE))
		fill_wp_params(img);
	
	//printf("frame picture %d %d %d\n",img->structure,img->ThisPOC,img->direct_type);
	
	
	// decode main slice information
	if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
		decode_one_slice(img,inp);//查img来源路线
    
	// setMB-Nr in case this slice was lost
	//  if(currSlice->ei_flag)  
	//    img->current_mb_nr = currSlice->last_mb_nr + 1;
	
}


/*!
************************************************************************
* \brief
*    Prepare field and frame buffer after frame decoding
************************************************************************
*/
void frame_postprocessing(struct img_par *img, struct inp_par_dec *inp)
{
}

/*!
************************************************************************
* \brief
*    Prepare field and frame buffer after field decoding
************************************************************************
*/
void field_postprocessing(struct img_par *img, struct inp_par_dec *inp)
{
	img->number /= 2;
}



void reset_wp_params(struct img_par *img)
{
	int i,comp;
	int log_weight_denom;
	
	for (i=0; i<MAX_REFERENCE_PICTURES; i++)
	{
		for (comp=0; comp<3; comp++)
		{
			log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
			img->wp_weight[0][i][comp] = 1<<log_weight_denom;
			img->wp_weight[1][i][comp] = 1<<log_weight_denom;
		}
	}
}


void fill_wp_params(struct img_par *img)
{
	int i, j, k;
	int comp;
	int log_weight_denom;
	int p0, pt;
	int bframe = (img->type==B_SLICE);
	int max_bwd_ref, max_fwd_ref;
	int x,z;
	
	max_fwd_ref = img->num_ref_idx_l0_active;
	max_bwd_ref = img->num_ref_idx_l1_active;
	
	if (active_pps->weighted_bipred_idc == 2 && bframe)
	{
		img->luma_log2_weight_denom = 5;
		img->chroma_log2_weight_denom = 5;
		img->wp_round_luma = 16;
		img->wp_round_chroma = 16;
		
		for (i=0; i<MAX_REFERENCE_PICTURES; i++)
		{
			for (comp=0; comp<3; comp++)
			{
				log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
				img->wp_weight[0][i][comp] = 1<<log_weight_denom;
				img->wp_weight[1][i][comp] = 1<<log_weight_denom;
				img->wp_offset[0][i][comp] = 0;
				img->wp_offset[1][i][comp] = 0;
			}
		}
	}
	
	if (bframe)
	{
		for (i=0; i<max_fwd_ref; i++)
		{
			for (j=0; j<max_bwd_ref; j++)
			{
				for (comp = 0; comp<3; comp++)
				{
					log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
					if (active_pps->weighted_bipred_idc == 1)
					{
						img->wbp_weight[0][i][j][comp] =  img->wp_weight[0][i][comp];
						img->wbp_weight[1][i][j][comp] =  img->wp_weight[1][j][comp];
					}
					else if (active_pps->weighted_bipred_idc == 2)
					{
						pt = listX[LIST_1][j]->poc - listX[LIST_0][i]->poc;
						if (pt == 0 || listX[LIST_1][j]->is_long_term || listX[LIST_0][i]->is_long_term)
						{
							img->wbp_weight[0][i][j][comp] =   32;
							img->wbp_weight[1][i][j][comp] =   32;
						}
						else
						{
							p0 = img->ThisPOC - listX[LIST_0][i]->poc;
							
							x = (16384 + (pt>>1))/pt;
							z = Clip3(-1024, 1023, (x*p0 + 32 )>>6);
							
							img->wbp_weight[1][i][j][comp] = z >> 2;
							img->wbp_weight[0][i][j][comp] = 64 - img->wbp_weight[1][i][j][comp];
							if (img->wbp_weight[1][i][j][comp] < -64 || img->wbp_weight[1][i][j][comp] > 128)
							{
								img->wbp_weight[0][i][j][comp] = 32;
								img->wbp_weight[1][i][j][comp] = 32;
								img->wp_offset[0][i][comp] = 0;
								img->wp_offset[1][i][comp] = 0;
							}
						}
					}
				}
			}
		}
	}
	
	if (bframe && img->MbaffFrameFlag)
	{
		for (i=0; i<2*max_fwd_ref; i++)
		{
			for (j=0; j<2*max_bwd_ref; j++)
			{
				for (comp = 0; comp<3; comp++)
				{
					for (k=2; k<6; k+=2)
					{
						img->wp_offset[k+0][i][comp] = img->wp_offset[0][i/2][comp];
						img->wp_offset[k+1][i][comp] = img->wp_offset[1][i/2][comp];
						
						log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
						if (active_pps->weighted_bipred_idc == 1)
						{
							img->wbp_weight[k+0][i][j][comp] =  img->wp_weight[0][i/2][comp];
							img->wbp_weight[k+1][i][j][comp] =  img->wp_weight[1][j/2][comp];
						}
						else if (active_pps->weighted_bipred_idc == 2)
						{
							pt = listX[k+LIST_1][j]->poc - listX[k+LIST_0][i]->poc;
							if (pt == 0 || listX[k+LIST_1][j]->is_long_term || listX[k+LIST_0][i]->is_long_term)
							{
								img->wbp_weight[k+0][i][j][comp] =   32;
								img->wbp_weight[k+1][i][j][comp] =   32;
							}
							else
							{
								p0 = ((k==2)?img->toppoc:img->bottompoc) - listX[k+LIST_0][i]->poc;
								
								x = (16384 + (pt>>1))/pt;
								z = Clip3(-1024, 1023, (x*p0 + 32 )>>6);
								
								img->wbp_weight[k+1][i][j][comp] = z >> 2;
								img->wbp_weight[k+0][i][j][comp] = 64 - img->wbp_weight[k+1][i][j][comp];
								if (img->wbp_weight[k+1][i][j][comp] < -64 || img->wbp_weight[k+1][i][j][comp] > 128)
								{
									img->wbp_weight[k+1][i][j][comp] = 32;
									img->wbp_weight[k+0][i][j][comp] = 32;
									img->wp_offset[k+0][i][comp] = 0;
									img->wp_offset[k+1][i][comp] = 0;
								}
							}
						}
					}
				}
			}
		}
	}
}

/////////////////////////////////////<
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////