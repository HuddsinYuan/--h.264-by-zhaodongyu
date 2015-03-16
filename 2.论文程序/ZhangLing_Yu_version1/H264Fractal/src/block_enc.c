#include "block_enc.h"
#include "compute.h"
#include "i_global.h"
#include "image.h"

extern int search_range;
int max_search_points;
typedef struct 
{
	char   x;
	char   y;
} MV;
// unsigned int S_P;
int WINDOW_SIZE;
//WIDOW_SIZE 7  cross-diamond search
extern short*  spiral_search_x;
extern short*  spiral_search_y;
extern short*  spiral_hpel_search_x;
extern short*  spiral_hpel_search_y;
void no_mem_exit1(char *where)
{
	printf("Could not allocate memory: %s",where);
	exit(0);
}
void initialcode(void)
{
	int k,l,i;
	max_search_points= imax(9, (2 * search_range + 1) * (2 * search_range + 1));
	if ((spiral_search_x = (short*)calloc(max_search_points, sizeof(short))) == NULL)
		no_mem_exit1("Init_Motion_Search_Module: spiral_search_x");
	if ((spiral_search_y = (short*)calloc(max_search_points, sizeof(short))) == NULL)
		no_mem_exit1("Init_Motion_Search_Module: spiral_search_y");
	if ((spiral_hpel_search_x = (short*)calloc(max_search_points, sizeof(short))) == NULL)
		no_mem_exit1("Init_Motion_Search_Module: spiral_hpel_search_x");
	if ((spiral_hpel_search_y = (short*)calloc(max_search_points, sizeof(short))) == NULL)
		no_mem_exit1("Init_Motion_Search_Module: spiral_hpel_search_y");
	spiral_search_x[0] = spiral_search_y[0] = 0;
	spiral_hpel_search_x[0] = spiral_hpel_search_y[0] = 0;
	
	for (k=1,l=1; l <= imax(1,search_range); l++)
	{
		for (i=-l+1; i< l; i++)
		{
			spiral_search_x[k] =  i;//整像素
			spiral_search_y[k] = -l;
			spiral_hpel_search_x[k] =  i<<1;//半像素
			spiral_hpel_search_y[k++] = -l<<1;
			
			spiral_search_x[k] =  i;
			spiral_search_y[k] =  l;
			spiral_hpel_search_x[k] =  i<<1;
			spiral_hpel_search_y[k++] =  l<<1;
		}
		for (i=-l;i<=l; i++)
		{
			spiral_search_x[k] = -l;
			spiral_search_y[k] =  i;
			spiral_hpel_search_x[k] = -l<<1;
			spiral_hpel_search_y[k++] = i<<1;
			
			spiral_search_x[k] =  l;
			spiral_search_y[k] =  i;
			spiral_hpel_search_x[k] =  l<<1;
			spiral_hpel_search_y[k++] = i<<1;
		}
}
}
void init_searchedpt(void)
{ 
	int i,j,k=WINDOW_SIZE*2+1;
	for(i=0;i<k;i++)
	{
		for(j=0;j<k;j++)
		{
			searchedpt[i][j]=0;
		}
	} 
}

void changeReferenceFrame(char video)
{
	switch(video)
	{
	case 'C':		
// 		imgY_ref_temp = imgY_ref;
// 		imgUV_ref_temp = imgUV_ref;

		imgY_ref_temp = imgY_ref/*_c*/;
		imgUV_ref_temp = imgUV_ref/*_c*/;
		plane_Y_domain_temp = plane_Y_domain;//plane_domain应该和img_ref同步更新，在decode_oneframe里
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
		imgY_ref_temp = imgY_ref_h;
		imgUV_ref_temp = imgUV_ref_h;
		
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
		imgY_ref_temp = imgY_ref_m;
		imgUV_ref_temp = imgUV_ref_m;
		
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
		imgY_ref_temp = imgY_ref_n;
		imgUV_ref_temp = imgUV_ref_n;
		
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
}

void SetRefAndMotionVectors_fract(int mode,int CurrMb,TRANS_NODE *trans,int con)//mode=0是初始化,mode=1: 16x16存储；mode=2: 16x18存储；mode=3: 8x16存储；mode=4: 8x8存储
{
	int block16_x_n,block16_y_n,block16_x,block16_y;
	int i,j,k,m,n;
	int Y,U,V;
	int yuv;
	Macroblock *currMB = &img->mb_data[img->current_mb_nr];
	
	block16_x_n=img->width/16;//图像一行包含多少个宏块
	block16_y_n=img->height/16;//图像一列包含多少个宏块
	
	block16_x=CurrMb%block16_x_n;//当前宏块在图像中的x坐标（以宏块为单位）
	block16_y=CurrMb/block16_x_n;//当前宏块在图像中的y坐标（以宏块为单位）
	Y=U=V=0;
	yuv=4;
	if (2==con)//U分量
	{
// 		con=1;
		U=img->width/4;//88
		yuv=2;
	} 
	if (3==con)//V分量
	{
// 		con=1;
		U=img->width/4;//88
		V=img->height/8;//36
		yuv=2;
	}//else
// 	{
// 		con=0;
// 	}
	if(0==mode)//初始化16*16宏块，以4*4为基本单位，因此要循环16次
	{
// 		currMB->mb_type = mode;

		for(j=0;j<yuv;j++)
			for(i=0;i<yuv;i++)//U,V分量存储在Y分量的后边
			{
				
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=0;//将MV初始化为（0,0）
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=0;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=0;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=0;//将list0，list1参考索引初始化为0
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=0;
			}
// 		printf("\nmode=%d",mode);
	}
	if(1==mode)//16*16，以4*4为基本单位，因此要循环16次，16次都一样
	{
// 		printf("\nmode=%d",mode);
		currMB->mb_type = mode;
		for(i=0;i<yuv;i++)
		{
			currMB->b8mode[i] = mode;//什么意思，为什么循环四次？
// 			currMB->b8pdir[i] = (bframe?direct_pdir[img->block_x+(i%2)*2][img->block_y+(i/2)*2]:0);
		}
		for(j=0;j<yuv;j++)
			for(i=0;i<yuv;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->offset/5;
// 				printf("\nscale=%ld",enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]);
// 				printf("offset=%ld",enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]);
			}

// 		printf("\nscale=%.0lf",trans->scale);
// 		printf("offset=%.0lf",trans->offset);

	}
	
	if(2==mode)//16*8
	{
// 		printf("\nmode=%d",mode);
		currMB->mb_type = mode;
		for(i=0;i<yuv;i++)
		{
			currMB->b8mode[i] = mode;
// 			currMB->b8pdir[i] = (bframe?direct_pdir[img->block_x+(i%2)*2][img->block_y+(i/2)*2]:0);
		}
		for(j=0;j<yuv;j++)//左边16*8块
			for(i=0;i<2;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[0].x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[0].y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].offset/5;
			}
// 		printf("\ntrans->next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[0].scale,trans->next[0].offset);


		for(j=0;j<yuv;j++)//右边16*8块
			for(i=2;i<yuv;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[1].x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[1].y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[1].reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[1].scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[1].offset/5;
			}
//		printf("\ntrans->next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[1].scale,trans->next[1].offset);

	}
	
	if(3==mode)//8*16
	{
// 		printf("\nmode=%d",mode);
		currMB->mb_type = mode;
		for(i=0;i<yuv;i++)
		{
			currMB->b8mode[i] = mode;
// 			currMB->b8pdir[i] = (bframe?direct_pdir[img->block_x+(i%2)*2][img->block_y+(i/2)*2]:0);
		}
		for(j=0;j<2;j++)//上边8*16块
			for(i=0;i<yuv;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[0].x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[0].y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].offset/5;
			}
// 		printf("\ntrans->next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[0].scale,trans->next[0].offset);

		for(j=2;j<yuv;j++)//下边8*16块
			for(i=0;i<yuv;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[1].x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[1].y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[1].reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[1].scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[1].offset/5;
			}
//		printf("\ntrans->next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[1].scale,trans->next[1].offset);

	}
	
	
	if(4==mode)//p8*8=8还要存每个8*8块的划分模式
	{
// 		printf("\nmode=%d",mode);
		currMB->mb_type = mode+4;

		for(m=0;m<yuv/2;m++)//列
			for(n=0;n<yuv/2;n++)//行
			{
				k=2*m+n;//k从0~3，分别表示4个8*8块
// 				printf("\nk=%d",k);
				if(0==trans->next[k].partition)//8*8不再分了
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k]:scale=.0lf,offset=.0lf\n",trans->next[k].scale,trans->next[k].offset);
					
					currMB->b8mode[k]   = 4;//8*8不再分了

					for(j=2*m;j<2*m+2;j++)//应该是以4*4为基本单位存储的吧，所以循环了4次？？？
						for(i=2*n;i<2*n+2;i++)
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].offset/5;
						}

				}//if
				
				if(1==trans->next[k].partition)//8*8被划分成两个8*4
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k].next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[0].scale,trans->next[k].next[0].offset);
// 					printf("\ntrans->next[k].next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[1].scale,trans->next[k].next[1].offset);

					currMB->b8mode[k]   = 5;//8*8被划分成两个8*4

// 						currMB->b8pdir[i]   = b8pdir[i];

					for(j=2*m;j<2*m+2;j++)
						for(i=2*n;i<2*n+1;i++)//存储第一个8*4块
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[0].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[0].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].offset/5;
						}
					for(j=2*m;j<2*m+2;j++)
						for(i=2*n+1;i<2*n+2;i++)//存储第二个8*4块
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[1].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[1].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].offset/5;
						}
				}//if
				
				
				if(2==trans->next[k].partition)//8*8被划分成两个4*8
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k].next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[0].scale,trans->next[k].next[0].offset);
// 					printf("\ntrans->next[k].next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[1].scale,trans->next[k].next[1].offset);

					
					currMB->b8mode[k]   = 6;//8*8被划分成两个4*8

					for(j=2*m;j<2*m+1;j++)
						for(i=2*n;i<2*n+2;i++)//存储第一个4*8块
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[0].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[0].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].offset/5;
						}
					for(j=2*m+1;j<2*m+2;j++)
						for(i=2*n;i<2*n+2;i++)//存储第二个4*8块
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[1].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[1].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].offset/5;
						}
				}//if
				
				if(3==trans->next[k].partition)//8*8被划分成四个4*4
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k].next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[0].scale,trans->next[k].next[0].offset);
// 					printf("\ntrans->next[k].next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[1].scale,trans->next[k].next[1].offset);
// 					printf("\ntrans->next[k].next[2]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[2].scale,trans->next[k].next[2].offset);
// 					printf("\ntrans->next[k].next[3]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[3].scale,trans->next[k].next[3].offset);
					
					currMB->b8mode[k]   = 7;//8*8被划分成四个4*4

					
					for(j=2*m;j<2*m+1;j++)
						for(i=2*n;i<2*n+1;i++)//存储第一个4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[0].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[0].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].offset/5;
						}
					for(j=2*m;j<2*m+1;j++)
						for(i=2*n+1;i<2*n+2;i++)//存储第二个4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[1].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[1].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].offset/5;
						}
					for(j=2*m+1;j<2*m+2;j++)
						for(i=2*n;i<2*n+1;i++)//存储第三个4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[2].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[2].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[2].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[2].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[2].offset/5;
						}
					for(j=2*m+1;j<2*m+2;j++)
						for(i=2*n+1;i<2*n+2;i++)//存储第四个4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[3].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[3].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[3].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[3].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[3].offset/5;
						}
									
				}
			}
	}
}




/**
   brief:对一个宏块进行编码
   param:CurMb 当前宏块标记
   param:trans 保存处理后的结果
   param:con 当前分量
*/
void encode_one_macroblock(int CurMb,TRANS_NODE **trans,int con)//con=1：Y分量；con=2:U分量；con=3:V分量
{
    int i,j,k,ii=0,q=2;        
	int mode;//块划分模式
	byte **plane;
	double tol=input->tol_16;//配置文件输入的tol_16
	int mb_X = CurMb%(img->frmWidthInMbs/min(2,con));//将宏块的光栅扫描坐标，转化成二维x坐标(以宏块为单位)
	int mb_Y = CurMb/(img->frmWidthInMbs/min(2,con));//将宏块的光栅扫描坐标，转化成二维y坐标(以宏块为单位)
	int b_x=mb_X*16;//以像素为单位的宏块左上角x坐标
	int b_y=mb_Y*16;//以像素为单位的宏块左上角y坐标
	double rms,chun,sumR=0,sumD=0,r=0,d=0,sR=0,sD=0,mr=0,md=0,MR=0,MD=0;
	int rectRms=0,i16mode;
	double R[256];
    double D[256];

	obj=region=0;
    

// 	if (con==1)
// 	{tol=input->tol_16;
// 	} 
//     else
// 	{tol=0;//input->tol_16/10000000000000;
// 	}


	if(num_regions==2)//no //if OB ，判断当前宏块属于哪一对象
	{
		switch(con)
		{
		case 1:
			plane=plane_Y;break;
		case 2:
			plane=plane_UV[0];break;
		case 3:
			plane=plane_UV[1];break;
		}
		for(j=b_y;j<b_y+16;j++)
			for(i=b_x;i<b_x+16;i++)
			{
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //判断当前块属于哪个region
			}	
		region--;	//0-region0;1--region1;2--边界块
		if(region==2)
		{
			trans[0][CurMb].use=trans[1][CurMb].use=1;
			trans[0][CurMb].region=trans[1][CurMb].region=1;
		}
		else
		{
			trans[region][CurMb].use=1;
			trans[1-region][CurMb].use=0;
			trans[0][CurMb].region=trans[1][CurMb].region=0;//是不是基于对象的
			obj=region;
		}
	}

	for(k=region;k<pow(2.0,region);k++) //region=0
	{
		if(region>1)
			obj=k-2;
		trans[obj][CurMb].partition=0;//初始化trans
		trans[obj][CurMb].block_type = 0;
		trans[obj][CurMb].x=0;
        trans[obj][CurMb].y=0;	
	    trans[obj][CurMb].reference=0;
        SetRefAndMotionVectors_fract(0,CurMb,&trans[obj][CurMb],con);//注意第一个参数为0,表示初始化参考索引和分形MV为0

		if (search_mode==0)//全搜索 16x16块
		{	
			rms =full_search(b_x,b_y,16,16,con,&trans[obj][CurMb]);  //利用前一帧作为参考计算出trans
		}
		if (search_mode==1)//新六边形搜索 16x16块
		{
			rms =new_hexagon_block_search(b_x,b_y,16,16,con,&trans[obj][CurMb]);  //利用前一帧作为参考计算出trans
			
		}
		if (search_mode==2)//UMHEX搜索 16x16块
		{
			rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans[obj][CurMb]);  //利用前一帧作为参考计算出trans
			UpCur.mv_x=trans[obj][CurMb].x;
			UpCur.mv_y=trans[obj][CurMb].y;
		}
		if (search_mode==3)//六边形搜索 16x16块
		{
			rms =hexagon_block_search(b_x,b_y,16,16,con,&trans[obj][CurMb]);
		}



		/////////////////////////////////////////////////////////////////////////////////////////////////////
		//预搜索限制条件，范数的判断，只有满足条件，才进行后面的搜索///////////////////////////////////////////////////////	 		
		//预搜索限制条件，范数的判断，只有满足条件，才进行后面的搜索///////////////////////////////////////////////////////	 		

			/*	for(ii=0;ii<256;ii++)
			{
			mr+=pow(fabs((R[ii]-r)/(sqrt(sR))),q/(q-1));
			//mr+=((R[ii]-r)/(sqrt(sR)))*((R[ii]-r)/(sqrt(sR)));
			//	chu+=((R[ii]-r)/(sqrt(sR)))*((D[ii]-d)/(sqrt(sD)));
			//	chu+=(((R[ii]-r)/(sqrt(sR)))*mr-((D[ii]-d)/(sqrt(sD))))*(((R[ii]-r)/(sqrt(sR)))*mr-((D[ii]-d)/(sqrt(sD))));
			md+=pow(fabs((D[ii]-d)/(sqrt(sD))),q);
			//	md+=((D[ii]-d)/(sqrt(sD)))*((D[ii]-d)/(sqrt(sD)));
			//	chun+=(R[ii]-r)/(sqrt(sR));
			//mrd+=(R[ii]-r)*(R[ii]-r);
			}
			MR=pow(mr,2*(q-1)/q);
			MD=pow(md,2/q);
			
			  //	chun=MR*MD;
			//chu=fabs(1-MR*MD);*/
			//开始判断是否进行下面的分形处理///////////////////////////////////////////////////////////////
			//if((fabs(1-MR*MD))<14.4)
			//  if(fabs(MR*MD)>10)    
			//	if(chun<0.9999999)
			//对这里进行不同的尝试，得出最佳的结果，也有让chun>=0的地方
    	if (currentVideo=='C')////如果是C目，则利用分数像素计算RMS
// 		if (currentVideo!='C')////撤销掉分数像素ZZZZZZZZZZZZZZZZZZZ
		{	
			double rms_rl;
			TRANS_NODE trans_rl;
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('H');//为什么参考H目了？？？
			if (search_mode==0)//全搜索 16x16块
			{	
				rms_rl =full_search(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出当前range块的RMS，分形参数s,o,x,y
			}
			if (search_mode==1)//新六边形搜索 16x16块
			{
				rms_rl =new_hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
				
			}
			if (search_mode==2)//UMHEX搜索 16x16块
			{
				rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
				UpCur.mv_x=trans[obj][CurMb].x;
				UpCur.mv_y=trans[obj][CurMb].y;
			}
			if (search_mode==3)//六边形搜索 16x16块
			{
				rms_rl =hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);
			}

			changeReferenceFrame('C');
			if (rms_rl<rms)//更新RMS误差和分形参数trans
			{
				num_rl++;
				trans[obj][CurMb].reference=1;//参考块是H目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
}
			else
				trans[obj][CurMb].reference=0;//参考块是C目
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('M');//为什么参考M目了？？？
			if (search_mode==0)//全搜索 16x16块
			{	
				rms_rl =full_search(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
			}
			if (search_mode==1)//新六边形搜索 16x16块
			{
				rms_rl =new_hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
				
			}
			if (search_mode==2)//UMHEX搜索 16x16块
			{
				rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
				UpCur.mv_x=trans[obj][CurMb].x;
				UpCur.mv_y=trans[obj][CurMb].y;
			}
			if (search_mode==3)//六边形搜索 16x16块
			{
				rms_rl =hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);
			}

			changeReferenceFrame('C');
			if (rms_rl<rms)//更新trans
			{
				num_rl++;
				trans[obj][CurMb].reference=2;//参考块是M目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('N');//为什么参考N目了？??
			if (search_mode==0)//全搜索 16x16块
			{	
				rms_rl =full_search(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
			}
			if (search_mode==1)//新六边形搜索 16x16块
			{
				rms_rl =new_hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
				
			}
			if (search_mode==2)//UMHEX搜索 16x16块
			{
				rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans_rl);  //利用前一帧作为参考计算出trans
				UpCur.mv_x=trans[obj][CurMb].x;
				UpCur.mv_y=trans[obj][CurMb].y;
			}
			if (search_mode==3)//六边形搜索 16x16块
			{
				rms_rl =hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);
			}
			changeReferenceFrame('C');
			if (rms_rl<rms)//更新trans
			{
				num_rl++;
				trans[obj][CurMb].reference=3;//参考块是N目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
		}
		
		if (currentVideo!='C')////如果不是C目，则利用快速搜索计算RMS
			//  else                   //如果不是C目,则进行快速搜索
		{
			double rms_rl;
			TRANS_NODE trans_rl;
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('C');
			//rms_rl = full_search(b_x,b_y,16,16,con,&trans_rl);
			rms_rl = full_search_R(b_x,b_y,16,16,con,&trans_rl);//参考C目
			changeReferenceFrame(currentVideo);
			if (rms_rl<rms)
			{
				num_rl++;
				trans[obj][CurMb].reference=0;//参考块是C目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
			else
			{	trans[obj][CurMb].block_type=0;
			trans[obj][CurMb].reference=4;
			}
			
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('H');
			rms_rl = full_search_R(b_x,b_y,16,16,con,&trans_rl);
			changeReferenceFrame(currentVideo);
			if (rms_rl<rms)
			{
				num_rl++;
				trans[obj][CurMb].reference=1;//参考块是H目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
			
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('M');
			rms_rl = full_search_R(b_x,b_y,16,16,con,&trans_rl);
			changeReferenceFrame(currentVideo);
			if (rms_rl<rms)
			{
				num_rl++;
				trans[obj][CurMb].reference=2;//参考块是M目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
			
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('N');
			//rms_rl = full_search(b_x,b_y,16,16,con,&trans_rl);
			rms_rl = full_search_R(b_x,b_y,16,16,con,&trans_rl);
			changeReferenceFrame(currentVideo);
			if (rms_rl<rms)
			{
				num_rl++;
				trans[obj][CurMb].reference=3;//参考块是N目
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}

		}


  		if(k==0||k==1)//对象或者背景//yes
		{
			for(j=b_x;j<b_x+16;j++)
				for(i=b_y;i<b_y+16;i++)     //b_x,b_y分别表示当前16*16 range块对应的坐标
				{ 
					switch(con)
					{
					case 1://Y分量
						R[ii]=imgY_org[i][j];
						D[ii]=imgY_ref_temp[i][j];
						break;
					case 2://U分量
						R[ii]=imgUV_org[0][i][j];
						D[ii]=imgUV_ref_temp[0][i][j];
						break;
					case 3://V分量
						R[ii]=imgUV_org[1][i][j];
						D[ii]=imgUV_ref_temp[1][i][j];
						//	D[ii]=imgY_ref[i][j];
						break;
					}
					sumR+=R[ii];//计算R块的像素和
					sumD+=D[ii];//计算D块的像素和
					ii++;
				}
			r=sumR/256;  //计算R块像素平均值
			d=sumD/256; //计算D块像素平均值
			for(ii=0;ii<256;ii++)
			{
				sR+=(R[ii]-r)*(R[ii]-r);//计算R块 各像素减去均值的差 的平方
				sD+=(D[ii]-d)*(D[ii]-d);//计算D块 各像素减去均值的差 的平方
			}
			for(ii=0;ii<256;ii++)
			{
				mr+=((R[ii]-r)/(sqrt(sR)))*((D[ii]-d)/(sqrt(sD)));
			}
			chun=mr*mr;//上面的计算都是为了chun


       		if(chun<=1&&chun>=0.9&&rms> (tol*tol*no))//如果16x16没有匹配，继续向下划分,这里用到了chun和RMS两个判据
			{
 				for(mode=1;mode<3;mode++)//16x8划分和8x16划分
				{
					rectRms=0;
					trans[obj][CurMb].partition=mode;
					for(i=0;i<2;i++)//每种划分都会将一个宏块划分为两个小块，i=0是第一个小块，i=1是第二个小块
					{				    
						trans[obj][CurMb].next[i].x=0;//next是指什么呢？？？
		        	    trans[obj][CurMb].next[i].y=0;
						if(encode_block_rect(b_x,b_y,i,&trans[obj][CurMb].next[i],con,mode,1)==0)//如果矩形匹配（16x8或8x16）失败
						{
							i=4;
						}
						else
						{
							rectRms++;
						}
					}
					if(rectRms==2)      //说明每个分块都匹配，停止mode循环，注意必须每个分块都匹配才行
					{
						partition[obj][partition_length[obj]]=mode;
						partition_length[obj]++;

						if (1==mode)//记录16*8
						{
						    SetRefAndMotionVectors_fract(2,CurMb,&trans[obj][CurMb],con);//16*8存储
						}

						if (2==mode)//记录8*16
						{
							SetRefAndMotionVectors_fract(3,CurMb,&trans[obj][CurMb],con);//8*16存储
						}

//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 						xy_p[obj][trans[obj][CurMb].next[0].x+input->search_range]++;
// 						xy_p[obj][trans[obj][CurMb].next[1].x+input->search_range]++;
// 						xy_p[obj][trans[obj][CurMb].next[0].y+input->search_range]++;
// 						xy_p[obj][trans[obj][CurMb].next[1].y+input->search_range]++;
// 
// 						s_p[obj][(int)((trans[obj][CurMb].next[0].scale-MIN_ALPHA)*100/5)]++;
// 						s_p[obj][(int)((trans[obj][CurMb].next[1].scale-MIN_ALPHA)*100/5)]++;
// 
// 						o_p[obj][(int)(trans[obj][CurMb].next[0].offset-MIN_BETA)/5]++;
// 						o_p[obj][(int)(trans[obj][CurMb].next[1].offset-MIN_BETA)/5]++;
// 
// 						x_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[0].x+search_range;
// 						y_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[0].y+search_range;
// 						s_trans[obj][trans_count[obj]]=(int)((trans[obj][CurMb].next[0].scale-MIN_ALPHA)*100)/5;
// 						o_trans[obj][trans_count[obj]]=(int)(trans[obj][CurMb].next[0].offset-MIN_BETA)/5;
// 						trans_reference[obj][trans_count[obj]] = trans[obj][CurMb].next[0].reference;
// 						trans_count[obj]++;
// 
// 						x_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[1].x+search_range;
// 						y_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[1].y+search_range;
// 						s_trans[obj][trans_count[obj]]=(int)((trans[obj][CurMb].next[1].scale-MIN_ALPHA)*100)/5;
// 						o_trans[obj][trans_count[obj]]=(int)(trans[obj][CurMb].next[1].offset-MIN_BETA)/5;
// 						trans_reference[obj][trans_count[obj]] = trans[obj][CurMb].next[1].reference;
// 						trans_count[obj]++;
// 					
// 						mode=4;
// 						if(region!=1)
// 							num_rect1+=2;
					}
				}
				if(mode<4)//16*16，16*8，8*16方式划分
				{
					trans[obj][CurMb].partition=3;//8*8划分
					partition[obj][partition_length[obj]]=trans[obj][CurMb].partition;
		    	    partition_length[obj]++;				
					for(i=0;i<2;i++)
						for(j=0;j<2;j++)
						{
							encode_block_8(b_x+j*8,b_y+i*8,i*2+j,&trans[obj][CurMb].next[i*2+j],con);//循环编码4个8*8块
						}

					if (3==mode)//记录8*8
					{
						SetRefAndMotionVectors_fract(4,CurMb,&trans[obj][CurMb],con);//8*8存储(含8*4，4*8，4*4)
					}

				
				
				}
			}
    		else    //16x16块就满足要求了
			{
 				if (region!=1)
 					num_16++;
				partition[obj][partition_length[obj]]=trans[obj][CurMb].partition;
				partition_length[obj]++;


				SetRefAndMotionVectors_fract(1,CurMb,&trans[obj][CurMb],con);//16*16存储

//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 				xy_p[obj][trans[obj][CurMb].x+input->search_range]++;//通过+7，使得肯定是正数？
// 				xy_p[obj][trans[obj][CurMb].y+input->search_range]++;
// 				s_p[obj][(int)((trans[obj][CurMb].scale-MIN_ALPHA)*100)/5]++;
// 				o_p[obj][(int)(trans[obj][CurMb].offset-MIN_BETA)/5]++;
// 
// 				x_trans[obj][trans_count[obj]]=trans[obj][CurMb].x+search_range;
// 				y_trans[obj][trans_count[obj]]=trans[obj][CurMb].y+search_range;
// 				s_trans[obj][trans_count[obj]]=(int)((trans[obj][CurMb].scale-MIN_ALPHA)*100)/5;
// 				o_trans[obj][trans_count[obj]]=(int)(trans[obj][CurMb].offset-MIN_BETA)/5;
// 
// 				trans_reference[obj][trans_count[obj]] = trans[obj][CurMb].reference;
// 				trans_count[obj]++;
			    
			}
		}
  		else//边界？
		{
      		if(rms> (tol*tol*no))//没有匹配，继续向下划分
			{
				for(mode=1;mode<3;mode++)
				{
					rectRms=0;
					trans[obj][CurMb].partition=mode;
					for(i=0;i<2;i++)
					{				    
						trans[obj][CurMb].next[i].x=0;
			            trans[obj][CurMb].next[i].y=0;
						if(encode_block_rect(b_x,b_y,i,&trans[obj][CurMb].next[i],con,mode,1)==0)//如果矩形匹配失败
						{
							i=4;
						}
						else
						{
						rectRms++;
						}
					}
					if(rectRms==2)      //说明每个分块都匹配，停止mode循环
					{
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 						partition[obj][partition_length[obj]]=mode;
// 						partition_length[obj]++;
// 						xy_p[obj][trans[obj][CurMb].next[0].x+input->search_range]++;
// 						xy_p[obj][trans[obj][CurMb].next[1].x+input->search_range]++;
// 						xy_p[obj][trans[obj][CurMb].next[0].y+input->search_range]++;
// 						xy_p[obj][trans[obj][CurMb].next[1].y+input->search_range]++;
// 
// 						s_p[obj][(int)((trans[obj][CurMb].next[0].scale-MIN_ALPHA)*100/5)]++;
// 						s_p[obj][(int)((trans[obj][CurMb].next[1].scale-MIN_ALPHA)*100/5)]++;
// 
// 						o_p[obj][(int)(trans[obj][CurMb].next[0].offset-MIN_BETA)/5]++;
// 						o_p[obj][(int)(trans[obj][CurMb].next[1].offset-MIN_BETA)/5]++;
// 
// 						x_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[0].x+search_range;
// 						y_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[0].y+search_range;
// 						s_trans[obj][trans_count[obj]]=(int)((trans[obj][CurMb].next[0].scale-MIN_ALPHA)*100)/5;
// 						o_trans[obj][trans_count[obj]]=(int)(trans[obj][CurMb].next[0].offset-MIN_BETA)/5;
// 						trans_reference[obj][trans_count[obj]] = trans[obj][CurMb].next[0].reference;
// 						trans_count[obj]++;
// 
// 						x_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[1].x+search_range;
// 						y_trans[obj][trans_count[obj]]=trans[obj][CurMb].next[1].y+search_range;
// 						s_trans[obj][trans_count[obj]]=(int)((trans[obj][CurMb].next[1].scale-MIN_ALPHA)*100)/5;
// 						o_trans[obj][trans_count[obj]]=(int)(trans[obj][CurMb].next[1].offset-MIN_BETA)/5;
// 						trans_reference[obj][trans_count[obj]] = trans[obj][CurMb].next[1].reference;
// 						trans_count[obj]++;
// 					
// 						mode=4;
// 						if(region!=1)
// 							num_rect1+=2;
					}
				}
				if(mode<4)
				{
					trans[obj][CurMb].partition=3;
					partition[obj][partition_length[obj]]=trans[obj][CurMb].partition;
					partition_length[obj]++;				
			 	    for(i=0;i<2;i++)
						for(j=0;j<2;j++)
						{
							encode_block_8(b_x+j*8,b_y+i*8,i*2+j,&trans[obj][CurMb].next[i*2+j],con);
						}
				}
			}
            else//16就满足（边界块）
			{
//huffman//////ZZZZZZZZZZZZZZZZZZZ
//  	    		if (region!=1)
// 		    		num_16++;
// 		    	partition[obj][partition_length[obj]]=trans[obj][CurMb].partition;
// 			    partition_length[obj]++;
// 				xy_p[obj][trans[obj][CurMb].x+input->search_range]++;/////
// 				xy_p[obj][trans[obj][CurMb].y+input->search_range]++;/////
// 				s_p[obj][(int)((trans[obj][CurMb].scale-MIN_ALPHA)*100)/5]++;
// 				o_p[obj][(int)(trans[obj][CurMb].offset-MIN_BETA)/5]++;
// 				
// 				x_trans[obj][trans_count[obj]]=trans[obj][CurMb].x+search_range;
// 				y_trans[obj][trans_count[obj]]=trans[obj][CurMb].y+search_range;
// 				s_trans[obj][trans_count[obj]]=(int)((trans[obj][CurMb].scale-MIN_ALPHA)*100)/5;
// 				o_trans[obj][trans_count[obj]]=(int)(trans[obj][CurMb].offset-MIN_BETA)/5;
// 				
// 				trans_reference[obj][trans_count[obj]] = trans[obj][CurMb].reference;
// 				trans_count[obj]++;
			}
		}
	}
}

// void p_encode_one_macroblock(int CurMb,TRANS_NODE **trans,int con)
// {
// // 	int block16_x_n,block16_y_n;
// // 	int block16_x,block16_y;
// // 
// // 	block16_x_n = img->width/16;//22
// // 	block16_y_n = img->height/16;//18
// // 	block16_x   = CurMb%block16_x_n;
// // 	block16_y   = CurMb/block16_x_n;
// 
// 	encode_one_macroblock(CurMb,trans,con);
// // 	decode_one_macroblock(CurMb,trans,con);
// // 
// // 	img->current_mb_nr=block16_y*block16_x_n+block16_x;
// // 	printf("\n\ncurmb=%3d,16x=%2d,16y=%2d",CurMb,block16_x,block16_y);
// // 
// // 	LumaResidualCoding_fract(block16_x,block16_y);
// }

int encode_block_rect(int b_x,int b_y,int Curblock_8,TRANS_NODE *trans,int con,int mode,int depth)
{
	int i,j;
	double tol=input->tol_8;
	double rms;
	byte **plane;
	int block_size_x,block_size_y;
	int region=0;

// 	if (con==1)
// 	{tol=input->tol_8;
// 	} 
//     else
// 	{tol=0;//input->tol_8/10000000000000;
// 	}
	switch(mode)
	{
		case 1://16x8划分（8x4划分）
			block_size_x=16/depth;block_size_y=8/depth;
			b_y=b_y+block_size_y*Curblock_8;break;//16x8（8x4）块的y坐标
		case 2://8x16划分（4x8划分）
			block_size_x=8/depth;block_size_y=16/depth;
			b_x=b_x+block_size_x*Curblock_8;break;//8x16（4x8）块的y坐标
	}

	if(num_regions==2)//no
	{
		switch(con)
		{
		case 1:
			plane=plane_Y;break;
		case 2:
			plane=plane_UV[0];break;
		case 3:
			plane=plane_UV[1];break;
		}
		for(j=b_y;j<b_y+block_size_y;j++)
			for(i=b_x;i<b_x+block_size_x;i++)
			{
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //判断当前块属于哪个region
			}	
		region--;	//0-0region;1--1region;2--边界块
		if(region==2||region==obj)
		{
			trans->use=1;
			if(region==2)
				trans->region=1; //边界块
			else
				trans->region=0;//对象内部块
		}
		else
		{
			trans->use=0;
			return 0;
		}
	}
	
	if (search_mode==0)
	{
		rms = full_search(b_x,b_y,block_size_x,block_size_y,con,trans);//全搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
	}
	if (search_mode==1)
	{
		rms = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,trans);//新六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
	}
	if (search_mode==2)
	{
		rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,trans);  //UMHEX搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
	}
	if (search_mode==3)
	   {
        rms =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,trans);//六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
	   }


	if (currentVideo=='C')//如果是C目，则利用分数像素计算RMS

	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));

		changeReferenceFrame('H');//参考H目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//全搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//新六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);  //UMHEX搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
	   }

		changeReferenceFrame('C');
		if (rms_rl<rms)//更新trans
		{
			num_rl++;
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->reference=1;//参考H目
			trans->block_type=0;
		}
		else
			trans->reference=0;//参考C目
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');//参考M目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//全搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//新六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);   //UMHEX搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}

		changeReferenceFrame('C');
		if (rms_rl<rms)//更新trans
		{
			num_rl++;
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->reference=2;//参考M目
			trans->block_type=0;
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');//参考N目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//全搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//新六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);  //UMHEX搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//六边形搜索 16x8块（或8x16\8x4\4x8块） 计算出trans
		}

		changeReferenceFrame('C');
		if (rms_rl<rms)//更新trans
		{
			num_rl++;
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->reference=3;//参考N目
			trans->block_type=0;
		}
	}

else//如果当前目不是C目

{
	double rms_rl;
	TRANS_NODE trans_rl;
	memset(&trans_rl,0,sizeof(TRANS_NODE));
	changeReferenceFrame('C');
	rms_rl = full_search_R(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);
	changeReferenceFrame(currentVideo);

	if (rms_rl<rms)
	{
		num_rl++;
		trans->reference=0;//参考块是C目
		rms = rms_rl;
		trans->x=trans_rl.x;
		trans->y=trans_rl.y;
		trans->offset=trans_rl.offset;
		trans->scale=trans_rl.scale;
		trans->block_type=0;
	}
	else
	{	trans->block_type=0;
    	trans->reference=4;
	}


	memset(&trans_rl,0,sizeof(TRANS_NODE));
	changeReferenceFrame('H');
	rms_rl = full_search_R(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);
	changeReferenceFrame(currentVideo);
	if (rms_rl<rms)
	{
		num_rl++;
		//	trans->reference=1;//参考块是C目
		rms = rms_rl;
		trans->x=trans_rl.x;
		trans->y=trans_rl.y;
		trans->offset=trans_rl.offset;
		trans->scale=trans_rl.scale;
		trans->reference=1;//参考H目
		trans->block_type=0;
	}
	memset(&trans_rl,0,sizeof(TRANS_NODE));
	changeReferenceFrame('M');
	rms_rl = full_search_R(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);
	changeReferenceFrame(currentVideo);
	if (rms_rl<rms)
	{
		num_rl++;
		rms = rms_rl;
		trans->x=trans_rl.x;
		trans->y=trans_rl.y;
		trans->offset=trans_rl.offset;
		trans->scale=trans_rl.scale;
		trans->reference=2;//参考M目
		trans->block_type=0;
	}
	memset(&trans_rl,0,sizeof(TRANS_NODE));
	changeReferenceFrame('N');
	rms_rl = full_search_R(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);
	changeReferenceFrame(currentVideo);
	if (rms_rl<rms)
	{
		num_rl++;
		//	trans->reference=1;//参考块是C目
		rms = rms_rl;
		trans->x=trans_rl.x;
		trans->y=trans_rl.y;
		trans->offset=trans_rl.offset;
		trans->scale=trans_rl.scale;
		trans->reference=3;//参考N目
		trans->block_type=0;
		}

}


	if(rms> (tol*tol*no))//16x8与8x16块的阈值为tol_8 * tol_8 * 128, 8x4与4x8块的阈值为tol_8 * tol_8 * 32
		return 0;//16x8与8x16块(8x4与4x8)不能够匹配
	else
	{
		return 1;//16x8与8x16块(8x4与4x8)足够匹配
	}
}


void encode_block_8(int b_x,int b_y,int Curblock_8,TRANS_NODE *trans,int con)//编码8*8划分的块
{
	int i,j;
	double tol=input->tol_8;
	double rms;
	byte **plane;
	int mode,rectRms;

	int region=0;

// 	if (con==1)
// 	{
// 		tol=input->tol_8;
// 	} 
//     else
// 	{
// 		tol=0;//input->tol_8/10000000000000;
// 	}


	if(num_regions==2)//no
	{
		switch(con)
		{
		case 1:
			plane=plane_Y;break;
		case 2:
			plane=plane_UV[0];break;
		case 3:
			plane=plane_UV[1];break;
		}
		for(j=b_y;j<b_y+8;j++)
			for(i=b_x;i<b_x+8;i++)
			{
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //判断当前块属于哪个region
			}	
		region--;	//0-0region;1--1region;2--边界块
		if(region==2||region==obj)
		{
			trans->use=1;
			if(region==2)
				trans->region=1;
			else
				trans->region=0;
		}
		else
		{
			trans->use=0;
			return;
		}
	}

	trans->partition=0;
	trans->reference = 0;
	trans->x=0;
	trans->y=0;

	if (search_mode==0)
	{
		rms = full_search(b_x,b_y,8,8,con,trans);//全搜索 8x8块 计算trans
	}
	if (search_mode==1)
	{
		rms =new_hexagon_block_search(b_x,b_y,8,8,con,trans);//新六边形搜索 8x8块 计算trans
	}
	if (search_mode==2)
	{
		rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,trans);  //UMHEX搜索 8x8块 计算trans
	}
	if (search_mode==3)
	{
        rms =hexagon_block_search(b_x,b_y,8,8,con,trans);//六边形搜索 8x8块 计算trans
	}

	if (currentVideo=='C')//如果是C目，则利用分数像素计算rms
// 	if (currentVideo!='C')////撤销掉分数像素ZZZZZZZZZZZZZZZZZZZ
	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');//参考H目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,8,8,con,&trans_rl);//全搜索 8x8块 计算trans
		}
		if (search_mode==1)
		{
			rms_rl =new_hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//新六边形搜索 8x8块 计算trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,&trans_rl);  //UMHEX搜索 8x8块 计算trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//六边形搜索 8x8块 计算trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)//更新trans
		{
			num_rl++;
			trans->reference=1;//参考块是H目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
			trans->reference=0;//参考块是C目
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');//参考M目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,8,8,con,&trans_rl);//全搜索 8x8块 计算trans
		}
		if (search_mode==1)
		{
			rms_rl =new_hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//新六边形搜索 8x8块 计算trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,&trans_rl);  //UMHEX搜索 8x8块 计算trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,8,8,con,&trans_rl); //六边形搜索 8x8块 计算trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)//更新trans
		{
			num_rl++;
			trans->reference=2;//参考块是M目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');//参考N目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,8,8,con,&trans_rl);//全搜索 8x8块 计算trans
		}
		if (search_mode==1)
		{
			rms_rl =new_hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//新六边形搜索 8x8块 计算trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,&trans_rl);  //UMHEX搜索 8x8块 计算trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//六边形搜索 8x8块 计算trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)//更新trans
		{
			num_rl++;
			trans->reference=3;//参考块是N目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
	}
    else//当前目非C目
	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('C');
		rms_rl = full_search_R(b_x,b_y,8,8,con,&trans_rl);
		changeReferenceFrame(currentVideo);
   
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=0;//参考块是C目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
		{	trans->block_type=0;
		    trans->reference=4;//参考当前目
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');
		rms_rl = full_search_R(b_x,b_y,8,8,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=1;//参考块是H目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');
		rms_rl = full_search_R(b_x,b_y,8,8,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=2;//参考块是M目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');
		rms_rl = full_search_R(b_x,b_y,8,8,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=3;//参考块是N目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}

	}
	
	if(rms> (tol*tol*no))//若8*8不能匹配，继续向下划分为8*4，4*8，4*4（匹配的阈值为tol_8*tol_8*64）
	{
		int skip_old=skip;
		for(mode=1;mode<3;mode++)//8*4，4*8划分
		{
			rectRms=0;
			trans->partition=mode;
			for(i=0;i<2;i++)//8*4或4*8划分都将1个8*8划分成了两个小块，所以，这里要循环两次
			{		
				trans->next[i].x=0;
				trans->next[i].y=0;
				if(encode_block_rect(b_x,b_y,i,&trans->next[i],con,mode,2)==0)//8*4（或4*8）划分不能匹配
				{
					i=4;//如果两个小块中的一个不能匹配就跳出for(i=0;i<2;i++)循环
				}
				else
					rectRms++;
			}
			if(rectRms==2)//两个8*4（或4*8）小块都能匹配
			{
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 				partition[obj][partition_length[obj]]=mode;
// 		        partition_length[obj]++;
// 
// 				xy_p[obj][trans->next[0].x+input->search_range]++;
// 				xy_p[obj][trans->next[1].x+input->search_range]++;
// 				xy_p[obj][trans->next[0].y+input->search_range]++;
// 				xy_p[obj][trans->next[1].y+input->search_range]++;
// 
// 				s_p[obj][(int)((trans->next[0].scale-MIN_ALPHA)*100/5)]++;
// 				s_p[obj][(int)((trans->next[1].scale-MIN_ALPHA)*100/5)]++;
// 				
// 				o_p[obj][(int)(trans->next[0].offset-MIN_BETA)/5]++;
// 				o_p[obj][(int)(trans->next[1].offset-MIN_BETA)/5]++;
// 
// 				x_trans[obj][trans_count[obj]]=trans->next[0].x+search_range;
// 				y_trans[obj][trans_count[obj]]=trans->next[0].y+search_range;
// 				s_trans[obj][trans_count[obj]]=(int)((trans->next[0].scale-MIN_ALPHA)*100)/5;
// 				o_trans[obj][trans_count[obj]]=(int)(trans->next[0].offset-MIN_BETA)/5;
// 				trans_reference[obj][trans_count[obj]] = trans->next[0].reference; 
// 				trans_count[obj]++;
// 
// 				x_trans[obj][trans_count[obj]]=trans->next[1].x+search_range;
// 				y_trans[obj][trans_count[obj]]=trans->next[1].y+search_range;
// 				s_trans[obj][trans_count[obj]]=(int)((trans->next[1].scale-MIN_ALPHA)*100)/5;
// 				o_trans[obj][trans_count[obj]]=(int)(trans->next[1].offset-MIN_BETA)/5;
// 				trans_reference[obj][trans_count[obj]] = trans->next[1].reference;
// 				trans_count[obj]++;
				
				mode=4;//两个小块（8*4或4*8）都能匹配，mode=4
// 				if(region!=1)
// 				num_rect1+=2;
// 				num_rect2+=2;
			}
		}
		if(mode<4)//否则4*4划分
		{
			trans->partition=3;
			partition[obj][partition_length[obj]]=3;
		    partition_length[obj]++;			
			for(i=0;i<2;i++)//循环4个4*4小块
				for(j=0;j<2;j++)
				{
					encode_block_4(b_x+j*4,b_y+i*4,&trans->next[i*2+j],con);//4*4划分
				}
		}
	}
	else//若8*8就能匹配
	{
		if (region!=1)
		num_8++;
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 		partition[obj][partition_length[obj]]=trans->partition;
// 		partition_length[obj]++;
// 		{
// 			xy_p[obj][trans->x+input->search_range]++;
// 			xy_p[obj][trans->y+input->search_range]++;
// 			s_p[obj][(int)((trans->scale-MIN_ALPHA)*100/5)]++;
// 			o_p[obj][(int)(trans->offset-MIN_BETA)/5]++;
// 
// 			x_trans[obj][trans_count[obj]]=trans->x+search_range;
// 			y_trans[obj][trans_count[obj]]=trans->y+search_range;
// 			s_trans[obj][trans_count[obj]]=(int)((trans->scale-MIN_ALPHA)*100)/5;
// 			o_trans[obj][trans_count[obj]]=(int)(trans->offset-MIN_BETA)/5;
// 
// 			trans_reference[obj][trans_count[obj]] = trans->reference;
// 			trans_count[obj]++;
// 		}
	}
}


void encode_block_4(int b_x,int b_y,TRANS_NODE *trans,int con)
{
	//当前trans是4x4块的trans
	int i,j;
	byte **plane;
	double tol=input->tol_4;
	double rms;
	int region=0;

// 	if (con==1)
// 	{tol=input->tol_4;
// 	} 
//     else
// 	{tol=input->tol_4/10000000000000;
// 	}
	if(num_regions==2)//no
	{
		switch(con)
		{
		case 1:
			plane=plane_Y;break;
		case 2:
			plane=plane_UV[0];break;
		case 3:
			plane=plane_UV[1];break;
		}
		for(j=b_y;j<b_y+4;j++)
			for(i=b_x;i<b_x+4;i++)
			{
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //判断当前块属于哪个region
			}	
		region--;	//0-0region;1--1region;2--边界块
		if(region==2||region==obj)
		{
			trans->use=1;
			if(region==2)
				trans->region=1;
			else
				trans->region=0;
		}
		else
		{
			trans->use=0;
			return;
		}
	}
	if (region!=1)//region=0
	num_4++;
	trans->x=trans->y=0;

	if(con==3&&b_x==0&&b_y==0)
		printf("con==3&&b_x==0&&b_y==0");
	if (search_mode==0)
	{
		rms = full_search(b_x,b_y,4,4,con,trans);//全搜索 4x4块 计算trans
	}
	if (search_mode==1)
	{
		rms =new_hexagon_block_search(b_x,b_y,4,4,con,trans);//新六边形搜索 4x4块 计算trans
	}
	if (search_mode==2)
	{
		rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,trans);  //UMHEX搜索 4x4块 计算trans
	}
	if (search_mode==3)
	   {
        rms =hexagon_block_search(b_x,b_y,4,4,con,trans);//六边形搜索 4x4块 计算trans
	   }

	if (currentVideo=='C')//如果是C目，则利用分数像素 并计算rms
// 	if (currentVideo!='C')////撤销掉分数像素ZZZZZZZZZZZZZZZZZZZ
	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');//参考H目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,4,4,con,&trans_rl);//全搜索 4x4块 计算trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//新六边形搜索 4x4块 计算trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,&trans_rl);   //UMHEX搜索 4x4块 计算trans
		}
		if (search_mode==3)
		{
			rms_rl = hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//六边形搜索 4x4块 计算trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)
		{
			num_rl++;
			trans->partition=trans->reference=1;//参考块是H目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
			trans->reference=0;//参考块是C目
		
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');//参考M目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,4,4,con,&trans_rl);//全搜索 4x4块 计算trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//新六边形搜索 4x4块 计算trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,&trans_rl);  //UMHEX搜索 4x4块 计算trans
		}
		if (search_mode==3)
		{
			rms_rl = hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//六边形搜索 4x4块 计算trans
		}
		changeReferenceFrame('C');
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=2;//参考块是M目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');//参考N目
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,4,4,con,&trans_rl);//全搜索 4x4块 计算trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//新六边形搜索 4x4块 计算trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,&trans_rl);   //UMHEX搜索 4x4块 计算trans
		}
		if (search_mode==3)
		{
			rms_rl = hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//六边形搜索 4x4块 计算trans
		}
		changeReferenceFrame('C');
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=3;//参考块是N目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
	}
    else//不是C目 就不用分数像素？？？？应该是哈
	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('C');
		rms_rl = full_search_R(b_x,b_y,4,4,con,&trans_rl);
		changeReferenceFrame(currentVideo);
	
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=0;//参考块是C目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
		{	trans->block_type=0;
	    	trans->reference=4;//参考块是当前目
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');
		rms_rl = full_search_R(b_x,b_y,4,4,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=1;//参考块是H目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');
		rms_rl = full_search_R(b_x,b_y,4,4,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=2;//参考块是M目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');
		rms_rl = full_search_R(b_x,b_y,4,4,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=3;//参考块是N目
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}

	}
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 	xy_p[obj][trans->x+input->search_range]++;
// 	xy_p[obj][trans->y+input->search_range]++;
// 	s_p[obj][(int)((trans->scale-MIN_ALPHA)*100/5)]++;
// 	o_p[obj][(int)(trans->offset-MIN_BETA)/5]++;
// 	x_trans[obj][trans_count[obj]]=trans->x+search_range;
// 	y_trans[obj][trans_count[obj]]=trans->y+search_range;
// 	s_trans[obj][trans_count[obj]]=(int)((trans->scale-MIN_ALPHA)*100)/5;
// 	o_trans[obj][trans_count[obj]]=(int)(trans->offset-MIN_BETA)/5;
// 	trans_reference[obj][trans_count[obj]] = trans->reference;
// 	trans_count[obj]++;
}


double full_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
	double best_rms=1e30;
    double rms;
	int l,k,i,j,m,n;
	double alpha,beta;
    int window_size=input->search_range;//输入的搜索范围
	
	best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);//计算mv=（0,0）位置的domain块和range块的RMS匹配误差，con=1:Y分量;con=2：U分量;con=3：V分量
	trans->scale=alpha;//将分形参数alpha传入trans
	trans->offset=beta;//将分形参数beta传入trans
	for(l=1;l<=window_size;l++)//在搜索窗口内,window_size等于配置文件输入的搜索范围
	{
		i=j=-l;//第一个搜索点MV是(-1,-1)
		for(k=0;k<8*l;k++)//在每一个螺旋圈内搜索点的个数
		{
			m=block_x+i;//domain块的x坐标
			n=block_y+j;//domain块的y坐标

			if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))//检查domain块是否超出图像边界，(m,n)是domain块的坐标，（block_x,block_y）是range块的坐标
			{
				rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);//计算domain块（m,n）和range块(block_x,block_y)的RMS匹配误差，分形参数alpha，beta
				if(rms<best_rms)//更新RMS误差及分形参数
				{
					best_rms=rms;
					trans->x=m-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
			}
			// spt++;
			
			if(k<2*l) //下面是更改搜索点坐标
				i++;
			else if(k<4*l)
				j++;
			else if(k<6*l)
				i--;
			else 
				j--;
		}
	}
	return best_rms;
}

// double full_search_R(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
// {
// 	double best_rms=1e30;
//     double rms;
// 	int l,k,i,j,m,n;
// 	double alpha,beta;
//     int window_size=50;
// 	
// 	best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
// 	trans->scale=alpha;
// 	trans->offset=beta;
// 	for(l=1;l<=window_size;l++)//在搜索窗口内
// 	{
// 		i=j=-l;
// 		for(k=0;k<8*l;k++)//在每一个螺旋圈内搜索点的个数
// 		{
// 			m=block_x+i;
// 			n=block_y+j;
// 			if((-2<=(n-block_y))&&((n-block_y)<=2)&&(0<=(m-block_x)))
// 			{
// 			
// 			if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 			{   
// 				S_P++;
// 				rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 				if(rms<best_rms)
// 				{
// 					best_rms=rms;
// 					trans->x=m-block_x;
// 					trans->y=n-block_y;
// 					trans->scale=alpha;
// 					trans->offset=beta;
// 				}
// 			}
// 			}
// 			// spt++;
// 			
// 			if(k<2*l) 
// 				i++;
// 			else if(k<4*l)
// 				j++;
// 			else if(k<6*l)
// 				i--;
// 			else 
// 				j--;
// 		}
// 	}
// 	return best_rms;
// }



double full_search_R(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
	double best_rms=1e30;
    double rms;//,r1,r2
	int i,m,n;//,ii,jj,jl,k,
	double alpha,beta;
	
	best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
	trans->scale=alpha;
	trans->offset=beta;
	
	for(i=1;i<3;i++)
	{
		m=disparity.mv_x+i*3;
		n=disparity.mv_y;
		if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{   
			S_P++;
			rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
			}
			/* else
			{
			S_P++;
			r1=compute_rms(block_x,block_y,m,n+2,&alpha,&beta,block_size_x,block_size_y,con);
			r2=compute_rms(block_x,block_y,m,n-2,&alpha,&beta,block_size_x,block_size_y,con);
			if(r1>r2)
			{
			for (ii=m-3;ii<=m;ii++)
			for (jj=n;jj<=n+2;jj++)
			{
			S_P++;
			rms=compute_rms(block_x,block_y,ii,jj,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
			//disparity.mv_x=ii;
			// disparity.mv_y=jj;
			best_rms=rms;
			trans->x=ii-block_x;
			trans->y=jj-block_y;
			trans->scale=alpha;
			trans->offset=beta;
			}
			} 
			}
			else
			{
			for (ii=m-3;ii<=m;ii++)
			for (jj=n-2;jj<=n;jj++)
			{
			S_P++;
			rms=compute_rms(block_x,block_y,ii,jj,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
			//	disparity.mv_x=ii;
			//  disparity.mv_y=jj;
			best_rms=rms;
			trans->x=ii-block_x;
			trans->y=jj-block_y;
			trans->scale=alpha;
			trans->offset=beta;
			}						  
			} 
			}
			break;
		}*/
		}
	}
	
    disparity.mv_x=trans->x;
	disparity.mv_y=trans->y;
	return best_rms;
}


/*double full_search_R(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
	double best_rms=1e30;
    double rms;
	int l,k,i,j,m,n;
	double alpha,beta;
    int window_size=50;
	
	best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
	trans->scale=alpha;
	trans->offset=beta;
	for(i=1;i<20;i++)
	{
      m=block_x+i*3;
	  n=block_y;
      if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
	  {   
		  S_P++;
		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
		  if(rms<best_rms)
		  {
			  best_rms=rms;
			  trans->x=m-block_x;
			  trans->y=n-block_y;
			  trans->scale=alpha;
			  trans->offset=beta;
		  }
		  else
		  {
			  S_P++;
			  rms=compute_rms(block_x,block_y,m-1,n,&alpha,&beta,block_size_x,block_size_y,con);
			  if(rms<best_rms)
			  {
				  best_rms=rms;
				  trans->x=m-1-block_x;
				  trans->y=n-block_y;
				  trans->scale=alpha;
				  trans->offset=beta;
			  }
              S_P++;
			  rms=compute_rms(block_x,block_y,m-2,n,&alpha,&beta,block_size_x,block_size_y,con);
			  if(rms<best_rms)
			  {
				  best_rms=rms;
				  trans->x=m-2-block_x;
				  trans->y=n-block_y;
				  trans->scale=alpha;
				  trans->offset=beta;
			  }
			  break;
		  }
	  }
	}
	for(i=1;i<20;i++)
	{
		m=block_x+i*3;
		n=block_y+1;
		if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{   
			S_P++;
			rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
			}
			else
			{
				S_P++;
				rms=compute_rms(block_x,block_y,m-1,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-1-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				S_P++;
				rms=compute_rms(block_x,block_y,m-2,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-2-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				break;
			}
		}
	}
	for(i=1;i<20;i++)
	{
		m=block_x+i*3;
		n=block_y+2;
		if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{   
			S_P++;
			rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
			}
			else
			{
				S_P++;
				rms=compute_rms(block_x,block_y,m-1,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-1-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				S_P++;
				rms=compute_rms(block_x,block_y,m-2,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-2-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				break;
			}
		}
	}
	for(i=1;i<20;i++)
	{
		m=block_x+i*3;
		n=block_y-1;
		if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{   
			S_P++;
			rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
			}
			else
			{
				S_P++;
				rms=compute_rms(block_x,block_y,m-1,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-1-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				S_P++;
				rms=compute_rms(block_x,block_y,m-2,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-2-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				break;
			}
		}
	}
	for(i=1;i<20;i++)
	{
		m=block_x+i*3;
		n=block_y-2;
		if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{   
			S_P++;
			rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
			}
			else
			{
				S_P++;
				rms=compute_rms(block_x,block_y,m-1,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-1-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				S_P++;
				rms=compute_rms(block_x,block_y,m-2,n,&alpha,&beta,block_size_x,block_size_y,con);
				if(rms<best_rms)
				{
					best_rms=rms;
					trans->x=m-2-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
				break;
			}
		}
	}
	//}
// 	for(l=1;l<=window_size;l++)//在搜索窗口内
// 	{
// 		i=j=-l;
// 		for(k=0;k<8*l;k++)//在每一个螺旋圈内搜索点的个数
// 		{
// 			m=block_x+i;
// 			n=block_y+j;
// 			if((-2<=(n-block_y))&&((n-block_y)<=2)&&(0<=(m-block_x)))
// 			{
// 				
// 				if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 				{   
// 					S_P++;
// 					rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 					if(rms<best_rms)
// 					{
// 						best_rms=rms;
// 						trans->x=m-block_x;
// 						trans->y=n-block_y;
// 						trans->scale=alpha;
// 						trans->offset=beta;
// 					}
// 				}
// 			}
// 			// spt++;
// 			
// 			if(k<2*l) 
// 				i++;
// 			else if(k<4*l)
// 				j++;
// 			else if(k<6*l)
// 				i--;
// 			else 
// 				j--;
// 		}
// 	}
	return best_rms;
}*/


double new_hexagon_block_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
   int i,j,k,m,n,c=9,t,b,ox=block_size_x,oy=block_size_y,step;//l,,err0
   MV mv,tmv;

   double best_rms=1e30;
   double rms=0.0;
   //int l,k,i,j,m,n;
   double alpha,beta;
  

   static MV crx[8]={{-2,0},{-1,0},{0,-2},{0,-1},{2,0},{1,0},{0,2},{0,1}};//8 point of cross
   static MV hex[6]={{-2,0},{-1,-2},{1,-2},{2,0},{1,2},{-1,2}};//六边形的六个搜索点
   static MV lcrx[4]={{-2,0},{0,-2},{2,0},{0,2}};
   static MV scrx[4]={{-1,0},{0,-1},{1,0},{0,1}};
   static MV dia[8]={{-2,0},{-1,-1},{0,-2},{1,-1},{2,0},{1,1},{0,2},{-1,1}};
   WINDOW_SIZE=input->search_range;
//    init_searchedpt();
//    
//   
//    
//    best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
//    trans->scale=alpha;
// 	trans->offset=beta;
//    //*err=mean_abs_err(cur,prev,x,y,x,y,BLOCK_WIDTH,BLOCK_HEIGHT);
//    searchedpt[WINDOW_SIZE][WINDOW_SIZE]=1;

  // init_searchedpt();
   best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
   trans->scale=alpha;
	trans->offset=beta;
   //*err=mean_abs_err(cur,prev,x,y,x,y,BLOCK_WIDTH,BLOCK_HEIGHT);
  // searchedpt[WINDOW_SIZE][WINDOW_SIZE]=1;

   tmv.x=0;
   tmv.y=0;
   k=5;
   step=1;   
   for(i=0;i<4;i++)
   {  
      m=block_x+scrx[i].x;
      n=block_y+scrx[i].y;
      	if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{
		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
		  
		 //searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;

	//	 S_P+=1;

         if(rms<best_rms)
		 {
			 best_rms=rms;
			 trans->x=m-block_x;
			 trans->y=n-block_y;
			 trans->scale=alpha;
			 trans->offset=beta;
			 k=i;
			 tmv=scrx[i];
			}
      }
   }

  mv=tmv;

   if(k==5) 
   {
     return best_rms;  //1_nd stop  
   }

   else
   {
	   c=k; k=5;
	   tmv.x=0; tmv.y=0; 
	   if(c==0)     j=3;
	   else if(c==1) j=0;
	   else if(c==2) j=1;
	   else if(c==3) j=2;
	   for(i=0;i<3;i++)
	   {
		   m=block_x+mv.x+scrx[j].x;
           n=block_y+mv.y+scrx[j].y;

		   if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		   {
			   rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			   //err0=novel_mean_abs_err(cur,prev,x,y,m,n,BLOCK_WIDTH,BLOCK_HEIGHT);
			   //searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;

			   S_P+=1;
			   
			   if(rms<best_rms)
			   {
				   best_rms=rms;
				   trans->x=m-block_x;
				   trans->y=n-block_y;
				   trans->scale=alpha;
				   trans->offset=beta;
				   k=j;
				   tmv=scrx[j];
			}
			   
			   
			   j=(j+1)<4?(j+1):(j+1-4);
		   }
	   }
   }
   
    mv.x=mv.x+tmv.x;
    mv.y=mv.y+tmv.y;
 
   if(k==5) 
   {
      return best_rms;  //2nd-step-stop
   }


   k=5;
   if(c==0)     j=1;
   else if(c==1) j=2;
   else if(c==2) j=3;
   else if(c==3) j=0;
   tmv.x=0; tmv.y=0;
   
   for(i=0;i<3;i++)
   {
	   m=block_x+lcrx[j].x;
	   n=block_y+lcrx[j].y;
	   
	   if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
	   {
          rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
           //err0=novel_mean_abs_err(cur,prev,x,y,m,n,BLOCK_WIDTH,BLOCK_HEIGHT);
		  // searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;

		   S_P+=1;
		   if(rms<best_rms)
		   {
			   best_rms=rms;
			   trans->x=m-block_x;
			   trans->y=n-block_y;
			   trans->scale=alpha;
			   trans->offset=beta;
			   tmv=lcrx[j];
			   k=j;
			}
          
		   j=(j+1)<4?(j+1):(j+1-4);

	   }
   }
  

   if(k!=5) 
   {
      mv = tmv ;
   }
  
   k = 0;

   
   //search 6 points once
   for(i=0;i<6;i++)
   {
         m=block_x+mv.x+hex[i].x;
         n=block_y+mv.y+hex[i].y;
         if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		 {
            rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
            //searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;
//#ifdef OPC
          //  op.c++;
//#endif

			S_P+=1;
            if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
				tmv=hex[i];
               k=i;
			}
           
         }
   }
   mv.x=mv.x+tmv.x;
   mv.y=mv.y+tmv.y;

   if(k==9)
	   return best_rms;

   while(1)
   {
      if(k==9)
         break;
      else
	  {
         j=k;k=9;
         tmv.x=0;tmv.y=0;
      }

      t=3;b=-1;
      j=(j+b)<0?(6+j+b):(j+b);/////////////////////
      for(i=0;i<t;i++)
	  {
         m=block_x+mv.x+ hex[j].x;
         n=block_y+mv.y+ hex[j].y;
         if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		 {
           rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			 //err0=novel_mean_abs_err(cur,prev,x,y,m,n,BLOCK_WIDTH,BLOCK_HEIGHT);
            //searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;
//#ifdef OPC
          //  op.c++;
//#endif

			S_P+=1;
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
			 tmv=dia[j];
               k=j;
			}
            
         }
         j=(j+1)>5?0:(j+1);//*******************************************************
      }
      mv.x=mv.x+tmv.x;
      mv.y=mv.y+tmv.y;
   }

   tmv.x=0;tmv.y=0;
   for(i=0;i<8;i+=2)
   {
      m=block_x+mv.x+(dia[i].x>>1);
      n=block_y+mv.y+(dia[i].y>>1);
      if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
	  {
        rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
       // searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;
//#ifdef OPC
       // op.c++;
//#endif
		S_P+=1;
		if(rms<best_rms)
		{
			best_rms=rms;
			trans->x=m-block_x;
			trans->y=n-block_y;
			trans->scale=alpha;
			trans->offset=beta;
			tmv.x=(dia[i].x>>1);
            tmv.y=(dia[i].y>>1);
			}
       
      }
   }
   mv.x=mv.x+tmv.x;
   mv.y=mv.y+tmv.y;
//   spt += count_searchedpt();
   return best_rms;
} 


//基于中心偏置的菱形搜索新方法////////////////////////////////////////////////////////////
//double horizon_bias_diamond_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
double hexagon_block_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
	int i,j,k,m,n,c=9;//b,err0,t,,l
	MV mv,tmv;	
	double best_rms=1e30;
	double rms=0.0;
	double alpha,beta;
	//static MV bias[4]={{-2,0},{0,-1},{2,0},{0,1}};//偏置菱形的4个点/////////////////
	static MV fbx[5]={{-2,0},{-1,-1},{1,-1},{2,0},{0,1}};//五边形的5个点/////////////////
	//static MV scr[4]={{-1,0},{0,-1},{1,0},{0,1}};//小十字形的4个点////////////
	static MV bcr[4]={{-1,0},{0,-2},{1,0},{0,2}};//垂直偏置菱形的4个点////////////
	
	WINDOW_SIZE=input->search_range;
	best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
	trans->scale=alpha;
	trans->offset=beta;
	tmv.x=0;
	tmv.y=0;
	k=5;

    for(i=0;i<5;i++)
   {  
		m=motion.mv_x+fbx[i].x;
		n=motion.mv_y+fbx[i].y;
    //  m=block_x+fbx[i].x;              /////////////////////////////////////////
     // n=block_y+fbx[i].y;                 //////// 五边形的5个点
      	if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		{
		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);		  
		 //searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;
	//	 S_P+=1;
         if(rms<best_rms)
		 {
			 //motion.mv_x=m;
			 //motion.mv_y=n;
			 best_rms=rms;
			 trans->x=m-block_x;
			 trans->y=n-block_y;
			 trans->scale=alpha;
			 trans->offset=beta;
			 k=i;
			 tmv=fbx[i];
			}
      }
   }

  mv=tmv;
  
 while(1)
  {
      if(k==5)
		  break;
      else
	  {
		  c=k;k=5;
		  tmv.x=0;tmv.y=0;
      }
	  
	  if(c==0)     j=4;
	  else if(c==1) j=0;
	  else if(c==2) j=2;
	   else if(c==3) j=3;
	   else if(c==4) j=4;
	   for(i=0;i<3;i++)    ///////////////////////没用，我删掉的，呵呵/////////////////
		                                 //错了，有用，用于其他三个点判断
	   {
		   // m=block_x+mv.x+fbx[j].x;
          // n=block_y+mv.y+fbx[j].y;     //4 小菱形///////////////////////////////
		  m=motion.mv_x+mv.x+fbx[j].x;
		   n=motion.mv_y+mv.y+fbx[j].y;
		   if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		   {
			   rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			   S_P+=1;			   
			   if(rms<best_rms)
			   {
				   //motion.mv_x=m;
			     //  motion.mv_y=n;
			       best_rms=rms;
				   trans->x=m-block_x;
				   trans->y=n-block_y;
				   trans->scale=alpha;
				   trans->offset=beta;
				   k=j;
				   tmv=fbx[j];
			   }			   			   
			   j=(j+1)<5?(j+1):(j+1-5);
		   }
	   }
   }

//垂直偏置菱形的4个点////////////////////////////////////////
 
  mv.x=mv.x+tmv.x;
  mv.y=mv.y+tmv.y;
   
   tmv.x=0; tmv.y=0;
   
   for(j=0;j<4;j++)
   {
	   // m=block_x+mv.x+bcr[j].x;
	   //n=block_y+mv.y+bcr[j].y;     //小十字搜索的左右两个点//////////////////////////
	  m=motion.mv_x+mv.x+bcr[j].x;
	   n=motion.mv_y+mv.y+bcr[j].y;     //垂直偏置菱形的4个点//////////////////////////////////////	
	   
	   if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
	   {
		   rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
           //err0=novel_mean_abs_err(cur,prev,x,y,m,n,BLOCK_WIDTH,BLOCK_HEIGHT);
		   // searchedpt[WINDOW_SIZE+m-ox][WINDOW_SIZE+n-oy]=1;
		   S_P+=1;
		   if(rms<best_rms)
		   {
			   //motion.mv_x=m;
		      // motion.mv_y=n;
			   best_rms=rms;
			   trans->x=m-block_x;
			   trans->y=n-block_y;
			   trans->scale=alpha;
			   trans->offset=beta;
		   }        
	   }
   }
   motion.mv_x=trans->x;
   motion.mv_y=trans->y;
   return best_rms;
}


/*double hexagon_block_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
   int i,j,k,b,t,m,n,err0,c=9,ox=block_x,oy=block_y,step;
   MV mv,tmv;
   double best_rms=1e30;
   double rms=0.0;

   double alpha,beta;
   static MV dia[8]={{-2,0},{-1,-1},{0,-2},{1,-1},{2,0},{1,1},{0,2},{-1,1}};//应用到4 points of big_diamond

   static MV hex[6]={{-2,0},{-1,-2},{1,-2},{2,0},{1,2},{-1,2}};//六边形的六个搜索点


   WINDOW_SIZE=input->search_range;
   best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
   trans->scale=alpha;
   trans->offset=beta;
   tmv.x=0;
   tmv.y=0;

   k=7;
   step=1;
   for(i=0;i<6;i++)
   {
   
	   m=block_x+hex[i].x;
	   n=block_y+hex[i].y;
     if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
	  {
		 
		 rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
		 
        
         if(rms<best_rms)
		 {
			 best_rms=rms;
			 trans->x=m-block_x;
			 trans->y=n-block_y;
			 trans->scale=alpha;
			 trans->offset=beta;
			 k=i;
			 tmv=hex[i];
		}
      }
   }
   mv=tmv;
// 
// 
   while(1)
   {
      if(k==7)
         break;
      else
	  {
         j=k;k=7;
         tmv.x=0;tmv.y=0;
      }
      step+=1;

	  t=3;b=-1;
      j=(j+b)<0?(6+j+b):(j+b);//
      for(i=0;i<t;i++)
	  {
         m=block_x+mv.x+hex[j].x;
         n=block_y+mv.y+hex[j].y;
         if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
		 {

			rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
			
		
			
			if(rms<best_rms)
			{
				best_rms=rms;
				trans->x=m-block_x;
				trans->y=n-block_y;
				trans->scale=alpha;
				trans->offset=beta;
				k=j;
				tmv=hex[j];
			}

         }
         j=(j+1)>5?0:(j+1);//*******************************************************
      }
      mv.x=mv.x+tmv.x;
      mv.y=mv.y+tmv.y;
   }

   tmv.x=0;tmv.y=0;
   for(i=0;i<8;i+=2)
   {
      m=block_x+mv.x+(dia[i].x>>1);
      n=block_y+mv.y+(dia[i].y>>1);
      if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
	  {
		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
		  
	
		  
		  if(rms<best_rms)
		  {
			  best_rms=rms;
			  trans->x=mv.x+(dia[i].x>>1);
			  trans->y=mv.y+(dia[i].y>>1);
			  trans->scale=alpha;
			  trans->offset=beta;
			 /* k=j;*/

		/*	}

      }
    }
      return best_rms;
   }*/


int bound_chk(int m, int n, int cx, int cy,int block_size_x,int block_size_y,int con)
//当前块搜索坐标[m,n]不超出当前帧图像范围，才做下一步运算
{
	int ilow,ihigh,jlow,jhigh;
	int b=0;
	int window_size=input->search_range;//配置文件输入的搜索范围
    int width=input->imagewidth/min(2,con);//图像宽度
	int height=input->imageheight/min(2,con);//图像高度

	ilow=cx-window_size;//搜索范围最左边界坐标
	ihigh=cx+window_size;//搜索范围最右边界坐标
	if(ilow<0)//如果搜索范围最左边界超出图像了
		ilow=0;
	if(ihigh>width-block_size_x)//如果搜索范围最右边界超出图像了
		ihigh=width-block_size_x;
	
	jlow=cy-window_size;//搜索范围最上边界坐标
	jhigh=cy+window_size;//搜索范围最下边界坐标
	if(jlow<0)//如果搜索范围最上边界超出图像了
		jlow=0;
	if(jhigh>height-block_size_y)//如果搜索范围最下边界超出图像了
		jhigh=height-block_size_y;
 
	b=m<=ihigh&&m>=ilow&&n<=jhigh&&n>=jlow;//若搜索范围超过图像边界，b=0
	return(b);
}




// 新的运动估计，基于田隽师姐的最新EPZS(Enhanced Predictive Zonal Search)算法/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double UMHEXIntegerPelBlockMotionSearch  (int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
       int   iYMinNow, iXMinNow,m,n,i,k,k1;
	   int   best_x = 0, best_y = 0;
       double alpha,beta,tmv_x=0,tmv_y=0,t1mv_x=0,t1mv_y=0,lx_x[2],lx_y[2];
       double best_rms=1e30;
       double rms;
	   static const int dia_x[4] = {-1, 0, 1,0};
	   static const int dia_y[4] = {0, -1, 0,1};
	   static const int zf_x[4] = {-1, 1, 1,-1};
	   static const int zf_y[4] = {1, 1, -1,-1};

//上层模式预测/////////////////////////	   
	   if (UpCur.mv_x!=0||UpCur.mv_y!=0)
	   {
	   if(bound_chk(UpCur.mv_x+block_x,UpCur.mv_y+block_y,block_x,block_y,block_size_x,block_size_y,con))
	   {
	 	  rms=compute_rms(block_x,block_y,UpCur.mv_x+block_x,UpCur.mv_y+block_y,&alpha,&beta,block_size_x,block_size_y,con);
	 	  if(rms<best_rms)
	 	  {
	 		  best_rms=rms;
	 		  trans->x=UpCur.mv_x;
	 		  trans->y=UpCur.mv_y;
	 		  trans->scale=alpha;
	 		  trans->offset=beta;
	 		  best_x=UpCur.mv_x+block_x;
	 		  best_y=UpCur.mv_y+block_y;
	 		  iXMinNow = UpCur.mv_x+block_x;
	           iYMinNow = UpCur.mv_y+block_y;
	 	  }
	 	}
	   }
      //如果当前块上层模式的运动矢量在整像素位置，则定义这样的块为平坦块。
	   //对于平坦块，可以跳过分数像素运动矢量搜索。
    //  if(best_x%16==0||best_x%8==0||best_x%4==0)
	  {	  
		    	  //相邻块预测
		          iXMinNow = best_x;
			      iYMinNow = best_y;	
				  m = iXMinNow-1;
				  n = iYMinNow ;
				  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
				  {
					  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
					  if(rms<best_rms)
					  {
						  best_rms=rms;
						  trans->x=m-block_x;
						  trans->y=n-block_y;
						  trans->scale=alpha;
						  trans->offset=beta;
						  best_x=m;
						  best_y=n;
					  }
				  }

				  //原始搜索中心(0,0)预测
				  iXMinNow = best_x;
				  iYMinNow = best_y;	
				  m = block_x;
				  n = block_y;
				  if(bound_chk(m,n,block_x,block_x,block_size_x,block_size_y,con))
				  {
					  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
					  if(rms<best_rms)
					  {
						  best_rms=rms;
						  trans->x=m-block_x;
						  trans->y=n-block_y;
						  trans->scale=alpha;
						  trans->offset=beta;
						  best_x=m;
						  best_y=n;
					  }
				  }
				  //比较搜索起始点和其周围的4个菱形搜索点
				  iXMinNow = best_x;
				  iYMinNow = best_y;	
				  for(i=0;i<4;i++)
				  {			  
				     m=block_x+dia_x[i];
					  n=block_y+dia_y[i];
					  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
					  {					  
						  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);					 		  
						  if(rms<best_rms)
						  {
							  best_rms=rms;
							  trans->x=m-block_x;
							  trans->y=n-block_y;
							  trans->scale=alpha;
							  trans->offset=beta;
							  k=i;
							  best_x=m;
						      best_y=n;
							  tmv_x=dia_x[i];
							  tmv_y=dia_y[i];
						  }
					  }
				  }
			
				  //如果最小绝对误差之和MSAD位于搜索起始点，则停止
				  if (best_x == iXMinNow && best_y == iYMinNow)
				 {
					       return best_rms;
				  }

				  //最佳匹配点和次最佳匹配点相对，则选择最佳匹配点MV
				  for(i=0;i<4;i++)
				  {			  
					if(i!=k)
					{
					  m=block_x+dia_x[i];
					  n=block_y+dia_y[i];
					  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
					  {					  
						  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);					 		  
						  if(rms<best_rms)
						  {
							  best_rms=rms;
							  trans->x=m-block_x;
							  trans->y=n-block_y;
							  trans->scale=alpha;
							  trans->offset=beta;
							  best_x=m;
						      best_y=n;
							  k1=i;
							  t1mv_x=dia_x[i];
							  t1mv_y=dia_y[i];
						  }
					  }
					}
				  }

				  if(tmv_x==t1mv_x||tmv_y==t1mv_y)
					  return best_rms;
				  
			  //计算与其相邻的正方形模板上点的匹配误差
// 				  if (tmv_x==1||tmv_y==1||t1mv_x==1||t1mv_y==1)
// 				  {
// 					  zf_x[3]={1,1,2};
// 					  zf_y[3]={1,1,2};
// 				  }
// 
// 				  else
// 				  {
// 					  zf_x[3]={-1,-1,-2};
// 					  zf_y[3]={-1,-1,-2};
// 				  }
// 				  for(i=0;i<3;i++)
// 				  {			  
// 					  m=block_x+zf_x[i];
// 					  n=block_y+zf_y[i];
// 					  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 					  {					  
// 						  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);					 		  
// 						  if(rms<best_rms)
// 						  {
// 							  best_rms=rms;
// 							  trans->x=m-block_x;
// 							  trans->y=n-block_y;
// 							  trans->scale=alpha;
// 							  trans->offset=beta;
// 							  k=i;
// 							  tmv_x=zf_x[i];
// 							  tmv_y=zf_y[i];
// 						  }
// 					  }
// 				  }
// 				  if (abs(tmv_x)==1&&abs(tmv_y)==1)
// 				  {
// 					  return best_rms;
// 				  }

					  for(i=0;i<4;i++)
					  {			  
						  m=block_x+zf_x[i];
						  n=block_y+zf_y[i];
						  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
						  {					  
							  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);					 		  
							  if(rms<best_rms)
							  {
								  best_rms=rms;
								  trans->x=m-block_x;
								  trans->y=n-block_y;
								  trans->scale=alpha;
								  trans->offset=beta;
								  k=i;
								  tmv_x=zf_x[i];
								  tmv_y=zf_y[i];
							  }
						  }
					  }

					  if (abs(tmv_x)!=1||abs(tmv_y)!=1)
					  {
						 return best_rms;
					  }

					  //以第四步中正方形模板上的搜索点为中心，用菱形模板搜索其周围的
					  if (tmv_x==1)
					  {
						  lx_x[0]=1;
						  lx_x[1]=2;
					  }
					  
					  else
					  {
						  lx_x[0]=-1;
						  lx_x[1]=-2;
					  }

					  if (tmv_y==1)
					  {
						  lx_y[0]=1;
						  lx_y[1]=2;	  
					  }
					  
					  else
					  {
						  lx_y[0]=-1;
						  lx_y[1]=-2;
					  }

					  for(i=0;i<2;i++)
					  {			  
						  m=block_x+lx_x[i];
						  n=block_y+lx_y[i];
						  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
						  {					  
							  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);					 		  
							  if(rms<best_rms)
							  {
								  best_rms=rms;
								  trans->x=m-block_x;
								  trans->y=n-block_y;
								  trans->scale=alpha;
								  trans->offset=beta;
							  }
						  }
					  }
			
   }  
 return best_rms;	   
}





//以前的非对称六边形算法/////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                   //  ==> minimum motion cost after search
// double UMHEXIntegerPelBlockMotionSearch  (int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
// { 
// 
//      int   best_x = 0, best_y = 0;
//      int   search_step, iYMinNow, iXMinNow;
//      int   pos, cand_x, cand_y,  mcost;
// //   int   i,m,j;
// //   float betaFourth_1,betaFourth_2;
//      int  temp_Big_Hexagon_x[16];//  temp for Big_Hexagon_x;
//      int  temp_Big_Hexagon_y[16];//  temp for Big_Hexagon_y;
// 
//     double best_rms=1e30;
//     double rms;
//     int l,k,i,j,m,n;
//     double alpha,beta;
//     int window_size=input->search_range;
//   
// 	static const int Diamond_x[4] = {-1, 0, 1, 0};
// 	static const int Diamond_y[4] = {0, 1, 0, -1};
// 	static const int Hexagon_x[6] = {2, 1, -1, -2, -1, 1};
// 	static const int Hexagon_y[6] = {0, -2, -2, 0,  2, 2};
// 	static const int Big_Hexagon_x[16] = {0,-2, -4,-4,-4, -4, -4, -2,  0,  2,  4,  4, 4, 4, 4, 2};
// 	static const int Big_Hexagon_y[16] = {4, 3, 2,  1, 0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2, 3};
// ///////////////////五边形/////////////////////////////////////////////
// 	static const int wbx_x[25] = {-2,-1, 0,1,2, -2,-1, 0,1,2,-2,-1, 0,1,2,-2,-1, 0,1,2,-2,-1, 0,1,2};
// 	static const int wbx_y[25] = {2,2,2,2,2,1,1,1,1,1,0,0,0,0,0,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2};
// 
//  /*sec_step:*/ //Unsymmetrical-cross search
//   //对应点预测
//   iXMinNow = block_x;
//   iYMinNow = block_y;
//   best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
//   trans->scale=alpha;
//   trans->offset=beta;
//   best_x=block_x;
//   best_y=block_y;
//   //上层模式预测
//   
//   if (UpCur.mv_x!=0||UpCur.mv_y!=0)
//   {
//   if(bound_chk(UpCur.mv_x+block_x,UpCur.mv_y+block_y,block_x,block_y,block_size_x,block_size_y,con))
//   {
// 	  rms=compute_rms(block_x,block_y,UpCur.mv_x+block_x,UpCur.mv_y+block_y,&alpha,&beta,block_size_x,block_size_y,con);
// 	  if(rms<best_rms)
// 	  {
// 		  best_rms=rms;
// 		  trans->x=UpCur.mv_x;
// 		  trans->y=UpCur.mv_y;
// 		  trans->scale=alpha;
// 		  trans->offset=beta;
// 		  best_x=UpCur.mv_x+block_x;
// 		  best_y=UpCur.mv_y+block_y;
// 		  iXMinNow = UpCur.mv_x+block_x;
//           iYMinNow = UpCur.mv_y+block_y;
// 	  }
// 	}
//   }
//   //L预测模式
//   if(bound_chk(Lrefpre.mv_x+block_x,Lrefpre.mv_y+block_y,block_x,block_y,block_size_x,block_size_y,con))
//   {
// 	  rms=compute_rms(block_x,block_y,Lrefpre.mv_x+block_x,Lrefpre.mv_y+block_y,&alpha,&beta,block_size_x,block_size_y,con);
// 	  if(rms<best_rms)
// 	  {
// 		  best_rms=rms;
// 		  trans->x=Lrefpre.mv_x;
// 		  trans->y=UpCur.mv_y;
// 		  trans->scale=alpha;
// 		  trans->offset=beta;
// 		  best_x=Lrefpre.mv_x+block_x;
// 		  best_y=Lrefpre.mv_y+block_y;
// 		  iXMinNow = Lrefpre.mv_x+block_x;
//           iYMinNow = Lrefpre.mv_y+block_y;
// 	  }
//   }
// 
// 
//   //Unsymmetrical-cross search
//   for(i = 1; i < search_range; i+=2)
//   {
//     search_step = i;
//     m = iXMinNow + search_step;
//     n= iYMinNow ;
// 	if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	{
// 		rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		if(rms<best_rms)
// 		{
// 			best_rms=rms;
// 			trans->x=m-block_x;
// 			trans->y=n-block_y;
// 			trans->scale=alpha;
// 			trans->offset=beta;
// 			best_x=m;
// 			best_y=n;
// 		}
// 	}
//     m = iXMinNow - search_step;
//     n = iYMinNow ;
// 	if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	{
// 		rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		if(rms<best_rms)
// 		{
// 			best_rms=rms;
// 			trans->x=m-block_x;
// 			trans->y=n-block_y;
// 			trans->scale=alpha;
// 			trans->offset=beta;
// 			best_x=m;
// 			best_y=n;
// 		}
// 	}
//   }
//   for(i = 1; i < (search_range/2);i+=2)
//   {
//     search_step = i;
//     m = iXMinNow ;
//     n = iYMinNow + search_step;
// 	if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	{
// 		rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		if(rms<best_rms)
// 		{
// 			best_rms=rms;
// 			trans->x=m-block_x;
// 			trans->y=n-block_y;
// 			trans->scale=alpha;
// 			trans->offset=beta;
// 			best_x=m;
// 			best_y=n;
// 		}
// 	}
//     m = iXMinNow ;
//     n = iYMinNow - search_step;
// 	if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	{
// 		rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		if(rms<best_rms)
// 		{
// 			best_rms=rms;
// 			trans->x=m-block_x;
// 			trans->y=n-block_y;
// 			trans->scale=alpha;
// 			trans->offset=beta;
// 			best_x=m;
// 			best_y=n;
// 		}
// 	}
//   }
// 
// 
//   //early termination alogrithm, refer to JVT-G016
//  // EARLY_TERMINATION
// 
//   iXMinNow = best_x;
//   iYMinNow = best_y;
// 
//   //third_step:    // Uneven Multi-Hexagon-grid Search
//   //sub step 1: 5x5 squre search
//   for(pos=1;pos<25;pos++)
//   {
//    // m = iXMinNow + spiral_search_x[pos];
//   //  n = iYMinNow + spiral_search_y[pos];
// m = iXMinNow + wbx_x[pos];
//     n = iYMinNow + wbx_y[pos];
// 
//     if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	{
// 		rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		if(rms<best_rms)
// 		{
// 			best_rms=rms;
// 			trans->x=m-block_x;
// 			trans->y=n-block_y;
// 			trans->scale=alpha;
// 			trans->offset=beta;
// 			best_x=m;
// 			best_y=n;
// 		}
// 	}
//   }
//    
//   //early termination alogrithm, refer to JVT-G016
//   //EARLY_TERMINATION
// 
//   //sub step 2:  Multi-Hexagon-grid search
//   memcpy(temp_Big_Hexagon_x,Big_Hexagon_x,64);
//   memcpy(temp_Big_Hexagon_y,Big_Hexagon_y,64);
//   for(i=1;i<=(search_range/4); i++)
//   {
// 
//     for (l = 0; l < 16; l++)
//     {
//       m = iXMinNow + temp_Big_Hexagon_x[l];
//       n = iYMinNow + temp_Big_Hexagon_y[l];
//       temp_Big_Hexagon_x[l] += Big_Hexagon_x[l];
//       temp_Big_Hexagon_y[l] += Big_Hexagon_y[l];
// 
//       if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	  {
// 		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		  if(rms<best_rms)
// 		  {
// 			  best_rms=rms;
// 			  trans->x=m-block_x;
// 			  trans->y=n-block_y;
// 			  trans->scale=alpha;
// 			  trans->offset=beta;
// 			  best_x=m;
// 			  best_y=n;
// 		  }
// 	  }
//     }
//     // ET_Thd2: early termination Threshold for strong motion
//    // if(min_mcost < ET_Thred)
//     //{
//     //  goto terminate_step;
//     //}
//   }
// 
// 
// 
//   //fourth_step:  //Extended Hexagon-based Search
//   // the fourth step with a small search pattern
// //fourth_1_step:  //sub step 1: small Hexagon search
//   for(i = 0; i < search_range; i++)
//   {
//     iXMinNow = best_x;
//     iYMinNow = best_y;
//     for (l = 0; l < 6; l++)
//     {
//       m = iXMinNow + Hexagon_x[l];
//       n = iYMinNow + Hexagon_y[l];
// 	  if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	  {
// 		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		  if(rms<best_rms)
// 		  {
// 			  best_rms=rms;
// 			  trans->x=m-block_x;
// 			  trans->y=n-block_y;
// 			  trans->scale=alpha;
// 			  trans->offset=beta;
// 			  best_x=m;
// 			  best_y=n;
// 		  }
// 	  }
//     }
// 
//     if (best_x == iXMinNow && best_y == iYMinNow)
//     {
//       break;
//     }
//   }
// //fourth_2_step: //sub step 2: small Diamond search
// 
// 
//   for(i = 0; i < search_range; i++)
//   {
//     iXMinNow = best_x;
//     iYMinNow = best_y;
//     for (l = 0; l < 4; l++)
//     {
//       m = iXMinNow + Diamond_x[l];
//       n = iYMinNow + Diamond_y[l];
//       if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))
// 	  {
// 		  rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);
// 		  if(rms<best_rms)
// 		  {
// 			  best_rms=rms;
// 			  trans->x=m-block_x;
// 			  trans->y=n-block_y;
// 			  trans->scale=alpha;
// 			  trans->offset=beta;
// 			  best_x=m;
// 			  best_y=n;
// 		  }
// 	  }
//     }
//     if(best_x == iXMinNow && best_y == iYMinNow)
//       break;
//   }
// // 
// // terminate_step:
//   Lrefpre.mv_x=trans->x;
//   Lrefpre.mv_y=trans->y;
//   return best_rms;
// }