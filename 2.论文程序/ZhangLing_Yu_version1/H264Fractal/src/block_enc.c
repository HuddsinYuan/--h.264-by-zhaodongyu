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
			spiral_search_x[k] =  i;//������
			spiral_search_y[k] = -l;
			spiral_hpel_search_x[k] =  i<<1;//������
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
		plane_Y_domain_temp = plane_Y_domain;//plane_domainӦ�ú�img_refͬ�����£���decode_oneframe��
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

void SetRefAndMotionVectors_fract(int mode,int CurrMb,TRANS_NODE *trans,int con)//mode=0�ǳ�ʼ��,mode=1: 16x16�洢��mode=2: 16x18�洢��mode=3: 8x16�洢��mode=4: 8x8�洢
{
	int block16_x_n,block16_y_n,block16_x,block16_y;
	int i,j,k,m,n;
	int Y,U,V;
	int yuv;
	Macroblock *currMB = &img->mb_data[img->current_mb_nr];
	
	block16_x_n=img->width/16;//ͼ��һ�а������ٸ����
	block16_y_n=img->height/16;//ͼ��һ�а������ٸ����
	
	block16_x=CurrMb%block16_x_n;//��ǰ�����ͼ���е�x���꣨�Ժ��Ϊ��λ��
	block16_y=CurrMb/block16_x_n;//��ǰ�����ͼ���е�y���꣨�Ժ��Ϊ��λ��
	Y=U=V=0;
	yuv=4;
	if (2==con)//U����
	{
// 		con=1;
		U=img->width/4;//88
		yuv=2;
	} 
	if (3==con)//V����
	{
// 		con=1;
		U=img->width/4;//88
		V=img->height/8;//36
		yuv=2;
	}//else
// 	{
// 		con=0;
// 	}
	if(0==mode)//��ʼ��16*16��飬��4*4Ϊ������λ�����Ҫѭ��16��
	{
// 		currMB->mb_type = mode;

		for(j=0;j<yuv;j++)
			for(i=0;i<yuv;i++)//U,V�����洢��Y�����ĺ��
			{
				
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=0;//��MV��ʼ��Ϊ��0,0��
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=0;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=0;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=0;//��list0��list1�ο�������ʼ��Ϊ0
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=0;
			}
// 		printf("\nmode=%d",mode);
	}
	if(1==mode)//16*16����4*4Ϊ������λ�����Ҫѭ��16�Σ�16�ζ�һ��
	{
// 		printf("\nmode=%d",mode);
		currMB->mb_type = mode;
		for(i=0;i<yuv;i++)
		{
			currMB->b8mode[i] = mode;//ʲô��˼��Ϊʲôѭ���ĴΣ�
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
		for(j=0;j<yuv;j++)//���16*8��
			for(i=0;i<2;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[0].x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[0].y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].offset/5;
			}
// 		printf("\ntrans->next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[0].scale,trans->next[0].offset);


		for(j=0;j<yuv;j++)//�ұ�16*8��
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
		for(j=0;j<2;j++)//�ϱ�8*16��
			for(i=0;i<yuv;i++)
			{
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[0].x;
				enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[0].y;
// 				enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].reference;
				enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].scale*20;
				enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[0].offset/5;
			}
// 		printf("\ntrans->next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[0].scale,trans->next[0].offset);

		for(j=2;j<yuv;j++)//�±�8*16��
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
	
	
	if(4==mode)//p8*8=8��Ҫ��ÿ��8*8��Ļ���ģʽ
	{
// 		printf("\nmode=%d",mode);
		currMB->mb_type = mode+4;

		for(m=0;m<yuv/2;m++)//��
			for(n=0;n<yuv/2;n++)//��
			{
				k=2*m+n;//k��0~3���ֱ��ʾ4��8*8��
// 				printf("\nk=%d",k);
				if(0==trans->next[k].partition)//8*8���ٷ���
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k]:scale=.0lf,offset=.0lf\n",trans->next[k].scale,trans->next[k].offset);
					
					currMB->b8mode[k]   = 4;//8*8���ٷ���

					for(j=2*m;j<2*m+2;j++)//Ӧ������4*4Ϊ������λ�洢�İɣ�����ѭ����4�Σ�����
						for(i=2*n;i<2*n+2;i++)
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].offset/5;
						}

				}//if
				
				if(1==trans->next[k].partition)//8*8�����ֳ�����8*4
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k].next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[0].scale,trans->next[k].next[0].offset);
// 					printf("\ntrans->next[k].next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[1].scale,trans->next[k].next[1].offset);

					currMB->b8mode[k]   = 5;//8*8�����ֳ�����8*4

// 						currMB->b8pdir[i]   = b8pdir[i];

					for(j=2*m;j<2*m+2;j++)
						for(i=2*n;i<2*n+1;i++)//�洢��һ��8*4��
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[0].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[0].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].offset/5;
						}
					for(j=2*m;j<2*m+2;j++)
						for(i=2*n+1;i<2*n+2;i++)//�洢�ڶ���8*4��
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[1].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[1].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].offset/5;
						}
				}//if
				
				
				if(2==trans->next[k].partition)//8*8�����ֳ�����4*8
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k].next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[0].scale,trans->next[k].next[0].offset);
// 					printf("\ntrans->next[k].next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[1].scale,trans->next[k].next[1].offset);

					
					currMB->b8mode[k]   = 6;//8*8�����ֳ�����4*8

					for(j=2*m;j<2*m+1;j++)
						for(i=2*n;i<2*n+2;i++)//�洢��һ��4*8��
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[0].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[0].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].offset/5;
						}
					for(j=2*m+1;j<2*m+2;j++)
						for(i=2*n;i<2*n+2;i++)//�洢�ڶ���4*8��
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[1].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[1].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].offset/5;
						}
				}//if
				
				if(3==trans->next[k].partition)//8*8�����ֳ��ĸ�4*4
				{
// 					printf("\ntrans->next[k].partition=%d",trans->next[k].partition);
// 					printf("\ntrans->next[k].next[0]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[0].scale,trans->next[k].next[0].offset);
// 					printf("\ntrans->next[k].next[1]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[1].scale,trans->next[k].next[1].offset);
// 					printf("\ntrans->next[k].next[2]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[2].scale,trans->next[k].next[2].offset);
// 					printf("\ntrans->next[k].next[3]:scale=%.0lf,offset=%.0lf\n",trans->next[k].next[3].scale,trans->next[k].next[3].offset);
					
					currMB->b8mode[k]   = 7;//8*8�����ֳ��ĸ�4*4

					
					for(j=2*m;j<2*m+1;j++)
						for(i=2*n;i<2*n+1;i++)//�洢��һ��4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[0].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[0].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[0].offset/5;
						}
					for(j=2*m;j<2*m+1;j++)
						for(i=2*n+1;i<2*n+2;i++)//�洢�ڶ���4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[1].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[1].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[1].offset/5;
						}
					for(j=2*m+1;j<2*m+2;j++)
						for(i=2*n;i<2*n+1;i++)//�洢������4*4
						{
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][0]=trans->next[k].next[2].x;
							enc_picture->mv[0][U+block16_x*yuv+i][V+block16_y*yuv+j][1]=trans->next[k].next[2].y;
// 							enc_picture->ref_idx[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[2].reference;
							enc_picture->ref_pic_id[0][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[2].scale*20;
							enc_picture->ref_pic_id[1][U+block16_x*yuv+i][V+block16_y*yuv+j]=trans->next[k].next[2].offset/5;
						}
					for(j=2*m+1;j<2*m+2;j++)
						for(i=2*n+1;i<2*n+2;i++)//�洢���ĸ�4*4
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
   brief:��һ�������б���
   param:CurMb ��ǰ�����
   param:trans ���洦���Ľ��
   param:con ��ǰ����
*/
void encode_one_macroblock(int CurMb,TRANS_NODE **trans,int con)//con=1��Y������con=2:U������con=3:V����
{
    int i,j,k,ii=0,q=2;        
	int mode;//�黮��ģʽ
	byte **plane;
	double tol=input->tol_16;//�����ļ������tol_16
	int mb_X = CurMb%(img->frmWidthInMbs/min(2,con));//�����Ĺ�դɨ�����꣬ת���ɶ�άx����(�Ժ��Ϊ��λ)
	int mb_Y = CurMb/(img->frmWidthInMbs/min(2,con));//�����Ĺ�դɨ�����꣬ת���ɶ�άy����(�Ժ��Ϊ��λ)
	int b_x=mb_X*16;//������Ϊ��λ�ĺ�����Ͻ�x����
	int b_y=mb_Y*16;//������Ϊ��λ�ĺ�����Ͻ�y����
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


	if(num_regions==2)//no //if OB ���жϵ�ǰ���������һ����
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
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //�жϵ�ǰ�������ĸ�region
			}	
		region--;	//0-region0;1--region1;2--�߽��
		if(region==2)
		{
			trans[0][CurMb].use=trans[1][CurMb].use=1;
			trans[0][CurMb].region=trans[1][CurMb].region=1;
		}
		else
		{
			trans[region][CurMb].use=1;
			trans[1-region][CurMb].use=0;
			trans[0][CurMb].region=trans[1][CurMb].region=0;//�ǲ��ǻ��ڶ����
			obj=region;
		}
	}

	for(k=region;k<pow(2.0,region);k++) //region=0
	{
		if(region>1)
			obj=k-2;
		trans[obj][CurMb].partition=0;//��ʼ��trans
		trans[obj][CurMb].block_type = 0;
		trans[obj][CurMb].x=0;
        trans[obj][CurMb].y=0;	
	    trans[obj][CurMb].reference=0;
        SetRefAndMotionVectors_fract(0,CurMb,&trans[obj][CurMb],con);//ע���һ������Ϊ0,��ʾ��ʼ���ο������ͷ���MVΪ0

		if (search_mode==0)//ȫ���� 16x16��
		{	
			rms =full_search(b_x,b_y,16,16,con,&trans[obj][CurMb]);  //����ǰһ֡��Ϊ�ο������trans
		}
		if (search_mode==1)//������������ 16x16��
		{
			rms =new_hexagon_block_search(b_x,b_y,16,16,con,&trans[obj][CurMb]);  //����ǰһ֡��Ϊ�ο������trans
			
		}
		if (search_mode==2)//UMHEX���� 16x16��
		{
			rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans[obj][CurMb]);  //����ǰһ֡��Ϊ�ο������trans
			UpCur.mv_x=trans[obj][CurMb].x;
			UpCur.mv_y=trans[obj][CurMb].y;
		}
		if (search_mode==3)//���������� 16x16��
		{
			rms =hexagon_block_search(b_x,b_y,16,16,con,&trans[obj][CurMb]);
		}



		/////////////////////////////////////////////////////////////////////////////////////////////////////
		//Ԥ���������������������жϣ�ֻ�������������Ž��к��������///////////////////////////////////////////////////////	 		
		//Ԥ���������������������жϣ�ֻ�������������Ž��к��������///////////////////////////////////////////////////////	 		

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
			//��ʼ�ж��Ƿ��������ķ��δ���///////////////////////////////////////////////////////////////
			//if((fabs(1-MR*MD))<14.4)
			//  if(fabs(MR*MD)>10)    
			//	if(chun<0.9999999)
			//��������в�ͬ�ĳ��ԣ��ó���ѵĽ����Ҳ����chun>=0�ĵط�
    	if (currentVideo=='C')////�����CĿ�������÷������ؼ���RMS
// 		if (currentVideo!='C')////��������������ZZZZZZZZZZZZZZZZZZZ
		{	
			double rms_rl;
			TRANS_NODE trans_rl;
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('H');//Ϊʲô�ο�HĿ�ˣ�����
			if (search_mode==0)//ȫ���� 16x16��
			{	
				rms_rl =full_search(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο��������ǰrange���RMS�����β���s,o,x,y
			}
			if (search_mode==1)//������������ 16x16��
			{
				rms_rl =new_hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
				
			}
			if (search_mode==2)//UMHEX���� 16x16��
			{
				rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
				UpCur.mv_x=trans[obj][CurMb].x;
				UpCur.mv_y=trans[obj][CurMb].y;
			}
			if (search_mode==3)//���������� 16x16��
			{
				rms_rl =hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);
			}

			changeReferenceFrame('C');
			if (rms_rl<rms)//����RMS���ͷ��β���trans
			{
				num_rl++;
				trans[obj][CurMb].reference=1;//�ο�����HĿ
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
}
			else
				trans[obj][CurMb].reference=0;//�ο�����CĿ
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('M');//Ϊʲô�ο�MĿ�ˣ�����
			if (search_mode==0)//ȫ���� 16x16��
			{	
				rms_rl =full_search(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
			}
			if (search_mode==1)//������������ 16x16��
			{
				rms_rl =new_hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
				
			}
			if (search_mode==2)//UMHEX���� 16x16��
			{
				rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
				UpCur.mv_x=trans[obj][CurMb].x;
				UpCur.mv_y=trans[obj][CurMb].y;
			}
			if (search_mode==3)//���������� 16x16��
			{
				rms_rl =hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);
			}

			changeReferenceFrame('C');
			if (rms_rl<rms)//����trans
			{
				num_rl++;
				trans[obj][CurMb].reference=2;//�ο�����MĿ
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('N');//Ϊʲô�ο�NĿ�ˣ�??
			if (search_mode==0)//ȫ���� 16x16��
			{	
				rms_rl =full_search(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
			}
			if (search_mode==1)//������������ 16x16��
			{
				rms_rl =new_hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
				
			}
			if (search_mode==2)//UMHEX���� 16x16��
			{
				rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,16,16,con,&trans_rl);  //����ǰһ֡��Ϊ�ο������trans
				UpCur.mv_x=trans[obj][CurMb].x;
				UpCur.mv_y=trans[obj][CurMb].y;
			}
			if (search_mode==3)//���������� 16x16��
			{
				rms_rl =hexagon_block_search(b_x,b_y,16,16,con,&trans_rl);
			}
			changeReferenceFrame('C');
			if (rms_rl<rms)//����trans
			{
				num_rl++;
				trans[obj][CurMb].reference=3;//�ο�����NĿ
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}
		}
		
		if (currentVideo!='C')////�������CĿ�������ÿ�����������RMS
			//  else                   //�������CĿ,����п�������
		{
			double rms_rl;
			TRANS_NODE trans_rl;
			memset(&trans_rl,0,sizeof(TRANS_NODE));
			changeReferenceFrame('C');
			//rms_rl = full_search(b_x,b_y,16,16,con,&trans_rl);
			rms_rl = full_search_R(b_x,b_y,16,16,con,&trans_rl);//�ο�CĿ
			changeReferenceFrame(currentVideo);
			if (rms_rl<rms)
			{
				num_rl++;
				trans[obj][CurMb].reference=0;//�ο�����CĿ
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
				trans[obj][CurMb].reference=1;//�ο�����HĿ
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
				trans[obj][CurMb].reference=2;//�ο�����MĿ
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
				trans[obj][CurMb].reference=3;//�ο�����NĿ
				rms = rms_rl;
				trans[obj][CurMb].x=trans_rl.x;
				trans[obj][CurMb].y=trans_rl.y;
				trans[obj][CurMb].offset=trans_rl.offset;
				trans[obj][CurMb].scale=trans_rl.scale;
				trans[obj][CurMb].block_type=0;
			}

		}


  		if(k==0||k==1)//������߱���//yes
		{
			for(j=b_x;j<b_x+16;j++)
				for(i=b_y;i<b_y+16;i++)     //b_x,b_y�ֱ��ʾ��ǰ16*16 range���Ӧ������
				{ 
					switch(con)
					{
					case 1://Y����
						R[ii]=imgY_org[i][j];
						D[ii]=imgY_ref_temp[i][j];
						break;
					case 2://U����
						R[ii]=imgUV_org[0][i][j];
						D[ii]=imgUV_ref_temp[0][i][j];
						break;
					case 3://V����
						R[ii]=imgUV_org[1][i][j];
						D[ii]=imgUV_ref_temp[1][i][j];
						//	D[ii]=imgY_ref[i][j];
						break;
					}
					sumR+=R[ii];//����R������غ�
					sumD+=D[ii];//����D������غ�
					ii++;
				}
			r=sumR/256;  //����R������ƽ��ֵ
			d=sumD/256; //����D������ƽ��ֵ
			for(ii=0;ii<256;ii++)
			{
				sR+=(R[ii]-r)*(R[ii]-r);//����R�� �����ؼ�ȥ��ֵ�Ĳ� ��ƽ��
				sD+=(D[ii]-d)*(D[ii]-d);//����D�� �����ؼ�ȥ��ֵ�Ĳ� ��ƽ��
			}
			for(ii=0;ii<256;ii++)
			{
				mr+=((R[ii]-r)/(sqrt(sR)))*((D[ii]-d)/(sqrt(sD)));
			}
			chun=mr*mr;//����ļ��㶼��Ϊ��chun


       		if(chun<=1&&chun>=0.9&&rms> (tol*tol*no))//���16x16û��ƥ�䣬�������»���,�����õ���chun��RMS�����о�
			{
 				for(mode=1;mode<3;mode++)//16x8���ֺ�8x16����
				{
					rectRms=0;
					trans[obj][CurMb].partition=mode;
					for(i=0;i<2;i++)//ÿ�ֻ��ֶ��Ὣһ����黮��Ϊ����С�飬i=0�ǵ�һ��С�飬i=1�ǵڶ���С��
					{				    
						trans[obj][CurMb].next[i].x=0;//next��ָʲô�أ�����
		        	    trans[obj][CurMb].next[i].y=0;
						if(encode_block_rect(b_x,b_y,i,&trans[obj][CurMb].next[i],con,mode,1)==0)//�������ƥ�䣨16x8��8x16��ʧ��
						{
							i=4;
						}
						else
						{
							rectRms++;
						}
					}
					if(rectRms==2)      //˵��ÿ���ֿ鶼ƥ�䣬ֹͣmodeѭ����ע�����ÿ���ֿ鶼ƥ�����
					{
						partition[obj][partition_length[obj]]=mode;
						partition_length[obj]++;

						if (1==mode)//��¼16*8
						{
						    SetRefAndMotionVectors_fract(2,CurMb,&trans[obj][CurMb],con);//16*8�洢
						}

						if (2==mode)//��¼8*16
						{
							SetRefAndMotionVectors_fract(3,CurMb,&trans[obj][CurMb],con);//8*16�洢
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
				if(mode<4)//16*16��16*8��8*16��ʽ����
				{
					trans[obj][CurMb].partition=3;//8*8����
					partition[obj][partition_length[obj]]=trans[obj][CurMb].partition;
		    	    partition_length[obj]++;				
					for(i=0;i<2;i++)
						for(j=0;j<2;j++)
						{
							encode_block_8(b_x+j*8,b_y+i*8,i*2+j,&trans[obj][CurMb].next[i*2+j],con);//ѭ������4��8*8��
						}

					if (3==mode)//��¼8*8
					{
						SetRefAndMotionVectors_fract(4,CurMb,&trans[obj][CurMb],con);//8*8�洢(��8*4��4*8��4*4)
					}

				
				
				}
			}
    		else    //16x16�������Ҫ����
			{
 				if (region!=1)
 					num_16++;
				partition[obj][partition_length[obj]]=trans[obj][CurMb].partition;
				partition_length[obj]++;


				SetRefAndMotionVectors_fract(1,CurMb,&trans[obj][CurMb],con);//16*16�洢

//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 				xy_p[obj][trans[obj][CurMb].x+input->search_range]++;//ͨ��+7��ʹ�ÿ϶���������
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
  		else//�߽磿
		{
      		if(rms> (tol*tol*no))//û��ƥ�䣬�������»���
			{
				for(mode=1;mode<3;mode++)
				{
					rectRms=0;
					trans[obj][CurMb].partition=mode;
					for(i=0;i<2;i++)
					{				    
						trans[obj][CurMb].next[i].x=0;
			            trans[obj][CurMb].next[i].y=0;
						if(encode_block_rect(b_x,b_y,i,&trans[obj][CurMb].next[i],con,mode,1)==0)//�������ƥ��ʧ��
						{
							i=4;
						}
						else
						{
						rectRms++;
						}
					}
					if(rectRms==2)      //˵��ÿ���ֿ鶼ƥ�䣬ֹͣmodeѭ��
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
            else//16�����㣨�߽�飩
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
		case 1://16x8���֣�8x4���֣�
			block_size_x=16/depth;block_size_y=8/depth;
			b_y=b_y+block_size_y*Curblock_8;break;//16x8��8x4�����y����
		case 2://8x16���֣�4x8���֣�
			block_size_x=8/depth;block_size_y=16/depth;
			b_x=b_x+block_size_x*Curblock_8;break;//8x16��4x8�����y����
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
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //�жϵ�ǰ�������ĸ�region
			}	
		region--;	//0-0region;1--1region;2--�߽��
		if(region==2||region==obj)
		{
			trans->use=1;
			if(region==2)
				trans->region=1; //�߽��
			else
				trans->region=0;//�����ڲ���
		}
		else
		{
			trans->use=0;
			return 0;
		}
	}
	
	if (search_mode==0)
	{
		rms = full_search(b_x,b_y,block_size_x,block_size_y,con,trans);//ȫ���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
	}
	if (search_mode==1)
	{
		rms = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,trans);//������������ 16x8�飨��8x16\8x4\4x8�飩 �����trans
	}
	if (search_mode==2)
	{
		rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,trans);  //UMHEX���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
	}
	if (search_mode==3)
	   {
        rms =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,trans);//���������� 16x8�飨��8x16\8x4\4x8�飩 �����trans
	   }


	if (currentVideo=='C')//�����CĿ�������÷������ؼ���RMS

	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));

		changeReferenceFrame('H');//�ο�HĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//ȫ���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//������������ 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);  //UMHEX���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//���������� 16x8�飨��8x16\8x4\4x8�飩 �����trans
	   }

		changeReferenceFrame('C');
		if (rms_rl<rms)//����trans
		{
			num_rl++;
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->reference=1;//�ο�HĿ
			trans->block_type=0;
		}
		else
			trans->reference=0;//�ο�CĿ
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');//�ο�MĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//ȫ���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//������������ 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);   //UMHEX���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//���������� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}

		changeReferenceFrame('C');
		if (rms_rl<rms)//����trans
		{
			num_rl++;
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->reference=2;//�ο�MĿ
			trans->block_type=0;
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');//�ο�NĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//ȫ���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//������������ 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);  //UMHEX���� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);//���������� 16x8�飨��8x16\8x4\4x8�飩 �����trans
		}

		changeReferenceFrame('C');
		if (rms_rl<rms)//����trans
		{
			num_rl++;
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->reference=3;//�ο�NĿ
			trans->block_type=0;
		}
	}

else//�����ǰĿ����CĿ

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
		trans->reference=0;//�ο�����CĿ
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
		//	trans->reference=1;//�ο�����CĿ
		rms = rms_rl;
		trans->x=trans_rl.x;
		trans->y=trans_rl.y;
		trans->offset=trans_rl.offset;
		trans->scale=trans_rl.scale;
		trans->reference=1;//�ο�HĿ
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
		trans->reference=2;//�ο�MĿ
		trans->block_type=0;
	}
	memset(&trans_rl,0,sizeof(TRANS_NODE));
	changeReferenceFrame('N');
	rms_rl = full_search_R(b_x,b_y,block_size_x,block_size_y,con,&trans_rl);
	changeReferenceFrame(currentVideo);
	if (rms_rl<rms)
	{
		num_rl++;
		//	trans->reference=1;//�ο�����CĿ
		rms = rms_rl;
		trans->x=trans_rl.x;
		trans->y=trans_rl.y;
		trans->offset=trans_rl.offset;
		trans->scale=trans_rl.scale;
		trans->reference=3;//�ο�NĿ
		trans->block_type=0;
		}

}


	if(rms> (tol*tol*no))//16x8��8x16�����ֵΪtol_8 * tol_8 * 128, 8x4��4x8�����ֵΪtol_8 * tol_8 * 32
		return 0;//16x8��8x16��(8x4��4x8)���ܹ�ƥ��
	else
	{
		return 1;//16x8��8x16��(8x4��4x8)�㹻ƥ��
	}
}


void encode_block_8(int b_x,int b_y,int Curblock_8,TRANS_NODE *trans,int con)//����8*8���ֵĿ�
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
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //�жϵ�ǰ�������ĸ�region
			}	
		region--;	//0-0region;1--1region;2--�߽��
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
		rms = full_search(b_x,b_y,8,8,con,trans);//ȫ���� 8x8�� ����trans
	}
	if (search_mode==1)
	{
		rms =new_hexagon_block_search(b_x,b_y,8,8,con,trans);//������������ 8x8�� ����trans
	}
	if (search_mode==2)
	{
		rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,trans);  //UMHEX���� 8x8�� ����trans
	}
	if (search_mode==3)
	{
        rms =hexagon_block_search(b_x,b_y,8,8,con,trans);//���������� 8x8�� ����trans
	}

	if (currentVideo=='C')//�����CĿ�������÷������ؼ���rms
// 	if (currentVideo!='C')////��������������ZZZZZZZZZZZZZZZZZZZ
	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');//�ο�HĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,8,8,con,&trans_rl);//ȫ���� 8x8�� ����trans
		}
		if (search_mode==1)
		{
			rms_rl =new_hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//������������ 8x8�� ����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,&trans_rl);  //UMHEX���� 8x8�� ����trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//���������� 8x8�� ����trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)//����trans
		{
			num_rl++;
			trans->reference=1;//�ο�����HĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
			trans->reference=0;//�ο�����CĿ
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');//�ο�MĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,8,8,con,&trans_rl);//ȫ���� 8x8�� ����trans
		}
		if (search_mode==1)
		{
			rms_rl =new_hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//������������ 8x8�� ����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,&trans_rl);  //UMHEX���� 8x8�� ����trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,8,8,con,&trans_rl); //���������� 8x8�� ����trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)//����trans
		{
			num_rl++;
			trans->reference=2;//�ο�����MĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');//�ο�NĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,8,8,con,&trans_rl);//ȫ���� 8x8�� ����trans
		}
		if (search_mode==1)
		{
			rms_rl =new_hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//������������ 8x8�� ����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,8,8,con,&trans_rl);  //UMHEX���� 8x8�� ����trans
		}
		if (search_mode==3)
		{
			rms_rl =hexagon_block_search(b_x,b_y,8,8,con,&trans_rl);//���������� 8x8�� ����trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)//����trans
		{
			num_rl++;
			trans->reference=3;//�ο�����NĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
	}
    else//��ǰĿ��CĿ
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
			trans->reference=0;//�ο�����CĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
		{	trans->block_type=0;
		    trans->reference=4;//�ο���ǰĿ
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');
		rms_rl = full_search_R(b_x,b_y,8,8,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=1;//�ο�����HĿ
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
			trans->reference=2;//�ο�����MĿ
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
			trans->reference=3;//�ο�����NĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}

	}
	
	if(rms> (tol*tol*no))//��8*8����ƥ�䣬�������»���Ϊ8*4��4*8��4*4��ƥ�����ֵΪtol_8*tol_8*64��
	{
		int skip_old=skip;
		for(mode=1;mode<3;mode++)//8*4��4*8����
		{
			rectRms=0;
			trans->partition=mode;
			for(i=0;i<2;i++)//8*4��4*8���ֶ���1��8*8���ֳ�������С�飬���ԣ�����Ҫѭ������
			{		
				trans->next[i].x=0;
				trans->next[i].y=0;
				if(encode_block_rect(b_x,b_y,i,&trans->next[i],con,mode,2)==0)//8*4����4*8�����ֲ���ƥ��
				{
					i=4;//�������С���е�һ������ƥ�������for(i=0;i<2;i++)ѭ��
				}
				else
					rectRms++;
			}
			if(rectRms==2)//����8*4����4*8��С�鶼��ƥ��
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
				
				mode=4;//����С�飨8*4��4*8������ƥ�䣬mode=4
// 				if(region!=1)
// 				num_rect1+=2;
// 				num_rect2+=2;
			}
		}
		if(mode<4)//����4*4����
		{
			trans->partition=3;
			partition[obj][partition_length[obj]]=3;
		    partition_length[obj]++;			
			for(i=0;i<2;i++)//ѭ��4��4*4С��
				for(j=0;j<2;j++)
				{
					encode_block_4(b_x+j*4,b_y+i*4,&trans->next[i*2+j],con);//4*4����
				}
		}
	}
	else//��8*8����ƥ��
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
	//��ǰtrans��4x4���trans
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
				region|=1<<(int)((double)plane[j][i]/(double)GREY_LEVELS); //�жϵ�ǰ�������ĸ�region
			}	
		region--;	//0-0region;1--1region;2--�߽��
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
		rms = full_search(b_x,b_y,4,4,con,trans);//ȫ���� 4x4�� ����trans
	}
	if (search_mode==1)
	{
		rms =new_hexagon_block_search(b_x,b_y,4,4,con,trans);//������������ 4x4�� ����trans
	}
	if (search_mode==2)
	{
		rms =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,trans);  //UMHEX���� 4x4�� ����trans
	}
	if (search_mode==3)
	   {
        rms =hexagon_block_search(b_x,b_y,4,4,con,trans);//���������� 4x4�� ����trans
	   }

	if (currentVideo=='C')//�����CĿ�������÷������� ������rms
// 	if (currentVideo!='C')////��������������ZZZZZZZZZZZZZZZZZZZ
	{
		double rms_rl;
		TRANS_NODE trans_rl;
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');//�ο�HĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,4,4,con,&trans_rl);//ȫ���� 4x4�� ����trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//������������ 4x4�� ����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,&trans_rl);   //UMHEX���� 4x4�� ����trans
		}
		if (search_mode==3)
		{
			rms_rl = hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//���������� 4x4�� ����trans
		}
		
		changeReferenceFrame('C');
		if (rms_rl<rms)
		{
			num_rl++;
			trans->partition=trans->reference=1;//�ο�����HĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
			trans->reference=0;//�ο�����CĿ
		
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('M');//�ο�MĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,4,4,con,&trans_rl);//ȫ���� 4x4�� ����trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//������������ 4x4�� ����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,&trans_rl);  //UMHEX���� 4x4�� ����trans
		}
		if (search_mode==3)
		{
			rms_rl = hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//���������� 4x4�� ����trans
		}
		changeReferenceFrame('C');
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=2;//�ο�����MĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('N');//�ο�NĿ
		if (search_mode==0)
		{
			rms_rl = full_search(b_x,b_y,4,4,con,&trans_rl);//ȫ���� 4x4�� ����trans
		}
		if (search_mode==1)
		{
			rms_rl = new_hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//������������ 4x4�� ����trans
		}
		if (search_mode==2)
		{
			rms_rl =UMHEXIntegerPelBlockMotionSearch(b_x,b_y,4,4,con,&trans_rl);   //UMHEX���� 4x4�� ����trans
		}
		if (search_mode==3)
		{
			rms_rl = hexagon_block_search(b_x,b_y,4,4,con,&trans_rl);//���������� 4x4�� ����trans
		}
		changeReferenceFrame('C');
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=3;//�ο�����NĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
	}
    else//����CĿ �Ͳ��÷������أ�������Ӧ���ǹ�
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
			trans->reference=0;//�ο�����CĿ
			rms = rms_rl;
			trans->x=trans_rl.x;
			trans->y=trans_rl.y;
			trans->offset=trans_rl.offset;
			trans->scale=trans_rl.scale;
			trans->block_type=0;
		}
		else
		{	trans->block_type=0;
	    	trans->reference=4;//�ο����ǵ�ǰĿ
		}
		memset(&trans_rl,0,sizeof(TRANS_NODE));
		changeReferenceFrame('H');
		rms_rl = full_search_R(b_x,b_y,4,4,con,&trans_rl);
		changeReferenceFrame(currentVideo);
		if (rms_rl<rms)
		{
			num_rl++;
			trans->reference=1;//�ο�����HĿ
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
			trans->reference=2;//�ο�����MĿ
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
			trans->reference=3;//�ο�����NĿ
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
    int window_size=input->search_range;//�����������Χ
	
	best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);//����mv=��0,0��λ�õ�domain���range���RMSƥ����con=1:Y����;con=2��U����;con=3��V����
	trans->scale=alpha;//�����β���alpha����trans
	trans->offset=beta;//�����β���beta����trans
	for(l=1;l<=window_size;l++)//������������,window_size���������ļ������������Χ
	{
		i=j=-l;//��һ��������MV��(-1,-1)
		for(k=0;k<8*l;k++)//��ÿһ������Ȧ��������ĸ���
		{
			m=block_x+i;//domain���x����
			n=block_y+j;//domain���y����

			if(bound_chk(m,n,block_x,block_y,block_size_x,block_size_y,con))//���domain���Ƿ񳬳�ͼ��߽磬(m,n)��domain������꣬��block_x,block_y����range�������
			{
				rms=compute_rms(block_x,block_y,m,n,&alpha,&beta,block_size_x,block_size_y,con);//����domain�飨m,n����range��(block_x,block_y)��RMSƥ�������β���alpha��beta
				if(rms<best_rms)//����RMS�����β���
				{
					best_rms=rms;
					trans->x=m-block_x;
					trans->y=n-block_y;
					trans->scale=alpha;
					trans->offset=beta;
				}
			}
			// spt++;
			
			if(k<2*l) //�����Ǹ�������������
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
// 	for(l=1;l<=window_size;l++)//������������
// 	{
// 		i=j=-l;
// 		for(k=0;k<8*l;k++)//��ÿһ������Ȧ��������ĸ���
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
// 	for(l=1;l<=window_size;l++)//������������
// 	{
// 		i=j=-l;
// 		for(k=0;k<8*l;k++)//��ÿһ������Ȧ��������ĸ���
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
   static MV hex[6]={{-2,0},{-1,-2},{1,-2},{2,0},{1,2},{-1,2}};//�����ε�����������
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


//��������ƫ�õ����������·���////////////////////////////////////////////////////////////
//double horizon_bias_diamond_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
double hexagon_block_search(int block_x,int block_y,int block_size_x,int block_size_y,int con,TRANS_NODE *trans)
{
	int i,j,k,m,n,c=9;//b,err0,t,,l
	MV mv,tmv;	
	double best_rms=1e30;
	double rms=0.0;
	double alpha,beta;
	//static MV bias[4]={{-2,0},{0,-1},{2,0},{0,1}};//ƫ�����ε�4����/////////////////
	static MV fbx[5]={{-2,0},{-1,-1},{1,-1},{2,0},{0,1}};//����ε�5����/////////////////
	//static MV scr[4]={{-1,0},{0,-1},{1,0},{0,1}};//Сʮ���ε�4����////////////
	static MV bcr[4]={{-1,0},{0,-2},{1,0},{0,2}};//��ֱƫ�����ε�4����////////////
	
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
     // n=block_y+fbx[i].y;                 //////// ����ε�5����
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
	   for(i=0;i<3;i++)    ///////////////////////û�ã���ɾ���ģ��Ǻ�/////////////////
		                                 //���ˣ����ã����������������ж�
	   {
		   // m=block_x+mv.x+fbx[j].x;
          // n=block_y+mv.y+fbx[j].y;     //4 С����///////////////////////////////
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

//��ֱƫ�����ε�4����////////////////////////////////////////
 
  mv.x=mv.x+tmv.x;
  mv.y=mv.y+tmv.y;
   
   tmv.x=0; tmv.y=0;
   
   for(j=0;j<4;j++)
   {
	   // m=block_x+mv.x+bcr[j].x;
	   //n=block_y+mv.y+bcr[j].y;     //Сʮ������������������//////////////////////////
	  m=motion.mv_x+mv.x+bcr[j].x;
	   n=motion.mv_y+mv.y+bcr[j].y;     //��ֱƫ�����ε�4����//////////////////////////////////////	
	   
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
   static MV dia[8]={{-2,0},{-1,-1},{0,-2},{1,-1},{2,0},{1,1},{0,2},{-1,1}};//Ӧ�õ�4 points of big_diamond

   static MV hex[6]={{-2,0},{-1,-2},{1,-2},{2,0},{1,2},{-1,2}};//�����ε�����������


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
//��ǰ����������[m,n]��������ǰ֡ͼ��Χ��������һ������
{
	int ilow,ihigh,jlow,jhigh;
	int b=0;
	int window_size=input->search_range;//�����ļ������������Χ
    int width=input->imagewidth/min(2,con);//ͼ����
	int height=input->imageheight/min(2,con);//ͼ��߶�

	ilow=cx-window_size;//������Χ����߽�����
	ihigh=cx+window_size;//������Χ���ұ߽�����
	if(ilow<0)//���������Χ����߽糬��ͼ����
		ilow=0;
	if(ihigh>width-block_size_x)//���������Χ���ұ߽糬��ͼ����
		ihigh=width-block_size_x;
	
	jlow=cy-window_size;//������Χ���ϱ߽�����
	jhigh=cy+window_size;//������Χ���±߽�����
	if(jlow<0)//���������Χ���ϱ߽糬��ͼ����
		jlow=0;
	if(jhigh>height-block_size_y)//���������Χ���±߽糬��ͼ����
		jhigh=height-block_size_y;
 
	b=m<=ihigh&&m>=ilow&&n<=jhigh&&n>=jlow;//��������Χ����ͼ��߽磬b=0
	return(b);
}




// �µ��˶����ƣ���������ʦ�������EPZS(Enhanced Predictive Zonal Search)�㷨/////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

//�ϲ�ģʽԤ��/////////////////////////	   
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
      //�����ǰ���ϲ�ģʽ���˶�ʸ����������λ�ã����������Ŀ�Ϊƽ̹�顣
	   //����ƽ̹�飬�����������������˶�ʸ��������
    //  if(best_x%16==0||best_x%8==0||best_x%4==0)
	  {	  
		    	  //���ڿ�Ԥ��
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

				  //ԭʼ��������(0,0)Ԥ��
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
				  //�Ƚ�������ʼ�������Χ��4������������
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
			
				  //�����С�������֮��MSADλ��������ʼ�㣬��ֹͣ
				  if (best_x == iXMinNow && best_y == iYMinNow)
				 {
					       return best_rms;
				  }

				  //���ƥ���ʹ����ƥ�����ԣ���ѡ�����ƥ���MV
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
				  
			  //�����������ڵ�������ģ���ϵ��ƥ�����
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

					  //�Ե��Ĳ���������ģ���ϵ�������Ϊ���ģ�������ģ����������Χ��
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





//��ǰ�ķǶԳ��������㷨/////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
// ///////////////////�����/////////////////////////////////////////////
// 	static const int wbx_x[25] = {-2,-1, 0,1,2, -2,-1, 0,1,2,-2,-1, 0,1,2,-2,-1, 0,1,2,-2,-1, 0,1,2};
// 	static const int wbx_y[25] = {2,2,2,2,2,1,1,1,1,1,0,0,0,0,0,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2};
// 
//  /*sec_step:*/ //Unsymmetrical-cross search
//   //��Ӧ��Ԥ��
//   iXMinNow = block_x;
//   iYMinNow = block_y;
//   best_rms=compute_rms(block_x,block_y,block_x,block_y,&alpha,&beta,block_size_x,block_size_y,con);
//   trans->scale=alpha;
//   trans->offset=beta;
//   best_x=block_x;
//   best_y=block_y;
//   //�ϲ�ģʽԤ��
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
//   //LԤ��ģʽ
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