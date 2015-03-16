#include "compute.h"

#include "i_global.h"


double compute_rms(int block_x,int block_y,int m,int n,double *alpha,double *beta,int block_size_x,int block_size_y,int con)
{
	double rms=1e30;
	double det;
//          scale,    /* the scale factor */
//          offset;     /* The offset */
	int a;
	int aa;
	no=block_size_x*block_size_y;//一个宏块包含的像素数目
	dsum1=0,dsum2=0,rsum1=0,rsum2=0,rdsum=0;

	if(2==region)//no
	{
		classify(block_x,block_y,m,n,obj,block_size_x,block_size_y,con);
	}
	else//yes
	{
		if(block_size_x==16&&block_size_y==16)//16x16块
		{
			if(1==con)//Y分量
			{
				dsum1=sum_16_ref_temp[n][m];dsum2=sum2_16_ref_temp[n][m];//domain块像素和与像素平方和
				rsum1=sum_16_org[block_y/16][block_x/16];rsum2=sum2_16_org[block_y/16][block_x/16];//range块像素和与像素平方和
			}
			else if(2==con)//U分量
			{
				dsum1=sum_16_U_ref_temp[n][m];dsum2=sum2_16_U_ref_temp[n][m];
				rsum1=sum_16_U_org[block_y/16][block_x/16];rsum2=sum2_16_U_org[block_y/16][block_x/16];
			}
			else//V分量
			{
				dsum1=sum_16_V_ref_temp[n][m];dsum2=sum2_16_V_ref_temp[n][m];
				rsum1=sum_16_V_org[block_y/16][block_x/16];rsum2=sum2_16_V_org[block_y/16][block_x/16];
			}
		}
		
		if(block_size_x==8&&block_size_y==8)//8x8块
		{
			if(1==con)
			{
				dsum1=sum_8_ref_temp[n][m];dsum2=sum2_8_ref_temp[n][m];
				rsum1=sum_8_org[block_y/8][block_x/8];rsum2=sum2_8_org[block_y/8][block_x/8];
			}
			else if(2==con)
			{
				dsum1=sum_8_U_ref_temp[n][m];dsum2=sum2_8_U_ref_temp[n][m];
				rsum1=sum_8_U_org[block_y/8][block_x/8];rsum2=sum2_8_U_org[block_y/8][block_x/8];
			}
			else
			{
				dsum1=sum_8_V_ref_temp[n][m];dsum2=sum2_8_V_ref_temp[n][m];
				rsum1=sum_8_V_org[block_y/8][block_x/8];rsum2=sum2_8_V_org[block_y/8][block_x/8];
			}
		}
		
		if(block_size_x==4&&block_size_y==4)//4x4块
		{
			if(1==con)
			{
				dsum1=sum_4_ref_temp[n][m];dsum2=sum2_4_ref_temp[n][m];
				rsum1=sum_4_org[block_y/4][block_x/4];rsum2=sum2_4_org[block_y/4][block_x/4];
			}
			else if(2==con)
			{
				dsum1=sum_4_U_ref_temp[n][m];dsum2=sum2_4_U_ref_temp[n][m];
				rsum1=sum_4_U_org[block_y/4][block_x/4];rsum2=sum2_4_U_org[block_y/4][block_x/4];
			}
			else
			{
				dsum1=sum_4_V_ref_temp[n][m];dsum2=sum2_4_V_ref_temp[n][m];
				rsum1=sum_4_V_org[block_y/4][block_x/4];rsum2=sum2_4_V_org[block_y/4][block_x/4];
			}
		}

		if(block_size_x==16&&block_size_y==8)//16x16块
		{
			if(1==con)
			{
				dsum1=sum_16_8_ref_temp[n][m];dsum2=sum2_16_8_ref_temp[n][m];
				rsum1=sum_16_8_org[block_y/8][block_x/16];rsum2=sum2_16_8_org[block_y/8][block_x/16];
			}
			else if(2==con)
			{
				dsum1=sum_16_8_U_ref_temp[n][m];dsum2=sum2_16_8_U_ref_temp[n][m];
				rsum1=sum_16_8_U_org[block_y/8][block_x/16];rsum2=sum2_16_8_U_org[block_y/8][block_x/16];
			}
			else
			{
				dsum1=sum_16_8_V_ref_temp[n][m];dsum2=sum2_16_8_V_ref_temp[n][m];
				rsum1=sum_16_8_V_org[block_y/8][block_x/16];rsum2=sum2_16_8_V_org[block_y/8][block_x/16];
			}
		}

		if(block_size_x==8&&block_size_y==16)//8x16块
		{
			if(1==con)
			{
				dsum1=sum_8_16_ref_temp[n][m];dsum2=sum2_8_16_ref_temp[n][m];
				rsum1=sum_8_16_org[block_y/16][block_x/8];rsum2=sum2_8_16_org[block_y/16][block_x/8];
			}
			else if(2==con)
			{
				dsum1=sum_8_16_U_ref_temp[n][m];dsum2=sum2_8_16_U_ref_temp[n][m];
				rsum1=sum_8_16_U_org[block_y/16][block_x/8];rsum2=sum2_8_16_U_org[block_y/16][block_x/8];
			}
			else
			{
				dsum1=sum_8_16_V_ref_temp[n][m];dsum2=sum2_8_16_V_ref_temp[n][m];
				rsum1=sum_8_16_V_org[block_y/16][block_x/8];rsum2=sum2_8_16_V_org[block_y/16][block_x/8];
			}
		}
	
		if(block_size_x==8&&block_size_y==4)//8x4块
		{
			if(1==con)
			{
				dsum1=sum_8_4_ref_temp[n][m];dsum2=sum2_8_4_ref_temp[n][m];
				rsum1=sum_8_4_org[block_y/4][block_x/8];rsum2=sum2_8_4_org[block_y/4][block_x/8];
			}
			else if(2==con)
			{
				dsum1=sum_8_4_U_ref_temp[n][m];dsum2=sum2_8_4_U_ref_temp[n][m];
				rsum1=sum_8_4_U_org[block_y/4][block_x/8];rsum2=sum2_8_4_U_org[block_y/4][block_x/8];
			}
			else
			{
				dsum1=sum_8_4_V_ref_temp[n][m];dsum2=sum2_8_4_V_ref_temp[n][m];
				rsum1=sum_8_4_V_org[block_y/4][block_x/8];rsum2=sum2_8_4_V_org[block_y/4][block_x/8];
			}
		}

		if(block_size_x==4&&block_size_y==8)//4x8块
		{
			if(1==con)
			{
				dsum1=sum_4_8_ref_temp[n][m];dsum2=sum2_4_8_ref_temp[n][m];
				rsum1=sum_4_8_org[block_y/8][block_x/4];rsum2=sum2_4_8_org[block_y/8][block_x/4];
			}
			else if(2==con)
			{
				dsum1=sum_4_8_U_ref_temp[n][m];dsum2=sum2_4_8_U_ref_temp[n][m];
				rsum1=sum_4_8_U_org[block_y/8][block_x/4];rsum2=sum2_4_8_U_org[block_y/8][block_x/4];
			}
			else
			{
				dsum1=sum_4_8_V_ref_temp[n][m];dsum2=sum2_4_8_V_ref_temp[n][m];
				rsum1=sum_4_8_V_org[block_y/8][block_x/4];rsum2=sum2_4_8_V_org[block_y/8][block_x/4];
			}
		}

		rdsum=compute_rdSum(block_x,block_y,m,n,block_size_x,block_size_y,con); //计算domain块和range块的像素乘积之和
	}

	det = no*dsum2 - dsum1*dsum1;
	if (det == 0.0) /* variance is 0 and no scaling is needed (s = 0) */
	{	*alpha = 0.0;} 
	else 
	{	*alpha =(no*rdsum - rsum1*dsum1)/det;//计算分形参数alpha
///////////////修改系数
// 		aa =(no*rdsum - rsum1*dsum1)/det*20;
// 		*alpha =aa/20;

	}
	a=(int)(*alpha*100);//alpha乘以100
	//计算分形参数beta？？？为什么根公式不一样呢？？？少减了后边的项   公式为：(rsum1-(*alpha)*dsum1)/no
	*beta = rsum1/no;   /* orthogonalization:[0 255], Used to be (rsum1-alpha*dsum1)/no. */

 	QUAN_A(a);//将a的个位数字量化为0或5
	QUAN_A(*beta);//将beta的个位数字量化为0或5
    *alpha=(double)(a)/100;//a除以100，回复到原来alpha的数量级
	//*beta=b/100;
	if(*alpha<MIN_ALPHA||*alpha>MAX_ALPHA)
		return rms;
	if (*beta<MIN_BETA||*beta>MAX_BETA)
		return rms;
	rms= rsum2 + (*alpha)*((*alpha)*dsum2 - 2.0*rdsum + 2.0*((*beta)-(*alpha)*dsum1/no)*dsum1) 
		+ ((*beta)-(*alpha)*dsum1/no)*(((*beta)-(*alpha)*dsum1/no)*no - 2.0*rsum1);//计算domain块和range块的RMS匹配误差，（因为上边beta少减了后边一项，所以计算RMS时要加上这项）

	if (rms<0.0)
	{
		printf("rms<0.0");
	}
	return rms;
}


double compute_rdSum(int block_x,int block_y,int m,int n,int block_size_x,int block_size_y,int con)
{
	
	int i,j;
	double rdsum=0.0;
    byte **img_domain,**img_range;

	switch(con)
	{
	case 1://Y分量
		img_domain=imgY_ref_temp;img_range=imgY_org;break;
	case 2://U分量
		img_domain=imgUV_ref_temp[0];img_range=imgUV_org[0];break;
	case 3://V分量
		img_domain=imgUV_ref_temp[1];img_range=imgUV_org[1];break;
	}

	for(i=0;i<block_size_y;i++)       //计算domain块和range块的像素乘积之和
		for(j=0;j<block_size_x;j++)
		{
			rdsum+=img_range[block_y+i][block_x+j]*img_domain[n+i][m+j];
		}
	return rdsum;	
}


void classify(int block_x,int block_y,int m,int n,int obj,int block_size_x,int block_size_y,int con)
{
	int i1,j1,i2,j2,num_range=0,num_domain=0;
	byte **plane,**plane_domain;
    byte **img_domain,**img_range;
	unsigned char average_domain=0;

	switch(con)
	{
	case 1:
		img_domain=imgY_ref_temp;img_range=imgY_org;
		plane=plane_Y;plane_domain=plane_Y_domain_temp;break;
	case 2:
		img_domain=imgUV_ref_temp[0];img_range=imgUV_org[0];
		plane=plane_UV[0];plane_domain=plane_UV_domain_temp[0];break;
	case 3:
		img_domain=imgUV_ref_temp[1];img_range=imgUV_org[1];
		plane=plane_UV[1];plane_domain=plane_UV_domain_temp[1];break;
	}

	for (j1=block_y,j2=n; j1<block_y+block_size_y; ++j1,j2++)
		for (i1=block_x,i2=m; i1<block_x+block_size_x; ++i1,i2++)
		{
			if(plane[j1][i1]/GREY_LEVELS==obj)
			{
				rsum1+=img_range[j1][i1];
				rsum2+=img_range[j1][i1]*img_range[j1][i1];
				num_range++;

				if(plane_domain[j2][i2]/GREY_LEVELS==obj)
				{
					dsum1+=img_domain[j2][i2];
					dsum2+=img_domain[j2][i2]*img_domain[j2][i2];
					rdsum+=img_range[j1][i1]*img_domain[j2][i2];
					num_domain++;
				}
			}
		}
	no=num_range;
	if(num_domain>0)
		average_domain=(unsigned char)(dsum1/num_domain);
	else
		average_domain=0;

	for (j1=block_y,j2=n; j1<block_y+block_size_y; ++j1,j2++)
		for (i1=block_x,i2=m; i1<block_x+block_size_x; ++i1,i2++)
		{
			if(plane[j1][i1]/GREY_LEVELS==obj)
			{
				if(plane_domain[j2][i2]/GREY_LEVELS!=obj)
				{
					dsum1+=average_domain;
					dsum2+=average_domain*average_domain;
					rdsum+=img_range[j1][i1]*average_domain;
				}
			}
		}
}

void compute_domain_Sum()
{
	int i,j;
	int width=input->imagewidth-3,height=input->imageheight-3;
	///////////////////////////计算Y分量
	for(i=0;i<height;i++)               //计算4*4块sum，sum2   domain
		for(j=0;j<width;j++)
		{
			sum_4_ref_temp[i][j]=0.0;
			sum2_4_ref_temp[i][j]=0.0;

			sum_4_ref_temp[i][j]+=imgY_ref_temp[i][j];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i][j+1];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i][j+2];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i][j+3];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j+1];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j+2];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j+3];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j+1];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j+2];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j+3];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j+1];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j+2];
			sum_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j+3];

			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i][j]*imgY_ref_temp[i][j];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i][j+1]*imgY_ref_temp[i][j+1];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i][j+2]*imgY_ref_temp[i][j+2];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i][j+3]*imgY_ref_temp[i][j+3];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j]*imgY_ref_temp[i+1][j];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j+1]*imgY_ref_temp[i+1][j+1];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j+2]*imgY_ref_temp[i+1][j+2];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+1][j+3]*imgY_ref_temp[i+1][j+3];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j]*imgY_ref_temp[i+2][j];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j+1]*imgY_ref_temp[i+2][j+1];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j+2]*imgY_ref_temp[i+2][j+2];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+2][j+3]*imgY_ref_temp[i+2][j+3];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j]*imgY_ref_temp[i+3][j];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j+1]*imgY_ref_temp[i+3][j+1];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j+2]*imgY_ref_temp[i+3][j+2];
			sum2_4_ref_temp[i][j]+=imgY_ref_temp[i+3][j+3]*imgY_ref_temp[i+3][j+3];
		}

	width = input->imagewidth-7;//633
	height = input->imageheight-3;//477
	for(i=0;i<height;i++)        //计算8*4块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_4_ref_temp[i][j]=0.0;
			sum2_8_4_ref_temp[i][j]=0.0;

			sum_8_4_ref_temp[i][j]+=sum_4_ref_temp[i][j];
			sum_8_4_ref_temp[i][j]+=sum_4_ref_temp[i][j+4];

			sum2_8_4_ref_temp[i][j]+=sum2_4_ref_temp[i][j];
			sum2_8_4_ref_temp[i][j]+=sum2_4_ref_temp[i][j+4];
		}

	width = input->imagewidth-3;
	height = input->imageheight-7;
	for(i=0;i<height;i++)        //计算4*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_4_8_ref_temp[i][j]=0.0;
			sum2_4_8_ref_temp[i][j]=0.0;

			sum_4_8_ref_temp[i][j]+=sum_4_ref_temp[i][j];
			sum_4_8_ref_temp[i][j]+=sum_4_ref_temp[i+4][j];

			sum2_4_8_ref_temp[i][j]+=sum2_4_ref_temp[i][j];
			sum2_4_8_ref_temp[i][j]+=sum2_4_ref_temp[i+4][j];
		}

	width = input->imagewidth-7;
	height = input->imageheight-7;
	for(i=0;i<height;i++)        //计算8*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_ref_temp[i][j]=0.0;
			sum2_8_ref_temp[i][j]=0.0;

			sum_8_ref_temp[i][j]+=sum_8_4_ref_temp[i][j];
			sum_8_ref_temp[i][j]+=sum_8_4_ref_temp[i+4][j];

			sum2_8_ref_temp[i][j]+=sum2_8_4_ref_temp[i][j];
			sum2_8_ref_temp[i][j]+=sum2_8_4_ref_temp[i+4][j];
		}
   
	width = input->imagewidth-15;
	height = input->imageheight-7;
	for(i=0;i<height;i++)        //计算16*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_16_8_ref_temp[i][j]=0.0;
			sum2_16_8_ref_temp[i][j]=0.0;

			sum_16_8_ref_temp[i][j]+=sum_8_ref_temp[i][j];
			sum_16_8_ref_temp[i][j]+=sum_8_ref_temp[i][j+8];

			sum2_16_8_ref_temp[i][j]+=sum2_8_ref_temp[i][j];
			sum2_16_8_ref_temp[i][j]+=sum2_8_ref_temp[i][j+8];
		}

	width = input->imagewidth-7;
	height = input->imageheight-15;
	for(i=0;i<height;i++)        //计8*16块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_16_ref_temp[i][j]=0.0;
			sum2_8_16_ref_temp[i][j]=0.0;

			sum_8_16_ref_temp[i][j]+=sum_8_ref_temp[i][j];
			sum_8_16_ref_temp[i][j]+=sum_8_ref_temp[i+8][j];

			sum2_8_16_ref_temp[i][j]+=sum2_8_ref_temp[i][j];
			sum2_8_16_ref_temp[i][j]+=sum2_8_ref_temp[i+8][j];
		}

	width = input->imagewidth-15;
	height = input->imageheight-15;	
	for(i=0;i<height;i++)       //计算16*16块     domain
		for(j=0;j<width;j++)
		{
			sum_16_ref_temp[i][j]=0.0;
			sum2_16_ref_temp[i][j]=0.0;

			sum_16_ref_temp[i][j]+=sum_16_8_ref_temp[i][j];
			sum_16_ref_temp[i][j]+=sum_16_8_ref_temp[i+8][j];

			sum2_16_ref_temp[i][j]+=sum2_16_8_ref_temp[i][j];
			sum2_16_ref_temp[i][j]+=sum2_16_8_ref_temp[i+8][j];
		}

		//////////////////////////////////////////////////////////////////////////计算U分量	

	width=input->imagewidth/2-3,height=input->imageheight/2-3;
	
	for(i=0;i<height;i++)               //计算4*4块sum，sum2   domain
		for(j=0;j<width;j++)
		{
			sum_4_U_ref_temp[i][j]=0.0;
			sum2_4_U_ref_temp[i][j]=0.0;

			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j+1];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j+2];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j+3];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j+1];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j+2];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j+3];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j+1];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j+2];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j+3];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j+1];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j+2];
			sum_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j+3];

			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j]*imgUV_ref_temp[0][i][j];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j+1]*imgUV_ref_temp[0][i][j+1];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j+2]*imgUV_ref_temp[0][i][j+2];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i][j+3]*imgUV_ref_temp[0][i][j+3];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j]*imgUV_ref_temp[0][i+1][j];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j+1]*imgUV_ref_temp[0][i+1][j+1];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j+2]*imgUV_ref_temp[0][i+1][j+2];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+1][j+3]*imgUV_ref_temp[0][i+1][j+3];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j]*imgUV_ref_temp[0][i+2][j];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j+1]*imgUV_ref_temp[0][i+2][j+1];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j+2]*imgUV_ref_temp[0][i+2][j+2];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+2][j+3]*imgUV_ref_temp[0][i+2][j+3];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j]*imgUV_ref_temp[0][i+3][j];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j+1]*imgUV_ref_temp[0][i+3][j+1];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j+2]*imgUV_ref_temp[0][i+3][j+2];
			sum2_4_U_ref_temp[i][j]+=imgUV_ref_temp[0][i+3][j+3]*imgUV_ref_temp[0][i+3][j+3];
		}

	width = input->imagewidth/2-7;
	height = input->imageheight/2-3;
	for(i=0;i<height;i++)        //计算8*4块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_4_U_ref_temp[i][j]=0.0;
			sum2_8_4_U_ref_temp[i][j]=0.0;

			sum_8_4_U_ref_temp[i][j]+=sum_4_U_ref_temp[i][j];
			sum_8_4_U_ref_temp[i][j]+=sum_4_U_ref_temp[i][j+4];

			sum2_8_4_U_ref_temp[i][j]+=sum2_4_U_ref_temp[i][j];
			sum2_8_4_U_ref_temp[i][j]+=sum2_4_U_ref_temp[i][j+4];
		}

	width = input->imagewidth/2-3;
	height = input->imageheight/2-7;
	for(i=0;i<height;i++)        //计算4*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_4_8_U_ref_temp[i][j]=0.0;
			sum2_4_8_U_ref_temp[i][j]=0.0;

			sum_4_8_U_ref_temp[i][j]+=sum_4_U_ref_temp[i][j];
			sum_4_8_U_ref_temp[i][j]+=sum_4_U_ref_temp[i+4][j];

			sum2_4_8_U_ref_temp[i][j]+=sum2_4_U_ref_temp[i][j];
			sum2_4_8_U_ref_temp[i][j]+=sum2_4_U_ref_temp[i+4][j];
		}

	width = input->imagewidth/2-7;
	height = input->imageheight/2-7;
	for(i=0;i<height;i++)        //计算8*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_U_ref_temp[i][j]=0.0;
			sum2_8_U_ref_temp[i][j]=0.0;

			sum_8_U_ref_temp[i][j]+=sum_8_4_U_ref_temp[i][j];
			sum_8_U_ref_temp[i][j]+=sum_8_4_U_ref_temp[i+4][j];

			sum2_8_U_ref_temp[i][j]+=sum2_8_4_U_ref_temp[i][j];
			sum2_8_U_ref_temp[i][j]+=sum2_8_4_U_ref_temp[i+4][j];
		}
   
	width = input->imagewidth/2-15;
	height = input->imageheight/2-7;
	for(i=0;i<height;i++)        //计算16*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_16_8_U_ref_temp[i][j]=0.0;
			sum2_16_8_U_ref_temp[i][j]=0.0;

			sum_16_8_U_ref_temp[i][j]+=sum_8_U_ref_temp[i][j];
			sum_16_8_U_ref_temp[i][j]+=sum_8_U_ref_temp[i][j+8];

			sum2_16_8_U_ref_temp[i][j]+=sum2_8_U_ref_temp[i][j];
			sum2_16_8_U_ref_temp[i][j]+=sum2_8_U_ref_temp[i][j+8];
		}

	width = input->imagewidth/2-7;
	height = input->imageheight/2-15;
	for(i=0;i<height;i++)        //计8*16块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_16_U_ref_temp[i][j]=0.0;
			sum2_8_16_U_ref_temp[i][j]=0.0;

			sum_8_16_U_ref_temp[i][j]+=sum_8_U_ref_temp[i][j];
			sum_8_16_U_ref_temp[i][j]+=sum_8_U_ref_temp[i+8][j];

			sum2_8_16_U_ref_temp[i][j]+=sum2_8_U_ref_temp[i][j];
			sum2_8_16_U_ref_temp[i][j]+=sum2_8_U_ref_temp[i+8][j];
		}

	width = input->imagewidth/2-15;
	height = input->imageheight/2-15;	
	for(i=0;i<height;i++)       //计算16*16块     domain
		for(j=0;j<width;j++)
		{
			sum_16_U_ref_temp[i][j]=0.0;
			sum2_16_U_ref_temp[i][j]=0.0;

			sum_16_U_ref_temp[i][j]+=sum_16_8_U_ref_temp[i][j];
			sum_16_U_ref_temp[i][j]+=sum_16_8_U_ref_temp[i+8][j];

			sum2_16_U_ref_temp[i][j]+=sum2_16_8_U_ref_temp[i][j];
			sum2_16_U_ref_temp[i][j]+=sum2_16_8_U_ref_temp[i+8][j];
		}
		
			//////////////////////////////////////////////////////////////////////////计算V分量	

	width=input->imagewidth/2-3;
	height=input->imageheight/2-3;
	
	for(i=0;i<height;i++)               //计算4*4块sum，sum2   domain
		for(j=0;j<width;j++)
		{
			sum_4_V_ref_temp[i][j]=0.0;
			sum2_4_V_ref_temp[i][j]=0.0;

			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j+1];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j+2];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j+3];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j+1];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j+2];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j+3];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j+1];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j+2];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j+3];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j+1];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j+2];
			sum_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j+3];

			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j]*imgUV_ref_temp[1][i][j];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j+1]*imgUV_ref_temp[1][i][j+1];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j+2]*imgUV_ref_temp[1][i][j+2];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i][j+3]*imgUV_ref_temp[1][i][j+3];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j]*imgUV_ref_temp[1][i+1][j];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j+1]*imgUV_ref_temp[1][i+1][j+1];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j+2]*imgUV_ref_temp[1][i+1][j+2];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+1][j+3]*imgUV_ref_temp[1][i+1][j+3];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j]*imgUV_ref_temp[1][i+2][j];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j+1]*imgUV_ref_temp[1][i+2][j+1];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j+2]*imgUV_ref_temp[1][i+2][j+2];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+2][j+3]*imgUV_ref_temp[1][i+2][j+3];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j]*imgUV_ref_temp[1][i+3][j];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j+1]*imgUV_ref_temp[1][i+3][j+1];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j+2]*imgUV_ref_temp[1][i+3][j+2];
			sum2_4_V_ref_temp[i][j]+=imgUV_ref_temp[1][i+3][j+3]*imgUV_ref_temp[1][i+3][j+3];
		}

	width = input->imagewidth/2-7;
	height = input->imageheight/2-3;
	for(i=0;i<height;i++)        //计算8*4块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_4_V_ref_temp[i][j]=0.0;
			sum2_8_4_V_ref_temp[i][j]=0.0;

			sum_8_4_V_ref_temp[i][j]+=sum_4_V_ref_temp[i][j];
			sum_8_4_V_ref_temp[i][j]+=sum_4_V_ref_temp[i][j+4];

			sum2_8_4_V_ref_temp[i][j]+=sum2_4_V_ref_temp[i][j];
			sum2_8_4_V_ref_temp[i][j]+=sum2_4_V_ref_temp[i][j+4];
		}

	width = input->imagewidth/2-3;
	height = input->imageheight/2-7;
	for(i=0;i<height;i++)        //计算4*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_4_8_V_ref_temp[i][j]=0.0;
			sum2_4_8_V_ref_temp[i][j]=0.0;

			sum_4_8_V_ref_temp[i][j]+=sum_4_V_ref_temp[i][j];
			sum_4_8_V_ref_temp[i][j]+=sum_4_V_ref_temp[i+4][j];

			sum2_4_8_V_ref_temp[i][j]+=sum2_4_V_ref_temp[i][j];
			sum2_4_8_V_ref_temp[i][j]+=sum2_4_V_ref_temp[i+4][j];
		}

	width = input->imagewidth/2-7;
	height = input->imageheight/2-7;
	for(i=0;i<height;i++)        //计算8*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_V_ref_temp[i][j]=0.0;
			sum2_8_V_ref_temp[i][j]=0.0;

			sum_8_V_ref_temp[i][j]+=sum_8_4_V_ref_temp[i][j];
			sum_8_V_ref_temp[i][j]+=sum_8_4_V_ref_temp[i+4][j];

			sum2_8_V_ref_temp[i][j]+=sum2_8_4_V_ref_temp[i][j];
			sum2_8_V_ref_temp[i][j]+=sum2_8_4_V_ref_temp[i+4][j];
		}
   
	width = input->imagewidth/2-15;
	height = input->imageheight/2-7;
	for(i=0;i<height;i++)        //计算16*8块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_16_8_V_ref_temp[i][j]=0.0;
			sum2_16_8_V_ref_temp[i][j]=0.0;

			sum_16_8_V_ref_temp[i][j]+=sum_8_V_ref_temp[i][j];
			sum_16_8_V_ref_temp[i][j]+=sum_8_V_ref_temp[i][j+8];

			sum2_16_8_V_ref_temp[i][j]+=sum2_8_V_ref_temp[i][j];
			sum2_16_8_V_ref_temp[i][j]+=sum2_8_V_ref_temp[i][j+8];
		}

	width = input->imagewidth/2-7;
	height = input->imageheight/2-15;
	for(i=0;i<height;i++)        //计8*16块sum，sum2	domain
		for(j=0;j<width;j++)
		{
			sum_8_16_V_ref_temp[i][j]=0.0;
			sum2_8_16_V_ref_temp[i][j]=0.0;

			sum_8_16_V_ref_temp[i][j]+=sum_8_V_ref_temp[i][j];
			sum_8_16_V_ref_temp[i][j]+=sum_8_V_ref_temp[i+8][j];

			sum2_8_16_V_ref_temp[i][j]+=sum2_8_V_ref_temp[i][j];
			sum2_8_16_V_ref_temp[i][j]+=sum2_8_V_ref_temp[i+8][j];
		}

	width = input->imagewidth/2-15;
	height = input->imageheight/2-15;	
	for(i=0;i<height;i++)       //计算16*16块     domain
		for(j=0;j<width;j++)
		{
			sum_16_V_ref_temp[i][j]=0.0;
			sum2_16_V_ref_temp[i][j]=0.0;

			sum_16_V_ref_temp[i][j]+=sum_16_8_V_ref_temp[i][j];
			sum_16_V_ref_temp[i][j]+=sum_16_8_V_ref_temp[i+8][j];

			sum2_16_V_ref_temp[i][j]+=sum2_16_8_V_ref_temp[i][j];
			sum2_16_V_ref_temp[i][j]+=sum2_16_8_V_ref_temp[i+8][j];
		}

}

void compute_range_Sum()
{
	int i,j;
	int width,height;
	///////////////////////////计算Y分量

	width = img->frmWidthInMbs*4;
	height = img->frmHeightInMbs*4;		
	for(i=0;i<height;i++)       //计算4*4块        range
		for(j=0;j<width;j++)
		{
			sum_4_org[i][j]=0.0;
			sum2_4_org[i][j]=0.0;

			sum_4_org[i][j]+=imgY_org[i*4][j*4];
			sum_4_org[i][j]+=imgY_org[i*4][j*4+1];
			sum_4_org[i][j]+=imgY_org[i*4][j*4+2];
			sum_4_org[i][j]+=imgY_org[i*4][j*4+3];
			sum_4_org[i][j]+=imgY_org[i*4+1][j*4];
			sum_4_org[i][j]+=imgY_org[i*4+1][j*4+1];
			sum_4_org[i][j]+=imgY_org[i*4+1][j*4+2];
			sum_4_org[i][j]+=imgY_org[i*4+1][j*4+3];
			sum_4_org[i][j]+=imgY_org[i*4+2][j*4];
			sum_4_org[i][j]+=imgY_org[i*4+2][j*4+1];
			sum_4_org[i][j]+=imgY_org[i*4+2][j*4+2];
			sum_4_org[i][j]+=imgY_org[i*4+2][j*4+3];
			sum_4_org[i][j]+=imgY_org[i*4+3][j*4];
			sum_4_org[i][j]+=imgY_org[i*4+3][j*4+1];
			sum_4_org[i][j]+=imgY_org[i*4+3][j*4+2];
			sum_4_org[i][j]+=imgY_org[i*4+3][j*4+3];

			sum2_4_org[i][j]+=imgY_org[i*4][j*4]*imgY_org[i*4][j*4];
			sum2_4_org[i][j]+=imgY_org[i*4][j*4+1]*imgY_org[i*4][j*4+1];
			sum2_4_org[i][j]+=imgY_org[i*4][j*4+2]*imgY_org[i*4][j*4+2];
			sum2_4_org[i][j]+=imgY_org[i*4][j*4+3]*imgY_org[i*4][j*4+3];
			sum2_4_org[i][j]+=imgY_org[i*4+1][j*4]*imgY_org[i*4+1][j*4];
			sum2_4_org[i][j]+=imgY_org[i*4+1][j*4+1]*imgY_org[i*4+1][j*4+1];
			sum2_4_org[i][j]+=imgY_org[i*4+1][j*4+2]*imgY_org[i*4+1][j*4+2];
			sum2_4_org[i][j]+=imgY_org[i*4+1][j*4+3]*imgY_org[i*4+1][j*4+3];
			sum2_4_org[i][j]+=imgY_org[i*4+2][j*4]*imgY_org[i*4+2][j*4];
			sum2_4_org[i][j]+=imgY_org[i*4+2][j*4+1]*imgY_org[i*4+2][j*4+1];
			sum2_4_org[i][j]+=imgY_org[i*4+2][j*4+2]*imgY_org[i*4+2][j*4+2];
			sum2_4_org[i][j]+=imgY_org[i*4+2][j*4+3]*imgY_org[i*4+2][j*4+3];
			sum2_4_org[i][j]+=imgY_org[i*4+3][j*4]*imgY_org[i*4+3][j*4];
			sum2_4_org[i][j]+=imgY_org[i*4+3][j*4+1]*imgY_org[i*4+3][j*4+1];
			sum2_4_org[i][j]+=imgY_org[i*4+3][j*4+2]*imgY_org[i*4+3][j*4+2];
			sum2_4_org[i][j]+=imgY_org[i*4+3][j*4+3]*imgY_org[i*4+3][j*4+3];
		}

	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs*4;		
	for(i=0;i<height;i++)       //计算8*4块   range
		for(j=0;j<width;j++)
		{
			sum_8_4_org[i][j]=0.0;
			sum2_8_4_org[i][j]=0.0;

			sum_8_4_org[i][j]+=sum_4_org[i][j<<1];
			sum_8_4_org[i][j]+=sum_4_org[i][(j<<1)+1];

			sum2_8_4_org[i][j]+=sum2_4_org[i][j<<1];
			sum2_8_4_org[i][j]+=sum2_4_org[i][(j<<1)+1];
		}
	
	width = img->frmWidthInMbs*4;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算4*8块   range
		for(j=0;j<width;j++)
		{
			sum_4_8_org[i][j]=0.0;
			sum2_4_8_org[i][j]=0.0;

			sum_4_8_org[i][j]+=sum_4_org[i<<1][j];
			sum_4_8_org[i][j]+=sum_4_org[(i<<1)+1][j];

			sum2_4_8_org[i][j]+=sum2_4_org[i<<1][j];
			sum2_4_8_org[i][j]+=sum2_4_org[(i<<1)+1][j];
		}

	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算8*8块   range
		for(j=0;j<width;j++)
		{
			sum_8_org[i][j]=0.0;
			sum2_8_org[i][j]=0.0;

			sum_8_org[i][j]+=sum_8_4_org[i<<1][j];
			sum_8_org[i][j]+=sum_8_4_org[(i<<1)+1][j];

			sum2_8_org[i][j]+=sum2_8_4_org[i<<1][j];
			sum2_8_org[i][j]+=sum2_8_4_org[(i<<1)+1][j];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算16*8块   range
		for(j=0;j<width;j++)
		{
			sum_16_8_org[i][j]=0.0;
			sum2_16_8_org[i][j]=0.0;

			sum_16_8_org[i][j]+=sum_8_org[i][j<<1];
			sum_16_8_org[i][j]+=sum_8_org[i][(j<<1)+1];

			sum2_16_8_org[i][j]+=sum2_8_org[i][j<<1];
			sum2_16_8_org[i][j]+=sum2_8_org[i][(j<<1)+1];
		}

	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算8*16块   range
		for(j=0;j<width;j++)
		{
			sum_8_16_org[i][j]=0.0;
			sum2_8_16_org[i][j]=0.0;

			sum_8_16_org[i][j]+=sum_8_org[i<<1][j];
			sum_8_16_org[i][j]+=sum_8_org[(i<<1)+1][j];

			sum2_8_16_org[i][j]+=sum2_8_org[i<<1][j];
			sum2_8_16_org[i][j]+=sum2_8_org[(i<<1)+1][j];
		}
		
	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算16*16块   range
		for(j=0;j<width;j++)
		{
			sum_16_org[i][j]=0.0;
			sum2_16_org[i][j]=0.0;

			sum_16_org[i][j]+=sum_16_8_org[i<<1][j];
			sum_16_org[i][j]+=sum_16_8_org[(i<<1)+1][j];

			sum2_16_org[i][j]+=sum2_16_8_org[i<<1][j];
			sum2_16_org[i][j]+=sum2_16_8_org[(i<<1)+1][j];
		}

    ///////////////////////////////////////////////////计算U分量
	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算4*4块
		for(j=0;j<width;j++)
		{
			sum_4_U_org[i][j]=0.0;
			sum2_4_U_org[i][j]=0.0;

			sum_4_U_org[i][j]+=imgUV_org[0][i*4][j*4];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4][j*4+1];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4][j*4+2];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4][j*4+3];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4+1];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4+2];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4+3];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4+1];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4+2];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4+3];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4+1];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4+2];
			sum_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4+3];

			sum2_4_U_org[i][j]+=imgUV_org[0][i*4][j*4]*imgUV_org[0][i*4][j*4];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4][j*4+1]*imgUV_org[0][i*4][j*4+1];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4][j*4+2]*imgUV_org[0][i*4][j*4+2];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4][j*4+3]*imgUV_org[0][i*4][j*4+3];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4]*imgUV_org[0][i*4+1][j*4];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4+1]*imgUV_org[0][i*4+1][j*4+1];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4+2]*imgUV_org[0][i*4+1][j*4+2];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+1][j*4+3]*imgUV_org[0][i*4+1][j*4+3];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4]*imgUV_org[0][i*4+2][j*4];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4+1]*imgUV_org[0][i*4+2][j*4+1];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4+2]*imgUV_org[0][i*4+2][j*4+2];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+2][j*4+3]*imgUV_org[0][i*4+2][j*4+3];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4]*imgUV_org[0][i*4+3][j*4];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4+1]*imgUV_org[0][i*4+3][j*4+1];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4+2]*imgUV_org[0][i*4+3][j*4+2];
			sum2_4_U_org[i][j]+=imgUV_org[0][i*4+3][j*4+3]*imgUV_org[0][i*4+3][j*4+3];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算8*4块   range
		for(j=0;j<width;j++)
		{
			sum_8_4_U_org[i][j]=0.0;
			sum2_8_4_U_org[i][j]=0.0;

			sum_8_4_U_org[i][j]+=sum_4_U_org[i][j<<1];
			sum_8_4_U_org[i][j]+=sum_4_U_org[i][(j<<1)+1];

			sum2_8_4_U_org[i][j]+=sum2_4_U_org[i][j<<1];
			sum2_8_4_U_org[i][j]+=sum2_4_U_org[i][(j<<1)+1];
		}
	
	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算4*8块   range
		for(j=0;j<width;j++)
		{
			sum_4_8_U_org[i][j]=0.0;
			sum2_4_8_U_org[i][j]=0.0;

			sum_4_8_U_org[i][j]+=sum_4_U_org[i<<1][j];
			sum_4_8_U_org[i][j]+=sum_4_U_org[(i<<1)+1][j];

			sum2_4_8_U_org[i][j]+=sum2_4_U_org[i<<1][j];
			sum2_4_8_U_org[i][j]+=sum2_4_U_org[(i<<1)+1][j];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算8*8块   range
		for(j=0;j<width;j++)
		{
			sum_8_U_org[i][j]=0.0;
			sum2_8_U_org[i][j]=0.0;

			sum_8_U_org[i][j]+=sum_8_4_U_org[i<<1][j];
			sum_8_U_org[i][j]+=sum_8_4_U_org[(i<<1)+1][j];

			sum2_8_U_org[i][j]+=sum2_8_4_U_org[i<<1][j];
			sum2_8_U_org[i][j]+=sum2_8_4_U_org[(i<<1)+1][j];
		}

	width = img->frmWidthInMbs/2;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算16*8块   range
		for(j=0;j<width;j++)
		{
			sum_16_8_U_org[i][j]=0.0;
			sum2_16_8_U_org[i][j]=0.0;

			sum_16_8_U_org[i][j]+=sum_8_U_org[i][j<<1];
			sum_16_8_U_org[i][j]+=sum_8_U_org[i][(j<<1)+1];

			sum2_16_8_U_org[i][j]+=sum2_8_U_org[i][j<<1];
			sum2_16_8_U_org[i][j]+=sum2_8_U_org[i][(j<<1)+1];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs/2;		
	for(i=0;i<height;i++)       //计算8*16块   range
		for(j=0;j<width;j++)
		{
			sum_8_16_U_org[i][j]=0.0;
			sum2_8_16_U_org[i][j]=0.0;

			sum_8_16_U_org[i][j]+=sum_8_U_org[i<<1][j];
			sum_8_16_U_org[i][j]+=sum_8_U_org[(i<<1)+1][j];

			sum2_8_16_U_org[i][j]+=sum2_8_U_org[i<<1][j];
			sum2_8_16_U_org[i][j]+=sum2_8_U_org[(i<<1)+1][j];
		}
		
	width = img->frmWidthInMbs/2;
	height = img->frmHeightInMbs/2;		
	for(i=0;i<height;i++)       //计算16*16块   range
		for(j=0;j<width;j++)
		{
			sum_16_U_org[i][j]=0.0;
			sum2_16_U_org[i][j]=0.0;

			sum_16_U_org[i][j]+=sum_16_8_U_org[i<<1][j];
			sum_16_U_org[i][j]+=sum_16_8_U_org[(i<<1)+1][j];

			sum2_16_U_org[i][j]+=sum2_16_8_U_org[i<<1][j];
			sum2_16_U_org[i][j]+=sum2_16_8_U_org[(i<<1)+1][j];
		}
		//////////////////////////////////////////////////////////////////////////计算V分量	

	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算4*4块
		for(j=0;j<width;j++)
		{
			sum_4_V_org[i][j]=0.0;
			sum2_4_V_org[i][j]=0.0;

			sum_4_V_org[i][j]+=imgUV_org[1][i*4][j*4];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4][j*4+1];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4][j*4+2];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4][j*4+3];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4+1];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4+2];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4+3];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4+1];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4+2];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4+3];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4+1];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4+2];
			sum_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4+3];

			sum2_4_V_org[i][j]+=imgUV_org[1][i*4][j*4]*imgUV_org[1][i*4][j*4];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4][j*4+1]*imgUV_org[1][i*4][j*4+1];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4][j*4+2]*imgUV_org[1][i*4][j*4+2];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4][j*4+3]*imgUV_org[1][i*4][j*4+3];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4]*imgUV_org[1][i*4+1][j*4];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4+1]*imgUV_org[1][i*4+1][j*4+1];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4+2]*imgUV_org[1][i*4+1][j*4+2];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+1][j*4+3]*imgUV_org[1][i*4+1][j*4+3];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4]*imgUV_org[1][i*4+2][j*4];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4+1]*imgUV_org[1][i*4+2][j*4+1];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4+2]*imgUV_org[1][i*4+2][j*4+2];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+2][j*4+3]*imgUV_org[1][i*4+2][j*4+3];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4]*imgUV_org[1][i*4+3][j*4];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4+1]*imgUV_org[1][i*4+3][j*4+1];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4+2]*imgUV_org[1][i*4+3][j*4+2];
			sum2_4_V_org[i][j]+=imgUV_org[1][i*4+3][j*4+3]*imgUV_org[1][i*4+3][j*4+3];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs*2;		
	for(i=0;i<height;i++)       //计算8*4块   range
		for(j=0;j<width;j++)
		{
			sum_8_4_V_org[i][j]=0.0;
			sum2_8_4_V_org[i][j]=0.0;

			sum_8_4_V_org[i][j]+=sum_4_V_org[i][j<<1];
			sum_8_4_V_org[i][j]+=sum_4_V_org[i][(j<<1)+1];

			sum2_8_4_V_org[i][j]+=sum2_4_V_org[i][j<<1];
			sum2_8_4_V_org[i][j]+=sum2_4_V_org[i][(j<<1)+1];
		}
	
	width = img->frmWidthInMbs*2;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算4*8块   range
		for(j=0;j<width;j++)
		{
			sum_4_8_V_org[i][j]=0.0;
			sum2_4_8_V_org[i][j]=0.0;

			sum_4_8_V_org[i][j]+=sum_4_V_org[i<<1][j];
			sum_4_8_V_org[i][j]+=sum_4_V_org[(i<<1)+1][j];

			sum2_4_8_V_org[i][j]+=sum2_4_V_org[i<<1][j];
			sum2_4_8_V_org[i][j]+=sum2_4_V_org[(i<<1)+1][j];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算8*8块   range
		for(j=0;j<width;j++)
		{
			sum_8_V_org[i][j]=0.0;
			sum2_8_V_org[i][j]=0.0;

			sum_8_V_org[i][j]+=sum_8_4_V_org[i<<1][j];
			sum_8_V_org[i][j]+=sum_8_4_V_org[(i<<1)+1][j];

			sum2_8_V_org[i][j]+=sum2_8_4_V_org[i<<1][j];
			sum2_8_V_org[i][j]+=sum2_8_4_V_org[(i<<1)+1][j];
		}

	width = img->frmWidthInMbs/2;
	height = img->frmHeightInMbs;		
	for(i=0;i<height;i++)       //计算16*8块   range
		for(j=0;j<width;j++)
		{
			sum_16_8_V_org[i][j]=0.0;
			sum2_16_8_V_org[i][j]=0.0;

			sum_16_8_V_org[i][j]+=sum_8_V_org[i][j<<1];
			sum_16_8_V_org[i][j]+=sum_8_V_org[i][(j<<1)+1];

			sum2_16_8_V_org[i][j]+=sum2_8_V_org[i][j<<1];
			sum2_16_8_V_org[i][j]+=sum2_8_V_org[i][(j<<1)+1];
		}

	width = img->frmWidthInMbs;
	height = img->frmHeightInMbs/2;		
	for(i=0;i<height;i++)       //计算8*16块   range
		for(j=0;j<width;j++)
		{
			sum_8_16_V_org[i][j]=0.0;
			sum2_8_16_V_org[i][j]=0.0;

			sum_8_16_V_org[i][j]+=sum_8_V_org[i<<1][j];
			sum_8_16_V_org[i][j]+=sum_8_V_org[(i<<1)+1][j];

			sum2_8_16_V_org[i][j]+=sum2_8_V_org[i<<1][j];
			sum2_8_16_V_org[i][j]+=sum2_8_V_org[(i<<1)+1][j];
		}
		
	width = img->frmWidthInMbs/2;
	height = img->frmHeightInMbs/2;		
	for(i=0;i<height;i++)       //计算16*16块   range
		for(j=0;j<width;j++)
		{
			sum_16_V_org[i][j]=0.0;
			sum2_16_V_org[i][j]=0.0;

			sum_16_V_org[i][j]+=sum_16_8_V_org[i<<1][j];
			sum_16_V_org[i][j]+=sum_16_8_V_org[(i<<1)+1][j];

			sum2_16_V_org[i][j]+=sum2_16_8_V_org[i<<1][j];
			sum2_16_V_org[i][j]+=sum2_16_8_V_org[(i<<1)+1][j];
		}
}