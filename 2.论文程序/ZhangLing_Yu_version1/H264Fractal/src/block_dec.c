#include "block_dec.h"
#include "i_global.h"

/////////////////////////////////////>
#include <stdlib.h>
#include "global.h"

#include "block.h"
#include "image.h"
#include "mb_access.h"
#include <i_defines.h>

#define Q_BITS          15



/////////////////////////////////////<


void decode_one_macroblock(int CurMb,TRANS_NODE **trans,int con)//解码一个宏块，con=1:Y分量；con=2:U分量；con=3:V分量
{
	int i,i1,i2,j1,j2,num_domain,k,temp;
	double sum_16;
	byte **img_ref,**img_rec,**plane,**plane_domain;
	int mb_X = CurMb%(img->frmWidthInMbs/min(2,con));//将当前宏块的光栅扫描坐标变为二维x坐标(以宏块为单位)
	int mb_Y = CurMb/(img->frmWidthInMbs/min(2,con));//将当前宏块的光栅扫描坐标变为二维y坐标(以宏块为单位)
	int b_x=mb_X*16;//当前宏块最左上点像素的x坐标（以像素为单位）
	int b_y=mb_Y*16;//当前宏块最左上点像素的y坐标（以像素为单位）
	double scale,offset,average_domain;
	int mv_x,mv_y,mode;

	for(i=0;i<num_regions;i++)//number_regions=1
	{
		if((num_regions==2&&trans[i][CurMb].use==1)||num_regions==1)//yes
		{
			if(0==trans[i][CurMb].partition)//16x16
			{
				scale=trans[i][CurMb].scale;//从trans中读取分形参数
				offset=trans[i][CurMb].offset;
				mv_x=trans[i][CurMb].x;
				mv_y=trans[i][CurMb].y;
	
				if(trans[i][CurMb].region==1&&num_regions==2)   //是边界块  //no
				{
					switch(con)
					{
					case 1:
						if (currentVideo!='C'&&trans[i][CurMb].reference==1)//参考块是C目的块
						{
							plane_domain=plane_Y_domain;
							img_ref=imgY_ref;
						}
						else
						{
							plane_domain=plane_Y_domain_temp;
							img_ref=imgY_ref_temp;
						}
						plane=plane_Y;
						img_rec=imgY_rec_region[i];break;
					case 2:
						if (currentVideo!='C'&&trans[i][CurMb].reference==1)
						{
							plane_domain=plane_UV_domain[0];
							img_ref=imgUV_ref[0];
						}
						else
						{
							plane_domain=plane_UV_domain_temp[0];
							img_ref=imgUV_ref_temp[0];
						}
						plane=plane_UV[0];
						img_rec=imgUV_rec_region[0][i];break;
					case 3:
						if (currentVideo!='C'&&trans[i][CurMb].reference==1)
						{
							plane_domain=plane_UV_domain[1];
							img_ref=imgUV_ref[1];
						}
						else
						{
							plane_domain=plane_UV_domain_temp[1];
							img_ref=imgUV_ref_temp[1];
						}
						plane=plane_UV[1];
						img_rec=imgUV_rec_region[1][i];break;
					}

					num_domain=0;
					dsum1=0.0;
					for (j1=b_y,j2=mv_y+b_y; j1<b_y+16; ++j1,j2++)
						for (i1=b_x,i2=mv_x+b_x; i1<b_x+16; ++i1,i2++)
						{
							if(plane[j1][i1]/GREY_LEVELS==i)
							{				
								if(plane_domain[j2][i2]/GREY_LEVELS==i)
								{
									dsum1+=img_ref[j2][i2];
									num_domain++;
								}
							}
						}
						
					if(num_domain>0)
						average_domain=(unsigned char)(dsum1/num_domain);
					else
						average_domain=0;
					//分形解码图像	
					for (j1=b_y,j2=mv_y+b_y; j1<b_y+16; ++j1,j2++)
						for (i1=b_x,i2=mv_x+b_x; i1<b_x+16; ++i1,i2++)
						{
							if(plane[j1][i1]/GREY_LEVELS==i)
							{
								if(plane_domain[j2][i2]/GREY_LEVELS==i)
								    img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*img_ref[j2][i2] 
								              + offset - scale*average_domain);
 								else
 									img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*average_domain 
 								              + offset - scale*average_domain);
							}
						}
				}
				else//非边界块 //yes
				{
					switch(con)
					{
					case 1://Y分量
						if (trans[i][CurMb].reference==0)  //参考C目
						{
							img_ref=imgY_ref;
							sum_16=sum_16_ref[b_y+mv_y][b_x+mv_x];//参考的16*16 domain块的像素和
						}
						else if(trans[i][CurMb].reference==1)//参考H目
						{
							img_ref=imgY_ref_h;
							sum_16=sum_16_ref_H[b_y+mv_y][b_x+mv_x];//参考的16*16 domain块的像素和
						}
                        else if (trans[i][CurMb].reference==2)
                        {
							img_ref=imgY_ref_m;
							sum_16=sum_16_ref_M[b_y+mv_y][b_x+mv_x];//参考的16*16 domain块的像素和
                        }
						else if (trans[i][CurMb].reference==3)
						{
							img_ref=imgY_ref_n;
							sum_16=sum_16_ref_N[b_y+mv_y][b_x+mv_x];//参考的16*16 domain块的像素和
						}
						else
						{
							img_ref=imgY_ref_temp;
							sum_16=sum_16_ref_temp[b_y+mv_y][b_x+mv_x];//参考的16*16 domain块的像素和
						}
						if(num_regions>1)//no
							img_rec=imgY_rec_region[i];
						else//yes
							img_rec=imgY_rec;
						break;
					case 2://U分量
						if (trans[i][CurMb].reference==0)
						{
							img_ref=imgUV_ref[0];
							sum_16=sum_16_U_ref[b_y+mv_y][b_x+mv_x];
						}
						else if(trans[i][CurMb].reference==1)
						{
							img_ref=imgUV_ref_h[0];
							sum_16=sum_16_U_ref_H[b_y+mv_y][b_x+mv_x];
						} else if (trans[i][CurMb].reference==2)
						{
							img_ref=imgUV_ref_m[0];
							sum_16=sum_16_U_ref_M[b_y+mv_y][b_x+mv_x];
						} else if (trans[i][CurMb].reference==3)
						{
							img_ref=imgUV_ref_n[0];
							sum_16=sum_16_U_ref_N[b_y+mv_y][b_x+mv_x];
						}
						else
						{
							img_ref=imgUV_ref_temp[0];
							sum_16=sum_16_U_ref_temp[b_y+mv_y][b_x+mv_x];
						}

						if(num_regions>1)
							img_rec=imgUV_rec_region[0][i];
						else
							img_rec=imgUV_rec[0];
						break;
					case 3://V分量
						if (trans[i][CurMb].reference==0)
						{
							img_ref=imgUV_ref[1];
							sum_16=sum_16_V_ref[b_y+mv_y][b_x+mv_x];
						}
						else if(trans[i][CurMb].reference==1)
						{
							img_ref=imgUV_ref_h[1];
							sum_16=sum_16_V_ref_H[b_y+mv_y][b_x+mv_x];
						}else if (trans[i][CurMb].reference==2)
						{
							img_ref=imgUV_ref_m[1];
							sum_16=sum_16_V_ref_M[b_y+mv_y][b_x+mv_x];
						}else if (trans[i][CurMb].reference==3)
						{
							img_ref=imgUV_ref_n[1];
							sum_16=sum_16_V_ref_N[b_y+mv_y][b_x+mv_x];
						}
						else
						{
							img_ref=imgUV_ref_temp[1];
							sum_16=sum_16_V_ref_temp[b_y+mv_y][b_x+mv_x];
						}

						if(num_regions>1)
							img_rec=imgUV_rec_region[1][i];
						else
							img_rec=imgUV_rec[1];
						break;
					}
					average_domain = sum_16/(double)(16*16);//参考的16*16 domain块的像素均值
//////////////////////////解码公式
						for(j1=b_y,j2=b_y+mv_y;j1<b_y+16;j1++,j2++)//（i1,j1）是当前16*16 range块的像素坐标，（i2,j2）是相应的16*16 domain块的像素坐标
							for(i1=b_x,i2=b_x+mv_x;i1<b_x+16;i1++,i2++)
							{
// 								temp=img_ref[j2][i2];
// 								for (k=0;k<6;k++)
// 								{
// 								img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*temp// img_ref[j2][i2]
// 									+ offset - scale*average_domain);
// 								temp=img_rec[j1][i1];
// 								}

								img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*img_ref[j2][i2] //重建每一个像素值
									+ offset- scale*average_domain);//为什么还要减去scale*average_domain呢？？？

							}
					}
			}
			else//对块进行划分
				if(trans[i][CurMb].partition==1||trans[i][CurMb].partition==2)//16*8或8*16划分
				{
					mode=trans[i][CurMb].partition;
					for(i1=0;i1<2;i1++)//一个宏块包含两个16*8或8*16划分的块
					{
						decode_block_rect(b_x,b_y,i1,trans[i][CurMb].next[i1],con,mode,i,1);//解码一个16*8或8*16划分的块
					}
				}
				else//8*8划分

				{
					for(i1=0;i1<2;i1++)
						for(i2=0;i2<2;i2++)//一个宏块包含四个8*8划分的块
						{
							decode_block_8(b_x+i2*8,b_y+i1*8,i1*2+i2,trans[i][CurMb].next[i1*2+i2],con,i);//解码一个8*8划分的块
						}
				}
		}
		//不处理
		else//no
		{
			switch(con)
			{
			case 1:
				if(num_regions>1)
					img_rec=imgY_rec_region[i];
				else
					img_rec=imgY_rec;break;
			case 2:
				if(num_regions>1)
					img_rec=imgUV_rec_region[0][i];
				else
					img_rec=imgUV_rec[0];break;
			case 3:
				if(num_regions>1)
					img_rec=imgUV_rec_region[1][i];
				else
					img_rec=imgUV_rec[1];break;
			}

			for(j1=b_y;j1<b_y+16;j1++)
				for(i1=b_x;i1<b_x+16;i1++)
					img_rec[j1][i1]=0;
		}
	}
}

void decode_block_rect(int b_x,int b_y,int Curblock_8,TRANS_NODE trans,int con,int mode,int reg,int depth)
{
	int i1,i2,j1,j2,num_domain,k,temp;

	double sum_rect;
	byte **img_ref,**img_rec,**plane,**plane_domain;

	double scale,offset,average_domain;
	int mv_x,mv_y;
	int block_size_x,block_size_y;

	switch(mode)
	{
	case 1:
		block_size_x=16/depth;block_size_y=8/depth;
		b_y=b_y+block_size_y*Curblock_8;break;
	case 2:
		block_size_x=8/depth;block_size_y=16/depth;
		b_x=b_x+block_size_x*Curblock_8;break;
	}

	if((num_regions==2&&trans.use==1)||num_regions==1)
	{
		scale=trans.scale;
		offset=trans.offset;
		mv_x=trans.x;
		mv_y=trans.y;
		if(trans.region==1&&num_regions==2)   //是边界块
		{
			switch(con)
			{
			case 1:
				if (currentVideo!='C'&&trans.reference==1)//参考块是C目的块
				{
					plane_domain=plane_Y_domain;
					img_ref=imgY_ref;
				}
				else
				{
					plane_domain=plane_Y_domain_temp;
					img_ref=imgY_ref_temp;
				}
				plane=plane_Y;
				img_rec=imgY_rec_region[reg];break;
			case 2:
				if (currentVideo!='C'&&trans.reference==1)
				{
					plane_domain=plane_UV_domain[0];
					img_ref=imgUV_ref[0];
				}
				else
				{
					plane_domain=plane_UV_domain_temp[0];
					img_ref=imgUV_ref_temp[0];
				}
				plane=plane_UV[0];
				img_rec=imgUV_rec_region[0][reg];break;
			case 3:
				if (currentVideo!='C'&&trans.reference==1)
				{
					plane_domain=plane_UV_domain[1];
					img_ref=imgUV_ref[1];
				}
				else
				{
					plane_domain=plane_UV_domain_temp[1];
					img_ref=imgUV_ref_temp[1];
				}
				plane=plane_UV[1];
				img_rec=imgUV_rec_region[1][reg];break;
			}
				
			num_domain=0;
			dsum1=0.0;
			for (j1=b_y,j2=mv_y+b_y; j1<b_y+block_size_y; ++j1,j2++)
				for (i1=b_x,i2=mv_x+b_x; i1<b_x+block_size_x; ++i1,i2++)
				{
					if(plane[j1][i1]/GREY_LEVELS==reg)
					{				
						if(plane_domain[j2][i2]/GREY_LEVELS==reg)
						{
							dsum1+=img_ref[j2][i2];
							num_domain++;
						}
					}
				}

			if(num_domain>0)
				average_domain=(unsigned char)(dsum1/num_domain);
			else
				average_domain=0;
					
			for (j1=b_y,j2=mv_y+b_y; j1<b_y+block_size_y; ++j1,j2++)
				for (i1=b_x,i2=mv_x+b_x; i1<b_x+block_size_x; ++i1,i2++)
				{
					if(plane[j1][i1]/GREY_LEVELS==reg)
					{
						if(plane_domain[j2][i2]/GREY_LEVELS==reg)
							img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*img_ref[j2][i2] 
							+ offset - scale*average_domain);
						else
							img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*average_domain 
							+ offset - scale*average_domain);
					}
				}
		}

		else
		{
			if(block_size_x==16&&block_size_y==8)
			{
				switch(con)
				{
				case 1:
					if (trans.reference==0)
					{
						img_ref=imgY_ref;
						sum_rect=sum_16_8_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgY_ref_h;
						sum_rect=sum_16_8_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgY_ref_m;
						sum_rect=sum_16_8_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgY_ref_n;
						sum_rect=sum_16_8_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgY_rec_region[reg];
					else
						img_rec=imgY_rec;
					break;
				case 2:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[0];
						sum_rect=sum_16_8_U_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[0];
						sum_rect=sum_16_8_U_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgUV_ref_m[0];
						sum_rect=sum_16_8_U_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[0];
						sum_rect=sum_16_8_U_ref_M[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[0][reg];
					else
						img_rec=imgUV_rec[0];
					break;
				case 3:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[1];
						sum_rect=sum_16_8_V_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[1];
						sum_rect=sum_16_8_V_ref_H[b_y+mv_y][b_x+mv_x];
					} else if (trans.reference==2)
					{
						img_ref=imgUV_ref_n[1];
						sum_rect=sum_16_8_V_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[1][reg];
					else
						img_rec=imgUV_rec[1];
					break;
				}
			}

			if(block_size_x==8&&block_size_y==16)
			{
				switch(con)
				{
				case 1:
					if (trans.reference==0)
					{
						img_ref=imgY_ref;
						sum_rect=sum_8_16_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgY_ref_h;
						sum_rect=sum_8_16_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgY_ref_m;
						sum_rect=sum_8_16_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgY_ref_n;
						sum_rect=sum_8_16_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
					    img_rec=imgY_rec_region[reg];
					else
						img_rec=imgY_rec;
						break;
				case 2:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[0];
						sum_rect=sum_8_16_U_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[0];
						sum_rect=sum_8_16_U_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgUV_ref_m[0];
						sum_rect=sum_8_16_U_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[0];
						sum_rect=sum_8_16_U_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[0][reg];
					else
						img_rec=imgUV_rec[0];
					break;
				case 3:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[1];
						sum_rect=sum_8_16_V_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[1];
						sum_rect=sum_8_16_V_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgUV_ref_m[1];
						sum_rect=sum_8_16_V_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[1];
						sum_rect=sum_8_16_V_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[1][reg];
					else
						img_rec=imgUV_rec[1];
					break;
				}
			}

			if(block_size_x==8&&block_size_y==4)
			{
				switch(con)
				{
				case 1:
					if (trans.reference==0)
					{
						img_ref=imgY_ref;
						sum_rect=sum_8_4_ref[b_y+mv_y][b_x+mv_x];
					}
					else if (trans.reference==1)
					{
						img_ref=imgY_ref_h;
						sum_rect=sum_8_4_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgY_ref_m;
						sum_rect=sum_8_4_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgY_ref_n;
						sum_rect=sum_8_4_ref_N[b_y+mv_y][b_x+mv_x];
					}
					if(num_regions>1)
						img_rec=imgY_rec_region[reg];
					else
						img_rec=imgY_rec;
					break;
				case 2:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[0];
						sum_rect=sum_8_4_U_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[0];
						sum_rect=sum_8_4_U_ref_H[b_y+mv_y][b_x+mv_x];
					}else if(trans.reference==2)
					{
						img_ref=imgUV_ref_m[0];
						sum_rect=sum_8_4_U_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[0];
						sum_rect=sum_8_4_U_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[0][reg];
					else
						img_rec=imgUV_rec[0];break;
				case 3:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[1];
						sum_rect=sum_8_4_V_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[1];
						sum_rect=sum_8_4_V_ref_H[b_y+mv_y][b_x+mv_x];
					}else if(trans.reference==2)
					{
						img_ref=imgUV_ref_m[1];
						sum_rect=sum_8_4_V_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[1];
						sum_rect=sum_8_4_V_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[1][reg];
					else
						img_rec=imgUV_rec[1];break;
				}
			}

			if(block_size_x==4&&block_size_y==8)
			{
				switch(con)
				{
				case 1:
					if (trans.reference==0)
					{
						img_ref=imgY_ref;
						sum_rect=sum_4_8_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgY_ref_h;
						sum_rect=sum_4_8_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgY_ref_m;
						sum_rect=sum_4_8_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgY_ref_n;
						sum_rect=sum_4_8_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgY_rec_region[reg];
					else
						img_rec=imgY_rec;break;
				case 2:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[0];
						sum_rect=sum_4_8_U_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[0];
						sum_rect=sum_4_8_U_ref_H[b_y+mv_y][b_x+mv_x];
					}else if (trans.reference==2)
					{
						img_ref=imgUV_ref_m[0];
						sum_rect=sum_4_8_U_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[0];
						sum_rect=sum_4_8_U_ref_N[b_y+mv_y][b_x+mv_x];
					}
					if(num_regions>1)
						img_rec=imgUV_rec_region[0][reg];
					else
						img_rec=imgUV_rec[0];break;
				case 3:
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[1];
						sum_rect=sum_4_8_V_ref[b_y+mv_y][b_x+mv_x];
					}
					else if(trans.reference==1)
					{
						img_ref=imgUV_ref_h[1];
						sum_rect=sum_4_8_V_ref_H[b_y+mv_y][b_x+mv_x];
					} else if (trans.reference==2)
					{
						img_ref=imgUV_ref_m[1];
						sum_rect=sum_4_8_V_ref_M[b_y+mv_y][b_x+mv_x];
					}else
					{
						img_ref=imgUV_ref_n[1];
						sum_rect=sum_4_8_V_ref_N[b_y+mv_y][b_x+mv_x];
					}

					if(num_regions>1)
						img_rec=imgUV_rec_region[1][reg];
					else
						img_rec=imgUV_rec[1];break;
				}
			}

			average_domain = sum_rect/(double)(block_size_y*block_size_x);
//////////////////////////解码增加迭代次数

				for(i1=b_y,i2=b_y+mv_y;i1<b_y+block_size_y;i1++,i2++)
					for(j1=b_x,j2=b_x+mv_x;j1<b_x+block_size_x;j1++,j2++)
					{
// 						temp=img_ref[i2][j2];
// 						for (k=0;k<6;k++)
// 						{
// 				       		img_rec[i1][j1]=(unsigned char) bound(0.5 + scale*temp//img_ref[i2][j2] 
// 							+ offset- scale*average_domain); 
//                             temp=img_rec[i1][j1];
// 						}

						img_rec[i1][j1]=(unsigned char) bound(0.5 + scale*img_ref[i2][j2] 
 							+ offset - scale*average_domain);
					}
		}
	}
	
	else
	{
		switch(con)
		{
		case 1:
			if(num_regions>1)
				img_rec=imgY_rec_region[reg];
			else
				img_rec=imgY_rec;
			break;
		case 2:
			if(num_regions>1)
				img_rec=imgUV_rec_region[0][reg];
			else
				img_rec=imgUV_rec[0];
			break;
		case 3:
			if(num_regions>1)
				img_rec=imgUV_rec_region[1][reg];
			else
				img_rec=imgUV_rec[1];
			break;
		}
		for(j1=b_y;j1<b_y+block_size_y;j1++)
			for(i1=b_x;i1<b_x+block_size_x;i1++)
				img_rec[j1][i1]=0;
	}
}

void decode_block_8(int b_x,int b_y,int Curblock_8,TRANS_NODE trans,int con,int reg)
{
	int i1,i2,j1,j2,num_domain,k,temp;

	double sum_8;
	byte **img_ref,**img_rec,**plane,**plane_domain;

	double scale,offset,average_domain;
	int mv_x,mv_y,mode;
	if((num_regions==2&&trans.use==1)||num_regions==1)//yes
	{
		if(0==trans.partition)//8*8划分
		{
			scale=trans.scale;
			offset=trans.offset;
			mv_x=trans.x;
			mv_y=trans.y;

			if(trans.region==1&&num_regions==2)   //是边界块 //no
			{
				switch(con)
				{
				case 1:
					if (currentVideo!='C'&&trans.reference==1)//参考块是C目的块
					{
						plane_domain=plane_Y_domain;
						img_ref=imgY_ref;
					}
					else
					{
						plane_domain=plane_Y_domain_temp;
						img_ref=imgY_ref_temp;
					}
					plane=plane_Y;
					img_rec=imgY_rec_region[reg];break;
				case 2:
					if (currentVideo!='C'&&trans.reference==1)
					{
						plane_domain=plane_UV_domain[0];
						img_ref=imgUV_ref[0];
					}
					else
					{
						plane_domain=plane_UV_domain_temp[0];
						img_ref=imgUV_ref_temp[0];
					}
					plane=plane_UV[0];
					img_rec=imgUV_rec_region[0][reg];break;
				case 3:
					if (currentVideo!='C'&&trans.reference==1)
					{
						plane_domain=plane_UV_domain[1];
						img_ref=imgUV_ref[1];
					}
					else
					{
						plane_domain=plane_UV_domain_temp[1];
						img_ref=imgUV_ref_temp[1];
					}
					plane=plane_UV[1];
					img_rec=imgUV_rec_region[1][reg];break;
				}
				
				num_domain=0;
				dsum1=0.0;
				for (j1=b_y,j2=mv_y+b_y; j1<b_y+8; ++j1,j2++)
					for (i1=b_x,i2=mv_x+b_x; i1<b_x+8; ++i1,i2++)
					{
						if(plane[j1][i1]/GREY_LEVELS==reg)
						{				
							if(plane_domain[j2][i2]/GREY_LEVELS==reg)
							{
								dsum1+=img_ref[j2][i2];
								num_domain++;
							}
						}
					}

				if(num_domain>0)
					average_domain=(unsigned char)(dsum1/num_domain);
				else
					average_domain=0;
					
				for (j1=b_y,j2=mv_y+b_y; j1<b_y+8; ++j1,j2++)
					for (i1=b_x,i2=mv_x+b_x; i1<b_x+8; ++i1,i2++)
					{
						if(plane[j1][i1]/GREY_LEVELS==reg)
						{
							if(plane_domain[j2][i2]/GREY_LEVELS==reg)
						    	img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*img_ref[j2][i2] 
									+ offset - scale*average_domain);
 							else
 								img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*average_domain 
 									+ offset - scale*average_domain);
						}
					}
			}
			else//非边界块 //yes
			{
				switch(con)
				{
				case 1://Y分量
					if (trans.reference==0)//参考C目
					{
						img_ref=imgY_ref;
						sum_8=sum_8_ref[b_y+mv_y][b_x+mv_x];//参考的8*8 domain块的像素和
					}
					else//参考H目
					{
						img_ref=imgY_ref_h;
						sum_8=sum_8_ref_H[b_y+mv_y][b_x+mv_x];
					}
					
					if(num_regions>1)//no
						img_rec=imgY_rec_region[reg];
					else
						img_rec=imgY_rec;
					break;
				case 2://U分量
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[0];
						sum_8=sum_8_U_ref[b_y+mv_y][b_x+mv_x];
					}
					else
					{
						img_ref=imgUV_ref_h[0];
						sum_8=sum_8_U_ref_H[b_y+mv_y][b_x+mv_x];
					}
					
					if(num_regions>1)
						img_rec=imgUV_rec_region[0][reg];
					else
						img_rec=imgUV_rec[0];
					break;
				case 3://V分量
					if (trans.reference==0)
					{
						img_ref=imgUV_ref[1];
						sum_8=sum_8_V_ref[b_y+mv_y][b_x+mv_x];
					}
					else
					{
						img_ref=imgUV_ref_h[1];
						sum_8=sum_8_V_ref_H[b_y+mv_y][b_x+mv_x];
					}
					if(num_regions>1)
						img_rec=imgUV_rec_region[1][reg];
					else
						img_rec=imgUV_rec[1];
					break;
				}
				average_domain = sum_8/(double)(8*8);//参考的8*8 domain块的像素平均值
//////////////////////////解码一个8*8块
					for(i1=b_y,i2=b_y+mv_y;i1<b_y+8;i1++,i2++)
						for(j1=b_x,j2=b_x+mv_x;j1<b_x+8;j1++,j2++)//（j1,i1）是当前8*8 range块坐标，（j2,i2）是相应的8*8 domain块坐标
						{
// 							temp=img_ref[i2][j2];
// 							for (k=0;k<6;k++)
// 							{
// 							img_rec[i1][j1]=(unsigned char) bound(0.5 + scale*temp//img_ref[i2][j2] 
// 								+ offset - scale*average_domain);
// 							temp=img_rec[i1][j1];
// 							}

							img_rec[i1][j1]=(unsigned char) bound(0.5 + scale*img_ref[i2][j2] 
								+ offset - scale*average_domain);//分形解码每一个像素

						}
			}
		}
		else//8*8继续划分
			if(trans.partition==1||trans.partition==2)//8*4或4*8划分
			{
				mode=trans.partition;
				for(i1=0;i1<2;i1++)
					{
						decode_block_rect(b_x,b_y,i1,trans.next[i1],con,mode,reg,2);
					}
			}
			else//4*4划分
			{
				for(i1=0;i1<2;i1++)
					for(i2=0;i2<2;i2++)//一个8*8块包含4个4*4块
					{
						decode_block_4(b_x+i2*4,b_y+i1*4,trans.next[i1*2+i2],con,reg);
					}
			}
	}
	else//no
	{
		switch(con)
		{
		case 1:
			if(num_regions>1)
				img_rec=imgY_rec_region[reg];
			else
				img_rec=imgY_rec;
			break;
		case 2:
			if(num_regions>1)
				img_rec=imgUV_rec_region[0][reg];
			else
				img_rec=imgUV_rec[0];
			break;
		case 3:
			if(num_regions>1)
				img_rec=imgUV_rec_region[1][reg];
			else
				img_rec=imgUV_rec[1];
			break;
		}
		for(j1=b_y;j1<b_y+8;j1++)
			for(i1=b_x;i1<b_x+8;i1++)
				img_rec[j1][i1]=0;
	}
}

void decode_block_4(int b_x,int b_y,TRANS_NODE trans,int con,int reg)
{
	int i1,i2,j1,j2,num_domain=0,k,temp;

	double sum_4;
	byte **img_ref,**img_rec,**plane,**plane_domain;

	double scale,offset,average_domain;
	int mv_x,mv_y;
	
	if((num_regions==2&&trans.use==1)||num_regions==1)//yes
	{
		scale=trans.scale;
		offset=trans.offset;
		mv_x=trans.x;
		mv_y=trans.y;

		if(trans.region==1&&num_regions==2)   //是边界块 //no
		{
			switch(con)
			{
			case 1:
				if (currentVideo!='C'&&trans.reference==1)//参考块是C目的块
				{
					plane_domain=plane_Y_domain;
					img_ref=imgY_ref;
				}
				else
				{
					plane_domain=plane_Y_domain_temp;
					img_ref=imgY_ref_temp;
				}
				plane=plane_Y;
				img_rec=imgY_rec_region[reg];break;
			case 2:
				if (currentVideo!='C'&&trans.reference==1)
				{
					plane_domain=plane_UV_domain[0];
					img_ref=imgUV_ref[0];
				}
				else
				{
					plane_domain=plane_UV_domain_temp[0];
					img_ref=imgUV_ref_temp[0];
				}
				plane=plane_UV[0];
				img_rec=imgUV_rec_region[0][reg];break;
			case 3:
				if (currentVideo!='C'&&trans.reference==1)
				{
					plane_domain=plane_UV_domain[1];
					img_ref=imgUV_ref[1];
				}
				else
				{
					plane_domain=plane_UV_domain_temp[1];
					img_ref=imgUV_ref_temp[1];
				}
				plane=plane_UV[1];
				img_rec=imgUV_rec_region[1][reg];break;
			}
			
			num_domain=0;
			dsum1=0.0;
			for (j1=b_y,j2=mv_y+b_y; j1<b_y+4; ++j1,j2++)
				for (i1=b_x,i2=mv_x+b_x; i1<b_x+4; ++i1,i2++)
				{
					if(plane[j1][i1]/GREY_LEVELS==reg)
					{				
						if(plane_domain[j2][i2]/GREY_LEVELS==reg)
						{
							dsum1+=img_ref[j2][i2];
							num_domain++;
						}
					}
				}
					
			if(num_domain>0)
				average_domain=(unsigned char)(dsum1/num_domain);
			else
				average_domain=0;
					
			for (j1=b_y,j2=mv_y+b_y; j1<b_y+4; ++j1,j2++)
				for (i1=b_x,i2=mv_x+b_x; i1<b_x+4; ++i1,i2++)
				{
					if(plane[j1][i1]/GREY_LEVELS==reg)
					{
						if(plane_domain[j2][i2]/GREY_LEVELS==reg)
						   	img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*img_ref[j2][i2] 
								+ offset - scale*average_domain);
 						else
							img_rec[j1][i1]=(unsigned char) bound(0.5 + scale*average_domain 
 								+ offset - scale*average_domain);
					}
				}
		}
		else//非边界块 //yes
		{
			switch(con)//Y分量
			{
			case 1:
				if (trans.reference==0)//参考C目
				{
					img_ref=imgY_ref;
					sum_4=sum_4_ref[b_y+mv_y][b_x+mv_x];//参考的4*4 domain块的像素和
				}
				else if(trans.reference==1)//参考H目
				{
					img_ref=imgY_ref_h;
					sum_4=sum_4_ref_H[b_y+mv_y][b_x+mv_x];//参考的4*4 domain块的像素和
				}else if (trans.reference==2)//参考M目
				{
					
					img_ref=imgY_ref_m;
					sum_4=sum_4_ref_M[b_y+mv_y][b_x+mv_x];//参考的4*4 domain块的像素和
				}else//参考N目
				{
					
					img_ref=imgY_ref_n;
					sum_4=sum_4_ref_N[b_y+mv_y][b_x+mv_x];//参考的4*4 domain块的像素和
				}
				
				if(num_regions>1)//no
					img_rec=imgY_rec_region[reg];
				else//yes
					img_rec=imgY_rec;
				break;
			case 2:
				if (trans.reference==0)
				{
					img_ref=imgUV_ref[0];
					sum_4=sum_4_U_ref[b_y+mv_y][b_x+mv_x];
				}
				else if(trans.reference==1)
				{
					img_ref=imgUV_ref_h[0];
					sum_4=sum_4_U_ref_H[b_y+mv_y][b_x+mv_x];
				}else if(trans.reference==2)
				{
					img_ref=imgUV_ref_m[0];
					sum_4=sum_4_U_ref_M[b_y+mv_y][b_x+mv_x];
				}else
				{
					img_ref=imgUV_ref_n[0];
					sum_4=sum_4_U_ref_N[b_y+mv_y][b_x+mv_x];
				}
				if(num_regions>1)
					img_rec=imgUV_rec_region[0][reg];
				else
					img_rec=imgUV_rec[0];
				break;
			case 3:
				if (trans.reference==0)
				{
					img_ref=imgUV_ref[1];
					sum_4=sum_4_V_ref[b_y+mv_y][b_x+mv_x];
				}
				else if(trans.reference==0)
				{
					img_ref=imgUV_ref_h[1];
					sum_4=sum_4_V_ref_H[b_y+mv_y][b_x+mv_x];
				}else if(trans.reference==2)
				{
					img_ref=imgUV_ref_m[1];
					sum_4=sum_4_V_ref_M[b_y+mv_y][b_x+mv_x];
				}else
				{
					img_ref=imgUV_ref_n[1];
					sum_4=sum_4_V_ref_N[b_y+mv_y][b_x+mv_x];
				}
				if(num_regions>1)
					img_rec=imgUV_rec_region[1][reg];
				else
					img_rec=imgUV_rec[1];
				break;
			}
			average_domain = sum_4/(double)(4*4);//参考的4*4 domain块的像素平均值
//////////////////////////解码一个4*4块
				for(i1=b_y,i2=b_y+mv_y;i1<b_y+4;i1++,i2++)
					for(j1=b_x,j2=b_x+mv_x;j1<b_x+4;j1++,j2++)
					{
// 						temp=img_ref[i2][j2];
// 						for (k=0; k<6 ;k++)
// 						{
// 						img_rec[i1][j1]=(unsigned char) bound(0.5 + scale*temp//img_ref[i2][j2] 
// 							+ offset - scale*average_domain);
// 						temp=img_rec[i1][j1];
// 						}

						img_rec[i1][j1]=(unsigned char) bound(0.5 + scale*img_ref[i2][j2] 
							+ offset - scale*average_domain);//分形解码每个像素

					}
		}
	}
	else//no
	{
		switch(con)
		{
		case 1:
			if(num_regions>1)
				img_rec=imgY_rec_region[reg];
			else
				img_rec=imgY_rec;
			break;
		case 2:
			if(num_regions>1)
				img_rec=imgUV_rec_region[0][reg];
			else
				img_rec=imgUV_rec[0];
			break;
		case 3:
			if(num_regions>1)
				img_rec=imgUV_rec_region[1][reg];
			else
				img_rec=imgUV_rec[1];
			break;
		}
		for(j1=b_y;j1<b_y+4;j1++)
			for(i1=b_x;i1<b_x+4;i1++)
				img_rec[j1][i1]=0;
	}
}



/////////////////////////////////////>




//先讨论下quant_coef如何得到地：(见相关文档)

static const int quant_coef[6][4][4] = {
	{{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
	{{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
	{{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
	{{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
	{{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
	{{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};



static const int A[4][4] = {
	{ 16, 20, 16, 20},
	{ 20, 25, 20, 25},
	{ 16, 20, 16, 20},
	{ 20, 25, 20, 25}
};

// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..L, and from above left X, as follows:
//
//  X A B C D E F G H
//  I a b c d
//  J e f g h
//  K i j k l
//  L m n o p
//


// Predictor array index definitions
#define P_X (PredPel[0])
#define P_A (PredPel[1])
#define P_B (PredPel[2])
#define P_C (PredPel[3])
#define P_D (PredPel[4])
#define P_E (PredPel[5])
#define P_F (PredPel[6])
#define P_G (PredPel[7])
#define P_H (PredPel[8])
#define P_I (PredPel[9])
#define P_J (PredPel[10])
#define P_K (PredPel[11])
#define P_L (PredPel[12])




/*!
***********************************************************************
* \brief
*    makes and returns 4x4 blocks with all 5 intra prediction modes
*
* \return
*    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
*    SEARCH_SYNC   search next sync element as errors while decoding occured
***********************************************************************
*/

int intrapred(
			  struct img_par *img,  //!< image parameters
			  int ioff,             //!< pixel offset X within MB
			  int joff,             //!< pixel offset Y within MB
			  int img_block_x,      //!< location of block X, multiples of 4
			  int img_block_y)      //!< location of block Y, multiples of 4
{
	int i,j;
	int s0;
	int img_y,img_x;
	int PredPel[13];  // array of predictor pels
	
	byte **imgY = dec_picture->imgY;
	
	PixelPos pix_a[4];
	PixelPos pix_b, pix_c, pix_d;
	
	int block_available_up;
	int block_available_left;
	int block_available_up_left;
	int block_available_up_right;
	
	int mb_nr=img->current_mb_nr;
	
	byte predmode = img->ipredmode[img_block_x][img_block_y];
	
	img_x=img_block_x*4;
	img_y=img_block_y*4;
	
	for (i=0;i<4;i++)
	{
		getNeighbour(mb_nr, ioff -1 , joff +i , 1, &pix_a[i]);
	}
	
	
	getNeighbour(mb_nr, ioff    , joff -1 , 1, &pix_b);
	getNeighbour(mb_nr, ioff +4 , joff -1 , 1, &pix_c);
	getNeighbour(mb_nr, ioff -1 , joff -1 , 1, &pix_d);
	
	pix_c.available = pix_c.available && !(((ioff==4)||(ioff==12)) && ((joff==4)||(joff==12)));
	
	if (active_pps->constrained_intra_pred_flag)
	{
		for (i=0, block_available_left=1; i<4;i++)
			block_available_left  &= pix_a[i].available ? img->intra_block[pix_a[i].mb_addr]: 0;
		block_available_up       = pix_b.available ? img->intra_block [pix_b.mb_addr] : 0;
		block_available_up_right = pix_c.available ? img->intra_block [pix_c.mb_addr] : 0;
		block_available_up_left  = pix_d.available ? img->intra_block [pix_d.mb_addr] : 0;
	}
	else
	{
		block_available_left     = pix_a[0].available;
		block_available_up       = pix_b.available;
		block_available_up_right = pix_c.available;
		block_available_up_left  = pix_d.available;
	}
	
	// form predictor pels
	if (block_available_up)
	{
		P_A = imgY[pix_b.pos_y][pix_b.pos_x+0];
		P_B = imgY[pix_b.pos_y][pix_b.pos_x+1];
		P_C = imgY[pix_b.pos_y][pix_b.pos_x+2];
		P_D = imgY[pix_b.pos_y][pix_b.pos_x+3];
		
	}
	else
	{
		P_A = P_B = P_C = P_D = 128;
	}
	
	if (block_available_up_right)
	{
		P_E = imgY[pix_c.pos_y][pix_c.pos_x+0];
		P_F = imgY[pix_c.pos_y][pix_c.pos_x+1];
		P_G = imgY[pix_c.pos_y][pix_c.pos_x+2];
		P_H = imgY[pix_c.pos_y][pix_c.pos_x+3];
	}
	else
	{
		P_E = P_F = P_G = P_H = P_D;
	}
	
	if (block_available_left)
	{
		P_I = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
		P_J = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
		P_K = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
		P_L = imgY[pix_a[3].pos_y][pix_a[3].pos_x];
	}
	else
	{
		P_I = P_J = P_K = P_L = 128;
	}
	
	if (block_available_up_left)
	{
		P_X = imgY[pix_d.pos_y][pix_d.pos_x];
	}
	else
	{
		P_X = 128;
	}
	
	
	switch (predmode)
	{
	case DC_PRED:     //2                    /* DC prediction */
// 		printf("DC_PRED://2    ");
		s0 = 0;
		if (block_available_up && block_available_left)
		{   
			// no edge
			s0 = (P_A + P_B + P_C + P_D + P_I + P_J + P_K + P_L + 4)/(2*BLOCK_SIZE);
		}
		else if (!block_available_up && block_available_left)
		{
			// upper edge
			s0 = (P_I + P_J + P_K + P_L + 2)/BLOCK_SIZE;             
		}
		else if (block_available_up && !block_available_left)
		{
			// left edge
			s0 = (P_A + P_B + P_C + P_D + 2)/BLOCK_SIZE;             
		}
		else //if (!block_available_up && !block_available_left)
		{
			// top left corner, nothing to predict from
			s0 = 128;                           
		}
		
		for (j=0; j < BLOCK_SIZE; j++)
		{
			for (i=0; i < BLOCK_SIZE; i++)
			{
				// store DC prediction
				img->mpr[i+ioff][j+joff] = s0;
			}
		}
		break;
		
	case VERT_PRED:     //0                  /* vertical prediction from block above */
// 		printf("VERT_PRED://0    ");
		if (!block_available_up)
			printf ("warning: Intra_4x4_Vertical prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		for(j=0;j<BLOCK_SIZE;j++)
			for(i=0;i<BLOCK_SIZE;i++)
				img->mpr[i+ioff][j+joff]=imgY[pix_b.pos_y][pix_b.pos_x+i];/* store predicted 4x4 block */
			break;
			
	case HOR_PRED:            //1            /* horizontal prediction from left block */
// 		printf("HOR_PRED://1    ");
		if (!block_available_left)
			printf ("warning: Intra_4x4_Horizontal prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		for(j=0;j<BLOCK_SIZE;j++)
			for(i=0;i<BLOCK_SIZE;i++)
				img->mpr[i+ioff][j+joff]=imgY[pix_a[j].pos_y][pix_a[j].pos_x]; /* store predicted 4x4 block */
			break;
			
	case DIAG_DOWN_RIGHT_PRED://4
// 		printf("DIAG_DOWN_RIGHT_PRED://4    ");
		if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
			printf ("warning: Intra_4x4_Diagonal_Down_Right prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		img->mpr[0+ioff][3+joff] = (P_L + 2*P_K + P_J + 2) / 4; 
		img->mpr[0+ioff][2+joff] =
			img->mpr[1+ioff][3+joff] = (P_K + 2*P_J + P_I + 2) / 4; 
		img->mpr[0+ioff][1+joff] =
			img->mpr[1+ioff][2+joff] = 
			img->mpr[2+ioff][3+joff] = (P_J + 2*P_I + P_X + 2) / 4; 
		img->mpr[0+ioff][0+joff] =
			img->mpr[1+ioff][1+joff] =
			img->mpr[2+ioff][2+joff] =
			img->mpr[3+ioff][3+joff] = (P_I + 2*P_X + P_A + 2) / 4; 
		img->mpr[1+ioff][0+joff] =
			img->mpr[2+ioff][1+joff] =
			img->mpr[3+ioff][2+joff] = (P_X + 2*P_A + P_B + 2) / 4;
		img->mpr[2+ioff][0+joff] =
			img->mpr[3+ioff][1+joff] = (P_A + 2*P_B + P_C + 2) / 4;
		img->mpr[3+ioff][0+joff] = (P_B + 2*P_C + P_D + 2) / 4;
		break;
		
	case DIAG_DOWN_LEFT_PRED://3
// 		printf("DIAG_DOWN_LEFT_PRED://3    ");
		if (!block_available_up)
			printf ("warning: Intra_4x4_Diagonal_Down_Left prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		img->mpr[0+ioff][0+joff] = (P_A + P_C + 2*(P_B) + 2) / 4;
		img->mpr[1+ioff][0+joff] = 
			img->mpr[0+ioff][1+joff] = (P_B + P_D + 2*(P_C) + 2) / 4;
		img->mpr[2+ioff][0+joff] =
			img->mpr[1+ioff][1+joff] =
			img->mpr[0+ioff][2+joff] = (P_C + P_E + 2*(P_D) + 2) / 4;
		img->mpr[3+ioff][0+joff] = 
			img->mpr[2+ioff][1+joff] = 
			img->mpr[1+ioff][2+joff] = 
			img->mpr[0+ioff][3+joff] = (P_D + P_F + 2*(P_E) + 2) / 4;
		img->mpr[3+ioff][1+joff] = 
			img->mpr[2+ioff][2+joff] = 
			img->mpr[1+ioff][3+joff] = (P_E + P_G + 2*(P_F) + 2) / 4;
		img->mpr[3+ioff][2+joff] = 
			img->mpr[2+ioff][3+joff] = (P_F + P_H + 2*(P_G) + 2) / 4;
		img->mpr[3+ioff][3+joff] = (P_G + 3*(P_H) + 2) / 4;
		break;
		
	case  VERT_RIGHT_PRED://5  /* diagonal prediction -22.5 deg to horizontal plane */
// 		printf("VERT_RIGHT_PRED://5    ");
		if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
			printf ("warning: Intra_4x4_Vertical_Right prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		img->mpr[0+ioff][0+joff] = 
			img->mpr[1+ioff][2+joff] = (P_X + P_A + 1) / 2;
		img->mpr[1+ioff][0+joff] = 
			img->mpr[2+ioff][2+joff] = (P_A + P_B + 1) / 2;
		img->mpr[2+ioff][0+joff] = 
			img->mpr[3+ioff][2+joff] = (P_B + P_C + 1) / 2;
		img->mpr[3+ioff][0+joff] = (P_C + P_D + 1) / 2;
		img->mpr[0+ioff][1+joff] = 
			img->mpr[1+ioff][3+joff] = (P_I + 2*P_X + P_A + 2) / 4;
		img->mpr[1+ioff][1+joff] = 
			img->mpr[2+ioff][3+joff] = (P_X + 2*P_A + P_B + 2) / 4;
		img->mpr[2+ioff][1+joff] = 
			img->mpr[3+ioff][3+joff] = (P_A + 2*P_B + P_C + 2) / 4;
		img->mpr[3+ioff][1+joff] = (P_B + 2*P_C + P_D + 2) / 4;
		img->mpr[0+ioff][2+joff] = (P_X + 2*P_I + P_J + 2) / 4;
		img->mpr[0+ioff][3+joff] = (P_I + 2*P_J + P_K + 2) / 4;
		break;
		
	case  VERT_LEFT_PRED://7  /* diagonal prediction -22.5 deg to horizontal plane */
// 		printf("VERT_LEFT_PRED://7    ");
		if (!block_available_up)
			printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		img->mpr[0+ioff][0+joff] = (P_A + P_B + 1) / 2;
		img->mpr[1+ioff][0+joff] = 
			img->mpr[0+ioff][2+joff] = (P_B + P_C + 1) / 2;
		img->mpr[2+ioff][0+joff] = 
			img->mpr[1+ioff][2+joff] = (P_C + P_D + 1) / 2;
		img->mpr[3+ioff][0+joff] = 
			img->mpr[2+ioff][2+joff] = (P_D + P_E + 1) / 2;
		img->mpr[3+ioff][2+joff] = (P_E + P_F + 1) / 2;
		img->mpr[0+ioff][1+joff] = (P_A + 2*P_B + P_C + 2) / 4;
		img->mpr[1+ioff][1+joff] = 
			img->mpr[0+ioff][3+joff] = (P_B + 2*P_C + P_D + 2) / 4;
		img->mpr[2+ioff][1+joff] = 
			img->mpr[1+ioff][3+joff] = (P_C + 2*P_D + P_E + 2) / 4;
		img->mpr[3+ioff][1+joff] = 
			img->mpr[2+ioff][3+joff] = (P_D + 2*P_E + P_F + 2) / 4;
		img->mpr[3+ioff][3+joff] = (P_E + 2*P_F + P_G + 2) / 4;
		break;
		
	case  HOR_UP_PRED://8  /* diagonal prediction -22.5 deg to horizontal plane */
// 		printf("HOR_UP_PRED://8    ");
		if (!block_available_left)
			printf ("warning: Intra_4x4_Horizontal_Up prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		img->mpr[0+ioff][0+joff] = (P_I + P_J + 1) / 2;
		img->mpr[1+ioff][0+joff] = (P_I + 2*P_J + P_K + 2) / 4;
		img->mpr[2+ioff][0+joff] = 
			img->mpr[0+ioff][1+joff] = (P_J + P_K + 1) / 2;
		img->mpr[3+ioff][0+joff] = 
			img->mpr[1+ioff][1+joff] = (P_J + 2*P_K + P_L + 2) / 4;
		img->mpr[2+ioff][1+joff] = 
			img->mpr[0+ioff][2+joff] = (P_K + P_L + 1) / 2;
		img->mpr[3+ioff][1+joff] = 
			img->mpr[1+ioff][2+joff] = (P_K + 2*P_L + P_L + 2) / 4;
		img->mpr[3+ioff][2+joff] = 
			img->mpr[1+ioff][3+joff] = 
			img->mpr[0+ioff][3+joff] = 
			img->mpr[2+ioff][2+joff] = 
			img->mpr[2+ioff][3+joff] = 
			img->mpr[3+ioff][3+joff] = P_L;
		break;
		
	case  HOR_DOWN_PRED://6  /* diagonal prediction -22.5 deg to horizontal plane */
// 		printf("HOR_DOWN_PRED://6    ");
		if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
			printf ("warning: Intra_4x4_Horizontal_Down prediction mode not allowed at mb %d\n",img->current_mb_nr);
		
		img->mpr[0+ioff][0+joff] = 
			img->mpr[2+ioff][1+joff] = (P_X + P_I + 1) / 2;
		img->mpr[1+ioff][0+joff] = 
			img->mpr[3+ioff][1+joff] = (P_I + 2*P_X + P_A + 2) / 4;
		img->mpr[2+ioff][0+joff] = (P_X + 2*P_A + P_B + 2) / 4;
		img->mpr[3+ioff][0+joff] = (P_A + 2*P_B + P_C + 2) / 4;
		img->mpr[0+ioff][1+joff] = 
			img->mpr[2+ioff][2+joff] = (P_I + P_J + 1) / 2;
		img->mpr[1+ioff][1+joff] = 
			img->mpr[3+ioff][2+joff] = (P_X + 2*P_I + P_J + 2) / 4;
		img->mpr[0+ioff][2+joff] = 
			img->mpr[2+ioff][3+joff] = (P_J + P_K + 1) / 2;
		img->mpr[1+ioff][2+joff] = 
			img->mpr[3+ioff][3+joff] = (P_I + 2*P_J + P_K + 2) / 4;
		img->mpr[0+ioff][3+joff] = (P_K + P_L + 1) / 2;
		img->mpr[1+ioff][3+joff] = (P_J + 2*P_K + P_L + 2) / 4;
		break;
		
	default:
// 		printf("default    ");
		printf("Error: illegal intra_4x4 prediction mode: %d\n",predmode);
		return SEARCH_SYNC;
		break;
  }
  
  return DECODING_OK;
}


/*!
***********************************************************************
* \return
*    best SAD
***********************************************************************
*/
int intrapred_luma_16x16_dec(struct img_par *img, //!< image parameters
                         int predmode)        //!< prediction mode
{
	int s0=0,s1,s2;
	
	int i,j;
	
	int ih,iv;
	int ib,ic,iaa;
	
	byte **imgY=dec_picture->imgY;
	
	int mb_nr=img->current_mb_nr;
	
	PixelPos up;          //!< pixel position p(0,-1)
	PixelPos left[17];    //!< pixel positions p(-1, -1..15)
	
	int up_avail, left_avail, left_up_avail;
	
	s1=s2=0;
	
	for (i=0;i<17;i++)
	{
		getNeighbour(mb_nr, -1 ,  i-1 , 1, &left[i]);
	}
	
	getNeighbour(mb_nr, 0     ,  -1 , 1, &up);
	
	if (!active_pps->constrained_intra_pred_flag)
	{
		up_avail   = up.available;
		left_avail = left[1].available;
		left_up_avail = left[0].available;
	}
	else
	{
		up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
		for (i=1, left_avail=1; i<17;i++)
			left_avail  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
		left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
	}
	
	switch (predmode)
	{
	case VERT_PRED_16:                       // vertical prediction from block above
		if (!up_avail)
			error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);
		for(j=0;j<MB_BLOCK_SIZE;j++)
			for(i=0;i<MB_BLOCK_SIZE;i++)
				img->mpr[i][j]=imgY[up.pos_y][up.pos_x+i];// store predicted 16x16 block
			break;
			
	case HOR_PRED_16:                        // horisontal prediction from left block
		if (!left_avail)
			error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);
		for(j=0;j<MB_BLOCK_SIZE;j++)
			for(i=0;i<MB_BLOCK_SIZE;i++)
				img->mpr[i][j]=imgY[left[j+1].pos_y][left[j+1].pos_x]; // store predicted 16x16 block
			break;
			
	case DC_PRED_16:                         // DC prediction
		s1=s2=0;
		for (i=0; i < MB_BLOCK_SIZE; i++)
		{
			if (up_avail)
				s1 += imgY[up.pos_y][up.pos_x+i];    // sum hor pix
			if (left_avail)
				s2 += imgY[left[i+1].pos_y][left[i+1].pos_x];    // sum vert pix
		}
		if (up_avail && left_avail)
			s0=(s1+s2+16)>>5;       // no edge
		if (!up_avail && left_avail)
			s0=(s2+8)>>4;              // upper edge
		if (up_avail && !left_avail)
			s0=(s1+8)>>4;              // left edge
		if (!up_avail && !left_avail)
			s0=128;                            // top left corner, nothing to predict from
		for(i=0;i<MB_BLOCK_SIZE;i++)
			for(j=0;j<MB_BLOCK_SIZE;j++)
			{
				img->mpr[i][j]=s0;
			}
			break;
	case PLANE_16:// 16 bit integer plan pred
		if (!up_avail || !left_up_avail  || !left_avail)
			error ("invalid 16x16 intra pred Mode PLANE_16",500);
		
		ih=0;
		iv=0;
		for (i=1;i<9;i++)
		{
			if (i<8)
				ih += i*(imgY[up.pos_y][up.pos_x+7+i] - imgY[up.pos_y][up.pos_x+7-i]);
			else
				ih += i*(imgY[up.pos_y][up.pos_x+7+i] - imgY[left[0].pos_y][left[0].pos_x]);
			
			iv += i*(imgY[left[8+i].pos_y][left[8+i].pos_x] - imgY[left[8-i].pos_y][left[8-i].pos_x]);
		}
		ib=(5*ih+32)>>6;
		ic=(5*iv+32)>>6;
		
		iaa=16*(imgY[up.pos_y][up.pos_x+15]+imgY[left[16].pos_y][left[16].pos_x]);
		for (j=0;j< MB_BLOCK_SIZE;j++)
		{
			for (i=0;i< MB_BLOCK_SIZE;i++)
			{
				img->mpr[i][j]=max(0,min((iaa+(i-7)*ib +(j-7)*ic + 16)>>5,255));
			}
		}// store plane prediction
		break;
		
	default:
		{                                    // indication of fault in bitstream,exit
			printf("illegal 16x16 intra prediction mode input: %d\n",predmode);
			return SEARCH_SYNC;
		}
	}
	
	return DECODING_OK;
}


void intrapred_chroma_dec(struct img_par *img, int uv)
{
	int i,j, ii, jj, ioff, joff;
	
	byte ***imgUV = dec_picture->imgUV;
	
	int js0=0;
	int js1=0;
	int js2=0;
	int js3=0;
	
	int js[2][2];
	
	int pred;
	int ih, iv, ib, ic, iaa;
	
	int mb_nr=img->current_mb_nr;
	Macroblock *currMB = &img->mb_data[img->current_mb_nr];
	
	PixelPos up;       //!< pixel position p(0,-1)
	PixelPos left[9];  //!< pixel positions p(-1, -1..8)
	
	int up_avail, left_avail[2], left_up_avail;
	
	for (i=0;i<9;i++)
	{
		getNeighbour(mb_nr, -1 ,  i-1 , 0, &left[i]);
	}
	
	getNeighbour(mb_nr, 0     ,  -1 , 0, &up);
	
	if (!active_pps->constrained_intra_pred_flag)
	{
		up_avail   = up.available;
		left_avail[0] = left_avail[1] = left[1].available;
		left_up_avail = left[0].available;
	}
	else
	{
		up_avail      = up.available ? img->intra_block[up.mb_addr] : 0;
		for (i=1, left_avail[0]=1; i<5;i++)
			left_avail[0]  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
		for (i=5, left_avail[1]=1; i<9;i++)
			left_avail[1]  &= left[i].available ? img->intra_block[left[i].mb_addr]: 0;
		left_up_avail = left[0].available ? img->intra_block[left[0].mb_addr]: 0;
	}
	
	if (currMB->c_ipred_mode == DC_PRED_8)
	{
		for(i=0;i<4;i++)
		{
			if(up_avail)
			{
				js0=js0+imgUV[uv][up.pos_y][up.pos_x+i];
				js1=js1+imgUV[uv][up.pos_y][up.pos_x+i+4];
			}
			if(left_avail[0])
			{
				js2=js2+imgUV[uv][left[1+i].pos_y][left[1+i].pos_x];
			}
			if(left_avail[1])
			{
				js3=js3+imgUV[uv][left[1+i+4].pos_y][left[1+i+4].pos_x];
			}
		}
		if(up_avail && left_avail[0])
		{
			js[0][0]=(js0+js2+4)/8;
			js[1][0]=(js1+2)/4;
		}
		if(up_avail && left_avail[1])
		{
			js[0][1]=(js3+2)/4;
			js[1][1]=(js1+js3+4)/8;
		}
		if(up_avail && !left_avail[0])
		{
			js[0][0]=(js0+2)/4;
			js[1][0]=(js1+2)/4;
		}
		if(up_avail && !left_avail[1])
		{
			js[0][1]=(js0+2)/4;
			js[1][1]=(js1+2)/4;
		}
		if(left_avail[0] && !up_avail)
		{
			js[0][0]=(js2+2)/4;
			js[1][0]=(js2+2)/4;
		}
		if(left_avail[1] && !up_avail)
		{
			js[0][1]=(js3+2)/4;
			js[1][1]=(js3+2)/4;
		}
		if(!up_avail && !left_avail[0])
		{
			js[0][0]=128;
			js[1][0]=128;
		}
		if(!up_avail && !left_avail[1])
		{
			js[0][1]=128;
			js[1][1]=128;
		}
	}
	
	for (j=0;j<2;j++)
	{
		joff=j*4;
		
		for(i=0;i<2;i++)
		{
			ioff=i*4;
			
			switch (currMB->c_ipred_mode)
			{
			case DC_PRED_8://0
// 				printf("DC_PRED_8  ");
				for (ii=0; ii<4; ii++)
					for (jj=0; jj<4; jj++)
					{
						img->mpr[ii+ioff][jj+joff]=js[i][j];
					}
					break;
			case HOR_PRED_8://1
// 				printf("HOR_PRED_8  ");
				if ( !left_avail[0] || !left_avail[1] )
					error("unexpected HOR_PRED_8 chroma intra prediction mode",-1);
				
				for (jj=0; jj<4; jj++)
				{
					pred = imgUV[uv][left[1+jj+joff].pos_y][left[1+jj+joff].pos_x];
					for (ii=0; ii<4; ii++)
						img->mpr[ii+ioff][jj+joff]=pred;
				}
				break;
			case VERT_PRED_8://2
// 				printf("VERT_PRED_8  ");
				if (!up_avail)
					error("unexpected VERT_PRED_8 chroma intra prediction mode",-1);
				
				for (ii=0; ii<4; ii++)
				{
					pred = imgUV[uv][up.pos_y][up.pos_x+ii+ioff];
					for (jj=0; jj<4; jj++)
						img->mpr[ii+ioff][jj+joff]=pred;
				}
				break;
			case PLANE_8://3
// 				printf("PLANE_8  ");
				if (!left_up_avail || !left_avail[0] || !left_avail[1] || !up_avail)
					error("unexpected PLANE_8 chroma intra prediction mode",-1);
				
				ih=iv=0;
				
				for (ii=1;ii<5;ii++)
				{
					if (ii<4)
						ih += ii*(imgUV[uv][up.pos_y][up.pos_x+3+ii] - imgUV[uv][up.pos_y][up.pos_x+3-ii]);
					else
						ih += ii*(imgUV[uv][up.pos_y][up.pos_x+3+ii] - imgUV[uv][left[0].pos_y][left[0].pos_x]);
					
					iv += ii*(imgUV[uv][left[4+ii].pos_y][left[4+ii].pos_x] - imgUV[uv][left[4-ii].pos_y][left[4-ii].pos_x]);
				}
				ib=(17*ih+16)>>5;
				ic=(17*iv+16)>>5;
				iaa=16*(imgUV[uv][up.pos_y][up.pos_x+7]+imgUV[uv][left[8].pos_y][left[8].pos_x]);
				for (ii=0; ii<4; ii++)
					for (jj=0; jj<4; jj++)
						img->mpr[ii+ioff][jj+joff]=max(0,min(255,(iaa+(ii+ioff-3)*ib +(jj+joff-3)*ic + 16)>>5));
					break;
			default:
				error("illegal chroma intra prediction mode", 600);
				break;
			}
		}
	}
}


/*!
***********************************************************************
* \brief
*    Inverse 4x4 transformation, transforms cof to m7
***********************************************************************
*/
void itrans(struct img_par *img, //!< image parameters
            int ioff,            //!< index to 4x4 block
            int joff,            //!<
            int i0,              //!<
            int j0)              //!<
{
	int i,j,i1,j1;
	int m5[4];
	int m6[4];
	
	// horizontal
	for (j=0;j<BLOCK_SIZE;j++)
	{
		for (i=0;i<BLOCK_SIZE;i++)
		{
			m5[i]=img->cof[i0][j0][i][j];
		}
		m6[0]=(m5[0]+m5[2]);
		m6[1]=(m5[0]-m5[2]);
		m6[2]=(m5[1]>>1)-m5[3];
		m6[3]=m5[1]+(m5[3]>>1);
		
		for (i=0;i<2;i++)
		{
			i1=3-i;
			img->m7[i][j]=m6[i]+m6[i1];
			img->m7[i1][j]=m6[i]-m6[i1];
		}
	}
	// vertical
	for (i=0;i<BLOCK_SIZE;i++)
	{
		for (j=0;j<BLOCK_SIZE;j++)
			m5[j]=img->m7[i][j];
		
		m6[0]=(m5[0]+m5[2]);
		m6[1]=(m5[0]-m5[2]);
		m6[2]=(m5[1]>>1)-m5[3];
		m6[3]=m5[1]+(m5[3]>>1);
		
		for (j=0;j<2;j++)
		{
			j1=3-j;
			img->m7[i][j] =max(0,min(255,(m6[j]+m6[j1]+(img->mpr[i+ioff][j+joff] <<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
			img->m7[i][j1]=max(0,min(255,(m6[j]-m6[j1]+(img->mpr[i+ioff][j1+joff]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
		}
	}
	
}


/*!
 ***********************************************************************
 * \brief
 *    invers  transform
 ***********************************************************************
 */
void itrans_2(
   struct img_par *img) //!< image parameters
{
  int i,j,i1,j1;
  int M5[4];
  int M6[4];

  int qp_per = (img->qp-MIN_QP)/6;
  int qp_rem = (img->qp-MIN_QP)%6;

  // horizontal
  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
      M5[i]=img->cof[i][j][0][0];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->cof[i ][j][0][0]= M6[i]+M6[i1];
      img->cof[i1][j][0][0]=M6[i]-M6[i1];
    }
  }

  // vertical
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      M5[j]=img->cof[i][j][0][0];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->cof[i][j][0][0] = (((M6[j]+M6[j1])*dequant_coef[qp_rem][0][0]<<qp_per)+2)>>2;
      img->cof[i][j1][0][0]= (((M6[j]-M6[j1])*dequant_coef[qp_rem][0][0]<<qp_per)+2)>>2;
    }
  }
}



void itrans_sp(struct img_par *img,  //!< image parameters
               int ioff,             //!< index to 4x4 block
               int joff,             //!<
               int i0,               //!<
               int j0)               //!<
{
	int i,j,i1,j1;
	int m5[4];
	int m6[4];
	int predicted_block[BLOCK_SIZE][BLOCK_SIZE],ilev;
	
	int qp_per = (img->qp-MIN_QP)/6;
	int qp_rem = (img->qp-MIN_QP)%6;
	int q_bits    = Q_BITS+qp_per;
	
	int qp_per_sp = (img->qpsp-MIN_QP)/6;
	int qp_rem_sp = (img->qpsp-MIN_QP)%6;
	int q_bits_sp    = Q_BITS+qp_per_sp;
	int qp_const2=(1<<q_bits_sp)/2;  //sp_pred
	
	if (img->type == SI_SLICE) //ES  modified
	{
		qp_per = (img->qpsp-MIN_QP)/6;
		qp_rem = (img->qpsp-MIN_QP)%6;
		q_bits = Q_BITS+qp_per;
	}
	
	for (j=0; j< BLOCK_SIZE; j++)
		for (i=0; i< BLOCK_SIZE; i++)
			predicted_block[i][j]=img->mpr[i+ioff][j+joff];
		for (j=0; j < BLOCK_SIZE; j++)
		{
			for (i=0; i < 2; i++)
			{
				i1=3-i;
				m5[i]=predicted_block[i][j]+predicted_block[i1][j];
				m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
			}
			predicted_block[0][j]=(m5[0]+m5[1]);
			predicted_block[2][j]=(m5[0]-m5[1]);
			predicted_block[1][j]=m5[3]*2+m5[2];
			predicted_block[3][j]=m5[3]-m5[2]*2;
		}
		
		//  Vertival transform
		
		for (i=0; i < BLOCK_SIZE; i++)
		{
			for (j=0; j < 2; j++)
			{
				j1=3-j;
				m5[j]=predicted_block[i][j]+predicted_block[i][j1];
				m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
			}
			predicted_block[i][0]=(m5[0]+m5[1]);
			predicted_block[i][2]=(m5[0]-m5[1]);
			predicted_block[i][1]=m5[3]*2+m5[2];
			predicted_block[i][3]=m5[3]-m5[2]*2;
		}
		
		for (j=0;j<BLOCK_SIZE;j++)
			for (i=0;i<BLOCK_SIZE;i++)
			{
				// recovering coefficient since they are already dequantized earlier
				img->cof[i0][j0][i][j]=(img->cof[i0][j0][i][j] >> qp_per) / dequant_coef[qp_rem][i][j]; 
				if(img->sp_switch || img->type==SI_SLICE)  //M.W. patched for SI
				{
					ilev=(abs(predicted_block[i][j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp; //ES added
					ilev= sign(ilev,predicted_block[i][j])+ img->cof[i0][j0][i][j];                           //ES added
					img->cof[i0][j0][i][j] = sign(abs(ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp ,ilev) ; //ES added 
				}                                                                                             //ES added
				else
				{                                                                                          //ES added
					ilev=((img->cof[i0][j0][i][j]*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6)+predicted_block[i][j] ;
					img->cof[i0][j0][i][j]=sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp, ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
				}
			}
			// horizontal
			for (j=0;j<BLOCK_SIZE;j++)
			{
				for (i=0;i<BLOCK_SIZE;i++)
				{
					m5[i]=img->cof[i0][j0][i][j];
				}
				m6[0]=(m5[0]+m5[2]);
				m6[1]=(m5[0]-m5[2]);
				m6[2]=(m5[1]>>1)-m5[3];
				m6[3]=m5[1]+(m5[3]>>1);
				
				for (i=0;i<2;i++)
				{
					i1=3-i;
					img->m7[i][j]=m6[i]+m6[i1];
					img->m7[i1][j]=m6[i]-m6[i1];
				}
			}
			// vertical
			for (i=0;i<BLOCK_SIZE;i++)
			{
				for (j=0;j<BLOCK_SIZE;j++)
					m5[j]=img->m7[i][j];
				
				m6[0]=(m5[0]+m5[2]);
				m6[1]=(m5[0]-m5[2]);
				m6[2]=(m5[1]>>1)-m5[3];
				m6[3]=m5[1]+(m5[3]>>1);
				
				for (j=0;j<2;j++)
				{
					j1=3-j;
					img->m7[i][j] =max(0,min(255,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
					img->m7[i][j1]=max(0,min(255,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
				}
			}
}

/*!
***********************************************************************
* \brief
*    The routine performs transform,quantization,inverse transform, adds the diff.
*    to the prediction and writes the result to the decoded luma frame. Includes the
*    RD constrained quantization also.
*
* \par Input:
*    block_x,block_y: Block position inside a macro block (0,4,8,12).
*
* \par Output:
*    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels. \n
*    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
************************************************************************
*/
void copyblock_sp_dec(struct img_par *img,int block_x,int block_y)
{
	int sign(int a,int b);
	
	int i,j,i1,j1,m5[4],m6[4];
	
	int predicted_block[BLOCK_SIZE][BLOCK_SIZE];
	int qp_per = (img->qpsp-MIN_QP)/6;
	int qp_rem = (img->qpsp-MIN_QP)%6;
	int q_bits    = Q_BITS+qp_per;
	int qp_const2=(1<<q_bits)/2;  //sp_pred
	
	
	//  Horizontal transform
	for (j=0; j< BLOCK_SIZE; j++)
		for (i=0; i< BLOCK_SIZE; i++)
			predicted_block[i][j]=img->mpr[i+block_x][j+block_y];
		
		for (j=0; j < BLOCK_SIZE; j++)
		{
			for (i=0; i < 2; i++)
			{
				i1=3-i;
				m5[i]=predicted_block[i][j]+predicted_block[i1][j];
				m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
			}
			predicted_block[0][j]=(m5[0]+m5[1]);
			predicted_block[2][j]=(m5[0]-m5[1]);
			predicted_block[1][j]=m5[3]*2+m5[2];
			predicted_block[3][j]=m5[3]-m5[2]*2;
		}
		
		//  Vertival transform
		
		for (i=0; i < BLOCK_SIZE; i++)
		{
			for (j=0; j < 2; j++)
			{
				j1=3-j;
				m5[j]=predicted_block[i][j]+predicted_block[i][j1];
				m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
			}
			predicted_block[i][0]=(m5[0]+m5[1]);
			predicted_block[i][2]=(m5[0]-m5[1]);
			predicted_block[i][1]=m5[3]*2+m5[2];
			predicted_block[i][3]=m5[3]-m5[2]*2;
		}
		
		// Quant
		for (j=0;j < BLOCK_SIZE; j++)
			for (i=0; i < BLOCK_SIZE; i++)
				img->m7[i][j]=sign((abs(predicted_block[i][j])* quant_coef[qp_rem][i][j]+qp_const2)>> q_bits,predicted_block[i][j])*dequant_coef[qp_rem][i][j]<<qp_per;
			
			//     IDCT.
			//     horizontal
			
			for (j=0;j<BLOCK_SIZE;j++)
			{
				for (i=0;i<BLOCK_SIZE;i++)
				{
					m5[i]=img->m7[i][j];
				}
				m6[0]=(m5[0]+m5[2]);
				m6[1]=(m5[0]-m5[2]);
				m6[2]=(m5[1]>>1)-m5[3];
				m6[3]=m5[1]+(m5[3]>>1);
				
				for (i=0;i<2;i++)
				{
					i1=3-i;
					img->m7[i][j]=m6[i]+m6[i1];
					img->m7[i1][j]=m6[i]-m6[i1];
				}
			}
			// vertical
			for (i=0;i<BLOCK_SIZE;i++)
			{
				for (j=0;j<BLOCK_SIZE;j++)
					m5[j]=img->m7[i][j];
				
				m6[0]=(m5[0]+m5[2]);
				m6[1]=(m5[0]-m5[2]);
				m6[2]=(m5[1]>>1)-m5[3];
				m6[3]=m5[1]+(m5[3]>>1);
				
				for (j=0;j<2;j++)
				{
					j1=3-j;
					img->m7[i][j] =max(0,min(255,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
					img->m7[i][j1]=max(0,min(255,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
				}
			}
			
			//  Decoded block moved to frame memory
			
			for (j=0; j < BLOCK_SIZE; j++)
				for (i=0; i < BLOCK_SIZE; i++)
					dec_picture->imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];
				
}

void itrans_sp_chroma(struct img_par *img,int ll)
{
	int i,j,i1,j2,ilev,n2,n1,j1,mb_y;
	int m5[BLOCK_SIZE];
	int predicted_chroma_block[MB_BLOCK_SIZE/2][MB_BLOCK_SIZE/2],mp1[BLOCK_SIZE];
	int qp_per,qp_rem,q_bits;
	int qp_per_sp,qp_rem_sp,q_bits_sp,qp_const2;
	
	qp_per    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)/6;
	qp_rem    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)%6;
	q_bits    = Q_BITS+qp_per;
	
	qp_per_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)/6;
	qp_rem_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)%6;
	q_bits_sp    = Q_BITS+qp_per_sp;
	qp_const2=(1<<q_bits_sp)/2;  //sp_pred
	
	if (img->type == SI_SLICE)
	{
		qp_per    = ((img->qpsp < 0 ? img->qpsp : QP_SCALE_CR[img->qpsp]) - MIN_QP) / 6;
		qp_rem    = ((img->qpsp < 0 ? img->qpsp : QP_SCALE_CR[img->qpsp]) - MIN_QP) % 6;
		q_bits    = Q_BITS + qp_per;
	}
	
	for (j=0; j < MB_BLOCK_SIZE/2; j++)
		for (i=0; i < MB_BLOCK_SIZE/2; i++)
		{
			predicted_chroma_block[i][j]=img->mpr[i][j];
			img->mpr[i][j]=0;
		}
		for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
		{
			for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
			{
				//  Horizontal transform.
				for (j=0; j < BLOCK_SIZE; j++)
				{
					mb_y=n2+j;
					for (i=0; i < 2; i++)
					{
						i1=3-i;
						m5[i]=predicted_chroma_block[i+n1][mb_y]+predicted_chroma_block[i1+n1][mb_y];
						m5[i1]=predicted_chroma_block[i+n1][mb_y]-predicted_chroma_block[i1+n1][mb_y];
					}
					predicted_chroma_block[n1][mb_y]  =(m5[0]+m5[1]);
					predicted_chroma_block[n1+2][mb_y]=(m5[0]-m5[1]);
					predicted_chroma_block[n1+1][mb_y]=m5[3]*2+m5[2];
					predicted_chroma_block[n1+3][mb_y]=m5[3]-m5[2]*2;
				}
				
				//  Vertical transform.
				
				for (i=0; i < BLOCK_SIZE; i++)
				{
					j1=n1+i;
					for (j=0; j < 2; j++)
					{
						j2=3-j;
						m5[j]=predicted_chroma_block[j1][n2+j]+predicted_chroma_block[j1][n2+j2];
						m5[j2]=predicted_chroma_block[j1][n2+j]-predicted_chroma_block[j1][n2+j2];
					}
					predicted_chroma_block[j1][n2+0]=(m5[0]+m5[1]);
					predicted_chroma_block[j1][n2+2]=(m5[0]-m5[1]);
					predicted_chroma_block[j1][n2+1]=m5[3]*2+m5[2];
					predicted_chroma_block[j1][n2+3]=m5[3]-m5[2]*2;
				}
			}
		}
		
		//     2X2 transform of DC coeffs.
		mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
		mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
		mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
		mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
		
		for (n1=0; n1 < 2; n1 ++)
			for (n2=0; n2 < 2; n2 ++)
			{
				if (img->sp_switch || img->type==SI_SLICE)  //M.W. patched for SI
				{
					//quantization fo predicted block
					ilev=(abs (mp1[n1+n2*2]) * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp + 1); 
					//addition 	  
					ilev=img->cof[n1+ll][4+n2][0][0]+sign(ilev,mp1[n1+n2*2]);                                   
					//dequantization
					mp1[n1+n2*2] =ilev*dequant_coef[qp_rem_sp][0][0]<<qp_per_sp;                                
				}   
				else
				{
					ilev=((img->cof[n1+ll][4+n2][0][0]*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)+mp1[n1+n2*2] ;
					mp1[n1+n2*2]=sign((abs(ilev)* quant_coef[qp_rem_sp][0][0]+ 2 * qp_const2)>> (q_bits_sp+1),ilev)*dequant_coef[qp_rem_sp][0][0]<<qp_per_sp;
				}
			}
			
			
			for (n2=0; n2 < 2; n2 ++)
				for (n1=0; n1 < 2; n1 ++)
					for (i=0;i< BLOCK_SIZE; i++)
						for (j=0;j< BLOCK_SIZE; j++)
						{
							// recovering coefficient since they are already dequantized earlier
							img->cof[n1+ll][4+n2][i][j] = (img->cof[n1+ll][4+n2][i][j] >> qp_per) / dequant_coef[qp_rem][i][j];
							
							if (img->sp_switch || img->type==SI_SLICE)  //M.W. patched for SI
							{
								//quantization of the predicted block
								ilev =  (abs(predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp;
								//addition of the residual
								ilev = sign(ilev,predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j]) + img->cof[n1+ll][4+n2][i][j];
								// Inverse quantization 
								img->cof[n1+ll][4+n2][i][j] = ilev * dequant_coef[qp_rem_sp][i][j] << qp_per_sp  ;
							}
							else
							{
								//dequantization and addition of the predicted block
								ilev=((img->cof[n1+ll][4+n2][i][j]*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6)+predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j] ;
								//quantization and dequantization
								img->cof[n1+ll][4+n2][i][j] = sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2)>> q_bits_sp,ilev)*dequant_coef[qp_rem_sp][i][j]<<qp_per_sp;
							}
						}
						img->cof[0+ll][4][0][0]=(mp1[0]+mp1[1]+mp1[2]+mp1[3])>>1;
						img->cof[1+ll][4][0][0]=(mp1[0]-mp1[1]+mp1[2]-mp1[3])>>1;
						img->cof[0+ll][5][0][0]=(mp1[0]+mp1[1]-mp1[2]-mp1[3])>>1;
						img->cof[1+ll][5][0][0]=(mp1[0]-mp1[1]-mp1[2]+mp1[3])>>1;
}

int sign(int a , int b)
{
	int x;
	
	x=abs(a);
	if (b>0)
		return(x);
	else return(-x);
}

/////////////////////////////////////<