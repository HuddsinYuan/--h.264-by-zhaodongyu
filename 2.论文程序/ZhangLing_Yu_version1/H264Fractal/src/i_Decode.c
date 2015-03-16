#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "i_codec_types.h"
#include "i_decode_global.h"
#include "i_encode_globals.h"
#include <malloc.h>

extern unsigned char **imgY_ref_temp,***imgUV_ref_temp;
I_FRM_ENC_GLOBAL *T_I;

short current_region;
BYTE block_region;
int num_byte=0;

static int Zig_Zag[8][8]={{0,1,5,6,14,15,27,28},
        {2,4,7,13,16,26,29,42},
        {3,8,12,17,25,30,41,43},
        {9,11,18,24,37,40,44,53},
        {10,19,23,32,39,45,52,54},
        {20,22,33,38,46,51,55,60},
        {21,34,37,47,50,56,59,61},
        {35,36,48,49,57,58,62,63}
       };

BYTE *blockTemp_Y,*blockTemp_UV[2];
BYTE *blockTemp;
BYTE *pQt;
DWORD width,height;
void set_quant_table2(BYTE *basic_table,BYTE scale_factor,BYTE *newtable)
{
	
	int tmpVal = 0;                  //临时变量
	DWORD i    = 0; 
	
	if (scale_factor < 1) 
		scale_factor = 1;               //限制质量系数
	if (scale_factor > 100) 
		scale_factor = 100;
	
	//非线性映射 1->5000, 10->500, 25->200, 50->100, 75->50, 100->0
	if (scale_factor < 50)
		scale_factor = 5000 / scale_factor;
	else
		scale_factor = 200 - scale_factor * 2;
	
	for (i = 0; i < 64; ++i)
	{
		tmpVal = (basic_table[i] * scale_factor + 50L) / 100L;
		
		if (tmpVal < 1)                 //数值范围限定
		{
			tmpVal = 1L;
		}
		if (tmpVal > 255)
		{
			tmpVal = 255L;
		}
		newtable[zigzag[i]] = (BYTE)(tmpVal);
	} 
}

void init_Table()
{
	short i,j,huftabindex;

	And[0]=0;And[1]=1;And[2]=3;And[3]=7;
	And[4]=0xf;And[5]=0x1f;And[6]=0x3f;And[7]=0x7f;And[8]=0xff;
	sizei=sizej=0;

	rrun=vvalue=0;
	BitPos=0;
	CurByte=0;
	restart=0;
    comp_num=3;
	HufTabIndex=0;
	
	set_quant_table2(std_luminance_qt,i_quality,qt_table[0]);
	set_quant_table2(std_chrominance_qt,i_quality,qt_table[1]);  //设置量化信息
	
	for(i=0;i<4;i++)
		for(j=0;j<16;j++)
		{
			code_pos_table[i][j]=0;
			huf_max_value[i][j]=0;
			huf_min_value[i][j]=0;
		}
///////////////////////////
	for(j=0;j<15;j++)
		code_len_table[0][j]=std_dc_luminance_nrcodes[j+1];
	code_len_table[0][j]=0;
	for(j=0;j<15;j++)
		code_len_table[1][j]=std_dc_chrominance_nrcodes[j+1];
	code_len_table[1][j]=0;
	for(j=0;j<15;j++)
		code_len_table[2][j]=std_ac_luminance_nrcodes[j+1];
	code_len_table[2][j]=std_ac_luminance_nrcodes[j+1];
	for(j=0;j<15;j++)
		code_len_table[3][j]=std_ac_chrominance_nrcodes[j+1];
	code_len_table[3][j]=std_ac_chrominance_nrcodes[j+1];
	////////////////////////////////////////////huffman
	for(j=0;j<12;j++)
		code_value_table[0][j]=std_dc_luminance_values[j];
	for (j=12;j<256;j++)
		code_value_table[0][j]=0;
	for(j=0;j<12;j++)
		code_value_table[1][j]=std_dc_chrominance_values[j];
	for (j=12;j<256;j++)
		code_value_table[1][j]=0;
	for(j=0;j<162;j++)
		code_value_table[2][j]=std_ac_luminance_values[j];
	for (j=162;j<256;j++)
		code_value_table[2][j]=0;
	for(j=0;j<162;j++)
		code_value_table[3][j]=std_ac_chrominance_values[j];
	for (j=162;j<256;j++)
		code_value_table[3][j]=0;
    //////////////////////////////////////////////////////
	for (huftabindex=0;huftabindex<4;huftabindex++)
	{
		i=0;
		while (code_len_table[huftabindex][i]==0)
			i++;
		for (j=0;j<i;j++)
		{
			huf_min_value[huftabindex][j]=0;
			huf_max_value[huftabindex][j]=0;
		}
		huf_min_value[huftabindex][i]=0;
		huf_max_value[huftabindex][i]=code_len_table[huftabindex][i]-1;
		for (j=i+1;j<16;j++)
		{
			huf_min_value[huftabindex][j]=(huf_max_value[huftabindex][j-1]+1)<<1;
			huf_max_value[huftabindex][j]=huf_min_value[huftabindex][j]+code_len_table[huftabindex][j]-1;
		}
		code_pos_table[huftabindex][0]=0;
		for (j=1;j<16;j++)
			code_pos_table[huftabindex][j]=code_len_table[huftabindex][j-1]+code_pos_table[huftabindex][j-1];
	}
	
/////////////////////////////////////////////////////////////
	for(i=0;i<64;i++)
	{
		Y[i]=0;
		U[i]=0;
		V[i]=0;
		BlockBuffer[i]=0;
	}
	coef=0;
}

BYTE  ReadByte()
{
	BYTE  i;
	num_byte++;
	i=*(lp++);
	if(i==0xff)
	{
		lp++;
		num_byte++;
	}
	BitPos=8;
	CurByte=i;
	return i;
}

int DecodeElement()
{
	int thiscode,tempcode;
	unsigned short temp,valueex;
	short codelen;
	BYTE hufexbyte,runsize,tempsize,sign;
	BYTE newbyte,lastbyte;
	
	if(BitPos>=1){
		BitPos--;
		thiscode=(BYTE)CurByte>>BitPos;
		CurByte=CurByte&And[BitPos];
	}
	else{
		lastbyte=ReadByte();
		BitPos--;
		newbyte=CurByte&And[BitPos];
		thiscode=lastbyte>>7;
		CurByte=newbyte;
	}
	codelen=1;
	while ((thiscode<huf_min_value[HufTabIndex][codelen-1])||
		(code_len_table[HufTabIndex][codelen-1]==0)||
		(thiscode>huf_max_value[HufTabIndex][codelen-1]))
	{
		if(BitPos>=1){
			BitPos--;
			tempcode=(BYTE)CurByte>>BitPos;
			CurByte=CurByte&And[BitPos];
		}
		else{
			lastbyte=ReadByte();
			BitPos--;
			newbyte=CurByte&And[BitPos];
			tempcode=(BYTE)lastbyte>>7;
			CurByte=newbyte;
		}
		thiscode=(thiscode<<1)+tempcode;
		codelen++;
		if(codelen>16)
			return FUNC_FORMAT_ERROR;
	}  //while
	temp=thiscode-huf_min_value[HufTabIndex][codelen-1]+code_pos_table[HufTabIndex][codelen-1];
	hufexbyte=(BYTE)code_value_table[HufTabIndex][temp];
	rrun=(short)(hufexbyte>>4);
	runsize=hufexbyte&0x0f;
	if(runsize==0){
		vvalue=0;
		return FUNC_OK;
	}
	tempsize=runsize;
	if(BitPos>=runsize){
		BitPos-=runsize;
		valueex=(BYTE)CurByte>>BitPos;
		CurByte=CurByte&And[BitPos];
	}
	else{
		valueex=CurByte;
		tempsize-=BitPos;
		while(tempsize>8){
			lastbyte=ReadByte();
			valueex=(valueex<<8)+(BYTE)lastbyte;
			tempsize-=8;
		}  //while
		lastbyte=ReadByte();
		BitPos-=tempsize;
		valueex=(valueex<<tempsize)+(lastbyte>>BitPos);
		CurByte=lastbyte&And[BitPos];
	}  //else
	sign=valueex>>(runsize-1);
	if(sign)
		vvalue=valueex;
	else{
		valueex=valueex^0xffff;
		temp=0xffff<<runsize;
		vvalue=-(short)(valueex^temp);
	}
	return FUNC_OK;
}

int HufBlock(BYTE dchufindex,BYTE achufindex)   //huffman 
{
	short count=0;
	short i;
	int funcret;
	
	//dc
	HufTabIndex=dchufindex;  //huffman 表的系数
	funcret=DecodeElement();
	if(funcret!=FUNC_OK)
		return funcret;
	
	BlockBuffer[count++]=vvalue;
	//ac
	HufTabIndex=achufindex;
	while (count<64){
		funcret=DecodeElement();
		if(funcret!=FUNC_OK)
			return funcret;
		if ((rrun==0)&&(vvalue==0)){
			for (i=count;i<64;i++)
				BlockBuffer[i]=0;
			count=64;
		}
		else{
			for (i=0;i<rrun;i++)
				BlockBuffer[count++]=0;
			BlockBuffer[count++]=vvalue;
		}
	}
	return FUNC_OK;
}

int DecodeBlock(BYTE dcHufindex,BYTE acHufindex)
{
	short i;
	int funcret;

	for (i=0;i<1;i++)  //Y
	{
		funcret=HufBlock(dcHufindex,acHufindex);
		if (funcret!=FUNC_OK)
			return funcret;
		BlockBuffer[0]=BlockBuffer[0]+coef;
		coef=BlockBuffer[0];
	}
	return FUNC_OK;
}

void Fast_IDCT(int * block)
{
 short i;

 for (i=0; i<8; i++)
  idctrow(block+8*i);

 for (i=0; i<8; i++)
  idctcol(block+i);
}

void idctrow(int * blk)
{
	int x0, x1, x2, x3, x4, x5, x6, x7, x8;
	//intcut
	if (!((x1 = blk[4]<<11) | (x2 = blk[6]) | (x3 = blk[2]) |
		(x4 = blk[1]) | (x5 = blk[7]) | (x6 = blk[5]) | (x7 = blk[3])))
	{
		blk[0]=blk[1]=blk[2]=blk[3]=blk[4]=blk[5]=blk[6]=blk[7]=blk[0]<<3;
		return;
	}
	x0 = (blk[0]<<11) + 128; // for proper rounding in the fourth stage 
	//first stage
	x8 = W7*(x4+x5);
	x4 = x8 + (W1-W7)*x4;
	x5 = x8 - (W1+W7)*x5;
	x8 = W3*(x6+x7);
	x6 = x8 - (W3-W5)*x6;
	x7 = x8 - (W3+W5)*x7;
	//second stage
	x8 = x0 + x1;
	x0 -= x1;
	x1 = W6*(x3+x2);
	x2 = x1 - (W2+W6)*x2;
	x3 = x1 + (W2-W6)*x3;
	x1 = x4 + x6;
	x4 -= x6;
	x6 = x5 + x7;
	x5 -= x7;
	//third stage
	x7 = x8 + x3;
	x8 -= x3;
	x3 = x0 + x2;
	x0 -= x2;
	x2 = (181*(x4+x5)+128)>>8;
	x4 = (181*(x4-x5)+128)>>8;
	//fourth stage
	blk[0] = (x7+x1)>>8;
	blk[1] = (x3+x2)>>8;
	blk[2] = (x0+x4)>>8;
	blk[3] = (x8+x6)>>8;
	blk[4] = (x8-x6)>>8;
	blk[5] = (x0-x4)>>8;
	blk[6] = (x3-x2)>>8;
	blk[7] = (x7-x1)>>8;
}
//////////////////////////////////////////////////////////////////////////////
void idctcol(int * blk)
{
	int x0, x1, x2, x3, x4, x5, x6, x7, x8;
	//intcut
	if (!((x1 = (blk[8*4]<<8)) | (x2 = blk[8*6]) | (x3 = blk[8*2]) |
		(x4 = blk[8*1]) | (x5 = blk[8*7]) | (x6 = blk[8*5]) | (x7 = blk[8*3])))
	{
		blk[8*0]=blk[8*1]=blk[8*2]=blk[8*3]=blk[8*4]=blk[8*5]
			=blk[8*6]=blk[8*7]=iclp[(blk[8*0]+32)>>6];
		return;
	}
	x0 = (blk[8*0]<<8) + 8192;
	//first stage
	x8 = W7*(x4+x5) + 4;
	x4 = (x8+(W1-W7)*x4)>>3;
	x5 = (x8-(W1+W7)*x5)>>3;
	x8 = W3*(x6+x7) + 4;
	x6 = (x8-(W3-W5)*x6)>>3;
	x7 = (x8-(W3+W5)*x7)>>3;
	//second stage
	x8 = x0 + x1;
	x0 -= x1;
	x1 = W6*(x3+x2) + 4;
	x2 = (x1-(W2+W6)*x2)>>3;
	x3 = (x1+(W2-W6)*x3)>>3;
	x1 = x4 + x6;
	x4 -= x6;
	x6 = x5 + x7;
	x5 -= x7;
	//third stage
	x7 = x8 + x3;
	x8 -= x3;
	x3 = x0 + x2;
	x0 -= x2;
	x2 = (181*(x4+x5)+128)>>8;
	x4 = (181*(x4-x5)+128)>>8;
	//fourth stage
	blk[8*0] = iclp[(x7+x1)>>14];
	blk[8*1] = iclp[(x3+x2)>>14];
	blk[8*2] = iclp[(x0+x4)>>14];
	blk[8*3] = iclp[(x8+x6)>>14];
	blk[8*4] = iclp[(x8-x6)>>14];
	blk[8*5] = iclp[(x0-x4)>>14];
	blk[8*6] = iclp[(x3-x2)>>14];
	blk[8*7] = iclp[(x7-x1)>>14];
}


void IQtIZzBlock()
{
	int i,j;
	short tag;
	int buffer2[8][8];
	int *buffer1;
	int offset;
	
	offset=128;
	for(i=0;i<8;i++)
		for(j=0;j<8;j++)
		{
			tag=Zig_Zag[i][j];
			buffer2[i][j]=(int)BlockBuffer[tag]*(int)pQt[tag];
		}
	
	buffer1=(int *)buffer2;
	Fast_IDCT(buffer1);
	if(num_regions>1)
	{
		for(i=0;i<8;i++)
			for(j=0;j<8;j++)
			{
				blockTemp[(sizei+i)*width+sizej+j]=(BYTE)(bound(buffer2[i][j]+offset))
					*(plane[(sizei+i)*width+sizej+j]/GREY_LEVELS==current_region);
			}
	}
	
	else
	{
		for(i=0;i<8;i++)
			for(j=0;j<8;j++)
			{
				blockTemp[(sizei+i)*width+sizej+j]=(BYTE)(bound(buffer2[i][j]+offset));
			}	
	}

}

void match_Block()
{
	WORD x,y;
	

	while(1)
	{
		block_region=0;
		for (y=0;y<8;y++)
		{
			for (x=0;x<8;x++)
			{
				block_region|=1<<(int)((double)plane[(sizei+y)*width+sizej+x]/(double)GREY_LEVELS); //1：对象0；2：对象1；3边界块
			}
		} 
		block_region--;
		
		if (block_region==current_region||block_region==2)
			break;
		else
		{
			sizej+=8;
			if(sizej>=width)
			{
				sizej=0;
				sizei+=8;
			}
			if ((sizej==0)&&(sizei>=height))
				break;
		}
	}
}

int Decode()
{
	BYTE i;
	int funcret;
	
	BYTE dcHufindex,acHufindex;
	
    num_byte=0;
	for (i=0;i<comp_num;i++)
	{
		
		sizei=sizej=0;
		coef=0; //pre_dc

		dcHufindex=(i==0?0:1);//huffman table
		acHufindex=(i==0?2:3);
		switch(i)
		{
		case 0:
			blockTemp=blockTemp_Y;
			pQt=qt_table[0];
			width=ImgWidth;
			height=ImgHeight;
			if(num_regions>1)
				plane=T_I->Y_Buffer_plane;
			break;
		case 1:
			blockTemp=blockTemp_UV[0];
			pQt=qt_table[1];
			width=ImgWidth/2;
			height=ImgHeight/2;	
			if(num_regions>1)
				plane=T_I->Cb_Buffer_plane;
			break;
		case 2:
			blockTemp=blockTemp_UV[1];
			pQt=qt_table[1];
			width=ImgWidth/2; 
			height=ImgHeight/2;
			if(num_regions>1)
				plane=T_I->Cr_Buffer_plane;
			break;
		}

		if (num_regions>1) //如果OB，先匹配出第一个块
			match_Block();
		while((funcret=DecodeBlock(dcHufindex,acHufindex))==FUNC_OK)
		{
			//把块匹配到对应正确位置
			IQtIZzBlock();
			
			sizej+=8;         ////此处有问题
			if(sizej>=width)
			{
				sizej=0;
				sizei+=8;
			}
			
			if(num_regions>1)
				match_Block();

			if ((sizej==0)&&(sizei>=height))
				break;
		}	
	}
	return funcret;	
}

void Initialize_Fast_IDCT()
{
	short i;
	
	iclp = iclip+512;
	for (i= -512; i<512; i++)
		iclp[i] = (i<-256) ? -256 : ((i>255) ? 255 : i);
}

void i_Frm_Decoder()
{
	int funcret,i,j;

	ImgHeight=imageHeight;
	ImgWidth=imageWidth;
	i_quality=i_quality;
	
	blockTemp_Y=(BYTE*)malloc(ImgHeight*ImgWidth*sizeof(BYTE));
	blockTemp_UV[0]=(BYTE*)malloc(ImgHeight*ImgWidth/4*sizeof(BYTE));
	blockTemp_UV[1]=(BYTE*)malloc(ImgHeight*ImgWidth/4*sizeof(BYTE));

	Initialize_Fast_IDCT();
	for(current_region=0;current_region<num_regions;current_region++){

		lp=T_I->fp_i_stream[current_region];
		init_Table();
		memset(blockTemp_Y,0,ImgHeight*ImgWidth);
	    memset(blockTemp_UV[0],0,ImgHeight*ImgWidth/4);
    	memset(blockTemp_UV[1],0,ImgHeight*ImgWidth/4);
		funcret=Decode();
		
		if(num_regions>1)
		{
			fwrite(blockTemp_Y,1,ImgWidth*ImgHeight,fp_out_region_rec[current_region]);
			fwrite(blockTemp_UV[0],1,ImgWidth*ImgHeight/4,fp_out_region_rec[current_region]);
			fwrite(blockTemp_UV[1],1,ImgWidth*ImgHeight/4,fp_out_region_rec[current_region]);
			for (i=0;i<ImgHeight;i++)
				for (j=0;j<ImgWidth;j++)
				{
					if(blockTemp_Y[i*ImgWidth+j]!=0)
						imgY_ref_temp[0][i*ImgWidth+j]=blockTemp_Y[i*ImgWidth+j];
				}
			
			for (i=0;i<ImgHeight/2;i++)
				for (j=0;j<ImgWidth/2;j++)
				{
					if(blockTemp_UV[0][i*ImgWidth/2+j]!=0)
						imgUV_ref_temp[0][0][i*ImgWidth/2+j]=blockTemp_UV[0][i*ImgWidth/2+j];
					if(blockTemp_UV[1][i*ImgWidth/2+j]!=0)
						imgUV_ref_temp[1][0][i*ImgWidth/2+j]=blockTemp_UV[1][i*ImgWidth/2+j];
				}
		}
		else
		{
			memcpy(imgY_ref_temp[0],blockTemp_Y,ImgWidth*ImgHeight);
			memcpy(imgUV_ref_temp[0][0],blockTemp_UV[0],ImgWidth*ImgHeight/4);
			memcpy(imgUV_ref_temp[1][0],blockTemp_UV[1],ImgWidth*ImgHeight/4);

			fwrite(blockTemp_Y,1,ImgWidth*ImgHeight,fp_out_all_rec);
			fwrite(blockTemp_UV[0],1,ImgWidth*ImgHeight/4,fp_out_all_rec);
			fwrite(blockTemp_UV[1],1,ImgWidth*ImgHeight/4,fp_out_all_rec);
		}
	}
	if (num_regions>1)
	{
		fwrite(imgY_ref_temp[0],1,ImgWidth*ImgHeight,fp_out_all_rec);
		fwrite(imgUV_ref_temp[0][0],1,ImgWidth*ImgHeight/4,fp_out_all_rec);
		fwrite(imgUV_ref_temp[1][0],1,ImgWidth*ImgHeight/4,fp_out_all_rec);
	}
}