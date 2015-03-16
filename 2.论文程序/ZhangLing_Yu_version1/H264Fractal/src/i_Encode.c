#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "i_codec_types.h"
#include "i_encode_globals.h"
#include <malloc.h>

// extern InputParameters *input;
// extern ImageParameters *image;
short current_region;
BYTE block_region;
FILE *fp_in,*fp_plane;

void writebyte(BYTE b,I_FRM_ENC_GLOBAL *TJE) 
{
	(*TJE->fp_i_stream[current_region]) = b;
	TJE->fp_i_stream[current_region]++;
	TJE->i_length[current_region] += 1;
}
void writeword(WORD w,I_FRM_ENC_GLOBAL *TJE) 
{
	writebyte((BYTE)(w/256),TJE);
	writebyte((BYTE)(w%256),TJE);
}


/*void write_DQTinfo(TJPEG_ENC_GLOBAL *TJE)
{
 BYTE i;
 writeword(TJE->DQTinfo.marker,TJE);
 writeword(TJE->DQTinfo.length,TJE);
 writebyte(TJE->DQTinfo.QTYinfo,TJE);for (i=0;i<64;i++) writebyte(TJE->DQTinfo.Ytable[i],TJE);
 writebyte(TJE->DQTinfo.QTCbinfo,TJE);for (i=0;i<64;i++) writebyte(TJE->DQTinfo.Cbtable[i],TJE);
}*/

/*
 *global valuable:	basic_table
 *argument:			scale_factor(控制质量的东东)
 *output:			newtable
 根据输入的scale_factor得到量化步长
 */
void set_quant_table(BYTE *basic_table,BYTE scale_factor,BYTE *newtable)// Set quantization table and zigzag reorder it
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

void set_DQTinfo(I_FRM_ENC_GLOBAL *t_i)
{
 BYTE scalefactor= t_i->quality;
 
 t_i->DQTinfo.QTYinfo=0;
 t_i->DQTinfo.QTCbinfo=1;
 set_quant_table(std_luminance_qt,scalefactor,t_i->DQTinfo.Ytable);
 set_quant_table(std_chrominance_qt,scalefactor,t_i->DQTinfo.Cbtable);
}

/*void write_DHTinfo(TJPEG_ENC_GLOBAL *TJE)
{
 BYTE i;
 writeword(TJE->DHTinfo.marker,TJE);
 writeword(TJE->DHTinfo.length,TJE);
 writebyte(TJE->DHTinfo.HTYDCinfo,TJE);
 for (i=0;i<16;i++)  writebyte(TJE->DHTinfo.YDC_nrcodes[i],TJE);
 for (i=0;i<=11;i++) writebyte(TJE->DHTinfo.YDC_values[i],TJE);
 writebyte(TJE->DHTinfo.HTYACinfo,TJE);
 for (i=0;i<16;i++)  writebyte(TJE->DHTinfo.YAC_nrcodes[i],TJE);
 for (i=0;i<=161;i++) writebyte(TJE->DHTinfo.YAC_values[i],TJE);
 writebyte(TJE->DHTinfo.HTCbDCinfo,TJE);
 for (i=0;i<16;i++)  writebyte(TJE->DHTinfo.CbDC_nrcodes[i],TJE);
 for (i=0;i<=11;i++)  writebyte(TJE->DHTinfo.CbDC_values[i],TJE);
 writebyte(TJE->DHTinfo.HTCbACinfo,TJE);
 for (i=0;i<16;i++)  writebyte(TJE->DHTinfo.CbAC_nrcodes[i],TJE);
 for (i=0;i<=161;i++) writebyte(TJE->DHTinfo.CbAC_values[i],TJE);
}*/

void set_DHTinfo(I_FRM_ENC_GLOBAL *t_i)
{
	BYTE i;
	
	t_i->DHTinfo.HTYDCinfo=0;
	for (i=0;i<16;i++)  
		t_i->DHTinfo.YDC_nrcodes[i]=std_dc_luminance_nrcodes[i+1];
	for (i=0;i<=11;i++)  
		t_i->DHTinfo.YDC_values[i]=std_dc_luminance_values[i];
	t_i->DHTinfo.HTYACinfo=0x10;
	for (i=0;i<16;i++)  
		t_i->DHTinfo.YAC_nrcodes[i]=std_ac_luminance_nrcodes[i+1];
	for (i=0;i<=161;i++) 
		t_i->DHTinfo.YAC_values[i]=std_ac_luminance_values[i];
	t_i->DHTinfo.HTCbDCinfo=1;
	for (i=0;i<16;i++)  
		t_i->DHTinfo.CbDC_nrcodes[i]=std_dc_chrominance_nrcodes[i+1];
	for (i=0;i<=11;i++)  
		t_i->DHTinfo.CbDC_values[i]=std_dc_chrominance_values[i];
	t_i->DHTinfo.HTCbACinfo=0x11;
	for (i=0;i<16;i++)  
		t_i->DHTinfo.CbAC_nrcodes[i]=std_ac_chrominance_nrcodes[i+1];
	for (i=0;i<=161;i++) 
		t_i->DHTinfo.CbAC_values[i]=std_ac_chrominance_values[i];
}


void writebits(bitstring bs,I_FRM_ENC_GLOBAL *TJE)// A portable version; it should be done in assembler
{
	WORD value;
	SBYTE posval;//bit position in the bitstring we read, should be<=15 and >=0
	value=bs.value;
	posval=bs.length-1;
	while (posval>=0)
	{
		if (value & mask[posval]) 
			TJE->bytenew|=mask[TJE->bytepos];
		posval--;TJE->bytepos--;
		if (TJE->bytepos<0) 
		{ 
			if (TJE->bytenew==0xFF) 
			{
				writebyte(0xFF,TJE);
				writebyte(0,TJE);
			}
			else 
				writebyte(TJE->bytenew,TJE);
			TJE->bytepos=7;TJE->bytenew=0;
		}
	}
}

void compute_Huffman_table(BYTE *nrcodes,BYTE *std_table,bitstring *HT)
{
 BYTE k,j;
 BYTE pos_in_table;
 WORD codevalue;
 codevalue=0; pos_in_table=0;
 for (k=1;k<=16;k++)
   {
     for (j=1;j<=nrcodes[k];j++) 
	 {
		HT[std_table[pos_in_table]].value=codevalue;
		HT[std_table[pos_in_table]].length=k;
		pos_in_table++;
		codevalue++;
	}
     codevalue*=2;
   }
}
void init_Huffman_tables(I_FRM_ENC_GLOBAL *t_i)
{
 compute_Huffman_table(std_dc_luminance_nrcodes,std_dc_luminance_values,t_i->YDC_HT);
 compute_Huffman_table(std_dc_chrominance_nrcodes,std_dc_chrominance_values,t_i->CbDC_HT);
 compute_Huffman_table(std_ac_luminance_nrcodes,std_ac_luminance_values,t_i->YAC_HT);
 compute_Huffman_table(std_ac_chrominance_nrcodes,std_ac_chrominance_values,t_i->CbAC_HT);
}

void set_numbers_category_and_bitcode(I_FRM_ENC_GLOBAL *t_i)
{
	SDWORD nr;
	SDWORD nrlower,nrupper;
	BYTE cat;
	
	t_i->category=t_i->category_alloc+32767; //allow negative subscripts  将地址偏移到中间
	
	t_i->bitcode=t_i->bitcode_alloc+32767;
	
	nrlower=1;nrupper=2;
	for (cat=1;cat<=15;cat++) 
	{//Positive numbers
		for (nr=nrlower;nr<nrupper;nr++)
		{
			t_i->category[nr]=cat;
			t_i->bitcode[nr].length=cat;
			t_i->bitcode[nr].value=(WORD)nr;
		}
		//Negative numbers
		for (nr=-(nrupper-1);nr<=-nrlower;nr++)
		{
			t_i->category[nr]=cat;
			t_i->bitcode[nr].length=cat;
			t_i->bitcode[nr].value=(WORD)(nrupper-1+nr);
		}
		nrlower<<=1;
		nrupper<<=1;
	}
}

void prepare_quant_tables(I_FRM_ENC_GLOBAL *t_i)
{
 BYTE row, col;
 BYTE i = 0;

 for (row = 0; row < 8; row++)
 {
   for (col = 0; col < 8; col++)
     {
       t_i->fdtbl_Y[i] = (float) (1.0 / ((double) t_i->DQTinfo.Ytable[zigzag[i]] *
			  aanscalefactor[row] * aanscalefactor[col] * 8.0));
	   
      t_i->fdtbl_Cb[i] = (float) (1.0 / ((double) t_i->DQTinfo.Cbtable[zigzag[i]] *
			  aanscalefactor[row] * aanscalefactor[col] * 8.0));
	   i++;
     }
 }
}

void fdct_and_quantization(SBYTE *data,float *fdtbl,SWORD *outdata,I_FRM_ENC_GLOBAL *t_i)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z1, z2, z3, z4, z5, z11, z13;
  float *dataptr;
  float datafloat[64];
  float temp;
  SBYTE ctr;
  BYTE i;
  for (i=0;i<64;i++) datafloat[i]=data[i];

  /* Pass 1: process rows. */
  dataptr=datafloat;
  for (ctr = 7; ctr >= 0; ctr--) {
	tmp0 = dataptr[0] + dataptr[7];
    tmp7 = dataptr[0] - dataptr[7];
    tmp1 = dataptr[1] + dataptr[6];
    tmp6 = dataptr[1] - dataptr[6];
    tmp2 = dataptr[2] + dataptr[5];
    tmp5 = dataptr[2] - dataptr[5];
    tmp3 = dataptr[3] + dataptr[4];
    tmp4 = dataptr[3] - dataptr[4];

	/* Even part */
    tmp10 = tmp0 + tmp3;	/* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    dataptr[0] = tmp10 + tmp11; /* phase 3 */
    dataptr[4] = tmp10 - tmp11;

    z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
	dataptr[2] = tmp13 + z1;	/* phase 5 */
    dataptr[6] = tmp13 - z1;

    /* Odd part */
    tmp10 = tmp4 + tmp5;	/* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

	/* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
    z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
    z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * ((float) 0.707106781); /* c4 */

    z11 = tmp7 + z3;		/* phase 5 */
    z13 = tmp7 - z3;

    dataptr[5] = z13 + z2;	/* phase 6 */
    dataptr[3] = z13 - z2;
	dataptr[1] = z11 + z4;
    dataptr[7] = z11 - z4;

    dataptr += 8;		/* advance pointer to next row */
  }

  /* Pass 2: process columns. */
  dataptr = datafloat;
  for (ctr = 7; ctr >= 0; ctr--) {
    tmp0 = dataptr[0] + dataptr[56];
    tmp7 = dataptr[0] - dataptr[56];
    tmp1 = dataptr[8] + dataptr[48];
    tmp6 = dataptr[8] - dataptr[48];
    tmp2 = dataptr[16] + dataptr[40];
    tmp5 = dataptr[16] - dataptr[40];
    tmp3 = dataptr[24] + dataptr[32];
    tmp4 = dataptr[24] - dataptr[32];

    /* Even part */
    tmp10 = tmp0 + tmp3;	/* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    dataptr[0] = tmp10 + tmp11; /* phase 3 */
    dataptr[32] = tmp10 - tmp11;

	z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
    dataptr[16] = tmp13 + z1; /* phase 5 */
    dataptr[48] = tmp13 - z1;

    /* Odd part */
    tmp10 = tmp4 + tmp5;	/* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
	z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
    z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
    z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * ((float) 0.707106781); /* c4 */

    z11 = tmp7 + z3;		/* phase 5 */
    z13 = tmp7 - z3;

    dataptr[40] = z13 + z2; /* phase 6 */
	dataptr[24] = z13 - z2;
    dataptr[8] = z11 + z4;
    dataptr[56] = z11 - z4;

    dataptr++;			/* advance pointer to next column */
  }

 for (i = 0; i < 64; i++) 
 {
   temp = datafloat[i] * fdtbl[i];
   outdata[i] = (SWORD) ((SWORD)(temp + 16384.5) - 16384);
 }
}


void process_YCBCR(SBYTE *ComponentDU,float *fdtbl,SWORD *DC,
				   bitstring *HTDC,bitstring *HTAC,I_FRM_ENC_GLOBAL *t_i)
{
	bitstring EOB=HTAC[0x00];
	bitstring M16zeroes=HTAC[0xF0];
	BYTE i;
	BYTE startpos;
	BYTE end0pos;
	BYTE nrzeroes;
	BYTE nrmarker;
	SWORD Diff;
	
	if ((input->num_regions>1)&&(block_region!=current_region&&block_region!=2))
		return;
	fdct_and_quantization(ComponentDU,fdtbl,t_i->DU_DCT,t_i);
	//zigzag reorder
	for (i=0;i<=63;i++) 
		t_i->DU[zigzag[i]]=t_i->DU_DCT[i];
	Diff=t_i->DU[0]-*DC;
	*DC=t_i->DU[0];
	///////////
	///////////////////
	//Encode DC
	if (Diff==0)
		writebits(HTDC[0],t_i); //Diff might be 0
	else {
		writebits(HTDC[t_i->category[Diff]],t_i);
		writebits(t_i->bitcode[Diff],t_i);
	}
	//Encode ACs
	for (end0pos=63;(end0pos>0)&&(t_i->DU[end0pos]==0);end0pos--)  //判断是否全为0，并找到最后不是0的数
		;
	//end0pos = first element in reverse order !=0
	if (end0pos==0) {writebits(EOB,t_i);return;}
	
	i=1;
	while (i<=end0pos)
	{
		startpos=i;
		for (; (t_i->DU[i]==0)&&(i<=end0pos);i++) 
			;
		nrzeroes=i-startpos;
		if (nrzeroes>=16) {
			for (nrmarker=1;nrmarker<=nrzeroes/16;nrmarker++) writebits(M16zeroes,t_i);
			nrzeroes=nrzeroes%16;
		}
		writebits(HTAC[nrzeroes*16+t_i->category[t_i->DU[i]]],t_i);
		writebits(t_i->bitcode[t_i->DU[i]],t_i);
		i++;
	}
	if (end0pos!=63) writebits(EOB,t_i);
}

void load_data_units_from_YCbCr_buffer(WORD xpos,WORD ypos,int Width,SBYTE *block,BYTE *Buffer)
{
	WORD x,y;
	WORD pos=0;
	DWORD location;
	
	location=ypos*Width+xpos;
	
	for (y=0;y<8;y++)
	{
		for (x=0;x<8;x++)
		{
			block[pos] = (SBYTE)Buffer[location]-128;
			//灰度图像，色度分量全为128，因此这里可简化	
			location++;
			pos++;
		}
		location+=Width-8;
	} 
}

void load_data_units_from_YCbCr_buffer_plane(WORD xpos,WORD ypos,int Width,BYTE *block,BYTE *Buffer)
{
	WORD x,y;
	WORD pos=0;
	DWORD location;
	
	location=ypos*Width+xpos;
	block_region=0;
	for (y=0;y<8;y++)
	{
		for (x=0;x<8;x++)
		{
			block[pos] = Buffer[location];
			block_region|=1<<(int)((double)block[pos]/(double)GREY_LEVELS); //1：对象0；2：对象1；3边界块
			location++;
			pos++;
		}
		location+=Width-8;
	} 
	block_region--;
	
}
void main_encoder(I_FRM_ENC_GLOBAL *t_i)
{
	SWORD DCY=0,DCCb=0,DCCr=0; //DC coefficients used for differential encoding
	WORD xpos,ypos;
	
	for (ypos=0;ypos<t_i->Height;ypos+=8)
		for (xpos=0;xpos<t_i->Width;xpos+=8)
		{
			load_data_units_from_YCbCr_buffer(xpos,ypos,t_i->Width,t_i->Y_block,t_i->Y_Buffer);
			if (input->num_regions>1)
				load_data_units_from_YCbCr_buffer_plane(xpos,ypos,t_i->Width,t_i->Y_block_plane,t_i->Y_Buffer_plane);
			process_YCBCR(t_i->Y_block,t_i->fdtbl_Y,&DCY,t_i->YDC_HT,t_i->YAC_HT,t_i);
		}	
	for (ypos=0;ypos<t_i->Height/2;ypos+=8)
		for (xpos=0;xpos<t_i->Width/2;xpos+=8)
		{
			load_data_units_from_YCbCr_buffer(xpos,ypos,t_i->Width/2,t_i->Cb_block,t_i->Cb_Buffer);
			if (input->num_regions>1)
				load_data_units_from_YCbCr_buffer_plane(xpos,ypos,t_i->Width/2,t_i->Cb_block_plane,t_i->Cb_Buffer_plane);
			process_YCBCR(t_i->Cb_block,t_i->fdtbl_Cb,&DCCb,t_i->CbDC_HT,t_i->CbAC_HT,t_i);
		}
			
	for (ypos=0;ypos<t_i->Height/2;ypos+=8)
		for (xpos=0;xpos<t_i->Width/2;xpos+=8)
		{
			load_data_units_from_YCbCr_buffer(xpos,ypos,t_i->Width/2,t_i->Cr_block,t_i->Cr_Buffer);
			if (input->num_regions>1)
				load_data_units_from_YCbCr_buffer_plane(xpos,ypos,t_i->Width/2,t_i->Cr_block_plane,t_i->Cr_Buffer_plane);
			process_YCBCR(t_i->Cr_block,t_i->fdtbl_Cb,&DCCr,t_i->CbDC_HT,t_i->CbAC_HT,t_i);
		}
}

void init_I_all(I_FRM_ENC_GLOBAL *t_i)
{
	BYTE i;
	
	t_i->bytenew=0; 		
	t_i->bytepos=7; //should be<=7 and >=0
	t_i->Width = imageWidth;//320;
	t_i->Height = imageHeight;//256;
	t_i->Y_Buffer  = (BYTE*)(malloc(t_i->Width*t_i->Height));
	t_i->Cb_Buffer = (BYTE*)(malloc(t_i->Width*t_i->Height/4));
	t_i->Cr_Buffer = (BYTE*)(malloc(t_i->Width*t_i->Height/4));

	t_i->Y_Buffer_plane  = (BYTE*)(malloc(t_i->Width*t_i->Height));
	t_i->Cb_Buffer_plane = (BYTE*)(malloc(t_i->Width*t_i->Height/4));
	t_i->Cr_Buffer_plane = (BYTE*)(malloc(t_i->Width*t_i->Height/4));
	
	t_i->quality = input->i_quality;
	//output
	t_i->fp_i_stream = (BYTE **)(malloc(input->num_regions*sizeof(BYTE*)));
	t_i->i_length = (WORD *)(malloc(input->num_regions*sizeof(WORD)));
	for(i=0;i<input->num_regions;i++){
		t_i->fp_i_stream[i] = (BYTE *)(malloc(t_i->Width*t_i->Height*3/2));
	    t_i->i_length[i] = 0;
	}
	//malloc info
	t_i->category_alloc=(BYTE *)malloc(65535*sizeof(BYTE));  //
    t_i->bitcode_alloc=(bitstring *)malloc(65535*sizeof(bitstring));
	set_DQTinfo(t_i);  //设置量化信息
	set_DHTinfo(t_i);  //设置huffman表信息
	init_Huffman_tables(t_i); //计算得到huffman表
	set_numbers_category_and_bitcode(t_i);
	prepare_quant_tables(t_i);
}

void free_I_all(I_FRM_ENC_GLOBAL* t_i)
{
	int i;
	free(t_i->Y_Buffer);
	free(t_i->Cb_Buffer);
	free(t_i->Cr_Buffer);
	
	free(t_i->Y_Buffer_plane);
	free(t_i->Cb_Buffer_plane);
	free(t_i->Cr_Buffer_plane);

	free(t_i->i_length);

	for(i=0;i<input->num_regions;i++){
		free(t_i->fp_i_stream[i]);
	}
	free(t_i->fp_i_stream);

	free(t_i->category_alloc);
	free(t_i->bitcode_alloc);
}

BYTE i_Frm_Encoder(I_FRM_ENC_GLOBAL *t_i)
{	
	bitstring fillbits; 
	BYTE *p_out;

	fseek(fp_in,img->current_frame*imageWidth*imageHeight*3/2,0);
	fread(t_i->Y_Buffer,1,t_i->Width*t_i->Height,fp_in);
	fread(t_i->Cb_Buffer,1,t_i->Width*t_i->Height/4,fp_in);
	fread(t_i->Cr_Buffer,1,t_i->Width*t_i->Height/4,fp_in);

	if (input->num_regions>1)
	{
		fseek(fp_plane,img->current_frame*imageWidth*imageHeight*3/2,0);
		fread(t_i->Y_Buffer_plane,1,t_i->Width*t_i->Height,fp_plane);
		fread(t_i->Cb_Buffer_plane,1,t_i->Width*t_i->Height/4,fp_plane);
		fread(t_i->Cr_Buffer_plane,1,t_i->Width*t_i->Height/4,fp_plane);	
	}

	for(current_region=0;current_region<input->num_regions;current_region++){
		p_out =t_i->fp_i_stream[current_region];  //p_out保存码流的起始地址
		
		main_encoder(t_i); 
		
		if (t_i->bytepos>=0) 
		{
			fillbits.length=t_i->bytepos+1;
			fillbits.value=(1<<(t_i->bytepos+1))-1;
			writebits(fillbits,t_i);
		}
		//	writeword(0xFFD9,t_i);  //结束标志
		t_i->fp_i_stream[current_region] = p_out;
	}
	return 1;
}
