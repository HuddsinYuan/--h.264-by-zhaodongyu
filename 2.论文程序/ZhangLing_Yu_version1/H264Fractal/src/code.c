#include "windows.h"//专门为了计算时间而定的
#include "code.h"
#include "i_encode_globals.h"
#include "i_decode_global.h"
#include "global.h"
#include "huffman.h"
#include "compute.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#if defined WIN32
#include <conio.h>
#endif
#include <assert.h>
#include "i_global.h"
#include "configfile.h"
#include "leaky_bucket.h"
#include "memalloc.h"
#include "mbuffer.h"
#include "intrarefresh.h"
#include "fmo.h"
#include "sei.h"
#include "parset.h"
#include "image.h"
#include "output.h"
#include "fast_me.h"
#include "ratectl.h"
/////////////////////////////////////>
#include "output.h"
/////////////////////////////////////<
#define JM      "8"
#define VERSION "8.6"
#include "defines_enc.h"
#include "cabac.h"
#include <annexb.h>
#include <rtp_.h>
InputParameters inputs, *input = &inputs;
ImageParameters images, *img   = &images;
StatParameters  stats,  *stat  = &stats;
SNRParameters   snrs,   *snr   = &snrs;
Decoders decoders, *decs=&decoders;
/////////////////////////////////////>
#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"
#define TRACEFILE   "trace_dec.txt"
struct inp_par_dec    *input_dec; 
extern StorablePicture* dec_picture;
int global_init_done = 0;
/////////////////////////////////////<
Boolean In2ndIGOP = FALSE1;
int    start_frame_no_in_this_IGOP = 0;
int    start_tr_in_this_IGOP = 0;
int    FirstFrameIn2ndIGOP=0;
int    cabac_encoding = 0;
extern ColocatedParams *Co_located;
I_FRM_ENC_GLOBAL I_FRM_ENC;
I_FRM_ENC_GLOBAL *T_I;
short*  spiral_search_x;
short*  spiral_search_y;
short*  spiral_hpel_search_x;
short*  spiral_hpel_search_y;
#ifdef _ADAPT_LAST_GROUP_
int initial_Bframes = 0;
#endif
short current_region;
void Init_Motion_Search_Module ();
void Clear_Motion_Search_Module ();

void write_picture(StorablePicture *p, FILE *p_out);
extern void DeblockFrame(ImageParameters *img, byte **, byte ***);

void PSNR();
void PSNR_object();
int main(int argc,char **argv)
{
	int i,j,y,x;//,k
	int M,N,n,np,nb;           //Rate control		
	DWORD start_I=0,end_I=0;
	time_t start_time, end_time;
	int currentframe;//当前帧的临时变量
	double time_diff;
	input = &inputs;
	img = &images;
    Configure(argc, argv);//得到配置文件（encoder.cfg的参数）
	search_mode=0;//选择搜索模式0：全搜索，1：新六边形搜索, 2:基于田隽师姐的最新EPZS ,3:六边形搜索
	img->frmHeightInMbs = input->imageheight/16;//  一帧的高度（以宏块为单位）
	img->frmWidthInMbs = input->imagewidth/16;//  一帧的宽度（以宏块为单位）
	img->frmSizeInMbs = img->frmWidthInMbs*img->frmHeightInMbs;//  一帧包含的宏块总数
	imageWidth=input->imagewidth;imageHeight=input->imageheight;//  图像宽度高度
	search_range = input->search_range;//  =7，配置文件中定义
	i_quality=input->i_quality;//  =70，配置文件中定义                        
	num_regions = input->num_regions;//  =1，配置文件中定义，与基于对象有关  
    Memoalloc();         //分配内存
	T_I=&I_FRM_ENC;//不懂啥意思？
	psnr_avg1=0;
	psnr_avg2=0;
	psnr_avg3=0;
	total_bit=0;
    ///////////////////////////out
	if ((snr =  (struct snr_par *)calloc(1, sizeof(struct snr_par)))==NULL) no_mem_exit("main: snr");
    ///////////////////////////////////////////////////////////////////////////////////////编码			
	p_dec = p_stat = p_log = p_trace = NULL;
	Configure_h264 (argc, argv);
	AllocNalPayloadBuffer();///////////////////////初始化NAL单元（分配NAL空间）
	init_poc();
	GenerateParameterSets();
	init_img();/////////////////////////////用合适的参数初试化图像
	frame_pic = malloc_picture();
	init_rdopt ();
	init_dpb(input);
	init_out_buffer();
	enc_picture = enc_frame_picture = enc_top_picture = enc_bottom_picture = NULL;
	init_global_buffers();//////////////////*按帧大小动态分配内存空间，分配的内存需要在free_global_buffers()函数中释放
	create_context_memory ();
	Init_Motion_Search_Module ();
	information_init();
	//Rate control 
	if(input->RCEnable)//no
		rc_init_seq();
	if(input->FMEnable)//是否使用快速运动估计，不太明白这个参数是怎么传过来的？
		DefineThreshold();
	// B pictures
	Bframe_ctr=0;
	tot_time=0;                 // time for total encoding session
    #ifdef _ADAPT_LAST_GROUP_
	if (input->last_frame > 0)//no
	input->no_frames_h264 = 1 + (input->last_frame + input->jumpd) / (input->jumpd + 1);//1
	initial_Bframes = input->successive_Bframe;//0
    #endif
 	PatchInputNoFrames();
	// Write sequence header (with parameter sets)
	stat->bit_ctr_parametersets = 0;
	stat->bit_slice = start_sequence();////returned 0///////////////////////////打开输出文件并产生合适的序列头 
	stat->bit_ctr_parametersets += stat->bit_ctr_parametersets_n;//176=0+176什么意思？？？
	start_frame_no_in_this_IGOP = 0;
///////////////////////////out
	

	img->frame_num =0;//全部设置为0
	for (img->current_frame=0;img->current_frame<input->no_frames;img->current_frame+=1)//循环每一帧编码
	{
		currentframe= img->current_frame;//当前帧号

		img->nal_reference_idc = 1;
		img->toppoc =(input->intra_period && input->idr_enable ? IMG_NUMBER % input->intra_period : IMG_NUMBER) * (2*(input->successive_Bframe+1)); //IMG_NUMBER=当前帧号，img->toppoc则等于当前帧号*2
		if ((input->PicInterlace==FRAME_CODING)&&(input->MbInterlace==FRAME_CODING))//yes//编码方式为帧编码
			img->bottompoc = img->toppoc;     //progressive//0
		else //no
			img->bottompoc = img->toppoc+1;   //hard coded
		img->framepoc = min (img->toppoc, img->bottompoc);//0
		img->delta_pic_order_cnt[0]=0;
		SetImgType();// 为I、P、SP图像设置图像类型
		if((img->current_frame%input->I_frame)==0)//---zl---I帧编码改进--->
		{   		
			UV_dct=0;			
///////////////////////////////Rate control/////////////////////////////////////////////////////
			// which layer the image belonged to?
			if ( IMG_NUMBER % (input->NumFramesInELSubSeq+1) == 0 )//yes
				img->layer = 0;
			else//no
				img->layer = 1;
			start_oneframe('C');
			encode_one_frame(); // 编码一个I帧
			write_picture(enc_picture,fp_out_all_c_rec);
			imgY_ref/*_c*/=enc_picture->imgY;//将enc_picture->imgY赋给参考图像
			imgUV_ref/*_c*/=enc_picture->imgUV;


			if (input->right==1)//右目
			{
// 				printf("R video processing...\n");
// 				init_dpb(input);///////ここです
				start_oneframe('C');
				compute_domain_Sum(); //C目的前一帧为参考帧domain
// 				start_oneframe('R');
// 				compute_domain_Sum(); //
				currentVideo = 'R';
				encode_oneframe(img->current_frame);
				write_picture(enc_picture,fp_out_all_r_rec);
				
				imgY_ref_r=enc_picture->imgY;
				imgUV_ref_r=enc_picture->imgUV;
				
			}
			
			if (input->left==1)//左目
			{
// 				printf("L video processing...\n");
// 				init_dpb(input);///////ここです
				start_oneframe('C');
				compute_domain_Sum(); //C目的前一帧为参考帧domain
// 				start_oneframe('L');
// 				compute_domain_Sum(); //
				currentVideo = 'L';
				encode_oneframe(img->current_frame);
				write_picture(enc_picture,fp_out_all_l_rec);
				imgY_ref_r=enc_picture->imgY;
				imgUV_ref_r=enc_picture->imgUV;
				
			}
			img->frame_num++;//帧号加1

		//<---zl---I帧编码改进---
               // terminate sequence
//          terminate_sequence();//////////////////////////////结束序列并关闭输出文件
//          flush_dpb();//dpb保存的是重建帧ZZZZZZZZZZZZZ

//		 write_picture(enc_picture,fp_out_264);

// 
////////////////////////////////////////////////////////////////////////////////////////////////////////编码
        }
		else//P帧分形编码
		{

			//the following is sent in the slice header
			img->delta_pic_order_cnt[0]=0;
/*
			memcpy(imgY_ref_h[0],imgY_ref[0],imageHeight*imageWidth);// imgY_ref_h 右目参考图像,水平像素差值
			memcpy(imgUV_ref_h[0][0],imgUV_ref[0][0],imageHeight*imageWidth/4);
			memcpy(imgUV_ref_h[1][0],imgUV_ref[1][0],imageHeight*imageWidth/4);
			memcpy(imgY_ref_m[0],imgY_ref[0],imageHeight*imageWidth);//右目参考图像，垂直像素差值
			memcpy(imgUV_ref_m[0][0],imgUV_ref[0][0],imageHeight*imageWidth/4);
			memcpy(imgUV_ref_m[1][0],imgUV_ref[1][0],imageHeight*imageWidth/4);
			memcpy(imgY_ref_n[0],imgY_ref[0],imageHeight*imageWidth);//右目参考图像，1/4像素差值
			memcpy(imgUV_ref_n[0][0],imgUV_ref[0][0],imageHeight*imageWidth/4);
			memcpy(imgUV_ref_n[1][0],imgUV_ref[1][0],imageHeight*imageWidth/4);
			for (i=0;i<imageHeight;i++)
				for (j=0;j<imageWidth;j++)
				{   
					if (j<imageWidth-1)
						imgY_ref_h[i][j]=(imgY_ref[i][j]+imgY_ref[i][j+1])/2;
					else
						imgY_ref_h[i][j]=imgY_ref[i][j]; 
				}
			for (i=0;i<imageHeight;i++)
				for (j=0;j<imageWidth;j++)
				{   
					if (i<imageHeight-1)
						imgY_ref_m[i][j]=(imgY_ref[i][j]+imgY_ref[i+1][j])/2;
					else
						imgY_ref_m[i][j]=imgY_ref[i][j]; 
				}
			for (i=0;i<imageHeight;i++)
				for(j=0;j<imageWidth;j++)
				{
					if ((i<imageHeight-1)&&(j<imageWidth-1))
						imgY_ref_n[i][j]=(imgY_ref[i][j]+imgY_ref[i+1][j]+imgY_ref[i][j+1]+imgY_ref[i+1][j+1])/4;
					else
						imgY_ref_n[i][j]=imgY_ref[i][j]; 
				}

 */
			start_oneframe('C');//设置参考domain
			compute_domain_Sum(); //计算domain的像素和，像素平方和，包括4*4，4*8，8*4，8*8，16*8，8*16，16*16
/* 
			start_oneframe('H');
			compute_domain_Sum(); //C目的前一帧为参考帧domain
			start_oneframe('M');
			compute_domain_Sum(); //C目的前一帧为参考帧domain
			start_oneframe('N');
			compute_domain_Sum(); //C目的前一帧为参考帧domain
 */
			currentVideo = 'C';
			encode_oneframe(img->current_frame);//分形编码一个P帧
			write_picture(enc_picture,fp_out_all_c_rec);
			imgY_ref/*_c*/=enc_picture->imgY;//将enc_picture->imgY赋给参考图像
			imgUV_ref/*_c*/=enc_picture->imgUV;

			if (input->right==1)//---zl---右目编码---
			{
 				printf("R video processing...\n");
// 				init_dpb(input);///////ここです
				start_oneframe('C');
				compute_domain_Sum(); //C目的前一帧为参考帧domain
//     			start_oneframe('R');
// 				compute_domain_Sum(); //
				currentVideo = 'R';
				encode_oneframe(img->current_frame);
// 				encode_one_frame();
				write_picture(enc_picture,fp_out_all_r_rec);

				imgY_ref_r=enc_picture->imgY;
				imgUV_ref_r=enc_picture->imgUV;

			}
	
			if (input->left==1)//---zl---左目编码---
			{
 				printf("L video processing...\n");
// 				init_dpb(input);///////ここです
				start_oneframe('C');
				compute_domain_Sum(); //C目的前一帧为参考帧domain
//      			start_oneframe('L');
// 				compute_domain_Sum(); //
				currentVideo = 'L';
				encode_oneframe(img->current_frame);
// 				encode_one_frame();
				write_picture(enc_picture,fp_out_all_l_rec);
				imgY_ref_r=enc_picture->imgY;
				imgUV_ref_r=enc_picture->imgUV;

			}

// 			printf("\n--------------------------------------------------------------------------------");


		}
        img->current_frame= currentframe;	
	}
    terminate_sequence();//////////////////////////////结束序列并关闭输出文件
//       flush_dpb();//dpb保存的是重建帧ZZZZZZZZZZZZZ

    fclose( fp_in );
 	if (p_dec)
	fclose( p_dec );
	if (p_trace)
	fclose( p_trace );
	Clear_Motion_Search_Module ();
	RandomIntraUninit();
	FmoUninit();
//  free structure for rd-opt. mode decision
	clear_rdopt ();
	#ifdef _LEAKYBUCKET_
	calc_buffer();
	#endif
//  report everything
	report();
// 	printf("P frame \n");


    //最后一个trans在写文件时，应添加pack(1,0,0,fp_out[i]);
//	printf("\naverage:Y %2.3f U %2.3f V %2.3f",psnr_avg1/input->no_frames,psnr_avg2/input->no_frames,psnr_avg3/input->no_frames);
  
    frame_rate = (float)img->framerate / ( (float) (input->jumpd + 1) );
    stat->bitrate= ((float) total_bit * frame_rate)/((float) input->no_frames );
    printf(" Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);
	printf("\ntotalbit: %d",total_bit);
	printf("\ntol_time: %d",tol_time);
//  	time(&end_time);
//  	time_diff = difftime(end_time, start_time);
// 	printf("\nTime elapsed : %d hours %d minutes %lf seconds\n\n\n",
// 				 (int)time_diff/3600, ((int)time_diff%3600)/60, fmod(time_diff,60));
	/*printf("\npsnr_ave= %f",psnr_ave/input->no_frames_h264);*/
	Memofree();
	return 0;
}

void write_CompInfo(compressionInfo info)
{
	info.frames = input->no_frames_h264;
	info.Width=input->imagewidth;
	info.Height=input->imageheight;
	info.stereo=0;
	info.stereo=input->right*1+input->left*2;//是否多目0：单目，1：右目，2：左目，3：3目
	info.regions = input->num_regions;
	info.i_quality=input->i_quality;
	info.i_interval=input->I_frame;
	info.search_range = input->search_range;
	if(input->num_regions>1)
	{
		fwrite(&info,sizeof(compressionInfo),1,fp_out_c[0]);
//		fwrite(&info,sizeof(compressionInfo),1,fp_out_c[1]);
	}
	else
	    fwrite(&info,sizeof(compressionInfo),1,fp_out_all_c);//所有参数都写在中目的bit中
}

void write_Istream()
{
	long filesizepre=0,filesizelas=0;
	if(input->num_regions==1)
	{
		filesizepre=pack(2,0,0,fp_out[0]);
		fwrite(&T_I->i_length[0],sizeof(WORD),1,fp_out[0]);
    	fwrite(T_I->fp_i_stream[0],1,T_I->i_length[0],fp_out[0]);
		filesizelas=pack(2,0,0,fp_out[0]);
		printf("  write_Istream:  compression ratio =  %0.3f  ",(float)(input->imagewidth*input->imageheight*3/2
			/(float)(filesizelas-filesizepre)));
		printf("\n\n");
	}
	else
	{
		BYTE reg;
		for (reg=0;reg<input->num_regions;reg++)
		{
			filesizepre=pack(2,0,0,fp_out[reg]);
			fwrite(&T_I->i_length[reg],sizeof(WORD),1,fp_out[reg]);
    	    fwrite(T_I->fp_i_stream[reg],1,T_I->i_length[reg],fp_out[reg]);
		    filesizelas=pack(2,0,0,fp_out[reg]);
			printf("\n  Object %d:\n",reg);
			printf("  write_Istream:  Compression ratio =  %0.3f  ",(float)(input->imagewidth*input->imageheight*3/2
			/(float)(filesizelas-filesizepre)));
		}
		printf("\n\n");
	}	
	
}

// void write_Codestream()
// {
// 	int i,j,n;
// 	long filesizepre=0,filesizelas=0;
// 	for(i=0;i<num_regions;i++)
// 	{
// 		filesizepre=pack(2,0,0,fp_out[i]);	//0，文件流指针
// 		pack(0, partition_length[i], LENGTH_QUAD_COUNT, fp_out[i]);//块划分的个数，=16
// 		for (j=0;j<partition_length[i];j++)    
// 		{
// 			pack(0,partition[i][j],2,fp_out[i]);
// 		}
// 		//////////////////////////////huffman 表
// 		for(j=0; j<(1<<huff_search_range); j++)
// 		{
// 			pack(0, Huff_ptr[i][0][j].bits, (int) ceil(log(MAX_CL)/log(2)),fp_out[i]);
// 			pack(0, Huff_ptr[i][0][j].code, Huff_ptr[i][0][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<(int)((MAX_ALPHA-MIN_ALPHA)*100)/5+1;j++)
// 		{
// 			pack(0, Huff_ptr[i][1][j].bits, (int) ceil(log(MAX_CL)/log(2)),fp_out[i]);
// 			pack(0, Huff_ptr[i][1][j].code, Huff_ptr[i][1][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<(MAX_BETA-MIN_BETA)/5+1;j++)
// 		{
// 			pack(0, Huff_ptr[i][2][j].bits, (int) ceil(log(MAX_CL)/log(2)),fp_out[i]);
// 			pack(0, Huff_ptr[i][2][j].code, Huff_ptr[i][2][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, Huff_ptr[i][3][j].bits, 1,fp_out[i]);
// 			pack(0, Huff_ptr[i][3][j].code, Huff_ptr[i][3][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, Huff_ptr[i][4][j].bits, 1,fp_out[i]);
// 			pack(0, Huff_ptr[i][4][j].code, Huff_ptr[i][4][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, Huff_ptr[i][5][j].bits, 1,fp_out[i]);
// 			pack(0, Huff_ptr[i][5][j].code, Huff_ptr[i][5][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, Huff_ptr[i][6][j].bits, 1,fp_out[i]);
// 			pack(0, Huff_ptr[i][6][j].code, Huff_ptr[i][6][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, Huff_ptr[i][7][j].bits, 1,fp_out[i]);
// 			pack(0, Huff_ptr[i][7][j].code, Huff_ptr[i][7][j].bits, fp_out[i]);
// 		}
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, Huff_ptr[i][8][j].bits, 1,fp_out[i]);
// 			pack(0, Huff_ptr[i][8][j].code, Huff_ptr[i][8][j].bits, fp_out[i]);
// 		}
// 
// 	////////////code stream
// 		for(j=0;j<10;j++)
// 		{
// 			pack(0, code_length[j][i], LENGTH_CODE_LENGTH, fp_out[i]);
// 			pack(0, trans_count[i], LENGTH_TRANS_COUNT, fp_out[i]);
// 			for(n=0; n<(code_length[j][i]-7); n+=8)
// 			{
// 				pack(0, code[j][i][0], 8, fp_out[i]);
// 				code[j][i]++;
// 			}
// 			if((code_length[j][i]&7)>0)//不够一个字节，只写相应的位数
// 				pack(0, code[j][i][0], (code_length[j][i]&7), fp_out[i]);
// 		}
// 		if (currentVideo!='C')
// 		{
// 			for (j=0;j<trans_count[i];j++)
// 			{
// 				pack(0,trans_reference[i][j],1,fp_out[i]);
// 			}
// 		}
// 		filesizelas=pack(2,0,0,fp_out[i]);
// 		pack(1,0,0,fp_out[i]);
// 		filesizelas=pack(2,0,0,fp_out[i]);
// 		printf("\n%2d th",img->current_frame);
// 		printf("       %0.3f ",(float)(input->imagewidth*input->imageheight*3/2/(float)(filesizelas-filesizepre)));
// 	}	
// }
// 
void PSNR_object()
{
	double psnr1;
	int i,j,k,l;
	psnr1=0.0;
	l=0;
	for (k=0;k<2;k++)
	{
	
	for(i=0;i<input->imageheight;i++)
		for(j=0;j<input->imagewidth;j++)
		{   
			if (plane_Y[i][j]/255==k)
		{
		    l++;
			psnr1+=((double)(imgY_rec[i][j])-(double)(imgY_org[i][j]))
			*((double)(imgY_rec[i][j])-(double)(imgY_org[i][j]));
		}
		}

	psnr1 = 10*log10((double)255*255*l/psnr1);
    printf("    PSNR_object: %0.3f",psnr1);
	l=0;
	}
}
void PSNR()
{
	double psnr1,psnr2,psnr3;
	int i,j,y,x;
	psnr1=0;
	for(i=0;i<input->imageheight;i++)
		for(j=0;j<input->imagewidth;j++)
			psnr1+=(enc_picture->imgY[i][j]-(imgY_org[i][j]))*(enc_picture->imgY[i][j]-(imgY_org[i][j]));
	psnr2=0;
	for(i=0;i<input->imageheight/2;i++)
		for(j=0;j<input->imagewidth/2;j++)
			psnr2+=((enc_picture->imgUV[0][i][j])-(imgUV_org[0][i][j]))*((enc_picture->imgUV[0][i][j])-(imgUV_org[0][i][j]));
	psnr3=0;
	for(i=0;i<input->imageheight/2;i++)
		for(j=0;j<input->imagewidth/2;j++)
			psnr3+=((enc_picture->imgUV[1][i][j])-(imgUV_org[1][i][j]))*((enc_picture->imgUV[1][i][j])-(imgUV_org[1][i][j]));

// 	snr->snr_y= psnr1 = 10*log10((double)255*255*input->imageheight*input->imagewidth/psnr1);
// 	snr->snr_u=	psnr2 = 10*log10((double)255*255*input->imageheight*input->imagewidth/psnr2);
// 	snr->snr_v=	psnr3 = 10*log10((double)255*255*input->imageheight*input->imagewidth/psnr3);

	frame_no=img->current_frame;
	psnr_avg1+=	psnr1;
	psnr_avg2+=	psnr2;
	psnr_avg3+=	psnr3;

// 	printf("\n%2d th",img->current_frame);
// 
//     printf("    %8d       %0.3f  %0.3f  %0.3f ",stat->bit_ctr - stat->bit_ctr_n,psnr1,psnr2,psnr3);
// 	printf("enc\n");
// 	for (y=0; y<4; y++)
// 	{
// 		for (x=0; x<4; x++)
// 		{
// 			if (x%2==0)
// 			{
// 				printf("%d(%d,%d)   ",enc_picture->imgY[y][x],enc_picture->imgUV[0][y/2][x/2],enc_picture->imgUV[1][y/2][x/2]);
// 			} 
// 			else if(x%4==3)
// 			{
// 				printf("%d\n",enc_picture->imgY[y][x]);
// 			}
// 			else
// 			{
// 				printf("%d   ",enc_picture->imgY[y][x]);
// 			}
// 		}
// 	}
//     printf("\n");
}
void Init_Motion_Search_Module ();
void Clear_Motion_Search_Module ();


/*!
***********************************************************************
* \brief
*    Terminates and reports statistics on error.
* 
***********************************************************************
*/
void report_stats_on_error()
{
	input->no_frames_h264=img->current_frame-1;
	terminate_sequence();
	
	flush_dpb();
	
	fclose(fp_in_c);
	if (p_dec)
		fclose(p_dec);
	if (p_trace)
		fclose(p_trace);
	
	Clear_Motion_Search_Module ();
	
	RandomIntraUninit();
	FmoUninit();
	
	// free structure for rd-opt. mode decision
	clear_rdopt ();
	
#ifdef _LEAKYBUCKET_
	calc_buffer();
#endif
	
	// report everything
	report();
	
	free_picture (frame_pic);
	if (top_pic)
		free_picture (top_pic);
	if (bottom_pic)
		free_picture (bottom_pic);
	
	free_dpb();
	free_collocated(Co_located);
	uninit_out_buffer();
	
	free_global_buffers();
	
	// free image mem
	free_img ();
	free_context_memory ();
	FreeNalPayloadBuffer();
	FreeParameterSets();
}

/*!
***********************************************************************
* \brief
*    Initializes the POC structure with appropriate parameters.
* 
***********************************************************************
*/
void init_poc()
{
	//  if(input->no_frames_h264 > (1<<15))error("too many frames",-998);
	
	//the following should probably go in sequence parameters
	// frame poc's increase by 2, field poc's by 1
	
	img->pic_order_cnt_type=input->pic_order_cnt_type;
	
	img->delta_pic_order_always_zero_flag=0;
	img->num_ref_frames_in_pic_order_cnt_cycle= 1;
	
	if (input->StoredBPictures)
	{
		img->offset_for_non_ref_pic  =  0;
		img->offset_for_ref_frame[0] =   2;
	}
	else
	{
		img->offset_for_non_ref_pic  =  -2*(input->successive_Bframe);
		img->offset_for_ref_frame[0] =   2*(input->successive_Bframe+1);
	}
	
	if ((input->PicInterlace==FRAME_CODING)&&(input->MbInterlace==FRAME_CODING))
		img->offset_for_top_to_bottom_field=0;
	else    
		img->offset_for_top_to_bottom_field=1;
	
	//the following should probably go in picture parameters
	//  img->pic_order_present_flag=0;    //img->delta_pic_order_cnt[1] not sent
	// POC200301
	if ((input->PicInterlace==FRAME_CODING)&&(input->MbInterlace==FRAME_CODING))
	{
		img->pic_order_present_flag=0;
		img->delta_pic_order_cnt_bottom = 0;
	}
	else    
	{
		img->pic_order_present_flag=1;
		img->delta_pic_order_cnt_bottom = 1;
	}
}


/*!
***********************************************************************
* \brief
*    Initializes the img->nz_coeff
* \par Input:
*    none
* \par  Output:
*    none
* \ side effects
*    sets omg->nz_coef[][][][] to -1
***********************************************************************
*/
void CAVLC_init()
{
	unsigned int i, k, l;
	
	for (i=0;i < img->PicSizeInMbs; i++)
		for (k=0;k<4;k++)
			for (l=0;l<6;l++)
				img->nz_coeff[i][k][l]=-1;
			
			
}


/*!
***********************************************************************
* \brief
*    Initializes the Image structure with appropriate parameters.
* \par Input:
*    Input Parameters struct inp_par_dec *inp
* \par  Output:
*    Image Parameters struct img_par *img
***********************************************************************
*/
void init_img()
{
	int i,j;
	
	img->num_reference_frames = active_sps->num_ref_frames;
	img->max_num_references   = active_sps->frame_mbs_only_flag ? active_sps->num_ref_frames : 2 * active_sps->num_ref_frames;
	
	img->buf_cycle = input->num_reference_frames;
	
	img->DeblockCall = 0;
    
	//  img->framerate=INIT_FRAME_RATE;   // The basic frame rate (of the original sequence)
	img->framerate=input->FrameRate;   // The basic frame rate (of the original sequence)
	
	
	get_mem_mv (&(img->pred_mv));//
	get_mem_mv (&(img->all_mv));//
	
	get_mem_ACcoeff (&(img->cofAC));
	get_mem_DCcoeff (&(img->cofDC));
	
	if(input->MbInterlace) 
	{
		get_mem_mv (&(rddata_top_frame_mb.pred_mv));
		get_mem_mv (&(rddata_top_frame_mb.all_mv));
		
		get_mem_mv (&(rddata_bot_frame_mb.pred_mv));
		get_mem_mv (&(rddata_bot_frame_mb.all_mv));
		
		get_mem_mv (&(rddata_top_field_mb.pred_mv));
		get_mem_mv (&(rddata_top_field_mb.all_mv));
		
		get_mem_mv (&(rddata_bot_field_mb.pred_mv));
		get_mem_mv (&(rddata_bot_field_mb.all_mv));
		
		get_mem_ACcoeff (&(rddata_top_frame_mb.cofAC));
		get_mem_DCcoeff (&(rddata_top_frame_mb.cofDC));
		
		get_mem_ACcoeff (&(rddata_bot_frame_mb.cofAC));
		get_mem_DCcoeff (&(rddata_bot_frame_mb.cofDC));
		
		get_mem_ACcoeff (&(rddata_top_field_mb.cofAC));
		get_mem_DCcoeff (&(rddata_top_field_mb.cofDC));
		
		get_mem_ACcoeff (&(rddata_bot_field_mb.cofAC));
		get_mem_DCcoeff (&(rddata_bot_field_mb.cofDC));
	}
	
	if ((img->quad = (int*)calloc (511, sizeof(int))) == NULL)
		no_mem_exit ("init_img: img->quad");
	img->quad+=255;
	for (i=0; i < 256; ++i)
	{
		img->quad[i]=img->quad[-i]=i*i;
	}
	
	img->width    = input->imagewidth;
	img->height   = input->imageheight;
	img->width_cr = input->imagewidth/2;
	img->height_cr= input->imageheight/2;
	
	img->PicWidthInMbs    = input->imagewidth/MB_BLOCK_SIZE;
	img->FrameHeightInMbs = input->imageheight/MB_BLOCK_SIZE;
	img->FrameSizeInMbs   = img->PicWidthInMbs * img->FrameHeightInMbs;
	
	img->PicHeightInMapUnits = ( active_sps->frame_mbs_only_flag ? img->FrameHeightInMbs : img->FrameHeightInMbs/2 );
	
	
	if(((img->mb_data) = (Macroblock *) calloc(img->FrameSizeInMbs,sizeof(Macroblock))) == NULL)
		no_mem_exit("init_img: img->mb_data");
	
	if(input->UseConstrainedIntraPred)
	{
		if(((img->intra_block) = (int*)calloc(img->FrameSizeInMbs,sizeof(int))) == NULL)
			no_mem_exit("init_img: img->intra_block");
	}
	
	get_mem2Dint(&(img->ipredmode), img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);        //need two extra rows at right and bottom
	
	if(input->MbInterlace) 
	{
		get_mem2Dint(&(rddata_top_frame_mb.ipredmode), img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);
		get_mem2Dint(&(rddata_bot_frame_mb.ipredmode), img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);
		get_mem2Dint(&(rddata_top_field_mb.ipredmode), img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);
		get_mem2Dint(&(rddata_bot_field_mb.ipredmode), img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);
	}
	// CAVLC mem
	get_mem3Dint(&(img->nz_coeff), img->FrameSizeInMbs, 4, 6);
	
	CAVLC_init();
	
	for (i=0; i < img->width/BLOCK_SIZE; i++)
		for (j=0; j < img->height/BLOCK_SIZE; j++)
		{
			img->ipredmode[i][j]=-1;
		}
		
		img->mb_y_upd=0;
		
		RandomIntraInit (img->width/16, img->height/16, input->RandomIntraMBRefresh);
		
// 		InitSEIMessages();  // Tian Dong (Sept 2002)
		
		// Initialize filtering parameters. If sending parameters, the offsets are 
		// multiplied by 2 since inputs are taken in "div 2" format.
		// If not sending paramters, all fields are cleared 
		if (input->LFSendParameters)
		{
			input->LFAlphaC0Offset <<= 1;
			input->LFBetaOffset <<= 1;
		}
		else
		{
			input->LFDisableIdc = 0;
			input->LFAlphaC0Offset = 0;
			input->LFBetaOffset = 0;
		}
}

/*!
***********************************************************************
* \brief
*    Free the Image structures
* \par Input:
*    Image Parameters struct img_par *img
***********************************************************************
*/
void free_img ()
{
// 	CloseSEIMessages(); // Tian Dong (Sept 2002)
	free_mem_mv (img->pred_mv);
	free_mem_mv (img->all_mv);
	
	free_mem_ACcoeff (img->cofAC);
	free_mem_DCcoeff (img->cofDC);
	
	free (img->quad-255);
}



/*!
************************************************************************
* \brief
*    Allocates the picture structure along with its dependent
*    data structures
* \return
*    Pointer to a Picture
************************************************************************
*/

Picture *malloc_picture()
{
	Picture *pic;
	if ((pic = calloc (1, sizeof (Picture))) == NULL) no_mem_exit ("malloc_picture: Picture structure");
	//! Note: slice structures are allocated as needed in code_a_picture
	return pic;
}

/*!
************************************************************************
* \brief
*    Frees a picture
* \param
*    pic: POinter to a Picture to be freed
************************************************************************
*/


void free_picture(Picture *pic)
{
	if (pic != NULL)
	{
		free_slice_list(pic);
		free (pic);
	}
}




/*!
************************************************************************
* \brief
*    Reports the gathered information to appropriate outputs
* \par Input:
*    struct inp_par_dec *inp,                                            \n
*    struct img_par *img,                                            \n
*    struct stat_par *stat,                                          \n
*    struct stat_par *stat                                           \n
*
* \par Output:
*    None
************************************************************************
*/
void report()
{
	int bit_use[NUM_PIC_TYPE][2] ;
	int i,j;
	char name[20];
	int total_bits=0;
	float frame_rate;
	float mean_motion_info_bit_use[2];
	
#ifndef WIN32
	time_t now;
	struct tm *l_time;
	char string[1000];
#else
	char timebuf[128];
#endif
	bit_use[I_SLICE][0] = 1;
	bit_use[P_SLICE][0] = max(1,input->no_frames_h264-1);
	bit_use[B_SLICE][0] = Bframe_ctr;
	
	//  Accumulate bit usage for inter and intra frames
	for (j=0;j<NUM_PIC_TYPE;j++)
	{
		bit_use[j][1] = 0;
	}
	
	for (j=0;j<NUM_PIC_TYPE;j++)
	{
		for(i=0; i<MAXMODE; i++)
			bit_use[j][1] += stat->bit_use_mode    [j][i]; 
		
		bit_use[j][1]+=stat->bit_use_header      [j];
		bit_use[j][1]+=stat->bit_use_mb_type     [j];
		bit_use[j][1]+=stat->tmp_bit_use_cbp     [j];
		bit_use[j][1]+=stat->bit_use_coeffY      [j];
		bit_use[j][1]+=stat->bit_use_coeffC      [j];
		bit_use[j][1]+=stat->bit_use_delta_quant [j];
		bit_use[j][1]+=stat->bit_use_stuffingBits[j];
	}
	
	// B pictures
	if(Bframe_ctr!=0)
	{
		frame_rate = (float)(img->framerate *(input->successive_Bframe + 1)) / (float) (input->jumpd+1);
		
		//! Currently adding NVB bits on P rate. Maybe additional stat info should be created instead and added in log file
		//    stat->bitrate_P=(stat->bit_ctr_0+stat->bit_ctr_P)*(float)(frame_rate)/(float) (input->no_frames_h264 + Bframe_ctr);
		stat->bitrate_P=(stat->bit_ctr_0+stat->bit_ctr_P + stat->bit_ctr_parametersets)*(float)(frame_rate)/(float) (input->no_frames_h264 + Bframe_ctr);
		
#ifdef _ADAPT_LAST_GROUP_
		stat->bitrate_B=(stat->bit_ctr_B)*(float)(frame_rate)/(float) (input->no_frames_h264 + Bframe_ctr);    
#else
		stat->bitrate_B=(stat->bit_ctr_B)*(float)(frame_rate)/(float) (input->no_frames_h264 + Bframe_ctr);    
#endif
		
	}
	else
	{
		if (input->no_frames_h264 > 1)
		{
			stat->bitrate=(bit_use[I_SLICE][1]+bit_use[P_SLICE][1])*(float)img->framerate/(input->no_frames_h264*(input->jumpd+1));
		}
	}
	
// 	fprintf(stdout,"-------------------------------------------------------------------------------\n");ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 	fprintf(stdout,   " Freq. for encoded bitstream       : %1.0f\n",(float)img->framerate/(float)(input->jumpd+1));
// 	if(input->hadamard)
// 		fprintf(stdout," Hadamard transform                : Used\n");
// 	else
// 		fprintf(stdout," Hadamard transform                : Not used\n");
// 	
// 	fprintf(stdout," Image format                      : %dx%d\n",input->imagewidth,input->imageheight);
// 	
// 	if(input->intra_upd)
// 		fprintf(stdout," Error robustness                  : On\n");
// 	else
// 		fprintf(stdout," Error robustness                  : Off\n");
// 	fprintf(stdout,    " Search range                      : %d\n",input->search_range);
// 	
// 	
// 	fprintf(stdout,   " No of ref. frames used in P pred  : %d\n",input->num_reference_frames);
// 	if(input->successive_Bframe != 0)
// 		fprintf(stdout,   " No of ref. frames used in B pred  : %d\n",input->num_reference_frames);
// 	fprintf(stdout,   " Total encoding time for the seq.  : %.3f sec \n",tot_time*0.001);
// 	fprintf(stdout,   " Total ME time for sequence        : %.3f sec \n",me_tot_time*0.001);
// 	
// 	// B pictures
// 	fprintf(stdout, " Sequence type                     :" );
// 	if(input->successive_Bframe==0 && input->StoredBPictures > 0)
// 		fprintf(stdout, " I-P-BS-BS (QP: I %d, P BS %d) \n",input->qp0, input->qpN);
// 	else if(input->successive_Bframe==1 && input->StoredBPictures > 0)
// 		fprintf(stdout, " I-B-P-BS-B-BS (QP: I %d, P BS %d, B %d) \n",input->qp0, input->qpN, input->qpB);
// 	else if(input->successive_Bframe==2 && input->StoredBPictures > 0)
// 		fprintf(stdout, " I-B-B-P-B-B-BS (QP: I %d, P BS %d, B %d) \n",input->qp0, input->qpN, input->qpB);
// 	else if(input->successive_Bframe==1)   fprintf(stdout, " IBPBP (QP: I %d, P %d, B %d) \n",
// 		input->qp0, input->qpN, input->qpB);
// 	else if(input->successive_Bframe==2) fprintf(stdout, " IBBPBBP (QP: I %d, P %d, B %d) \n",
// 		input->qp0, input->qpN, input->qpB);
// 	else if(input->successive_Bframe==0 && input->sp_periodicity==0) fprintf(stdout, " IPPP (QP: I %d, P %d) \n",   input->qp0, input->qpN);
// 	
// 	else fprintf(stdout, " I-P-P-SP-P (QP: I %d, P %d, SP (%d, %d)) \n",  input->qp0, input->qpN, input->qpsp, input->qpsp_pred);
// 	
// 	// report on entropy coding  method
// 	if (input->symbol_mode == UVLC)
// 		fprintf(stdout," Entropy coding method             : CAVLC\n");
// 	else
// 		fprintf(stdout," Entropy coding method             : CABAC\n");
// 	
// 	fprintf(stdout,  " Profile/Level IDC                 : (%d,%d)\n",input->ProfileIDC,input->LevelIDC);
// #ifdef _FULL_SEARCH_RANGE_
// 	if (input->full_search == 2)
// 		fprintf(stdout," Search range restrictions         : none\n");
// 	else if (input->full_search == 1)
// 		fprintf(stdout," Search range restrictions         : older reference frames\n");
// 	else
// 		fprintf(stdout," Search range restrictions         : smaller blocks and older reference frames\n");
// #endif
// 	
// 	if (input->rdopt)
// 		fprintf(stdout," RD-optimized mode decision        : used\n");
// 	else
// 		fprintf(stdout," RD-optimized mode decision        : not used\n");
// 	
// 	switch(input->partition_mode)
// 	{
// 	case PAR_DP_1:
// 		fprintf(stdout," Data Partitioning Mode            : 1 partition \n");
// 		break;
// 	case PAR_DP_3:
// 		fprintf(stdout," Data Partitioning Mode            : 3 partitions \n");
// 		break;
// 	default:
// 		fprintf(stdout," Data Partitioning Mode            : not supported\n");
// 		break;
// 	}
// 	
// 	switch(input->of_mode)//"OutFileMode"
// 	{
// 	case PAR_OF_ANNEXB:
// 		fprintf(stdout," Output File Format                : H.264 Bit Stream File Format \n");
// 		break;
// 	case PAR_OF_RTP:
// 		fprintf(stdout," Output File Format                : RTP Packet File Format \n");
// 		break;
// 	default:
// 		fprintf(stdout," Output File Format                : not supported\n");
// 		break;
// 	}
// 	
// 	
// 	
	fprintf(stdout,"------------------ Average data all frames  -----------------------------------\n");//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	fprintf(stdout," SNR Y(dB)                         : %5.2f\n",snr->snr_ya);
	fprintf(stdout," SNR U(dB)                         : %5.2f\n",snr->snr_ua);
	fprintf(stdout," SNR V(dB)                         : %5.2f\n",snr->snr_va);
	
	if(Bframe_ctr!=0)
	{
		
		//      fprintf(stdout, " Total bits                        : %d (I %5d, P %5d, B %d) \n",
		//            total_bits=stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_B, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_B);
// 		fprintf(stdout, " Total bits                        : %d (I %5d, P %5d, B %d NVB %d) \n",ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
             total_bits=stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_B + stat->bit_ctr_parametersets;//, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_B,stat->bit_ctr_parametersets);
		
		frame_rate = (float)(img->framerate *(input->successive_Bframe + 1)) / (float) (input->jumpd+1);
		stat->bitrate= ((float) total_bits * frame_rate)/((float) (input->no_frames_h264 + Bframe_ctr));
		
// 		fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
		
	}
	else if (input->sp_periodicity==0)
	{
		//    fprintf(stdout, " Total bits                        : %d (I %5d, P %5d) \n",
		//    total_bits=stat->bit_ctr_P + stat->bit_ctr_0 , stat->bit_ctr_0, stat->bit_ctr_P);
// 		fprintf(stdout, " Total bits                        : %d (I %5d, P %5d, NVB %d) \n",ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 			total_bits=stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_parametersets, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_parametersets);
		
		
		frame_rate = (float)img->framerate / ( (float) (input->jumpd + 1) );
		stat->bitrate= ((float) total_bits * frame_rate)/((float) input->no_frames_h264 );
		
// 		fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	}else
	{
		//fprintf(stdout, " Total bits                        : %d (I %5d, P %5d) \n",
		//total_bits=stat->bit_ctr_P + stat->bit_ctr_0 , stat->bit_ctr_0, stat->bit_ctr_P);
// 		fprintf(stdout, " Total bits                        : %d (I %5d, P %5d, NVB %d) \n",ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 			total_bits=stat->bit_ctr_P + stat->bit_ctr_0 + stat->bit_ctr_parametersets, stat->bit_ctr_0, stat->bit_ctr_P, stat->bit_ctr_parametersets);
		
		
		frame_rate = (float)img->framerate / ( (float) (input->jumpd + 1) );
		stat->bitrate= ((float) total_bits * frame_rate)/((float) input->no_frames_h264 );
		
// 		fprintf(stdout, " Bit rate (kbit/s)  @ %2.2f Hz     : %5.2f\n", frame_rate, stat->bitrate/1000);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	}
	
// 	fprintf(stdout, " Bits to avoid Startcode Emulation : %d \n", stat->bit_ctr_emulationprevention);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 	fprintf(stdout, " Bits for parameter sets           : %d \n", stat->bit_ctr_parametersets);
// 	
// 	fprintf(stdout,"-------------------------------------------------------------------------------\n");
// 	fprintf(stdout,"Exit JM %s encoder ver %s ", JM, VERSION);
// 	fprintf(stdout,"\n");
	
	// status file
	if ((p_stat=fopen("stat.dat","wt"))==0)
	{
		printf(errortext, ET_SIZE, "Error open file %s", "stat.dat");
		error(errortext, 500);
	}
	fprintf(p_stat," -------------------------------------------------------------- \n");
	fprintf(p_stat,"  This file contains statistics for the last encoded sequence   \n");
	fprintf(p_stat," -------------------------------------------------------------- \n");
	fprintf(p_stat,   " Sequence                     : %s\n",input->infile);
	fprintf(p_stat,   " No.of coded pictures         : %4d\n",input->no_frames_h264+Bframe_ctr);
	fprintf(p_stat,   " Freq. for encoded bitstream  : %4.0f\n",frame_rate);
	
	// B pictures
	if(input->successive_Bframe != 0)
	{
		fprintf(p_stat,   " BaseLayer Bitrate(kb/s)      : %6.2f\n", stat->bitrate_P/1000);
		fprintf(p_stat,   " EnhancedLayer Bitrate(kb/s)  : %6.2f\n", stat->bitrate_B/1000);
	}
	else
		fprintf(p_stat,   " Bitrate(kb/s)                : %6.2f\n", stat->bitrate/1000);
	
	if(input->hadamard)
		fprintf(p_stat," Hadamard transform           : Used\n");
	else
		fprintf(p_stat," Hadamard transform           : Not used\n");
	
	fprintf(p_stat,  " Image format                 : %dx%d\n",input->imagewidth,input->imageheight);
	
	if(input->intra_upd)
		fprintf(p_stat," Error robustness             : On\n");
	else
		fprintf(p_stat," Error robustness             : Off\n");
	
	fprintf(p_stat,  " Search range                 : %d\n",input->search_range);
	
	fprintf(p_stat,   " No of frame used in P pred   : %d\n",input->num_reference_frames);
	if(input->successive_Bframe != 0)
		fprintf(p_stat, " No of frame used in B pred   : %d\n",input->num_reference_frames);
	
	if (input->symbol_mode == UVLC)
		fprintf(p_stat,   " Entropy coding method        : CAVLC\n");
	else
		fprintf(p_stat,   " Entropy coding method        : CABAC\n");
	
    fprintf(p_stat,   " Profile/Level IDC            : (%d,%d)\n",input->ProfileIDC,input->LevelIDC);
	if(input->MbInterlace)
		fprintf(p_stat, " MB Field Coding : On \n");
	
#ifdef _FULL_SEARCH_RANGE_
	if (input->full_search == 2)
		fprintf(p_stat," Search range restrictions    : none\n");
	else if (input->full_search == 1)
		fprintf(p_stat," Search range restrictions    : older reference frames\n");
	else
		fprintf(p_stat," Search range restrictions    : smaller blocks and older reference frames\n");
#endif
	if (input->rdopt)
		fprintf(p_stat," RD-optimized mode decision   : used\n");
	else
		fprintf(p_stat," RD-optimized mode decision   : not used\n");
	
	fprintf(p_stat," ---------------------|----------------|---------------|\n");
	fprintf(p_stat,"     Item             |      Intra     |   All frames  |\n");
	fprintf(p_stat," ---------------------|----------------|---------------|\n");
	fprintf(p_stat," SNR Y(dB)            |");
	fprintf(p_stat,"  %5.2f         |",snr->snr_y1);
	fprintf(p_stat," %5.2f         |\n",snr->snr_ya);
	fprintf(p_stat," SNR U/V (dB)         |");
	fprintf(p_stat,"  %5.2f/%5.2f   |",snr->snr_u1,snr->snr_v1);
	fprintf(p_stat," %5.2f/%5.2f   |\n",snr->snr_ua,snr->snr_va);
	
	// QUANT.
	fprintf(p_stat," Average quant        |");
	fprintf(p_stat,"  %5d         |",absm(input->qp0));
	fprintf(p_stat," %5.2f         |\n",(float)stat->quant1/max(1.0,(float)stat->quant0));
	
	// MODE
	fprintf(p_stat,"\n ---------------------|----------------|\n");
	fprintf(p_stat,"   Intra              |    Mode used   |\n");
	fprintf(p_stat," ---------------------|----------------|\n");
	
	fprintf(p_stat," Mode 0  intra 4x4    |  %5d         |\n",stat->mode_use[I_SLICE][I4MB]);
	fprintf(p_stat," Mode 1+ intra 16x16  |  %5d         |\n",stat->mode_use[I_SLICE][I16MB]);
	
	fprintf(p_stat,"\n ---------------------|----------------|-----------------|\n");
	fprintf(p_stat,"   Inter              |    Mode used   | MotionInfo bits |\n");
	fprintf(p_stat," ---------------------|----------------|-----------------|");
	fprintf(p_stat,"\n Mode  0  (copy)      |  %5d         |    %8.2f     |",stat->mode_use[P_SLICE][0   ],(float)stat->bit_use_mode[P_SLICE][0   ]/(float)bit_use[P_SLICE][0]);
	fprintf(p_stat,"\n Mode  1  (16x16)     |  %5d         |    %8.2f     |",stat->mode_use[P_SLICE][1   ],(float)stat->bit_use_mode[P_SLICE][1   ]/(float)bit_use[P_SLICE][0]);
	fprintf(p_stat,"\n Mode  2  (16x8)      |  %5d         |    %8.2f     |",stat->mode_use[P_SLICE][2   ],(float)stat->bit_use_mode[P_SLICE][2   ]/(float)bit_use[P_SLICE][0]);
	fprintf(p_stat,"\n Mode  3  (8x16)      |  %5d         |    %8.2f     |",stat->mode_use[P_SLICE][3   ],(float)stat->bit_use_mode[P_SLICE][3   ]/(float)bit_use[P_SLICE][0]);
	fprintf(p_stat,"\n Mode  4  (8x8)       |  %5d         |    %8.2f     |",stat->mode_use[P_SLICE][P8x8],(float)stat->bit_use_mode[P_SLICE][P8x8]/(float)bit_use[P_SLICE][0]);
	fprintf(p_stat,"\n Mode  5  intra 4x4   |  %5d         |-----------------|",stat->mode_use[P_SLICE][I4MB]);
	fprintf(p_stat,"\n Mode  6+ intra 16x16 |  %5d         |",stat->mode_use[P_SLICE][I16MB]);
	mean_motion_info_bit_use[0] = (float)(stat->bit_use_mode[P_SLICE][0] + stat->bit_use_mode[P_SLICE][1] + stat->bit_use_mode[P_SLICE][2] 
		+ stat->bit_use_mode[P_SLICE][3] + stat->bit_use_mode[P_SLICE][P8x8])/(float) bit_use[P_SLICE][0]; 
	
	// B pictures
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
	{
		
		fprintf(p_stat,"\n\n ---------------------|----------------|-----------------|\n");
		fprintf(p_stat,"   B frame            |    Mode used   | MotionInfo bits |\n");
		fprintf(p_stat," ---------------------|----------------|-----------------|");
		fprintf(p_stat,"\n Mode  0  (copy)      |  %5d         |    %8.2f     |",stat->mode_use[B_SLICE][0   ],(float)stat->bit_use_mode[B_SLICE][0   ]/(float)Bframe_ctr);
		fprintf(p_stat,"\n Mode  1  (16x16)     |  %5d         |    %8.2f     |",stat->mode_use[B_SLICE][1   ],(float)stat->bit_use_mode[B_SLICE][1   ]/(float)Bframe_ctr);
		fprintf(p_stat,"\n Mode  2  (16x8)      |  %5d         |    %8.2f     |",stat->mode_use[B_SLICE][2   ],(float)stat->bit_use_mode[B_SLICE][2   ]/(float)Bframe_ctr);
		fprintf(p_stat,"\n Mode  3  (8x16)      |  %5d         |    %8.2f     |",stat->mode_use[B_SLICE][3   ],(float)stat->bit_use_mode[B_SLICE][3   ]/(float)Bframe_ctr);
		fprintf(p_stat,"\n Mode  4  (8x8)       |  %5d         |    %8.2f     |",stat->mode_use[B_SLICE][P8x8],(float)stat->bit_use_mode[B_SLICE][P8x8]/(float)Bframe_ctr);
		fprintf(p_stat,"\n Mode  5  intra 4x4   |  %5d         |-----------------|",stat->mode_use[B_SLICE][I4MB]);
		fprintf(p_stat,"\n Mode  6+ intra 16x16 |  %5d         |",stat->mode_use[B_SLICE][I16MB]);
		mean_motion_info_bit_use[1] = (float)(stat->bit_use_mode[B_SLICE][0] + stat->bit_use_mode[B_SLICE][1] + stat->bit_use_mode[B_SLICE][2] 
			+ stat->bit_use_mode[B_SLICE][3] + stat->bit_use_mode[B_SLICE][P8x8])/(float) Bframe_ctr; 
		
	}
	
	fprintf(p_stat,"\n\n ---------------------|----------------|----------------|----------------|\n");
	fprintf(p_stat,"  Bit usage:          |      Intra     |      Inter     |    B frame     |\n");
	fprintf(p_stat," ---------------------|----------------|----------------|----------------|\n");
	
	fprintf(p_stat," Header               |");
	fprintf(p_stat," %10.2f     |",(float) stat->bit_use_header[I_SLICE]/bit_use[I_SLICE][0]);
	fprintf(p_stat," %10.2f     |",(float) stat->bit_use_header[P_SLICE]/bit_use[P_SLICE][0]);
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," %10.2f      |",(float) stat->bit_use_header[B_SLICE]/Bframe_ctr);
	else fprintf(p_stat," %10.2f      |", 0.);
	fprintf(p_stat,"\n");
	
	fprintf(p_stat," Mode                 |");
	fprintf(p_stat," %10.2f     |",(float)stat->bit_use_mb_type[I_SLICE]/bit_use[I_SLICE][0]);
	fprintf(p_stat," %10.2f     |",(float)stat->bit_use_mb_type[P_SLICE]/bit_use[P_SLICE][0]);
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," %10.2f     |",(float)stat->bit_use_mb_type[B_SLICE]/Bframe_ctr);
	else fprintf(p_stat," %10.2f     |", 0.);
	fprintf(p_stat,"\n");
	
	fprintf(p_stat," Motion Info          |");
	fprintf(p_stat,"        ./.     |");
	fprintf(p_stat," %10.2f     |",mean_motion_info_bit_use[0]);
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," %10.2f     |",mean_motion_info_bit_use[1]);
	else fprintf(p_stat," %10.2f     |", 0.);
	fprintf(p_stat,"\n");
	
	fprintf(p_stat," CBP Y/C              |");
	fprintf(p_stat," %10.2f     |", (float)stat->tmp_bit_use_cbp[I_SLICE]/bit_use[I_SLICE][0]);
	fprintf(p_stat," %10.2f     |", (float)stat->tmp_bit_use_cbp[P_SLICE]/bit_use[P_SLICE][0]);
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," %10.2f     |", (float)stat->tmp_bit_use_cbp[B_SLICE]/Bframe_ctr);
	else fprintf(p_stat," %10.2f     |", 0.);
	fprintf(p_stat,"\n");
	
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," Coeffs. Y            | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_coeffY[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_coeffY[P_SLICE]/bit_use[P_SLICE][0], (float)stat->bit_use_coeffY[B_SLICE]/Bframe_ctr);
	else
		fprintf(p_stat," Coeffs. Y            | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_coeffY[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_coeffY[P_SLICE]/(float)bit_use[P_SLICE][0], 0.);
	
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," Coeffs. C            | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_coeffC[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_coeffC[P_SLICE]/bit_use[P_SLICE][0], (float)stat->bit_use_coeffC[B_SLICE]/Bframe_ctr);
	else
		fprintf(p_stat," Coeffs. C            | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_coeffC[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_coeffC[P_SLICE]/bit_use[P_SLICE][0], 0.);
	
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," Delta quant          | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_delta_quant[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_delta_quant[P_SLICE]/bit_use[P_SLICE][0], (float)stat->bit_use_delta_quant[B_SLICE]/Bframe_ctr);
	else
		fprintf(p_stat," Delta quant          | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_delta_quant[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_delta_quant[P_SLICE]/bit_use[P_SLICE][0], 0.);
	
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," Stuffing Bits        | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_stuffingBits[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_stuffingBits[P_SLICE]/bit_use[P_SLICE][0], (float)stat->bit_use_stuffingBits[B_SLICE]/Bframe_ctr);
	else
		fprintf(p_stat," Stuffing Bits        | %10.2f     | %10.2f     | %10.2f     |\n",
		(float)stat->bit_use_stuffingBits[I_SLICE]/bit_use[I_SLICE][0], (float)stat->bit_use_stuffingBits[P_SLICE]/bit_use[P_SLICE][0], 0.);
	
	
	
	fprintf(p_stat," ---------------------|----------------|----------------|----------------|\n");
	
	fprintf(p_stat," average bits/frame   |");
	
	fprintf(p_stat," %10.2f     |", (float) bit_use[I_SLICE][1]/(float) bit_use[I_SLICE][0] );
	fprintf(p_stat," %10.2f     |", (float) bit_use[P_SLICE][1]/(float) bit_use[P_SLICE][0] );
	
	if(input->successive_Bframe!=0 && Bframe_ctr!=0)
		fprintf(p_stat," %10.2f     |", (float) bit_use[B_SLICE][1]/ (float) Bframe_ctr );
	else fprintf(p_stat," %10.2f     |", 0.);
	
	fprintf(p_stat,"\n");
	fprintf(p_stat," ---------------------|----------------|----------------|----------------|\n");
	
	fclose(p_stat);
	
	// write to log file
	if ((p_log=fopen("log.dat","r"))==0)                      // check if file exist
	{
		if ((p_log=fopen("log.dat","a"))==NULL)            // append new statistic at the end
		{
			printf(errortext, ET_SIZE, "Error open file %s  \n","log.dat");
			error(errortext, 500);
		}
		else                                            // Create header for new log file
		{
			fprintf(p_log," ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ \n");
			fprintf(p_log,"|            Encoder statistics. This file is generated during first encoding session, new sessions will be appended                                                                                          |\n");
			fprintf(p_log," ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ \n");
			fprintf(p_log,"| Date  | Time  |    Sequence        | #Img | QPI| QPP| QPB| Format  | #B | Hdmd | S.R |#Ref | Freq |Intra upd|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|#Bitr P|#Bitr B|     Total Time   |      Me Time     |\n");
			fprintf(p_log," ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ \n");
			/*
			fprintf(p_log," ----------------------------------------------------------------------------------------------------------------------------------------------------------------- \n");
			fprintf(p_log,"|            Encoder statistics. This file is generated during first encoding session, new sessions will be appended                                              |\n");
			fprintf(p_log," ----------------------------------------------------------------------------------------------------------------------------------------------------------------- \n");
			fprintf(p_log,"| Date  | Time  |    Sequence        |#Img|Quant1|QuantN|Format|Hadamard|Search r|#Ref | Freq |Intra upd|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|#Bitr P|#Bitr B|\n");
			fprintf(p_log," ----------------------------------------------------------------------------------------------------------------------------------------------------------------- \n");
			*/
		}
	}
	else
	{
		fclose (p_log);
		if ((p_log=fopen("log.dat","a"))==NULL)            // File exist,just open for appending
		{
			printf(errortext, ET_SIZE, "Error open file %s  \n","log.dat");
			error(errortext, 500);
		}
	}
	
#ifdef WIN32
	_strdate( timebuf );
	fprintf(p_log,"| %1.5s |",timebuf );
	
	_strtime( timebuf);
	fprintf(p_log," % 1.5s |",timebuf);
#else
	now = time ((time_t *) NULL); // Get the system time and put it into 'now' as 'calender time'
	time (&now);
	l_time = localtime (&now);
	strftime (string, sizeof string, "%d-%b-%Y", l_time);
	fprintf(p_log,"| %1.5s |",string );
	
	strftime (string, sizeof string, "%H:%M:%S", l_time);
	fprintf(p_log," %1.5s |",string );
#endif
	
	for (i=0;i<20;i++)
		name[i]=input->infile[i+max(0,((int)strlen(input->infile))-20)]; // write last part of path, max 20 chars
	fprintf(p_log,"%20.20s|",name);
	
	fprintf(p_log,"%5d |",input->no_frames_h264);
	fprintf(p_log,"  %2d  |",input->qp0);
	fprintf(p_log,"  %2d  |",input->qpN);
	fprintf(p_log," %2d |",input->qpB);
	
	fprintf(p_log,"%4dx%-4d|",input->imagewidth,input->imageheight);
	
	fprintf(p_log,"%3d |",input->successive_Bframe); 
	
	
	if (input->hadamard==1)
		fprintf(p_log,"   ON   |");
	else
		fprintf(p_log,"   OFF  |");
	
	fprintf(p_log,"   %3d   |",input->search_range );
	
	fprintf(p_log," %2d  |",input->num_reference_frames);
	
	//  fprintf(p_log," %3d  |",img->framerate/(input->jumpd+1));
    fprintf(p_log," %3d  |",(img->framerate *(input->successive_Bframe + 1)) / (input->jumpd+1));
	
	
	if (input->intra_upd==1)
		fprintf(p_log,"   ON    |");
	else
		fprintf(p_log,"   OFF   |");
	
	fprintf(p_log,"%-5.3f|",snr->snr_y1);
	fprintf(p_log,"%-5.3f|",snr->snr_u1);
	fprintf(p_log,"%-5.3f|",snr->snr_v1);
	fprintf(p_log,"%-5.3f|",snr->snr_ya);
	fprintf(p_log,"%-5.3f|",snr->snr_ua);
	fprintf(p_log,"%-5.3f|",snr->snr_va);
	if(input->successive_Bframe != 0)
	{
		fprintf(p_log,"%7.0f|",stat->bitrate_P);
		fprintf(p_log,"%7.0f|",stat->bitrate_B);
	}
	else
	{
		fprintf(p_log,"%7.0f|",stat->bitrate);
		fprintf(p_log,"%7.0f|",0.0);
	}
	
	
	fprintf(p_log,"   %12d   |", tot_time);
	fprintf(p_log,"   %12d   |\n", me_tot_time);
	
	/*
	if(input->successive_Bframe != 0)
	{
    fprintf(p_log,"%7.0f|",stat->bitrate_P);
    fprintf(p_log,"%7.0f|\n",stat->bitrate_B);
	}
	else
	{
    fprintf(p_log,"%7.0f|",stat->bitrate);
    fprintf(p_log,"%7.0f|\n",0.0);
	}
	*/
	
	fclose(p_log);
	
	p_log=fopen("data.txt","a");
	
	if(input->successive_Bframe != 0 && Bframe_ctr != 0) // B picture used
	{
		fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d %5d %.3f\n",
			input->no_frames_h264, input->qp0, input->qpN,
			snr->snr_y1,
			snr->snr_u1,
			snr->snr_v1,
			stat->bit_ctr_0,
			0.0,
			0.0,
			0.0,
			0,
			snr->snr_ya,
			snr->snr_ua,
			snr->snr_va,
			(stat->bit_ctr_0+stat->bit_ctr)/(input->no_frames_h264+Bframe_ctr),
			stat->bit_ctr_B/Bframe_ctr,
			(double)0.001*tot_time/(input->no_frames_h264+Bframe_ctr));
	}
	else
	{
		fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d %5d %.3f\n",
			input->no_frames_h264, input->qp0, input->qpN,
			snr->snr_y1,
			snr->snr_u1,
			snr->snr_v1,
			stat->bit_ctr_0,
			0.0,
			0.0,
			0.0,
			0,
			snr->snr_ya,
			snr->snr_ua,
			snr->snr_va,
			(stat->bit_ctr_0+stat->bit_ctr)/input->no_frames_h264,
			0,
			(double)0.001*tot_time/input->no_frames_h264);
	}
	
	fclose(p_log);
	
	
}


/*!
************************************************************************
* \brief
*    Prints the header of the protocol.
* \par Input:
*    struct inp_par_dec *inp
* \par Output:
*    none
************************************************************************
*/
void information_init()
{
	
	printf("-------------------------------------------------------------------------------\n");
	printf(" Input YUV file                    : %s \n",input->infile_c);
	printf(" Output H.264 bitstream            : %s \n",input->outfile);
	if (p_dec != NULL)
		printf(" Output YUV file                   : %s \n",input->ReconFile);
	printf(" Output log file                   : log.dat \n");
	printf(" Output statistics file            : stat.dat \n");
	printf("-------------------------------------------------------------------------------\n");
	printf("  Frame  Bit/pic WP QP   SnrY    SnrU    SnrV    Time(ms) MET(ms) Frm/Fld   I D\n");
	//  printf(" Frame   Bit/pic   QP   SnrY    SnrU    SnrV    Time(ms) Frm/Fld IntraMBs\n");
	printf("-------------------------------------------------------------------------------\n");
}


/*!
************************************************************************
* \brief
*    Free allocated memory of frame size related global buffers
*    buffers are defined in global.h, allocated memory is allocated in
*    int get_mem4global_buffers()
* \par Input:
*    Input Parameters struct inp_par_dec *inp,                             \n
*    Image Parameters struct img_par *img
* \par Output:
*    none
************************************************************************
*/
void free_global_buffers()
{
	int  i,j;
	
#ifdef _ADAPT_LAST_GROUP_
	extern int *last_P_no_frm;
	extern int *last_P_no_fld;
	free (last_P_no_frm);
	free (last_P_no_fld);
#endif
	
	free_mem2D(imgY_org_frm);      // free ref frame buffers
	free_mem3D(imgUV_org_frm,2);
	
	// free multiple ref frame buffers
	// number of reference frames increased by one for next P-frame
	
	if (input->WeightedPrediction || input->WeightedBiprediction)
	{
		free_mem3Dint(wp_weight,6);
		free_mem3Dint(wp_offset,6);
		free_mem4Dint(wbp_weight,6,MAX_REFERENCE_PICTURES);
	}
	
	if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
	{
		free_mem3Dint(direct_ref_idx,2);//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ 
		free_mem2Dint(direct_pdir);
	} // end if B frame
	
	
	free_mem2Dint(img4Y_tmp);    // free temp quarter pel frame buffer
	
	// free mem, allocated in init_img()
	// free intra pred mode buffer for blocks
	free_mem2Dint(img->ipredmode);
	free(img->mb_data);
	
	if(input->UseConstrainedIntraPred)
	{
		free (img->intra_block);
	}
	
	if (input->rdopt==2)
	{
		free(decs->resY[0]);
		free(decs->resY);
		free(decs->RefBlock[0]);
		free(decs->RefBlock);
		for (j=0; j<input->NoOfDecoders; j++)
		{
			free(decs->decY[j][0]);
			free(decs->decY[j]);
			free(decs->decY_best[j][0]);
			free(decs->decY_best[j]);
			for (i=0; i<img->max_num_references+1; i++)
			{
				free(decs->decref[j][i][0]);
				free(decs->decref[j][i]);
			}
			free(decs->decref[j]);
		}
		free(decs->decY);
		free(decs->decY_best);
		free(decs->decref);
		free(decs->status_map[0]);
		free(decs->status_map);
		free(decs->dec_mb_mode[0]);
		free(decs->dec_mb_mode);
	}
	if (input->RestrictRef)
	{
		free(pixel_map[0]);
		free(pixel_map);
		free(refresh_map[0]);
		free(refresh_map);
	}
	
	if(!active_sps->frame_mbs_only_flag)
	{
		free_mem2D(imgY_com);
		free_mem3D(imgUV_com,2);
		free_mem2D(imgY_org_top);      // free ref frame buffers
		free_mem3D(imgUV_org_top,2);
		free_mem2D(imgY_org_bot);      // free ref frame buffers
		free_mem3D(imgUV_org_bot,2);
		
	}
	
	free_mem3Dint(img->nz_coeff, img->FrameSizeInMbs);
	
	if(input->FMEnable)
		free_mem_FME();
}

/*!
************************************************************************
* \brief
*    Allocate memory for mv
* \par Input:
*    Image Parameters struct img_par *img                             \n
*    int****** mv
* \return memory size in bytes
************************************************************************
*/
int get_mem_mv (int******* mv)
{
	int i, j, k, l, m;
	
	if ((*mv = (int******)calloc(4,sizeof(int*****))) == NULL)
		no_mem_exit ("get_mem_mv: mv");
	for (i=0; i<4; i++)
	{
		if (((*mv)[i] = (int*****)calloc(4,sizeof(int****))) == NULL)
			no_mem_exit ("get_mem_mv: mv");
		for (j=0; j<4; j++)
		{
			if (((*mv)[i][j] = (int****)calloc(2,sizeof(int***))) == NULL)
				no_mem_exit ("get_mem_mv: mv");
			for (k=0; k<2; k++)
			{
				if (((*mv)[i][j][k] = (int***)calloc(img->max_num_references,sizeof(int**))) == NULL)
					no_mem_exit ("get_mem_mv: mv");
				for (l=0; l<img->max_num_references; l++)
				{
					if (((*mv)[i][j][k][l] = (int**)calloc(9,sizeof(int*))) == NULL)
						no_mem_exit ("get_mem_mv: mv");
					for (m=0; m<9; m++)
						if (((*mv)[i][j][k][l][m] = (int*)calloc(2,sizeof(int))) == NULL)
							no_mem_exit ("get_mem_mv: mv");
				}
			}
		}
	}
	return 4*4*img->max_num_references*9*2*sizeof(int);
}


/*!
************************************************************************
* \brief
*    Free memory from mv
* \par Input:
*    int****** mv
************************************************************************
*/
void free_mem_mv (int****** mv)
{
	int i, j, k, l, m;
	
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			for (k=0; k<2; k++)
			{
				for (l=0; l<img->max_num_references; l++)
				{
					for (m=0; m<9; m++)
					{
						free(mv[i][j][k][l][m]);
					}
					free (mv[i][j][k][l]);
				}
				free (mv[i][j][k]);
			}
			free (mv[i][j]);
		}
		free (mv[i]);
	}
	free (mv);
}





/*!
************************************************************************
* \brief
*    Allocate memory for AC coefficients
************************************************************************
*/
int get_mem_ACcoeff (int***** cofAC)
{
	int i, j, k;
	//开辟内存时候 cofAC 第一维是开的几。如果是 6，那前四个就是亮度，后两个就是色度；
	if ((*cofAC = (int****)calloc (6, sizeof(int***))) == NULL)              no_mem_exit ("get_mem_ACcoeff: cofAC");
	for (k=0; k<6; k++)
	{
		if (((*cofAC)[k] = (int***)calloc (4, sizeof(int**))) == NULL)         no_mem_exit ("get_mem_ACcoeff: cofAC");
		for (j=0; j<4; j++)
		{
			if (((*cofAC)[k][j] = (int**)calloc (2, sizeof(int*))) == NULL)      no_mem_exit ("get_mem_ACcoeff: cofAC");
			for (i=0; i<2; i++)
			{
				if (((*cofAC)[k][j][i] = (int*)calloc (65, sizeof(int))) == NULL)  no_mem_exit ("get_mem_ACcoeff: cofAC"); // 18->65 for ABT
			}
		}
	}
	return 6*4*2*65*sizeof(int);// 18->65 for ABT
}

/*!
************************************************************************
* \brief
*    Allocate memory for DC coefficients
************************************************************************
*/
int get_mem_DCcoeff (int**** cofDC)
{
	int j, k;
	
	if ((*cofDC = (int***)calloc (3, sizeof(int**))) == NULL)           no_mem_exit ("get_mem_DCcoeff: cofDC");
	for (k=0; k<3; k++)
	{
		if (((*cofDC)[k] = (int**)calloc (2, sizeof(int*))) == NULL)      no_mem_exit ("get_mem_DCcoeff: cofDC");
		for (j=0; j<2; j++)
		{
			if (((*cofDC)[k][j] = (int*)calloc (65, sizeof(int))) == NULL)  no_mem_exit ("get_mem_DCcoeff: cofDC"); // 18->65 for ABT
		}
	}
	return 3*2*65*sizeof(int); // 18->65 for ABT
}


/*!
************************************************************************
* \brief
*    Free memory of AC coefficients
************************************************************************
*/
void free_mem_ACcoeff (int**** cofAC)
{
	int i, j, k;
	
	for (k=0; k<6; k++)
	{
		for (i=0; i<4; i++)
		{
			for (j=0; j<2; j++)
			{
				free (cofAC[k][i][j]);
			}
			free (cofAC[k][i]);
		}
		free (cofAC[k]);
	}
	free (cofAC);
}

/*!
************************************************************************
* \brief
*    Free memory of DC coefficients
************************************************************************
*/
void free_mem_DCcoeff (int*** cofDC)
{
	int i, j;
	
	for (j=0; j<3; j++)
	{
		for (i=0; i<2; i++)
		{
			free (cofDC[j][i]);
		}
		free (cofDC[j]);
	}
	free (cofDC);
}


/*!
************************************************************************
* \brief
*    form frame picture from two field pictures 
************************************************************************
*/
void combine_field()
{
	int i;
	
	for (i=0; i<img->height / 2; i++)
	{
		memcpy(imgY_com[i*2], enc_top_picture->imgY[i], img->width);     // top field
		memcpy(imgY_com[i*2 + 1], enc_bottom_picture->imgY[i], img->width); // bottom field
	}
	
	for (i=0; i<img->height_cr / 2; i++)
	{
		memcpy(imgUV_com[0][i*2], enc_top_picture->imgUV[0][i], img->width_cr);
		memcpy(imgUV_com[0][i*2 + 1], enc_bottom_picture->imgUV[0][i], img->width_cr);
		memcpy(imgUV_com[1][i*2], enc_top_picture->imgUV[1][i], img->width_cr);
		memcpy(imgUV_com[1][i*2 + 1], enc_bottom_picture->imgUV[1][i], img->width_cr);
	}
}

/*!
************************************************************************
* \brief
*    RD decision of frame and field coding 
************************************************************************
*/
int decide_fld_frame(float snr_frame_Y, float snr_field_Y, int bit_field, int bit_frame, double lambda_picture)
{
	double cost_frame, cost_field;
	
	cost_frame = bit_frame * lambda_picture + snr_frame_Y;
	cost_field = bit_field * lambda_picture + snr_field_Y;
	
	if (cost_field > cost_frame)
		return (0);
	else
		return (1);
}

/*!
************************************************************************
* \brief
*    Do some initializaiton work for encoding the 2nd IGOP
************************************************************************
*/
void process_2nd_IGOP()
{
	Boolean FirstIGOPFinished = FALSE;
	if ( img->current_frame == input->no_frames_h264-1 )
		FirstIGOPFinished = TRUE;
	if (input->NumFrameIn2ndIGOP==0) return;
	if (!FirstIGOPFinished || In2ndIGOP) return;
	In2ndIGOP = TRUE;
	
	//  img->current_frame = -1;
	start_frame_no_in_this_IGOP = input->no_frames_h264;
	start_tr_in_this_IGOP = (input->no_frames_h264-1)*(input->jumpd+1) +1;
	input->no_frames_h264 = input->no_frames_h264 + input->NumFrameIn2ndIGOP;
	
	/*  reset_buffers();
	
	  frm->picbuf_short[0]->used=0;
	  frm->picbuf_short[0]->picID=-1;
	  frm->picbuf_short[0]->lt_picID=-1;
	frm->short_used = 0; */
	img->nb_references = 0;
}


/*!
************************************************************************
* \brief
*    Set the image type for I,P and SP pictures (not B!)
************************************************************************
*/
void SetImgType()
{
	if ((img->current_frame%input->I_frame)==0)
	{
// 		if (IMG_NUMBER == 0)//( img->current_frame - start_frame_no_in_this_IGOP )
// 		{
			img->type = I_SLICE;        // set image type for first image to I-frame
// 		}
// 		else
// 		{
// 			img->type = P_SLICE;        // P-frame
// 			
// 			if (input->sp_periodicity)
// 			{
// 				if ((IMG_NUMBER % input->sp_periodicity) ==0)
// 				{
// 					img->type=SP_SLICE;
// 				}
// 			}
// 		}
	}
	else
	{
// 		if ((IMG_NUMBER%input->intra_period) == 0)
// 		{
			img->type = P_SLICE;
// 		}
// 		else
// 		{
// 			img->type = P_SLICE;        // P-frame
// 			if (input->sp_periodicity)
// 			{
// 				if ((IMG_NUMBER % input->sp_periodicity) ==0)
// 					img->type=SP_SLICE;
// 			}
// 		}
	}
}
/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 * \par Input:
 *    Input Parameters struct inp_par *inp,                            \n
 *    Image Parameters struct img_par *img
 * \return Number of allocated bytes
 ************************************************************************
 */
int init_global_buffers()
{
  int j,memory_size=0;
  int height_field = img->height/2;
#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no_frm;
  extern int *last_P_no_fld;

  if ((last_P_no_frm = (int*)malloc(2*img->max_num_references*sizeof(int))) == NULL)
    no_mem_exit("init_global_buffers: last_P_no");
  if(!active_sps->frame_mbs_only_flag)
    if ((last_P_no_fld = (int*)malloc(4*img->max_num_references*sizeof(int))) == NULL)
      no_mem_exit("init_global_buffers: last_P_no");
#endif

  // allocate memory for encoding frame buffers: imgY, imgUV

  // allocate memory for reference frame buffers: imgY_org, imgUV_org
  // byte imgY_org[288][352];
  // byte imgUV_org[2][144][176];
  memory_size += get_mem2D(&imgY_org_frm, img->height, img->width);
  memory_size += get_mem3D(&imgUV_org_frm, 2, img->height_cr, img->width_cr);


    if (input->WeightedPrediction || input->WeightedBiprediction)
    {
      // Currently only use up to 20 references. Need to use different indicator such as maximum num of references in list
      memory_size += get_mem3Dint(&wp_weight,6,MAX_REFERENCE_PICTURES,3);
      memory_size += get_mem3Dint(&wp_offset,6,MAX_REFERENCE_PICTURES,3);

      memory_size += get_mem4Dint(&wbp_weight, 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);
    }

  // allocate memory for reference frames of each block: refFrArr
  // int  refFrArr[72][88];

  if(input->successive_Bframe!=0 || input->StoredBPictures > 0)
  {    
    memory_size += get_mem3Dint(&direct_ref_idx, 2, img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);
    memory_size += get_mem2Dint(&direct_pdir, img->width/BLOCK_SIZE, img->height/BLOCK_SIZE);
  }

  // allocate memory for temp quarter pel luma frame buffer: img4Y_tmp
  // int img4Y_tmp[576][704];  (previously int imgY_tmp in global.h)
  memory_size += get_mem2Dint(&img4Y_tmp, img->height+2*IMG_PAD_SIZE, (img->width+2*IMG_PAD_SIZE)*4);

  if (input->rdopt==2)
  {
    memory_size += get_mem2Dint(&decs->resY, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    if ((decs->decref = (byte****) calloc(input->NoOfDecoders,sizeof(byte***))) == NULL) 
      no_mem_exit("init_global_buffers: decref");
    for (j=0 ; j<input->NoOfDecoders; j++)
    {
      memory_size += get_mem3D(&decs->decref[j], img->max_num_references+1, img->height, img->width);
    }
    memory_size += get_mem2D(&decs->RefBlock, BLOCK_SIZE,BLOCK_SIZE);
    memory_size += get_mem3D(&decs->decY, input->NoOfDecoders, img->height, img->width);
    memory_size += get_mem3D(&decs->decY_best, input->NoOfDecoders, img->height, img->width);
    memory_size += get_mem2D(&decs->status_map, img->height/MB_BLOCK_SIZE,img->width/MB_BLOCK_SIZE);
    memory_size += get_mem2D(&decs->dec_mb_mode, img->width/MB_BLOCK_SIZE,img->height/MB_BLOCK_SIZE);
  }
  if (input->RestrictRef)
  {
    memory_size += get_mem2D(&pixel_map, img->height,img->width);
    memory_size += get_mem2D(&refresh_map, img->height/8,img->width/8);
  }

  if(!active_sps->frame_mbs_only_flag)
  {
    // allocate memory for encoding frame buffers: imgY, imgUV
    memory_size += get_mem2D(&imgY_com, img->height, img->width);
    memory_size += get_mem3D(&imgUV_com, 2, img->height/2, img->width_cr);

    // allocate memory for reference frame buffers: imgY_org, imgUV_org
    memory_size += get_mem2D(&imgY_org_top, height_field, img->width);
    memory_size += get_mem3D(&imgUV_org_top, 2, height_field/2, img->width_cr);
    memory_size += get_mem2D(&imgY_org_bot, height_field, img->width);
    memory_size += get_mem3D(&imgUV_org_bot, 2, height_field/2, img->width_cr);

  }

  if(input->FMEnable)
    memory_size += get_mem_FME();

  return (memory_size);
}

////////////////////////////////////////////>

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 *
 *  \par Input:
 *    Input Parameters struct inp_par_dec *inp, Image Parameters struct img_par *img
 *
 *  \par Output:
 *     Number of allocated bytes
 ***********************************************************************
 */
int init_global_buffers_dec()
{
  int memory_size=0;

  if (global_init_done)
  {
    free_global_buffers();
  }

  // allocate memory for reference frame in find_snr
  memory_size += get_mem2D(&imgY_ref, img->height, img->width);
  memory_size += get_mem3D(&imgUV_ref, 2, img->height_cr, img->width_cr);
  /////////////////////////////////////>解决CopyFrameToOldImgOrgVariables中imgY_org_frm[y][x]分配不足内存问题
  
  
  memory_size += get_mem2D(&imgY_org_frm, img->height, img->width);
  memory_size += get_mem3D(&imgUV_org_frm, 2, img->height_cr, img->width_cr);
//解决内存分配不足
  memory_size += get_mem2Dint(&img4Y_tmp, img->height+2*IMG_PAD_SIZE, (img->width+2*IMG_PAD_SIZE)*4);

  
/////////////////////////////////////<

  // allocate memory in structure img
  if(((img->mb_data) = (Macroblock *) calloc(img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
    no_mem_exit("init_global_buffers: img->mb_data");

  if(((img->intra_block) = (int*)calloc(img->FrameSizeInMbs, sizeof(int))) == NULL)
    no_mem_exit("init_global_buffers: img->intra_block");

  memory_size += get_mem2Dint(&(img->ipredmode), 4*img->PicWidthInMbs , 4*img->FrameHeightInMbs);

  memory_size += get_mem2Dint(&(img->field_anchor),4*img->FrameHeightInMbs, 4*img->PicWidthInMbs);

  memory_size += get_mem3Dint(&(img->wp_weight), 2, MAX_REFERENCE_PICTURES, 3);
  memory_size += get_mem3Dint(&(img->wp_offset), 6, MAX_REFERENCE_PICTURES, 3);
  memory_size += get_mem4Dint(&(img->wbp_weight), 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);

  // CAVLC mem
  memory_size += get_mem3Dint(&(img->nz_coeff), img->FrameSizeInMbs, 4, 6);

  memory_size += get_mem2Dint(&(img->siblock),img->PicWidthInMbs  , img->FrameHeightInMbs);

  global_init_done = 1;

  img->oldFrameSizeInMbs = img->FrameSizeInMbs;

  return (memory_size);
}



/*!
***********************************************************************
* \brief
*    Initilize some arrays
***********************************************************************
*/
void init_dec(struct img_par *img)  //!< image parameters
{
	int i;
	// initilize quad matrix used in snr routine
	for (i=0; i <  256; i++)
	{
		img->quad_dec[i]=i*i; // fix from TML 1, truncation removed
	}
	
	img->oldFrameSizeInMbs = -1;
}


/*!
************************************************************************
* \brief
*    Read input from configuration file
*
* \par Input:
*    Name of configuration filename
*
* \par Output
*    none
************************************************************************
*/
void init_conf_dec(struct inp_par_dec *inp,
               char *config_filename)
{
	FILE *fd;
	int NAL_mode;
	
	// read the decoder configuration file
	if((fd=fopen("decoder.cfg","r")) == NULL)
	{
		printf(errortext, ET_SIZE, "Error: Control file %s not found\n",config_filename);
		error(errortext, 300);
	}


	fscanf(fd,"%s",inp->infile);                // H.264 compressed input bitsream
	fscanf(fd,"%*[^\n]");
	
	fscanf(fd,"%s",inp->outfile);               // YUV 4:2:2 input format
	fscanf(fd,"%*[^\n]");
	
	fscanf(fd,"%s",inp->reffile);               // reference file
	fscanf(fd,"%*[^\n]");
	
	
	// Frame buffer size
	fscanf(fd,"%d,",&inp->dpb_size);   // may be overwritten in case of RTP NAL
	fscanf(fd,"%*[^\n]");
	if (inp->dpb_size < 1)
	{
		printf(errortext, ET_SIZE, "Decoded Picture Buffer Size is %d. It has to be at least 1",inp->dpb_size);
		error(errortext,1);
	}
	
	fscanf(fd,"%d",&(NAL_mode));                // NAL mode
    fscanf(fd,"%*[^\n]");
	
	switch(NAL_mode)
	{
	case 0:
		inp->FileFormat = PAR_OF_ANNEXB;
		break;
	case 1:
		inp->FileFormat = PAR_OF_RTP;
		break;
	default:
		printf(errortext, ET_SIZE, "NAL mode %i is not supported", NAL_mode);
		error(errortext,400);
	}
	
	fscanf(fd,"%d,",&inp->ref_offset);   // offset used for SNR computation
	fscanf(fd,"%*[^\n]");
	
	fscanf(fd,"%d,",&inp->poc_scale);   // offset used for SNR computation
	fscanf(fd,"%*[^\n]");
	
	
	if (inp->poc_scale < 1 || inp->poc_scale > 2)
	{
		printf(errortext, ET_SIZE, "Poc Scale is %d. It has to be 1 or 2",inp->poc_scale);
		error(errortext,1);
	}
	
// #ifdef _LEAKYBUCKET_ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ解码的
// 	fscanf(fd,"%ld,",&inp->R_decoder);             // Decoder rate
// 	fscanf(fd, "%*[^\n]");
// 	fscanf(fd,"%ld,",&inp->B_decoder);             // Decoder buffer size
// 	fscanf(fd, "%*[^\n]");
// 	fscanf(fd,"%ld,",&inp->F_decoder);             // Decoder initial delay
// 	fscanf(fd, "%*[^\n]"); 
// 	fscanf(fd,"%s",inp->LeakyBucketParamFile);    // file where Leaky Bucket params (computed by encoder) are stored
// 	fscanf(fd,"%*[^\n]");
// #endif
	
	fclose (fd);
	
	
#if TRACE
	if ((p_trace=fopen(TRACEFILE,"w"))==0)             // append new statistic at the end
	{
		printf(errortext, ET_SIZE, "Error open file %s!",TRACEFILE);
		error(errortext,500);
	}
#endif
	
	
	if ((p_out=fopen(inp->outfile,"wb"))==0)
	{
		printf(errortext, ET_SIZE, "Error open file %s ",inp->outfile);
		error(errortext,500);
	}
	/*  if ((p_out2=fopen("out.yuv","wb"))==0)
	{
    printf(errortext, ET_SIZE, "Error open file %s ",inp->outfile);
    error(errortext,500);
}*/
	
	fprintf(stdout,"--------------------------------------------------------------------------\n");
	fprintf(stdout," Decoder config file                    : %s \n",config_filename);
	fprintf(stdout,"--------------------------------------------------------------------------\n");
	fprintf(stdout," Input H.264 bitstream                  : %s \n",inp->infile);
	fprintf(stdout," Output decoded YUV 4:2:0               : %s \n",inp->outfile);
	fprintf(stdout," Output status file                     : %s \n",LOGFILE);
	if ((p_ref=fopen(inp->reffile,"rb"))==0)
	{
		fprintf(stdout," Input reference file                   : %s does not exist \n",inp->reffile);
		fprintf(stdout,"                                          SNR values are not available\n");
	}
	else
		fprintf(stdout," Input reference file                   : %s \n",inp->reffile);
	
	fprintf(stdout,"--------------------------------------------------------------------------\n");
// #ifdef _LEAKYBUCKET_
// 	fprintf(stdout," Rate_decoder        : %8ld \n",inp->R_decoder);
// 	fprintf(stdout," B_decoder           : %8ld \n",inp->B_decoder);
// 	fprintf(stdout," F_decoder           : %8ld \n",inp->F_decoder);
// 	fprintf(stdout," LeakyBucketParamFile: %s \n",inp->LeakyBucketParamFile); // Leaky Bucket Param file
// 	calc_buffer(inp);
// 	fprintf(stdout,"--------------------------------------------------------------------------\n");
// #endif
	fprintf(stdout,"POC must = frame# or field# for SNRs to be correct\n");
	fprintf(stdout,"Frame    POC   QP  SnrY    SnrU    SnrV   Time(ms)\n");
}

/*!
************************************************************************
* \brief
*    Reports the gathered information to appropriate outputs
*
* \par Input:
*    struct inp_par_dec *inp,
*    struct img_par *img,
*    struct snr_par *stat
*
* \par Output:
*    None
************************************************************************
*/
void report_dec(struct inp_par_dec *inp, struct img_par *img, struct snr_par *snr)
{
#define OUTSTRING_SIZE 255
	char string[OUTSTRING_SIZE];
	FILE *p_log;
	
#ifndef WIN32
	time_t  now;
	struct tm *l_time;
#else
	char timebuf[128];
#endif
	
	fprintf(stdout,"-------------------- Average SNR all frames ------------------------------\n");
	fprintf(stdout," SNR Y(dB)           : %5.2f\n",snr->snr_ya);
	fprintf(stdout," SNR U(dB)           : %5.2f\n",snr->snr_ua);
	fprintf(stdout," SNR V(dB)           : %5.2f\n",snr->snr_va);
	fprintf(stdout," Total decoding time : %.3f sec \n",tot_time*0.001);
	fprintf(stdout,"--------------------------------------------------------------------------\n");
	fprintf(stdout," Exit JM %s decoder, ver %s ",JM,VERSION);
	fprintf(stdout,"\n");
	// write to log file
	
	printf(string, OUTSTRING_SIZE, "%s", LOGFILE);
	if ((p_log=fopen(string,"r"))==0)                    // check if file exist
	{
		if ((p_log=fopen(string,"a"))==0)
		{
			printf(errortext, ET_SIZE, "Error open file %s for appending",string);
			error(errortext, 500);
		}
		else                                              // Create header to new file
		{
			fprintf(p_log," ------------------------------------------------------------------------------------------\n");
			fprintf(p_log,"|  Decoder statistics. This file is made first time, later runs are appended               |\n");
			fprintf(p_log," ------------------------------------------------------------------------------------------ \n");
			fprintf(p_log,"| Date  | Time  |    Sequence        |#Img|Format|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|\n");
			fprintf(p_log," ------------------------------------------------------------------------------------------\n");
		}
	}
	else
	{ 
		fclose(p_log);
		p_log=fopen(string,"a");                    // File exist,just open for appending
	}
	
#ifdef WIN32
	_strdate( timebuf );
	fprintf(p_log,"| %1.5s |",timebuf );
	
	_strtime( timebuf);
	fprintf(p_log," % 1.5s |",timebuf);
#else
	now = time ((time_t *) NULL); // Get the system time and put it into 'now' as 'calender time'
	time (&now);
	l_time = localtime (&now);
	strftime (string, sizeof string, "%d-%b-%Y", l_time);
	fprintf(p_log,"| %1.5s |",string );
	
	strftime (string, sizeof string, "%H:%M:%S", l_time);
	fprintf(p_log,"| %1.5s |",string );
#endif
	
	fprintf(p_log,"%20.20s|",inp->infile);
	
	fprintf(p_log,"%3d |",img->number);
	
	fprintf(p_log,"%6.3f|",snr->snr_y1);
	fprintf(p_log,"%6.3f|",snr->snr_u1);
	fprintf(p_log,"%6.3f|",snr->snr_v1);
	fprintf(p_log,"%6.3f|",snr->snr_ya);
	fprintf(p_log,"%6.3f|",snr->snr_ua);
	fprintf(p_log,"%6.3f|\n",snr->snr_va);
	
	fclose(p_log);
	
	printf(string, OUTSTRING_SIZE,"%s", DATADECFILE);
	p_log=fopen(string,"a");
	
	if(Bframe_ctr != 0) // B picture used
	{
		fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d %.3f\n",
			img->number, 0, img->qp,
			snr->snr_y1,
			snr->snr_u1,
			snr->snr_v1,
			0,
			0.0,
			0.0,
			0.0,
			0,
			snr->snr_ya,
			snr->snr_ua,
			snr->snr_va,
			0,
			(double)0.001*tot_time/(img->number+Bframe_ctr-1));
	}
	else
	{
		fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d "
			"%2.2f %2.2f %2.2f %5d %.3f\n",
			img->number, 0, img->qp,
			snr->snr_y1,
			snr->snr_u1,
			snr->snr_v1,
			0,
			0.0,
			0.0,
			0.0,
			0,
			snr->snr_ya,
			snr->snr_ua,
			snr->snr_va,
			0,
			(double)0.001*tot_time/img->number);
	}
	fclose(p_log);
}

/*!
************************************************************************
* \brief
*    Allocates a stand-alone partition structure.  Structure should
*    be freed by FreePartition();
*    data structures
*
* \par Input:
*    n: number of partitions in the array
* \par return
*    pointer to DataPartition Structure, zero-initialized
************************************************************************
*/

DataPartition *AllocPartition(int n)
{
	DataPartition *partArr, *dataPart;
	int i;
	
	partArr = (DataPartition *) calloc(n, sizeof(DataPartition));
	if (partArr == NULL)
	{
		printf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Data Partition failed");
		error(errortext, 100);
	}
	
	for (i=0; i<n; i++) // loop over all data partitions
	{
		dataPart = &(partArr[i]);
		dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
		if (dataPart->bitstream == NULL)
		{
			printf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Bitstream failed");
			error(errortext, 100);
		}
		dataPart->bitstream->streamBuffer = (byte *) calloc(MAX_CODED_FRAME_SIZE, sizeof(byte));
		if (dataPart->bitstream->streamBuffer == NULL)
		{
			printf(errortext, ET_SIZE, "AllocPartition: Memory allocation for streamBuffer failed");
			error(errortext, 100);
		}
	}
	return partArr;
}




/*!
************************************************************************
* \brief
*    Frees a partition structure (array).  
*
* \par Input:
*    Partition to be freed, size of partition Array (Number of Partitions)
*
* \par return
*    None
*
* \note
*    n must be the same as for the corresponding call of AllocPartition
************************************************************************
*/


void FreePartition (DataPartition *dp, int n)
{
	int i;
	
	assert (dp != NULL);
	assert (dp->bitstream != NULL);
	assert (dp->bitstream->streamBuffer != NULL);
	for (i=0; i<n; i++)
	{
		free (dp[i].bitstream->streamBuffer);
		free (dp[i].bitstream);
	}
	free (dp);
}


/*!
************************************************************************
* \brief
*    Allocates the slice structure along with its dependent
*    data structures
*
* \par Input:
*    Input Parameters struct inp_par_dec *inp,  struct img_par *img
************************************************************************
*/
void malloc_slice_dec(struct inp_par_dec *inp, struct img_par *img)
{
	Slice *currSlice;
	
	img->currentSlice = (Slice *) calloc(1, sizeof(Slice));
	if ( (currSlice = img->currentSlice) == NULL)
	{
		printf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", inp->FileFormat);
		error(errortext,100);
	}
	//  img->currentSlice->rmpni_buffer=NULL;
	//! you don't know whether we do CABAC hre, hence initialize CABAC anyway
	// if (inp->symbol_mode == CABAC)
	if (1)
	{
		// create all context models
		currSlice->mot_ctx = create_contexts_MotionInfo();
		currSlice->tex_ctx = create_contexts_TextureInfo();
	}
	currSlice->max_part_nr = 3;  //! assume data partitioning (worst case) for the following mallocs()
	currSlice->partArr = AllocPartition(currSlice->max_part_nr);
}


/*!
************************************************************************
* \brief
*    Memory frees of the Slice structure and of its dependent
*    data structures
*
* \par Input:
*    Input Parameters struct inp_par_dec *inp,  struct img_par *img
************************************************************************
*/
void free_slice_dec(struct inp_par_dec *inp, struct img_par *img)
{
	Slice *currSlice = img->currentSlice;
	
	FreePartition (currSlice->partArr, 3);
	//      if (inp->symbol_mode == CABAC)
	if (1)
	{
		// delete all context models
		delete_contexts_MotionInfo(currSlice->mot_ctx);
		delete_contexts_TextureInfo(currSlice->tex_ctx);
	}
	free(img->currentSlice);
	
	currSlice = NULL;
}

// /*!
// ************************************************************************
// * \brief
// *    Dynamic memory allocation of frame size related global buffers
// *    buffers are defined in global.h, allocated memory must be freed in
// *    void free_global_buffers()
// *
// *  \par Input:
// *    Input Parameters struct inp_par_dec *inp, Image Parameters struct img_par *img
// *
// *  \par Output:
// *     Number of allocated bytes
// ***********************************************************************
// */
// int init_global_buffers_dec()
// {
// 	int memory_size=0;
// 	
// 	if (global_init_done)
// 	{
// 		free_global_buffers_dec();
// 	}
// 	
// 	// allocate memory for reference frame in find_snr
// 	memory_size += get_mem2D(&imgY_ref, img->height, img->width);
// 	memory_size += get_mem3D(&imgUV_ref, 2, img->height_cr, img->width_cr);
// 	
// 	// allocate memory in structure img
// 	if(((img->mb_data) = (Macroblock *) calloc(img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
// 		no_mem_exit("init_global_buffers: img->mb_data");
// 	
// 	if(((img->intra_block) = (int*)calloc(img->FrameSizeInMbs, sizeof(int))) == NULL)
// 		no_mem_exit("init_global_buffers: img->intra_block");
// 	
// 	memory_size += get_mem2Dint(&(img->ipredmode), 4*img->PicWidthInMbs , 4*img->FrameHeightInMbs);
// 	
// 	memory_size += get_mem2Dint(&(img->field_anchor),4*img->FrameHeightInMbs, 4*img->PicWidthInMbs);
// 	
// 	memory_size += get_mem3Dint(&(img->wp_weight), 2, MAX_REFERENCE_PICTURES, 3);
// 	memory_size += get_mem3Dint(&(img->wp_offset), 6, MAX_REFERENCE_PICTURES, 3);
// 	memory_size += get_mem4Dint(&(img->wbp_weight), 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);
// 	
// 	// CAVLC mem
// 	memory_size += get_mem3Dint(&(img->nz_coeff), img->FrameSizeInMbs, 4, 6);
// 	
// 	memory_size += get_mem2Dint(&(img->siblock),img->PicWidthInMbs  , img->FrameHeightInMbs);
// 	
// 	global_init_done = 1;
// 	
// 	img->oldFrameSizeInMbs = img->FrameSizeInMbs;
// 	
// 	return (memory_size);
// }
// 
/*!
************************************************************************
* \brief
*    Free allocated memory of frame size related global buffers
*    buffers are defined in global.h, allocated memory is allocated in
*    int init_global_buffers()
*
* \par Input:
*    Input Parameters struct inp_par_dec *inp, Image Parameters struct img_par *img
*
* \par Output:
*    none
*
************************************************************************
*/
void free_global_buffers_dec()
{
	free_mem2D (imgY_ref);
	free_mem3D (imgUV_ref,2);
	
	// CAVLC free mem
// 	free_mem3Dint(img->nz_coeff, img->oldFrameSizeInMbs);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	
// 	free_mem2Dint(img->siblock);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	
	// free mem, allocated for structure img
	if (img->mb_data       != NULL) free(img->mb_data);
	
	free (img->intra_block);
	
// 	free_mem2Dint (img->ipredmode);//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

// 	free_mem2Dint(img->field_anchor);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	
// 	free_mem3Dint(img->wp_weight, 2);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 	free_mem3Dint(img->wp_offset, 6);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
// 	free_mem4Dint(img->wbp_weight, 6, MAX_REFERENCE_PICTURES);ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	
	global_init_done = 0;
	
}



