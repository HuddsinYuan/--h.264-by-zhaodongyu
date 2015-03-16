
/*!
 ***************************************************************************
 * \file ratectl.c
 *
 * \brief
 *    Rate Control algorithm
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *     - Siwei Ma <swma@jdl.ac.cn>
 *     - Zhengguo LI<ezgli@lit.a-star.edu.sg>
 *
 * \date
 *   16 Jan. 2003
 **************************************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "ratectl.h"


const double THETA=1.3636;//yyn注：在帧层码率控制前编码阶段计算P帧复杂度时用的常数。
const int Switch=0;

int Iprev_bits=0;
int Pprev_bits=0;


/* rate control variables */
int Xp, Xb;// 分别为 P 图像和 B 图像的图像复杂度（定义为图像编码比特数与图像所有基本单元平均量化参数之积）？？？？？
static int R,T_field;
static int Np, Nb, bits_topfield, Q;
long T,T1;
//HRD consideration
long UpperBound1, UpperBound2, LowerBound;
double InitialDelayOffset;
const double OMEGA=0.9;//yyn注：HRD的上限计算式中的常数 ，JVT-H017式(4)

double Wp,Wb;  // 分别为 P 图像和 B 图像的平均图像复杂度（定义为图像编码比特数与图像所有基本单元平均量化参数之积）
int TotalPFrame;
int DuantQp; 
int PDuantQp;
FILE *BitRate;
double DeltaP;//yyn注：基本单元层码率控制量化参数调整时使用的△，

// Initiate rate control parameters初始化码率控制参数
void rc_init_seq()
{
   double L1,L2,L3,bpp;//yyn注:GOP层码率控制时初始化QP0时使用，bpp为每图像比特数//
   int qp;
   int i;
  
   Xp=0;//用于复杂度Wp和Wb
   Xb=0;
   
   bit_rate=input->bit_rate;//输出比特率 //yyn注：bit_rate信道带宽
   frame_rate = (float)(img->framerate *(input->successive_Bframe + 1)) / (float) (input->jumpd + 1);///实际桢速率
                                                                        //successive_Bframe是后继B的数目.jumpd表示非参考帧可以丢弃的数目.
   PreviousBit_Rate=bit_rate;
   
   /*compute the total number of MBs in a frame*///计算1桢宏块数

   img->Frame_Total_Number_MB=img->height*img->width/256;//yyn注：图像中当前帧的宏块总数=图像的高x图像的宽/256
   if(input->basicunit>img->Frame_Total_Number_MB)//yyn注：如果输入文件的基本单元数>帧中宏块总数，则令基本单元数=帧中宏块总数
	   input->basicunit=img->Frame_Total_Number_MB;//基本单元宏块数目如果等于桢宏块数目 基本单元速率控制 用图像层速率控制
   if(input->basicunit<img->Frame_Total_Number_MB) //yyn注：如果输入文件的基本单元数<帧中宏块总数，则总基本单元数=帧中宏块
     TotalNumberofBasicUnit=img->Frame_Total_Number_MB/input->basicunit;//总数/输入文件的基本单元数。以上两句对应G012的公式
   
   MINVALUE=4.0;//在基本单元层 分配基本单元的位数R R= MAX(R, (int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));
               // TotalNumberofBasicUnit是1桢的BU总数 如果第1个P的纹理编码数目小于门限bit_rate/(MINVALUE*frame_rate考虑丢弃

   /*initialize the parameters of fluid flow traffic model流量往返模型的初始化参数*/
   
   BufferSize=bit_rate*2.56;//设定buffer的大小 //yyn注：缓冲器容量=比特率*2.56
   CurrentBufferFullness=0;//当前buffer满度0 //yyn注：当前缓冲器充盈度=0。此时初始化
   GOPTargetBufferLevel=CurrentBufferFullness; //yyn注：GOP目标缓冲级=当前缓冲器充盈度。
   /*HRD consideration*/
   InitialDelayOffset=BufferSize*0.8; //yyn注：初始化延迟偏移=缓冲容量*80%//用于UpperBound1=(long)(R+InitialDelayOffset);规定下hrd设想参考解码器的上限
   
   /*initialize the previous window size*/ //**yyn注：窗口尺寸用于在帧层码率控制后编码阶段时更新MAD的系数和R-D参数。//初始化之前窗尺寸  因为还未编码 所以都设置为0
   m_windowSize=0;
   MADm_windowSize=0;
   img->NumberofCodedBFrame=0;///已编码的B为0
   img->NumberofCodedPFrame=0;
   img->NumberofGOP=0;//gop序号为0
   /*remaining # of bits in GOP *///R表示gop中剩余的位数
   R = 0;//当然R<0说明流量小于了编码位数 初始为0 因为还没编码
   /*control parameter *///在图象层控制参数第i个gop中的第j个图象的目标位分配基于目标buffer等级 流出量 实际buffer占用率
   //T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5); 
   //有B图象的桢 gamma=0.25  无B图象的桢  gamma=0.5
   //beta是平滑细数 T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5)); T 图象目标位数
   //有B桢的时候beta=.5 无B桢的时候beta=.9

   if(input->successive_Bframe>0) //yyn注：如果连续B帧>0，即有B帧时
   {
     GAMMAP=0.25;//**yyn注：在帧层码率控制微观控制时，为了减小实际缓冲充盈度与目标缓存等级之间的差值，所使用的一个调整系数，没有B帧时为0.75，有B帧时为0.25。JVT-G012文档的式（16）
     BETAP=0.9;//yyn注：计算目标比特数时的加权模型的系数。即JVT-G012文档的式（19）
   }
   else
   {
     GAMMAP=0.5;
     BETAP=0.5;
   }
   
   /*quadratic rate-distortion model*///2次r-d模型参数初始值,Ti=X1xMAD1/Qi+X2xMAD2/Qi2    Ti用来编码当前第i帧图像的比特数，Qi当前第i帧的量化级别，MADi当前第i帧图像的MAD值。
   PPreHeader=0;
   
   Pm_X1=bit_rate*1.0; //yyn注：二次率失真模型系数//这里如何取初始值地 搞不清楚 如此粗燥的取流量 考虑量化参数和MAD值了不 难道都是无穷小相除极限是1
   Pm_X2=0.0;
   /* linear prediction model for P picture*///线性mad预测参数初始值
   PMADPictureC1=1.0;
   PMADPictureC2=0.0;
   //这里窗口尺寸20是为了平滑一个突然的场景改变，所以采用以前编码20个图像做参考
  // 窗口最大尺寸是20一般由MADC*20/MAD_prev决定
   for(i=0;i<20;i++)
   {
     Pm_rgQp[i]=0;
     Pm_rgRp[i]=0.0;
     PPictureMAD[i]=0.0;
   }
   PPictureMAD[20]=0.0;

   //yyn注：Pm_rgQp、Pm_rgRp、Pm_X1、Pm_X2 的含义是：
   //double Pm_rgQp[20];    //++ 参数值传递过程中的中间临时变量，可直接用 m_rgQp 替换
   //double Pm_rgRp[20];    //++ 参数值传递过程中的中间临时变量，可直接用 m_rgRp 替换
   //double Pm_X1;         //++ 参数值传递过程中的中间临时变量，可直接用 m_X1 替换
   //double Pm_X2;         //++ 参数值传递过程中的中间临时变量，可直接用 m_X2 替换
   
   //而 m_rgQp、m_rgRp、m_X1、m_X2 的含义是：
   //double  m_rgQp[21];     //++ FIFO 队列用来存储各个基本单元的量化步长
   //double m_rgRp[21];      //++ FIFO 队列用来存储各个基本单元编码完成后的二次方程左边项
   // double m_X1;          //++ 二次模型第一个系数
   // double m_X2;          //++ 二次模型第二个系数

   //Define the largest variation of quantization parameters得到DDquant的中间变量pDuantQp
   PDuantQp=2;//估计可用于对图像Qpc=min{Qpp+2,max(Qpp-2,Qpc)}进行clip，保证图像层Qpc在之前图像Qpp+-2范围内
   
   /*basic unit layer rate control*/
   //基本单元用到地参数 如果基本单元数目大于8 DDquant=2 小于取1，因为基本单元越多，希望量化步长越大 这样编码比特降低  
   //PAveHeaderBits头部的比特数，因为r-d模型是去除头部信息的模型
    PAveHeaderBits1=0;//yyn注：PAveHeaderBits1表示当前帧图像至前一个basic unit为止的平均头信息比特数。
    PAveHeaderBits3=0;      //yyn注：PAveHeaderBits3表示前一帧图像的平均头信息比特数。
//yyn注：PAveHeaderBits2表示当前basic unit的平均头信息比特数的估计值。

    if(TotalNumberofBasicUnit>=9)
      DDquant=1;
    else
      DDquant=2;
	//求出1行的宏块数
	MBPerRow=input->imagewidth/16;//为了第2 3种情况，第1种是画框第1个基本单元取前一图像所有基本单元平均值
	//第2 种在R<0时候 说明编码的位数没完全被流出 这样量化参数要大于前1量化单元，这样可以减少位数
		//基本单元宏块数目小于MBPerRow，m_Qc = MIN(m_Qc, PAveFrameQP+3);其他m_Qc = MIN(m_Qc, PAveFrameQP+6);
	//为何如此呢.降低方块效应.平滑图像.让QP增大
    
    /*adaptive field/frame coding*/
	img->FieldControl=0;//取0代表场控制无效 ，只有幀编码.qp范围（0-51）
    
    RC_MAX_QUANT = 51;  // clipping
    RC_MIN_QUANT = 0;//clipping
    
    /*compute thei initial QP*/
	//在第1个gop预先定义量化参数Qp 在I和第1张P就是用GOP预先定义的QP
    bpp = 1.0*bit_rate /(frame_rate*img->width*img->height);
	if (img->width == 176) //如果是Qcif格式门限L1 L2 L3取值如下
    {
      L1 = 0.1;
      L2 = 0.3;
      L3 = 0.6;
	}else if (img->width == 352)//如果是cif格式门限L1 L2 L3取值如下
    {
      L1 = 0.2;
      L2 = 0.6;
      L3 = 1.2;
	}else //其他格式取值
	{
      L1 = 0.6;
      L2 = 1.4;
      L3 = 2.4;
    }
	//根据bpp门限求出对应的初始Qp
    if (input->SeinitialQP==0)
    {
      if(bpp<= L1)
        qp = 35;
      else
        if(bpp<=L2)
          qp = 25;
        else
          if(bpp<=L3)
            qp  = 20;
          else
            qp =10;
          input->SeinitialQP = qp;
    }
}

// Initiate one GOP目的就是求出GOP的QP初始值
void rc_init_GOP(int np, int nb)//输入未编码的p数目np 未编码的b数目nb
{
  Boolean Overum=FALSE1;
  int OverBits;
  int OverDuantQp;
  int AllocatedBits;
  int GOPDquant;  //yyn注：GOP量化步长

  /*check if the last GOP over uses its budget. If yes, the initial QP of the I frame in 
 the coming  GOP will be increased.*/

  if(R<0)//当然R<0说明流量小于了编码位数
    Overum=TRUE1;
  OverBits=-R;

  /*initialize the lower bound and the upper bound for the target bits of each frame, HRD consideration*/
  //考虑到hrd（假设参考解码）来决定R中可用位数的上限和下限
  LowerBound=(long)(R+bit_rate/frame_rate);
  UpperBound1=(long)(R+InitialDelayOffset);

  /*compute the total number of bits for the current GOP*/ //计算当前gop的可用位
  AllocatedBits = (int) floor((1 + np + nb) * bit_rate / frame_rate + 0.5);//np nb未编码的p b桢 按当前流量流出的位数
  R +=AllocatedBits;//R表示当前gop可用位数 考虑gop的桢数 算出流量 
  //和以前R相加因为以前R<0可能存在，说明上1个gop还有未流出的数据 必须降低当前gop目标位数 
  Np  = np;
  Nb  = nb;

  OverDuantQp=(int)(8*OverBits/AllocatedBits+0.5);//原因就是R<0未流出的量 求出修正OverDuantQp 8*Tr(i-1,Ngop)/Tr(i,0)上次剩余的位数处以当前分配的位数
  GOPOverdue=FALSE1;//用来控制gop是否迟到，过期 可能由于网络原因有些gop阻塞了 后来才收到
  
/*field coding*/
  img->IFLAG=1;//场

/*Compute InitialQp for each GOP*/
  TotalPFrame=np;//P桢数目等于未编码的np数目，因为调用此函数 gop内搜有图像都未编码
  img->NumberofGOP++;//gop数目加1
  if(img->NumberofGOP==1)//第1个GOP取由bpp初始设定的QP值
  {
    MyInitialQp=input->SeinitialQP;
    PreviousQp2=MyInitialQp-1; //recent change -0;
    QPLastGOP=MyInitialQp;
  
  }
  else
  {//其他GOP的Qp计算公式如下QP=SumPqp/Np-8*Tr(i-1,Ngop)/Tr(i,0)-min{2,Ngop/15}
	//  =img->TotalQpforPPicture/img->NumberofPPicture-OverDuantQp-min{2,Ngop/15}

// 	  /*adaptive field/frame coding*///自适应桢场编码
    if((input->PicInterlace==ADAPTIVE_CODING)\
		||(input->MbInterlace))//如何是图像自适应桢场编码，或者宏块是场模式或者宏块是桢场自适应模式
    {
		if (img->FieldFrame == 1)//如果是幀编码  TotalQpforPPicture加上桢的FrameQPBuffer TotalQpforPPicture表示所有P图像的QP和
      {
        img->TotalQpforPPicture += FrameQPBuffer;
        QPLastPFrame = FrameQPBuffer;
      }
		else//场编码
      {
        img->TotalQpforPPicture += FieldQPBuffer;
        QPLastPFrame = FieldQPBuffer;
      }
      
    }
	/*compute the average QP of P frames in the previous GOP*///计算之前GOP所有P图像中的QP平均值 NumberofPPicture表示P图像数目
    PAverageQp=(int)(1.0*img->TotalQpforPPicture/img->NumberofPPicture+0.5);

    GOPDquant=(int)(0.5+1.0*(np+nb+1)/15);
    if(GOPDquant>2)
        GOPDquant=2;

    PAverageQp-=GOPDquant;

	if (PAverageQp > (QPLastPFrame - 2))//为了图像相关性考虑 如果平均的QP大于之前P图像的QP-2，减少PAverageQp
      PAverageQp--;
	PAverageQp = MAX(QPLastGOP-2,  PAverageQp);//限定范围在上个GOP的QP值+-2范围内
    PAverageQp = MIN(QPLastGOP+2, PAverageQp);
	PAverageQp = MIN(RC_MAX_QUANT, PAverageQp);//限定范围在0-51
    PAverageQp = MAX(RC_MIN_QUANT, PAverageQp);
  

    MyInitialQp=PAverageQp;
	QPLastGOP = MyInitialQp;//更新上一个gop的QP
    Pm_Qp=PAverageQp;
    PAveFrameQP=PAverageQp;
    PreviousQp1=PreviousQp2;
    PreviousQp2=MyInitialQp-1;  
  }

  img->TotalQpforPPicture=0;//清TotalQpforPPicture为当前GOP重新QP重新累加用
  img->NumberofPPicture=0;//当前GOP还没编码 编码的P图像数目 NumberofPPicture为0
  NumberofBFrames=0; //当前GOP还没编码 编码的B图像数目 NumberofPPicture为0
}


//初始化图像层输入fieldpic＝０场编码图像　１　桢编码图像　topfield １　是顶场
void rc_init_pict(int fieldpic,int topfield,int targetcomputation)//targetcomputation是否计算目标缓冲等级为1表示计算
{
  int i;

/*compute the total number of basic units in a frame*///计算１桢中BU（基本单元）的数目，用１桢中总宏块数目除BU中宏块数

 // 这里说下图像img->BasicUnit和输入参数　input->basicunit区别毕竟在
	//  在场模式和桢场自适应模式
	 // img->BasicUnit=input->basicunit*2;
  //只有桢模式
	 // img->BasicUnit=input->basicunit;
  //也就是说在场模式和桢场自适应模式中BU中宏块数目是２倍只有在桢模式中BU宏块数目

  if(input->MbInterlace)//在场模式和桢场自适应模式
    TotalNumberofBasicUnit=img->Frame_Total_Number_MB/img->BasicUnit;
  img->NumberofCodedMacroBlocks=0;//当前图像还未编码　所以编码的宏块数目为０

/*Normally, the bandwith for the VBR case is estimated by 
a congestion control algorithm. A bandwidth curve can be predefined if we only want to 
test the proposed algorithm*/
  //这里为了做VBR（变速率）逻辑测试　定义个带宽的曲线设置　当第５９个已编码P图像时候速率　上升到1.5倍　当第６０个P图像速率恢复以前
  if(input->channel_type==1)
  {
    if(img->NumberofCodedPFrame==58)
      bit_rate *=1.5;
    else if(img->NumberofCodedPFrame==59)
      PreviousBit_Rate=bit_rate;
  }

  /*predefine a target buffer level for each frame*///为每桢图像.预先设置目标缓冲水平

  if((fieldpic||topfield)&&targetcomputation)
  {
    switch (img->type)
    {
	case P_SLICE://P图像

		 //   因为带宽是实时变化地，所以R中可用位数必须１个图像１个图像地更新
		 //	  按照原理公式Tr(i,j)=Tr(i,j-1)+ (u(i,j)-u(i,j-1))*(Ngop-j)/Fr-b(i,j-1)　
		 //   i代表第i个gop　j代表第j个图像Tr是gop中可分配的位元，u是带宽,Fr桢速率，b是已编码产生的位数
		 //   上面公式这样理解，当前图像可用的位数，等于上一个图像可用位数减去带宽变化和上个图像产生的编码位数
       /*Since the available bandwidth may vary at any time, the total number of 
       bits is updated picture by picture*/
        if(PreviousBit_Rate!=bit_rate)
			R +=(int) floor((bit_rate-PreviousBit_Rate)*(Np+Nb)/frame_rate+0.5);//是 (u(i,j)-u(i,j-1))*(Ngop-j)/Fr
              
       /* predefine the  target buffer level for each picture.
       frame layer rate control*///设定桢层速率控制为每个图像设定缓冲水平TargetBufferLevel
	   //	Tbl(i,2)=Bc(i,2)      Bc(i,2)是第i个gop第1个P图像编码完毕后实际缓冲占有率
	   //   TBL(i,j+1)=TBL(i,j)-Tbl(i,2)/Np-1+AWp(i,j)*(L+1)*u(i,j)/Fr*AWp(i,j)*AWb(i,j)*L-u(i,j)/Fr AWp AWb为了计算复杂度的权

		if(img->BasicUnit==img->Frame_Total_Number_MB)//如果BU内宏块数目等于桢中宏块数目
        {
			if(img->NumberofPPicture==1)//已经编码一个P图像 目标缓冲水平等于当前缓冲占有率 Tbl(i,2)=Bc(i,2)
          {
			//    因为起始CurrentBufferFullness=0 GOPTargetBufferLevel=0 更新  CurrentBufferFullness += nbits - bit_rate/frame_rate;编码位数减去流出
			//	  其实也可以CurrentBufferFullness=BufferSize/8 GOPTargetBufferLevel=BufferSize/8 在带宽波动的情况下留点裕度
			//	  这样可以保证流体流动通信模型（fluid flow traffic model）HRD缓存是不会向上或者向下溢出的

            TargetBufferLevel=CurrentBufferFullness;
			DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/(TotalPFrame-1);//这里DeltaP=-Tbl(i,2)/Np-1 GOPTargetBufferLevel=0不考虑了
			TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
          }
          else if(img->NumberofPPicture>1)
			  TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
        }
        /*basic unit layer rate control*/
        else
        {
			if(img->NumberofCodedPFrame>0)//已编码P图像的数目大于0
          {
            /*adaptive frame/filed coding*/
            if(((input->PicInterlace==ADAPTIVE_CODING)||(input->MbInterlace))\
				&&(img->FieldControl==1))//在图像桢场自适应 和宏块场编码或宏块桢场自适应且FieldControl=1场编码情况 一句话宏块场编码情况
			{//TotalNumberofBasicUnit是图像中BU数目注意这里的图像即可以指场图像也可以指桢图像
				for(i=0;i<TotalNumberofBasicUnit;i++)//FCBUPFMAD是之前场控制编码mad,FCBUCFMAD是当前场控制编码mad
					FCBUPFMAD[i]=FCBUCFMAD[i];//更新之前场控制编码mad,因为该图像已经编码完毕
            }
            else
            {
				for(i=0;i<TotalNumberofBasicUnit;i++)//桢编码模式 BUPFMAD是之前mad,BUCFMAD是当前mad
                BUPFMAD[i]=BUCFMAD[i];
            }     
          }

			if(img->NumberofGOP==1)//第1个gop第2个p图像开始计算tbl
          {
			  if(img->NumberofPPicture==1)//已经编码一个P图像 目标缓冲水平等于当前缓冲占有率 Tbl(i,2)=Bc(i,2)
            {
			//	    因为起始CurrentBufferFullness=0 GOPTargetBufferLevel=0 更新  CurrentBufferFullness += nbits - bit_rate/frame_rate;编码位数减去流出
			//		其实也可以CurrentBufferFullness=BufferSize/8 GOPTargetBufferLevel=BufferSize/8 在带宽波动的情况下留点裕度
			//		这样可以保证流体流动通信模型（fluid flow traffic model）HRD缓存是不会向上或者向下溢出的   
              TargetBufferLevel=CurrentBufferFullness;
			  DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/(TotalPFrame-1);//这里DeltaP=-Tbl(i,2)/Np-1 GOPTargetBufferLevel=0不考虑了
			  TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
            }
			  else if(img->NumberofPPicture>1)//已当前GOP编码P图像数目大于1个
				TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
          }
			else if(img->NumberofGOP>1)//第2 3....个GOP第1个P开始就计算tbl用来自GOP的I桢编码后的当前缓冲占有率
          {
            if(img->NumberofPPicture==0)
            {
              TargetBufferLevel=CurrentBufferFullness;
              DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/TotalPFrame;
              TargetBufferLevel -=DeltaP;
            }
            else if(img->NumberofPPicture>0)
              TargetBufferLevel -=DeltaP;
          }
        }
		//AWp和AWb的更新公式AWP=AWp/8+7*Wp/8  AWb=AWb/8+7*Wb/8  Wb=b(i,j)*QPb(i,j)/1.3636 Wp=b(i,j)*QPp(i,j) 
		if(img->NumberofCodedPFrame==1)//如果第1个P图像已编码AWp=Wp
          AWp=Wp;
        //如果以编码P图像小于8且大于1 用已编码图像1/NumberofCodedPFrame做加权 毕竟已编码P数目小于8,用1/8加权平滑不合适 
        if((img->NumberofCodedPFrame<8)&&(img->NumberofCodedPFrame>1))
            AWp=Wp*(img->NumberofCodedPFrame-1)/img->NumberofCodedPFrame+\
              AWp/img->NumberofCodedPFrame;
		else if(img->NumberofCodedPFrame>1)//已编码P图像大于8的话用已编码图像1/8做加权AWP=AWp/8+7*Wp/8
            AWp=Wp/8+7*AWp/8;
          
        //compute the average complexity of B frames
		//这里tbl+=AWp(i,j)*(L+1)*u(i,j)/Fr*AWp(i,j)*AWb(i,j)*L  L表示后续B桢数目 理解成2个P图像之间B的数目就好
		//input->successive_Bframe就是L
        if(input->successive_Bframe>0)
        {
          //compute the target buffer level
          TargetBufferLevel +=(AWp*(input->successive_Bframe+1)*bit_rate\
            /(frame_rate*(AWp+AWb*input->successive_Bframe))-bit_rate/frame_rate);
        }
        
        break;

         case B_SLICE:
         /* update the total number of bits if the bandwidth is changed*/
			 //按照原理公式Tr(i,j)=Tr(i,j-1)+ (u(i,j)-u(i,j-1))*(Ngop-j)/Fr-b(i,j-1)
			 if(PreviousBit_Rate!=bit_rate)//带宽变动情况下Gop中可用位数R必须修正,修正值u(i,j)-u(i,j-1))*(Ngop-j)/Fr
             R +=(int) floor((bit_rate-PreviousBit_Rate)*(Np+Nb)/frame_rate+0.5);
			 if((img->NumberofCodedPFrame==1)&&(img->NumberofCodedBFrame==1))//编码第1个P和编码第1个B图像后如下初始化复杂度参数AWp

          {
            AWp=Wp;
            AWb=Wb;
          }
			 else if(img->NumberofCodedBFrame>1)//如果已编码的B数目大于1
          {
            //compute the average weight
			  //AWb的更新公式AWb=AWb/8+7*Wb/8 在已编码的B数目大小于8 加权系数取1/NumberofCodedBFrame
            if(img->NumberofCodedBFrame<8)
              AWb=Wb*(img->NumberofCodedBFrame-1)/img->NumberofCodedBFrame+\
                AWb/img->NumberofCodedBFrame;
            else
              AWb=Wb/8+7*AWb/8;
          }

            break;
    }
     /*Compute the target bit for each frame*/
	// 计算每个桢的目标位数,为什么要这一步呢,理解如下
    //	如果实际缓冲的等级等同于我们的目标缓冲等级,当然可以保证每个GOP都精确分配到位数,
	//	但是R-D模型和mad预测模型是不精确地,导致实际缓冲水平和目标位元等级的误差较大,
	//	那么我们需要计算每张图像的目标位数用来降低实际缓冲水平和目标位元之间的误差,这就是细节控制
	//	如同往杯子里面倒水,如果要正好倒1升水, 当然1滴1滴倒最好了!
	if(img->type==P_SLICE)//P图像模式
    {
      /*frame layer rate control*/
		if(img->BasicUnit==img->Frame_Total_Number_MB)//BU中宏块数目等于桢中宏块数目
      {
		  if(img->NumberofCodedPFrame>0)//已编码的P数目大于1
        {

		//	    通过Gop中R可用位数计算T是图像可用位数,公式如下
		//		f(i,j)=Wp(i,j-1)*Tr/(Wp(i,j-1)*Np+Wb(i,j-1)*Nb) 
		//		也就是T = Wp*R/(Np*Wp+Nb*Wb)

			T = (long) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);//T是图像可用位数,当然是Gop中R可用位数
		//	    下面考虑到第i个Gop的第J个图像的目标分配位数由  
		//		目标缓冲等级Tbl=TargetBufferLevel 实际缓冲区占有率Bc(i,j)=CurrentBufferFullness 频宽u(i,j)=bit_rate    
		//		GAMMAP在有B图像的情况下GAMMAP=0.25;  无B图像的情况下GAMMAP=0.5;
		//	    为什么有B图像比无B图像GAMMAP要小,我这样理解,B编码产生位数少,当然需要减去的实际和目标的误差系数小
		//		公式f1=u(i,j)/Fr+gamma*(Tbl(i,j)-Bc(i,j))

                
          T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1=MAX(0,T1);
		//  目标可用位数f(i,j)由f(i,j)和f1加权得到
		//	  f(i,j)=beta*f(i,j)+(1-beta)*f1
		//	  beta在没有B图像情况下beta=0.5          有B图像的情况下是0.9这说明有B图像情况主要由f(i,j)自身决定
          T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
       }
		/*basic unit layer rate control*///基本单元层控制
      else
      {
		  if((img->NumberofGOP==1)&&(img->NumberofCodedPFrame>0))//第1个GOP并且已编码P数目大于0
        {
			//  通过Gop中R可用位数计算T是图像可用位数,公式如下
			//	:f(i,j)=Wp(i,j-1)*Tr/(Wp(i,j-1)*Np+Wb(i,j-1)*Nb) 

          T = (int) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (int) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1=MAX(0,T1);
          T = (int)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
		  else if(img->NumberofGOP>1)//第2  3  .....个GOP
        {
          T = (long) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1 = MAX(0,T1);
          T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
      }

		/*reserve some bits for smoothing*///保留一些位数为了平滑
		//	很荒唐吧0.0*input->successive_Bframe估计0.0是个可调系数 调整T中可用位数

      T=(long)((1.0-0.0*input->successive_Bframe)*T);
	  /*HRD consideration*///保证图像可用位数T永远在[LowerBound,UpperBound2]范围内
	  //这里的LowerBound和UpperBound2在图像编码结束后会调用rc_update_pict更新地
	  //通过流出减去编码位数,如果流出大于流入(图像编码位数),瓶子空的更多了 表示有更多可用的位数 那么当前图像可用位数上下限应该增加
	  //如果流出大于流入(图像编码位数),表示有可用的位数减少 瓶子空间更少了 那么当前图像可用位数上下限应该减少       
	  //LowerBound  +=(long)(bit_rate/frame_rate-nbits);
	  //UpperBound1 +=(long)(bit_rate/frame_rate-nbits); UpperBound2 = (long)(OMEGA*UpperBound1);         OMEGA=0.9保证裕度
      T = MAX(T, (long) LowerBound);
        T = MIN(T, (long) UpperBound2);

      if((topfield)||(fieldpic&&((input->PicInterlace==ADAPTIVE_CODING)\
		  ||(input->MbInterlace))))//场编码情况
        T_field=T;
    }
  }

  if(fieldpic||topfield)//桢编码或者顶场编码
  {
	  /*frame layer rate control*///图像层速率控制 也叫桢层速率控制
	  img->NumberofHeaderBits=0;//还未编码 先把头信息位数NumberofHeaderBits和纹理信息位数NumberofTextureBits清0
    img->NumberofTextureBits=0;

	/*basic unit layer rate control*///只有在BU中宏块小于桢中宏块数目才存在BU基本单元控制
    if(img->BasicUnit<img->Frame_Total_Number_MB)
    {
		TotalFrameQP=0;//图像内所有BU的QP和
		img->NumberofBasicUnitHeaderBits=0;//头信息位数NumberofHeaderBits
		img->NumberofBasicUnitTextureBits=0;//纹理信息位数NumberofTextureBits
		img->TotalMADBasicUnit=0;//图像内所有BU的MAD和TotalMADBasicUnit
	  if(img->FieldControl==0)//桢编码
        NumberofBasicUnit=TotalNumberofBasicUnit;
      else//场编码
        NumberofBasicUnit=TotalNumberofBasicUnit/2;
    }
  }
    
  if((img->type==P_SLICE)&&(img->BasicUnit<img->Frame_Total_Number_MB)\
    &&(img->FieldControl==1))
  {
  /*top filed at basic unit layer rate control*///顶场编码在BU层速率控制
    if(topfield)
    {
      bits_topfield=0;
      T=(long)(T_field*0.6);//顶场可用位数//顶场和底场的一共的可用位数T_field 0.6系数说明顶场编码重要性 可能作参考场
    }
  /*bottom filed at basic unit layer rate control*/
    else
    {
      T=T_field-bits_topfield;//场可用的位数T_field减去顶场编码产生的位数bits_topfield得到底场可用的位数
      img->NumberofBasicUnitHeaderBits=0;//头信息位数NumberofHeaderBits
      img->NumberofBasicUnitTextureBits=0;//纹理信息位数NumberofTextureBits
      img->TotalMADBasicUnit=0;//图像内所有BU的MAD和TotalMADBasicUnit
      NumberofBasicUnit=TotalNumberofBasicUnit/2;//场编码每个场中BU数目自然为桢BU数目的1半
    }
  }
}

//calculate MAD for the current macroblock 计算当前宏块的MAD,diff是像素Y原值减去预测值
// MAD=sum(abs(diff))/256
double calc_MAD()
{
  int k,l;
    int s = 0;
  double MAD;

  for (k = 0; k < 16; k++)
    for (l = 0; l < 16; l++)
      s+= abs(diffy[k][l]);
  
  MAD=s*1.0/256;
  return MAD;
}

// update one picture after frame/field encoding在桢场编码后更新一个图像参数
void rc_update_pict(int nbits)
{
  R-= nbits; /* remaining # of bits in GOP *///GOP中剩余可用位数R,必须减去上个图像编码产生位数
  CurrentBufferFullness += nbits - bit_rate/frame_rate;//当前缓冲满度用流入减去流出修正

  /*update the lower bound and the upper bound for the target bits of each frame, HRD consideration*/
  //通过流出减去编码位数,如果流出大于流入(图像编码位数),瓶子空的更多了 表示有更多可用的位数 那么当前图像可用位数上下限应该增加
  //如果流出大于流入(图像编码位数),表示有可用的位数减少 瓶子空间更少了 那么当前图像可用位数上下限应该减少
  LowerBound  +=(long)(bit_rate/frame_rate-nbits);
  UpperBound1 +=(long)(bit_rate/frame_rate-nbits);
  UpperBound2 = (long)(OMEGA*UpperBound1);
  
  return;
}

// update after frame encoding在桢编码后更新RC
void rc_update_pict_frame(int nbits)
{

/*update the
complexity weight of I, P, B frame*///更新I P  B的复杂度权
  int Avem_Qc;// Avem_Qc是所有BU的QP平均值
  int X;
    
/*frame layer rate control*///图象层速率控制
// I P  B的复杂度权Wb=b(i,j)*QPb(i,j)/1.3636 Wp=b(i,j)*QPp(i,j)  nbits=b(i,j)

  if(img->BasicUnit==img->Frame_Total_Number_MB)//BU的宏块数目等于桢内宏块数目
    X = (int) floor(nbits*m_Qc+ 0.5);//注意是m_Qc nbits编码产生位数
/*basic unit layer rate control*///基本单元层速率控制
  else
  {
    if(img->type==P_SLICE)//P图像模式
    {
      if(((img->IFLAG==0)&&(img->FieldControl==1))         //\场编码且不是第1个GOP的第1个I场图像           IFLAG=1可能代表seq的第1个场
        ||(img->FieldControl==0))//桢编码
      {
        Avem_Qc=TotalFrameQP/TotalNumberofBasicUnit;//求出平均基本单元Qc
        X=(int)floor(nbits*Avem_Qc+0.5);
      }
    }
    else if(img->type==B_SLICE)//B图像模式
      X = (int) floor(nbits*m_Qc+ 0.5);
  }

//I P  B的复杂度权Wb=b(i,j)*QPb(i,j)/1.3636 Wp=b(i,j)*QPp(i,j)  nbits=b(i,j)
  switch (img->type)
  {
    case P_SLICE:
 /*filed coding*/
      if(((img->IFLAG==0)&&(img->FieldControl==1))      //\ 场编码且不是第1个GOP的第1个I场图像
        ||(img->FieldControl==0))//桢图像
      {
        Xp = X;
        Np--;//已编码P图像了所以未编码P图像减1
        Wp=Xp;
        Pm_Hp=img->NumberofHeaderBits;//头信息位数
        img->NumberofCodedPFrame++;//已编码P图像数目NumberofCodedPFrame
        img->NumberofPPicture++;//已编码P图像数目NumberofPPicture

      }
      else if((img->IFLAG!=0)&&(img->FieldControl==1))//第1个GOP的第1个I场图像已经编码,那么此标志清0
        img->IFLAG=0;
        break;
        case B_SLICE:
        Xb = X;
      Nb--;//已编码B图像了所以未编码P图像减1
        Wb=Xb/THETA; //来自Wb=b(i,j)*QPb(i,j)/1.3636 THETA=1.3636
      img->NumberofCodedBFrame++;//已编码B图像数目NumberofCodedBFrame
      NumberofBFrames++;//已编码P图像数目NumberofCodedBFrame

        break;
    }
}

// coded bits for top field顶场编码产生的位数bits_topfield
void setbitscount(int nbits)
{
  bits_topfield = nbits;
}

//compute a  quantization parameter for each frame计算每桢量化参数

int updateQuantizationParameter(int topfield)
{
  double dtmp;//记得初中学过的二次方程式吧 其实就是delta=sqrt(b*b-4*a*c)

  int m_Bits;
  int BFrameNumber;
  int StepSize;
  int PAverageQP;
  int SumofBasicUnit;
  int i;
  
/*frame layer rate control*///桢层速率控制

  if(img->BasicUnit==img->Frame_Total_Number_MB)//BU内宏块数目等于桢内宏块数目 自然不存在BU控制问题啦

  {
/*fixed quantization parameter is used to coded I frame, the first P frame and the first B frame
  the quantization parameter is adjusted according the available channel bandwidth and 
  the type of vide*/  
/*top field*/
    if((topfield)||(img->FieldControl==0))//是顶场或者场控制信号为0

    {
      if(img->type==I_SLICE)//I条带情况 m_qc取初始值

      {
        m_Qc=MyInitialQp;
        return m_Qc;
      }
      else if(img->type==B_SLICE)//B条带

      {
        if(input->successive_Bframe==1)//P与P之间只有1个B情况

        {//successive_Bframe=1对1求余任何整数都是0阿下面的就可笑了,毫无意义

          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
          if(BFrameNumber==0)//在这里又修正BFrameNumber=1说明BFrameNumber表示在2个P图像中的第几个Ｂ图像　真垃圾代码阿

            BFrameNumber=input->successive_Bframe;//因为successive_Bframe＝１说明在2个P图像中只有１个Ｂ图像BFrameNumber只可能为１

          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {//当L=1表示2个P图像之间有1个B图像 公式如下:
        //  QB1=(QP1+QP2+2)/2  当QP1!=QP2
         // QB1=QP1+2          当QP1 =QP2

            if((input->PicInterlace==ADAPTIVE_CODING)\
              ||(input->MbInterlace))
            {
              if(img->FieldControl==0)//场控制关闭

              {                   
                /*previous choice is frame coding*/
                if(img->FieldFrame==1)//之前选择是桢编码

                {
                  PreviousQp1=PreviousQp2;//例如P1B1P2B2P3当前是B1则PreviousQp1代表P1图像,PreviousQp2=P2

                  PreviousQp2=FrameQPBuffer;
                }           
                /*previous choice is field coding*/
                else//之前选择是场编码

                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }
//          m_Qc=QB1=(QP1+QP2+2)/2  当QP1!=QP2
          //          m_Qc=QB1=QP1+2          当QP1 =QP2

          if(PreviousQp1==PreviousQp2)
            m_Qc=PreviousQp1+2;
          else
            m_Qc=(PreviousQp1+PreviousQp2)/2+1;
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping钳位0-51

          m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        else
        {
//L=successive_Bframe>1 下面求余就有意义了 说明BFrameNumber表示在2个P图像中的第几个Ｂ图像
    //      这样举例把L=2 PB1B2PB3B4P那么NumberofBFrames=2已编码的B数目是2那么当前是第3个B实际B3,BFrameNumber=1在实际P与P的位置是1
       //    已编码的B数目是3那么当前是第3个B实际B3,BFrameNumber=2
       //   通过对后续B的数目求余 求得当前B在P 与P的位置


          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
          if(BFrameNumber==0)//上面再NumberofBFrames=3时候求于地BFrameNumber＝０必须修正为２

            BFrameNumber=input->successive_Bframe;
          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {//更新之前两个P图像的QP值因为是第１个B所以需要更新，其他B与第１个B采用的参考P是一致的作考虑

            if((input->PicInterlace==ADAPTIVE_CODING)\
              ||(input->MbInterlace))//图像自适应编码或帧编码

            {
              if(img->FieldControl==0)//场控制信号关闭当前图像可以使桢编码也可以是场编码

              {
                /*previous choice is frame coding*///更新之前两个P图像的QP值

                if(img->FieldFrame==1)//图像选择是桢编码
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FrameQPBuffer;
                }
                /*previous choice is field coding*/
                else//图像选择是场编码

                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }
// 当L>1表示两张P图像插入多张B图像 公式如下:
//         QBi=QP1+alpha+max{min{(QP2-QP1)/(L-1),2*(i-1)},-2*(i-1)}
//         alpha取值如下:
//         alpha=-3 if QP2-QP1<=-2*L-3
//         alpha=-2 if QP2-QP1=-2*L-2
//         alpha=-1 if QP2-QP1=-2*L-1
//         alpha=0  if QP2-QP1=-2*L
//         alpha=1  if QP2-QP1=-2*L+1
//         alpha=2   其他情况
//         为什么需要alpha=StepSize修正QB,估计是按照２个P图像的Qp差deltaQ和L的关系，做个平滑处理，避免图像质量差别过大

          if((PreviousQp2-PreviousQp1)<=(-2*input->successive_Bframe-3))
            StepSize=-3;
          else  if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-2))
            StepSize=-2;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-1))
            StepSize=-1;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe))
            StepSize=0;
          else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe+1))
            StepSize=1;
          else
            StepSize=2;
          
          m_Qc=PreviousQp1+StepSize;
// 下面按照max{min{(QP2-QP1)/(L-1),2*(i-1)},-2*(i-1)}修正　估计保证播放顺序靠近第１个P的修正更加接近第１个P的Qp,
//           接近第２个P图像的修正更加接近第２个P的Qp，也是为了做个平滑处理，避免图像质量差别过大

          m_Qc +=MIN(2*(BFrameNumber-1),MAX(-2*(BFrameNumber-1), \
            (BFrameNumber-1)*(PreviousQp2-PreviousQp1)/(input->successive_Bframe-1)));
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        return m_Qc;
      }
      else if((img->type==P_SLICE)&&(img->NumberofPPicture==0))//当前是P图像并且已编码的P图像数目是０

      {
        m_Qc=MyInitialQp;//MyInitialQp是GOP的QP

        
        if(img->FieldControl==0)//场控制关闭

        {
          if(active_sps->frame_mbs_only_flag)//只有桢模式桢宏块编码

          {
            img->TotalQpforPPicture +=m_Qc;//TotalQpforPPicture表示GOP内已编码P图像Qp和

            PreviousQp1=PreviousQp2;//更新2个P图像Qc

            PreviousQp2=m_Qc;
            Pm_Qp=m_Qc;
          }
          /*adaptive field/frame coding*///自适应桢场模式

          else
            FrameQPBuffer=m_Qc;//FrameQPBuffer桢QP缓冲

        }
        
        return m_Qc;  
      }
      else
      {
        /*adaptive field/frame coding*/
        if(((input->PicInterlace==ADAPTIVE_CODING)//\图像自适应桢场模式

          ||(input->MbInterlace))//\宏块场模式或自适应桢场模式

          &&(img->FieldControl==0))//场控制关闭

        {
          /*previous choice is frame coding*/
          if(img->FieldFrame==1)//桢模式
          {
            img->TotalQpforPPicture +=FrameQPBuffer;//TotalQpforPPicture表示GOP内已编码P图像Qp和

            Pm_Qp=FrameQPBuffer;
          }
          /*previous choice is field coding*/
          else//场模式

          {
            img->TotalQpforPPicture +=FieldQPBuffer;
            Pm_Qp=FieldQPBuffer;
          }
        }
// 下面是2次r-d模型系数 具体模型推导 由taylor级数推导 以后有时间相信写个RD模型推导
//         (R-H)/MAD=C1/Qstep+C2/(Qstep*Qstep) 这里m_X1 m_X2就是C1 C2
//         这里是(T-H)/MAD=MADPictureC1/Qstep+MADPictureC2/(Qstep*Qstep) 

        m_X1=Pm_X1;
        m_X2=Pm_X2;
        m_Hp=PPreHeader;//之前的头信息位数

        m_Qp=Pm_Qp;//m_QP是m_可能代表宏块，因为BU宏块等于桢宏块数目，所以m_qp等于P图像Pm_Qp

        DuantQp=PDuantQp;
        MADPictureC1=PMADPictureC1;//PMADPictureC1是预测MAD

        MADPictureC2=PMADPictureC2;//PMADPictureC2是预测MAD

        PreviousPictureMAD=PPictureMAD[0];//PreviousPictureMAD之前图像MAD

        
        /* predict the MAD of current picture*///预测当前图像MAD
// MADc=a1*MADp+a2 MADPictureC1=a1 MADPictureC2=a2 MADp是前1张图像对应位置的MAD,预测待编码基本单元的MAD
//         因为BU宏块等于桢宏块数目，所以直接用图像的MAD值PreviousPictureMAD

        CurrentFrameMAD=MADPictureC1*PreviousPictureMAD+MADPictureC2;
        
        /*compute the number of bits for the texture*/      
        
        if(T<0)//T是图像可用位数，如果小于0，需要调大当前图像QP，减少位数

        {
          m_Qc=m_Qp+DuantQp;//调大当前图像QP

          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping保证不大于最大量化步长51

        }
        else
        {//公式(T-H)/MAD=MADPictureC1/Qstep+MADPictureC2/(Qstep*Qstep)

          m_Bits =T-m_Hp;//去除头信息

//MINVALUE=4.0;     在基本单元层 分配基本单元的位数R R= MAX(R, (int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));
          //这里TotalNumberofBasicUnit=1所以不存在基本单元层 直接如下拉


          m_Bits = MAX(m_Bits, (int)(bit_rate/(MINVALUE*frame_rate)));
          dtmp = CurrentFrameMAD * m_X1 * CurrentFrameMAD * m_X1 //\dtmp是2次方程式delta如果这都不懂 回去练初中

            + 4 * m_X2 * CurrentFrameMAD * m_Bits;
//下面是2次方程式求解 初中数学 不详细推了

          if ((m_X2 == 0.0) || (dtmp < 0) || ((sqrt (dtmp) - m_X1 * CurrentFrameMAD) <= 0.0)) // fall back 1st order mode
            m_Qstep = (float) (m_X1 * CurrentFrameMAD / (double) m_Bits);
          else // 2nd order mode
            m_Qstep = (float) ((2 * m_X2 * CurrentFrameMAD) / (sqrt (dtmp) - m_X1 * CurrentFrameMAD));
//           量化步长m_Qstep和量化索引QP关系：量化索引+6那么量化步长翻倍！量化索引+1那么量化索引增加12.5%
//           和我们的uniform均一量化有区别吧 具体为什么可能由于概率密度函数决定地量化概率分布和失真的关系吧
//           求出0-51的量化索引m_qc吧

          m_Qc=Qstep2QP(m_Qstep);
          //下面是保证+-2和1-51量化钳位 保证图像质量前后差别不大，量化参数不超过上下限

          m_Qc = MIN(m_Qp+DuantQp,  m_Qc);  // control variation
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = MAX(m_Qp-DuantQp, m_Qc); // control variation
          m_Qc = MAX(RC_MIN_QUANT, m_Qc);
        }
        
        if(img->FieldControl==0)//场信号控制关闭

        {
          /*frame coding*/
          if(active_sps->frame_mbs_only_flag)//只有宏块桢模式

          {
            img->TotalQpforPPicture +=m_Qc;//TotalQpforPPicture表示GOP内已编码P图像Qp和

            PreviousQp1=PreviousQp2;//2个P图像的QP中间参数更新

            PreviousQp2=m_Qc;
            Pm_Qp=m_Qc;
          }
          /*adaptive field/frame coding*///自适应桢场模式

          else
            FrameQPBuffer=m_Qc;
        }
        
        return m_Qc;
      }
   }
   /*bottom field*///底场不说了 写到手麻

   else
   {
     if((img->type==P_SLICE)&&(img->IFLAG==0))
     {
       /*field coding*/
       if(input->PicInterlace==FIELD_CODING)
       {
         img->TotalQpforPPicture +=m_Qc;
         PreviousQp1=PreviousQp2+1; 
         PreviousQp2=m_Qc;//+0 Recent change 13/1/2003
         Pm_Qp=m_Qc;
       }
       /*adaptive field/frame coding*/
       else
         FieldQPBuffer=m_Qc;     
     }
     return m_Qc;
   }
  }
  /*basic unit layer rate control*///基本单元速率控制 下面和前面基本相似不多少说了 
  else
  {
    /*top filed of I frame*///顶场

    if(img->type==I_SLICE)
    {
      m_Qc=MyInitialQp;//直接用GOP的QP

      return m_Qc;
    }
    /*bottom field of I frame*///是P条带场控制开启//还有待修改注释

    else if((img->type==P_SLICE)&&(img->IFLAG==1)&&(img->FieldControl==1))
    {
      m_Qc=MyInitialQp;
      return m_Qc;
    }
    else if(img->type==B_SLICE)//B条带

    {
      /*top filed of B frame*///顶场
      if((topfield)||(img->FieldControl==0))//场控制关闭

      {
        if(input->successive_Bframe==1)
        {//BFrameNumber只能为1前面已经分析过了

          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
        if(BFrameNumber==0)
          BFrameNumber=input->successive_Bframe;
        /*adaptive field/frame coding*/
        else if(BFrameNumber==1)
        {
          if((input->PicInterlace==ADAPTIVE_CODING)//\图像自适应桢场编码

            ||(input->MbInterlace))//宏块自适应或场编码

          {
            if(img->FieldControl==0)
            {       // 更新之前两个P图像的QP值
     
        /*previous choice is frame coding*/
              if(img->FieldFrame==1)//当前图像是桢编码

              {
                PreviousQp1=PreviousQp2;
                PreviousQp2=FrameQPBuffer;
              }
        /*previous choice is field coding*/
              else//当前图像是场编码

              {
                PreviousQp1=PreviousQp2;
                PreviousQp2=FieldQPBuffer;
              }
            }
          }
        }
//          m_Qc=QB1=(QP1+QP2+2)/2  当QP1!=QP2
          //          m_Qc=QB1=QP1+2          当QP1 =QP2

          if(PreviousQp1==PreviousQp2)
            m_Qc=PreviousQp1+2;
          else
            m_Qc=(PreviousQp1+PreviousQp2)/2+1;
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping钳位0-51

          m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        else//后续B桢大于1的情况

        {////BFrameNumber只能为1-successive_Bframe前面已经分析过了,代表2个P之间第几个B图像

          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
        if(BFrameNumber==0)
          BFrameNumber=input->successive_Bframe;
        /*adaptive field/frame coding*/
        else if(BFrameNumber==1)
        {
          if((input->PicInterlace==ADAPTIVE_CODING)\

            ||(input->MbInterlace))//图像自适应编码或帧编码
          {//更新之前两个P图像的QP值

            if(img->FieldControl==0)
            {
          /*previous choice is frame coding*/
              if(img->FieldFrame==1)
              {
                PreviousQp1=PreviousQp2;
                PreviousQp2=FrameQPBuffer;
              }
          /*previous choice is field coding*/
              else
              {
                PreviousQp1=PreviousQp2;
                PreviousQp2=FieldQPBuffer;
              }
            } 
          }
        }
// 参见前面写的alpha取值和为什么要这样取值

        if((PreviousQp2-PreviousQp1)<=(-2*input->successive_Bframe-3))
          StepSize=-3;
        else  if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-2))
          StepSize=-2;
        else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe-1))
          StepSize=-1;
        else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe))
          StepSize=0;//0
        else if((PreviousQp2-PreviousQp1)==(-2*input->successive_Bframe+1))
          StepSize=1;//1
        else
          StepSize=2;//2
        m_Qc=PreviousQp1+StepSize;
		//         //保证播放顺序靠近第１个P的修正更加接近第１个P的Qp,
//接近第２个P图像的修正更加接近第２个P的Qp，也是为了做个平滑处理，避免图像质量差别过大

        m_Qc +=MIN(2*(BFrameNumber-1),MAX(-2*(BFrameNumber-1), \
          (BFrameNumber-1)*(PreviousQp2-PreviousQp1)/(input->successive_Bframe-1)));
        m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping
        m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        return m_Qc;
      }
      /*bottom field of B frame*/
      else
        return m_Qc;
    }
    else if(img->type==P_SLICE)//P条带

    {
      if((img->NumberofGOP==1)&&(img->NumberofPPicture==0))//第1个GOP并且已经编码P数目为0,也就是GOP中所有P图像都还没编码时候

      {
        if((img->FieldControl==0)||((img->FieldControl==1)\
          &&(img->IFLAG==0)))
        {
        /*top field of the first P frame*///第1个P的顶场,取GOP的QP初始值,因为P未编码,
 //       索以把基本单元内头信息NumberofBasicUnitHeaderBits和纹理信息NumberofBasicUnitTextureBits都清零

          m_Qc=MyInitialQp;
          img->NumberofBasicUnitHeaderBits=0;
          img->NumberofBasicUnitTextureBits=0;
          NumberofBasicUnit--;//BU数目减1

        /*bottom field of the first P frame*///第1个P的底场

          if((!topfield)&&(NumberofBasicUnit==0))
          {
            /*frame coding or field coding*/
            if((active_sps->frame_mbs_only_flag)||(input->PicInterlace==FIELD_CODING))//只有宏块桢编码 或者图像场编码

            {
              img->TotalQpforPPicture +=m_Qc;//当前P的Qc总数和更新

              PreviousQp1=PreviousQp2;
              PreviousQp2=m_Qc;
              PAveFrameQP=m_Qc;
// 介绍下PAveHeaderBits2 PAveHeaderBits3 PAveHeaderBits1吧用于计算所有单元的头信息的位数
//               
//               Mhdrl第l个基本单元实际产生的头信息位数
//               AveMhdrl=AveMhdrl*(1-1/l)+Mhdrl/l 可见l=1时AveMhdrl=Mhdrl
//               PMhdrl是前一章图像所有基本单元预测得到地 PMhdrl=AveMhdrl/Nunit+PMhdrl*(1-1/Nunit)初始化PMhdrl=0
//               PAveHeaderBits2是表示预测BU的头信息位数也就是平均BU的头信息位数


              PAveHeaderBits3=PAveHeaderBits2;
            }
            /*adaptive frame/field coding*/
            else if((input->PicInterlace==ADAPTIVE_CODING)\
              ||(input->MbInterlace))
            {
              if(img->FieldControl==0)
              {
                FrameQPBuffer=m_Qc;
                FrameAveHeaderBits=PAveHeaderBits2;//平均BU的头信息位数

              }
              else
              {
                FieldQPBuffer=m_Qc;
                FieldAveHeaderBits=PAveHeaderBits2;
              }
            }
          }
          Pm_Qp=m_Qc;
          TotalFrameQP +=m_Qc;
          return m_Qc;
        }
      }
      else
      {
        m_X1=Pm_X1;
        m_X2=Pm_X2;
        m_Hp=PPreHeader;
        m_Qp=Pm_Qp;
        DuantQp=PDuantQp;
        MADPictureC1=PMADPictureC1;
        MADPictureC2=PMADPictureC2;

        if(img->FieldControl==0)//场控制关闭 基本单元数目不变

          SumofBasicUnit=TotalNumberofBasicUnit;
        else//场控制开启 基本单元数目减半

          SumofBasicUnit=TotalNumberofBasicUnit/2;

        /*the average QP of the previous frame is used to coded the first basic unit of the current frame or field*/
        if(NumberofBasicUnit==SumofBasicUnit)//图像内基本单元编码完毕计算之前图像的平均Qp用于第1个BU单元编码

        {

          /*adaptive field/frame coding*/
          if(((input->PicInterlace==ADAPTIVE_CODING)\
            ||(input->MbInterlace))\
            &&(img->FieldControl==0))
          {
            /*previous choice is frame coding*/
            if(img->FieldFrame==1)
            {
              if(img->NumberofPPicture>0)
                img->TotalQpforPPicture +=FrameQPBuffer;
              PAveFrameQP=FrameQPBuffer;
              PAveHeaderBits3=FrameAveHeaderBits;
            }       
            /*previous choice is field coding*/
            else
            {
              if(img->NumberofPPicture>0)
                img->TotalQpforPPicture +=FieldQPBuffer;
              PAveFrameQP=FieldQPBuffer;
              PAveHeaderBits3=FieldAveHeaderBits;
            }
          }

          if(T<=0)//Ｔ小于０　图像代表可用位数小于０

          {
            m_Qc=PAveFrameQP+2;
            if(m_Qc>RC_MAX_QUANT)
              m_Qc=RC_MAX_QUANT;
            if(topfield||(img->FieldControl==0))
              GOPOverdue=TRUE1;
          }
          else
          {
            m_Qc=PAveFrameQP; 
          }
          TotalFrameQP +=m_Qc;
          NumberofBasicUnit--;
          Pm_Qp=PAveFrameQP;
          return m_Qc;
        }else
         {
          /*compute the number of remaining bits*/
           TotalBasicUnitBits=img->NumberofBasicUnitHeaderBits+img->NumberofBasicUnitTextureBits;
           T -=TotalBasicUnitBits;
           img->NumberofBasicUnitHeaderBits=0;
           img->NumberofBasicUnitTextureBits=0;
           if(T<0)
           {
             if(GOPOverdue==TRUE1)
               m_Qc=m_Qp+2;
             else 
               m_Qc=m_Qp+DDquant;//2 
             m_Qc = MIN(m_Qc, RC_MAX_QUANT);  // clipping
             if(input->basicunit>=MBPerRow)
               m_Qc = MIN(m_Qc, PAveFrameQP+6); 
             else
               m_Qc = MIN(m_Qc, PAveFrameQP+3);

             TotalFrameQP +=m_Qc;
             NumberofBasicUnit--;
             if(NumberofBasicUnit==0)
             {
               if((!topfield)||(img->FieldControl==0))
               {
                 /*frame coding or field coding*/
                 if((active_sps->frame_mbs_only_flag)||(input->PicInterlace==FIELD_CODING))
                 {
                   PAverageQP=(int)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                   if (img->NumberofPPicture == (input->intra_period - 2))
                          QPLastPFrame = PAverageQP;

                   img->TotalQpforPPicture +=PAverageQP;
                   if(GOPOverdue==TRUE1)
                   {
                     PreviousQp1=PreviousQp2+1;
                     PreviousQp2=PAverageQP;                   
                   }
                   else
                   {
                     if((img->NumberofPPicture==0)&&(img->NumberofGOP>1))
                     {
                       PreviousQp1=PreviousQp2;
                       PreviousQp2=PAverageQP;
                     }
                     else if(img->NumberofPPicture>0)
                     {
                        PreviousQp1=PreviousQp2+1;
                        PreviousQp2=PAverageQP;
                     }
                   }
                   PAveFrameQP=PAverageQP;
                   PAveHeaderBits3=PAveHeaderBits2;
                 }
                 /*adaptive field/frame coding*/
                 else if((input->PicInterlace==ADAPTIVE_CODING)\
                   ||(input->MbInterlace))
                 {
                   if(img->FieldControl==0)
                   {
                     PAverageQP=(int)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                     FrameQPBuffer=PAverageQP;
                     FrameAveHeaderBits=PAveHeaderBits2;
                   }
                   else
                   {
                     PAverageQP=(int)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                     FieldQPBuffer=PAverageQP;
                     FieldAveHeaderBits=PAveHeaderBits2;
                   }
                 }
               }
             }
             if(GOPOverdue==TRUE1)
               Pm_Qp=PAveFrameQP;
             else
               Pm_Qp=m_Qc;
             return m_Qc;
           }
           else
           {//计算第l个BU的目标位数Bi=Tr*MADl*MADl/(sum(k=l-Nunit)MADk)其中sum(k=l-Nunit)MADk是已编码BU的MAD和

             /*predict the MAD of current picture*/
             if(((input->PicInterlace==ADAPTIVE_CODING)||(input->MbInterlace))\
               &&(img->FieldControl==1))
             {//预测当前MAD通过之前FCBUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]；FCBUPFMAD先对第TotalNumberofBasicUnit个到第1个存储地

               CurrentFrameMAD=MADPictureC1*FCBUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]+MADPictureC2;
               TotalBUMAD=0;
			   ////TotalNumberofBasicUnit-NumberofBasicUni表示已编码的BU数目

               for(i=TotalNumberofBasicUnit-1; i>=(TotalNumberofBasicUnit-NumberofBasicUnit);i--)
               {//根据前面公式更新当前基本单元MADc=a1*MADp+a2

                 CurrentBUMAD=MADPictureC1*FCBUPFMAD[i]+MADPictureC2;
                 TotalBUMAD +=CurrentBUMAD*CurrentBUMAD;
               }
             }
             else
             {
               CurrentFrameMAD=MADPictureC1*BUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]+MADPictureC2;
               TotalBUMAD=0;
               for(i=TotalNumberofBasicUnit-1; i>=(TotalNumberofBasicUnit-NumberofBasicUnit);i--)
               {
                 CurrentBUMAD=MADPictureC1*BUPFMAD[i]+MADPictureC2;
                 TotalBUMAD +=CurrentBUMAD*CurrentBUMAD;
               }
             }
                
             /*compute the total number of bits for the current basic unit*/
             m_Bits =(int)(T*CurrentFrameMAD*CurrentFrameMAD/TotalBUMAD);
             /*compute the number of texture bits*///减去头信息位数得到纹理信息可用位数
// Mhdrl第l个基本单元实际产生的头信息位数
//               AveMhdrl=AveMhdrl*(1-1/l)+Mhdrl/l 可见l=1时AveMhdrl=Mhdrl
//               PMhdrl是前一章图像所有基本单元预测得到地 PMhdrl=AveMhdrl/Nunit+PMhdrl*(1-1/Nunit)初始化PMhdrl=0
//               这里 PAveHeaderBits2=PMhdrl也就是表示预测的头信息位数    
             m_Bits -=PAveHeaderBits2;
             // 在基本单元层 分配基本单元的位数R R= MAX(R, (int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));

             m_Bits=MAX(m_Bits,(int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));
        //m_qc求法上面说过 不提了

             dtmp = CurrentFrameMAD * m_X1 * CurrentFrameMAD * m_X1 \
               + 4 * m_X2 * CurrentFrameMAD * m_Bits;
             if ((m_X2 == 0.0) || (dtmp < 0) || ((sqrt (dtmp) - m_X1 * CurrentFrameMAD) <= 0.0))  // fall back 1st order mode
               m_Qstep = (float)(m_X1 * CurrentFrameMAD / (double) m_Bits);
             else // 2nd order mode
               m_Qstep = (float) ((2 * m_X2 * CurrentFrameMAD) / (sqrt (dtmp) - m_X1 * CurrentFrameMAD));

             m_Qc=Qstep2QP(m_Qstep);
             m_Qc = MIN(m_Qp+DDquant,  m_Qc); // control variation

             if(input->basicunit>=MBPerRow)
               m_Qc = MIN(PAveFrameQP+6, m_Qc);
             else
               m_Qc = MIN(PAveFrameQP+3, m_Qc);

             m_Qc = MIN(m_Qc, RC_MAX_QUANT);  // clipping
             m_Qc = MAX(m_Qp-DDquant, m_Qc);  // control variation 
             if(input->basicunit>=MBPerRow)//前面说过基本单元内宏块数目大于１行中宏块数目，修正QP是为了平滑图像质量

               m_Qc = MAX(PAveFrameQP-6, m_Qc);
             else
               m_Qc = MAX(PAveFrameQP-3, m_Qc);

             m_Qc = MAX(RC_MIN_QUANT, m_Qc);
             TotalFrameQP +=m_Qc;
             Pm_Qp=m_Qc;
             NumberofBasicUnit--;
             if((NumberofBasicUnit==0)&&(img->type==P_SLICE))//NumberofBasicUnit表示最后1个BU

             {
               if((!topfield)||(img->FieldControl==0))
               {
                 /*frame coding or field coding*/
                 if((active_sps->frame_mbs_only_flag)||(input->PicInterlace==FIELD_CODING))
                 {
                   PAverageQP=(int)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                   if (img->NumberofPPicture == (input->intra_period - 2))
                          QPLastPFrame = PAverageQP;

                   img->TotalQpforPPicture +=PAverageQP;
                   PreviousQp1=PreviousQp2;
                   PreviousQp2=PAverageQP; 
                   PAveFrameQP=PAverageQP;
                   PAveHeaderBits3=PAveHeaderBits2;
                 }
                 else if((input->PicInterlace==ADAPTIVE_CODING)\
                   ||(input->MbInterlace))
                 {
                   if(img->FieldControl==0)
                   {
                     PAverageQP=(int)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                     FrameQPBuffer=PAverageQP;
                     FrameAveHeaderBits=PAveHeaderBits2;
                   }
                   else
                   {
                     PAverageQP=(int)(1.0*TotalFrameQP/TotalNumberofBasicUnit+0.5);
                     FieldQPBuffer=PAverageQP;
                     FieldAveHeaderBits=PAveHeaderBits2;
                   }
                 }
               }
             }
             return m_Qc;
           }
         }
      }
    } 
  }
  return m_Qc;
}

//update the parameters of quadratic R-D model
void updateRCModel ()
{

  int n_windowSize;
  int i;
  double error[20], std = 0.0, threshold;
  int m_Nc;
  Boolean MADModelFlag = FALSE1;
   
  if(img->type==P_SLICE)
  {
    /*frame layer rate control*/
    if(img->BasicUnit==img->Frame_Total_Number_MB)
    {
      CurrentFrameMAD=ComputeFrameMAD();
      m_Nc=img->NumberofCodedPFrame;
    }
    /*basic unit layer rate control*/
    else
    {
      /*compute the MAD of the current basic unit*/
      if((input->MbInterlace)&&(img->FieldControl==0))
        CurrentFrameMAD=img->TotalMADBasicUnit/img->BasicUnit/2;
      else
        CurrentFrameMAD=img->TotalMADBasicUnit/img->BasicUnit;

        
      img->TotalMADBasicUnit=0;
              
      /* compute the average number of header bits*/
      
        CodedBasicUnit=TotalNumberofBasicUnit-NumberofBasicUnit;
        if(CodedBasicUnit>0)
        {
          PAveHeaderBits1=(int)(1.0*(PAveHeaderBits1*(CodedBasicUnit-1)+\
            +img->NumberofBasicUnitHeaderBits)/CodedBasicUnit+0.5);
          if(PAveHeaderBits3==0)
            PAveHeaderBits2=PAveHeaderBits1;
          else
            PAveHeaderBits2=(int)(1.0*(PAveHeaderBits1*CodedBasicUnit+\
            +PAveHeaderBits3*NumberofBasicUnit)/TotalNumberofBasicUnit+0.5);
        }
          /*update the record of MADs for reference*/
          if(((input->PicInterlace==ADAPTIVE_CODING)||(input->MbInterlace))\
          &&(img->FieldControl==1))
            FCBUCFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit]=CurrentFrameMAD;
          else
            BUCFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit]=CurrentFrameMAD;

        if(NumberofBasicUnit!=0)
          m_Nc=img->NumberofCodedPFrame*TotalNumberofBasicUnit+CodedBasicUnit;
        else
          m_Nc=(img->NumberofCodedPFrame-1)*TotalNumberofBasicUnit+CodedBasicUnit;
      
    }
    
    
   
    if(m_Nc>1)
      MADModelFlag=TRUE1;
    PPreHeader=img->NumberofHeaderBits;
    for (i = 19; i > 0; i--) {// update the history
      Pm_rgQp[i] = Pm_rgQp[i - 1];
      m_rgQp[i]=Pm_rgQp[i];
        Pm_rgRp[i] = Pm_rgRp[i - 1];
      m_rgRp[i]=Pm_rgRp[i];
    }
    Pm_rgQp[0] = QP2Qstep(m_Qc); //*1.0/CurrentFrameMAD;
    /*frame layer rate control*/
    if(img->BasicUnit==img->Frame_Total_Number_MB)
      Pm_rgRp[0] = img->NumberofTextureBits*1.0/CurrentFrameMAD;
    /*basic unit layer rate control*/
    else
      Pm_rgRp[0]=img->NumberofBasicUnitTextureBits*1.0/CurrentFrameMAD;

    m_rgQp[0]=Pm_rgQp[0];
    m_rgRp[0]=Pm_rgRp[0];
    m_X1=Pm_X1;
      m_X2=Pm_X2;
  

/*compute the size of window*/


    n_windowSize = (CurrentFrameMAD>PreviousFrameMAD)?(int)(PreviousFrameMAD/CurrentFrameMAD*20)\
      :(int)(CurrentFrameMAD/PreviousFrameMAD*20);
    n_windowSize=MAX(n_windowSize, 1);
    n_windowSize=MIN(n_windowSize,m_Nc);
    n_windowSize=MIN(n_windowSize,m_windowSize+1);
    n_windowSize=MIN(n_windowSize,20);

      /*update the previous window size*/
  m_windowSize=n_windowSize;
  


  for (i = 0; i < 20; i++) {
    m_rgRejected[i] = FALSE1;
  }

  // initial RD model estimator
  RCModelEstimator (n_windowSize);
 
  n_windowSize = m_windowSize;
  // remove outlier 
  
  for (i = 0; i < (int) n_windowSize; i++) {
    error[i] = m_X1 / m_rgQp[i] + m_X2 / (m_rgQp[i] * m_rgQp[i]) - m_rgRp[i];
    std += error[i] * error[i]; 
  }
  threshold = (n_windowSize == 2) ? 0 : sqrt (std / n_windowSize);
  for (i = 0; i < (int) n_windowSize; i++) {
    if (fabs(error[i]) > threshold)
      m_rgRejected[i] = TRUE1;
  }
    // always include the last data point
  m_rgRejected[0] = FALSE1;

  // second RD model estimator
  RCModelEstimator (n_windowSize);

  if(MADModelFlag)
    updateMADModel();
  else if(img->type==P_SLICE)
    PPictureMAD[0]=CurrentFrameMAD;
  } 
}

/*Boolean skipThisFrame ()
{
  if (m_B > (int) ((RC_SKIP_MARGIN / 100.0) * m_Bs)) {  // buffer full!
    m_skipNextFrame = TRUE;           // set the status
    m_Nr--;                   // skip one frame
    m_B -= m_Rp;                // decrease current buffer level
  } else
    m_skipNextFrame = FALSE;
  return m_skipNextFrame;
}
*/

void RCModelEstimator (int n_windowSize)
{
  int n_realSize = n_windowSize;
  int i;
  double oneSampleQ;
  double a00 = 0.0, a01 = 0.0, a10 = 0.0, a11 = 0.0, b0 = 0.0, b1 = 0.0;
  double MatrixValue;
  Boolean estimateX2 = FALSE1;

  for (i = 0; i < n_windowSize; i++) {// find the number of samples which are not rejected
    if (m_rgRejected[i])
      n_realSize--;
  }

  // default RD model estimation results

  m_X1 = m_X2 = 0.0;

  for (i = 0; i < n_windowSize; i++)  {
    if (!m_rgRejected[i])
      oneSampleQ = m_rgQp[i];
  }
  for (i = 0; i < n_windowSize; i++)  {// if all non-rejected Q are the same, take 1st order model
    if ((m_rgQp[i] != oneSampleQ) && !m_rgRejected[i])
      estimateX2 = TRUE1;
    if (!m_rgRejected[i])
      m_X1 += (m_rgQp[i] * m_rgRp[i]) / n_realSize;
  }

  // take 2nd order model to estimate X1 and X2
  if ((n_realSize >= 1) && estimateX2) {
      for (i = 0; i < n_windowSize; i++) {
      if (!m_rgRejected[i]) {
        a00 = a00 + 1.0;
        a01 += 1.0 / m_rgQp[i];
        a10 = a01;
        a11 += 1.0 / (m_rgQp[i] * m_rgQp[i]);
        b0 += m_rgQp[i] * m_rgRp[i];
        b1 += m_rgRp[i];
      }
    }
    // solve the equation of AX = B
      MatrixValue=a00*a11-a01*a10;
      if(fabs(MatrixValue)>0.000001)
      {
        m_X1=(b0*a11-b1*a01)/MatrixValue;
        m_X2=(b1*a00-b0*a10)/MatrixValue;
      }
      else
      {
        m_X1=b0/a00;
        m_X2=0.0;
      }
  
  }
  if(img->type==P_SLICE)
  {
    Pm_X1=m_X1;
    Pm_X2=m_X2;
  }
}

double ComputeFrameMAD()
{
  double TotalMAD;
  int i;
  TotalMAD=0.0;
//  CurrentFrameMAD=0.0;
  for(i=0;i<img->Frame_Total_Number_MB;i++)
    TotalMAD +=img->MADofMB[i];
  TotalMAD /=img->Frame_Total_Number_MB;
  return TotalMAD;
}


//update the parameters of linear prediction model
void updateMADModel ()
{
  
  int n_windowSize;
  int i;
  double error[20], std = 0.0, threshold;
  int m_Nc;
  
  if(img->NumberofCodedPFrame>0)
  {
    
    if(img->type==P_SLICE)
    {
      /*frame layer rate control*/
      if(img->BasicUnit==img->Frame_Total_Number_MB)
        m_Nc=img->NumberofCodedPFrame;
      /*basic unit layer rate control*/
      else
        m_Nc=img->NumberofCodedPFrame*TotalNumberofBasicUnit+CodedBasicUnit;
      for (i = 19; i > 0; i--) {// update the history
        PPictureMAD[i] = PPictureMAD[i - 1];
        PictureMAD[i]=PPictureMAD[i];
        ReferenceMAD[i]= ReferenceMAD[i-1];
      }
      PPictureMAD[0] = CurrentFrameMAD;
      PictureMAD[0]=PPictureMAD[0];
      if(img->BasicUnit==img->Frame_Total_Number_MB)
        ReferenceMAD[0]=PictureMAD[1];
      else
      {
        if(((input->PicInterlace==ADAPTIVE_CODING)||(input->MbInterlace))\
          &&(img->FieldControl==1))
          ReferenceMAD[0]=FCBUPFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit];
        else
          ReferenceMAD[0]=BUPFMAD[TotalNumberofBasicUnit-1-NumberofBasicUnit];
      }
      MADPictureC1=PMADPictureC1;
      MADPictureC2=PMADPictureC2;
    }
    
    
    /*compute the size of window*/
    
    n_windowSize = (CurrentFrameMAD>PreviousFrameMAD)?(int)(PreviousFrameMAD/CurrentFrameMAD*20)\
      :(int)(CurrentFrameMAD/PreviousFrameMAD*20);
    n_windowSize=MIN(n_windowSize,(m_Nc-1));
    n_windowSize=MAX(n_windowSize, 1);
    n_windowSize=MIN(n_windowSize,MADm_windowSize+1);
    n_windowSize=MIN(20,n_windowSize);
    /*update the previous window size*/
    MADm_windowSize=n_windowSize;
    
    for (i = 0; i < 20; i++) {
      PictureRejected[i] = FALSE1;
    }
    //update the MAD for the previous frame
    if(img->type==P_SLICE)
      PreviousFrameMAD=CurrentFrameMAD;
    
    // initial MAD model estimator
    MADModelEstimator (n_windowSize);
    
    // remove outlier 
    
    for (i = 0; i < (int) n_windowSize; i++) {
      error[i] = MADPictureC1*ReferenceMAD[i]+MADPictureC2-PictureMAD[i];
      std += error[i] * error[i]; 
    }
    threshold = (n_windowSize == 2) ? 0 : sqrt (std / n_windowSize);
    for (i = 0; i < (int) n_windowSize; i++) {
      if (fabs(error[i]) > threshold)
        PictureRejected[i] = TRUE1;
    }
    // always include the last data point
    PictureRejected[0] = FALSE1;
    
    // second MAD model estimator
    MADModelEstimator (n_windowSize);
  }
}

void MADModelEstimator (int n_windowSize)
{
  int n_realSize = n_windowSize;
  int i;
  double oneSampleQ;
  double a00 = 0.0, a01 = 0.0, a10 = 0.0, a11 = 0.0, b0 = 0.0, b1 = 0.0;
  double MatrixValue;
  Boolean estimateX2 = FALSE1;
  
  for (i = 0; i < n_windowSize; i++) {// find the number of samples which are not rejected
    if (PictureRejected[i])
      n_realSize--;
  }
  
  // default MAD model estimation results
  
  MADPictureC1 = MADPictureC2 = 0.0;
  
  for (i = 0; i < n_windowSize; i++)  {
    if (!PictureRejected[i])
      oneSampleQ = PictureMAD[i];
  }
  for (i = 0; i < n_windowSize; i++)  {// if all non-rejected MAD are the same, take 1st order model
    if ((PictureMAD[i] != oneSampleQ) && !PictureRejected[i])
      estimateX2 = TRUE1;
    if (!PictureRejected[i])
      MADPictureC1 += PictureMAD[i] / (ReferenceMAD[i]*n_realSize);
  }
  
  // take 2nd order model to estimate X1 and X2
  if ((n_realSize >= 1) && estimateX2) {
    for (i = 0; i < n_windowSize; i++) {
      if (!PictureRejected[i]) {
        a00 = a00 + 1.0;
        a01 += ReferenceMAD[i];
        a10 = a01;
        a11 += ReferenceMAD[i]*ReferenceMAD[i];
        b0 += PictureMAD[i];
        b1 += PictureMAD[i]*ReferenceMAD[i];
      }
    }
    // solve the equation of AX = B
    MatrixValue=a00*a11-a01*a10;
    if(fabs(MatrixValue)>0.000001)
    {
      MADPictureC2=(b0*a11-b1*a01)/MatrixValue;
      MADPictureC1=(b1*a00-b0*a10)/MatrixValue;
    }
    else
    {
      MADPictureC1=b0/a01;
      MADPictureC2=0.0;
    }
    
  }
  if(img->type==P_SLICE)
  {
    PMADPictureC1=MADPictureC1;
    PMADPictureC2=MADPictureC2;
  }
}


double QP2Qstep( int QP )
{
  int i; 
  double Qstep;
  static const double QP2QSTEP[6] = { 0.625, 0.6875, 0.8125, 0.875, 1.0, 1.125 };
  
  Qstep = QP2QSTEP[QP % 6];
  for( i=0; i<(QP/6); i++)
    Qstep *= 2;
  
  return Qstep;
}

int Qstep2QP( double Qstep )
{
  int q_per = 0, q_rem = 0;
  
  //  assert( Qstep >= QP2Qstep(0) && Qstep <= QP2Qstep(51) );
  if( Qstep < QP2Qstep(0))
    return 0;
  else if (Qstep > QP2Qstep(51) )
    return 51;
  
  while( Qstep > QP2Qstep(5) )
  {
    Qstep /= 2;
    q_per += 1;
  }
  
  if (Qstep <= (0.625+0.6875)/2) 
  {
    Qstep = 0.625;
    q_rem = 0;
  }
  else if (Qstep <= (0.6875+0.8125)/2)
  {
    Qstep = 0.6875;
    q_rem = 1;
  }
  else if (Qstep <= (0.8125+0.875)/2)
  {
    Qstep = 0.8125;
    q_rem = 2;
  }
  else if (Qstep <= (0.875+1.0)/2)
  {
    Qstep = 0.875;
    q_rem = 3;
  }
  else if (Qstep <= (1.0+1.125)/2)
  {
    Qstep = 1.0;  
    q_rem = 4;
  }
  else 
  {
    Qstep = 1.125;
    q_rem = 5;
  }
  
  return (q_per * 6 + q_rem);
}
