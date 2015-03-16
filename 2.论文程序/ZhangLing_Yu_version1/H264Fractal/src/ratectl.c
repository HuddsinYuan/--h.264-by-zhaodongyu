
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


const double THETA=1.3636;//yynע����֡�����ʿ���ǰ����׶μ���P֡���Ӷ�ʱ�õĳ�����
const int Switch=0;

int Iprev_bits=0;
int Pprev_bits=0;


/* rate control variables */
int Xp, Xb;// �ֱ�Ϊ P ͼ��� B ͼ���ͼ���Ӷȣ�����Ϊͼ������������ͼ�����л�����Ԫƽ����������֮��������������
static int R,T_field;
static int Np, Nb, bits_topfield, Q;
long T,T1;
//HRD consideration
long UpperBound1, UpperBound2, LowerBound;
double InitialDelayOffset;
const double OMEGA=0.9;//yynע��HRD�����޼���ʽ�еĳ��� ��JVT-H017ʽ(4)

double Wp,Wb;  // �ֱ�Ϊ P ͼ��� B ͼ���ƽ��ͼ���Ӷȣ�����Ϊͼ������������ͼ�����л�����Ԫƽ����������֮����
int TotalPFrame;
int DuantQp; 
int PDuantQp;
FILE *BitRate;
double DeltaP;//yynע��������Ԫ�����ʿ���������������ʱʹ�õġ���

// Initiate rate control parameters��ʼ�����ʿ��Ʋ���
void rc_init_seq()
{
   double L1,L2,L3,bpp;//yynע:GOP�����ʿ���ʱ��ʼ��QP0ʱʹ�ã�bppΪÿͼ�������//
   int qp;
   int i;
  
   Xp=0;//���ڸ��Ӷ�Wp��Wb
   Xb=0;
   
   bit_rate=input->bit_rate;//��������� //yynע��bit_rate�ŵ�����
   frame_rate = (float)(img->framerate *(input->successive_Bframe + 1)) / (float) (input->jumpd + 1);///ʵ��������
                                                                        //successive_Bframe�Ǻ��B����Ŀ.jumpd��ʾ�ǲο�֡���Զ�������Ŀ.
   PreviousBit_Rate=bit_rate;
   
   /*compute the total number of MBs in a frame*///����1������

   img->Frame_Total_Number_MB=img->height*img->width/256;//yynע��ͼ���е�ǰ֡�ĺ������=ͼ��ĸ�xͼ��Ŀ�/256
   if(input->basicunit>img->Frame_Total_Number_MB)//yynע����������ļ��Ļ�����Ԫ��>֡�к�����������������Ԫ��=֡�к������
	   input->basicunit=img->Frame_Total_Number_MB;//������Ԫ�����Ŀ�������������Ŀ ������Ԫ���ʿ��� ��ͼ������ʿ���
   if(input->basicunit<img->Frame_Total_Number_MB) //yynע����������ļ��Ļ�����Ԫ��<֡�к�����������ܻ�����Ԫ��=֡�к��
     TotalNumberofBasicUnit=img->Frame_Total_Number_MB/input->basicunit;//����/�����ļ��Ļ�����Ԫ�������������ӦG012�Ĺ�ʽ
   
   MINVALUE=4.0;//�ڻ�����Ԫ�� ���������Ԫ��λ��R R= MAX(R, (int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));
               // TotalNumberofBasicUnit��1���BU���� �����1��P�����������ĿС������bit_rate/(MINVALUE*frame_rate���Ƕ���

   /*initialize the parameters of fluid flow traffic model��������ģ�͵ĳ�ʼ������*/
   
   BufferSize=bit_rate*2.56;//�趨buffer�Ĵ�С //yynע������������=������*2.56
   CurrentBufferFullness=0;//��ǰbuffer����0 //yynע����ǰ��������ӯ��=0����ʱ��ʼ��
   GOPTargetBufferLevel=CurrentBufferFullness; //yynע��GOPĿ�껺�弶=��ǰ��������ӯ�ȡ�
   /*HRD consideration*/
   InitialDelayOffset=BufferSize*0.8; //yynע����ʼ���ӳ�ƫ��=��������*80%//����UpperBound1=(long)(R+InitialDelayOffset);�涨��hrd����ο�������������
   
   /*initialize the previous window size*/ //**yynע�����ڳߴ�������֡�����ʿ��ƺ����׶�ʱ����MAD��ϵ����R-D������//��ʼ��֮ǰ���ߴ�  ��Ϊ��δ���� ���Զ�����Ϊ0
   m_windowSize=0;
   MADm_windowSize=0;
   img->NumberofCodedBFrame=0;///�ѱ����BΪ0
   img->NumberofCodedPFrame=0;
   img->NumberofGOP=0;//gop���Ϊ0
   /*remaining # of bits in GOP *///R��ʾgop��ʣ���λ��
   R = 0;//��ȻR<0˵������С���˱���λ�� ��ʼΪ0 ��Ϊ��û����
   /*control parameter *///��ͼ�����Ʋ�����i��gop�еĵ�j��ͼ���Ŀ��λ�������Ŀ��buffer�ȼ� ������ ʵ��bufferռ����
   //T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5); 
   //��Bͼ����� gamma=0.25  ��Bͼ�����  gamma=0.5
   //beta��ƽ��ϸ�� T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5)); T ͼ��Ŀ��λ��
   //��B���ʱ��beta=.5 ��B���ʱ��beta=.9

   if(input->successive_Bframe>0) //yynע���������B֡>0������B֡ʱ
   {
     GAMMAP=0.25;//**yynע����֡�����ʿ���΢�ۿ���ʱ��Ϊ�˼�Сʵ�ʻ����ӯ����Ŀ�껺��ȼ�֮��Ĳ�ֵ����ʹ�õ�һ������ϵ����û��B֡ʱΪ0.75����B֡ʱΪ0.25��JVT-G012�ĵ���ʽ��16��
     BETAP=0.9;//yynע������Ŀ�������ʱ�ļ�Ȩģ�͵�ϵ������JVT-G012�ĵ���ʽ��19��
   }
   else
   {
     GAMMAP=0.5;
     BETAP=0.5;
   }
   
   /*quadratic rate-distortion model*///2��r-dģ�Ͳ�����ʼֵ,Ti=X1xMAD1/Qi+X2xMAD2/Qi2    Ti�������뵱ǰ��i֡ͼ��ı�������Qi��ǰ��i֡����������MADi��ǰ��i֡ͼ���MADֵ��
   PPreHeader=0;
   
   Pm_X1=bit_rate*1.0; //yynע��������ʧ��ģ��ϵ��//�������ȡ��ʼֵ�� �㲻��� ��˴����ȡ���� ��������������MADֵ�˲� �ѵ���������С���������1
   Pm_X2=0.0;
   /* linear prediction model for P picture*///����madԤ�������ʼֵ
   PMADPictureC1=1.0;
   PMADPictureC2=0.0;
   //���ﴰ�ڳߴ�20��Ϊ��ƽ��һ��ͻȻ�ĳ����ı䣬���Բ�����ǰ����20��ͼ�����ο�
  // �������ߴ���20һ����MADC*20/MAD_prev����
   for(i=0;i<20;i++)
   {
     Pm_rgQp[i]=0;
     Pm_rgRp[i]=0.0;
     PPictureMAD[i]=0.0;
   }
   PPictureMAD[20]=0.0;

   //yynע��Pm_rgQp��Pm_rgRp��Pm_X1��Pm_X2 �ĺ����ǣ�
   //double Pm_rgQp[20];    //++ ����ֵ���ݹ����е��м���ʱ��������ֱ���� m_rgQp �滻
   //double Pm_rgRp[20];    //++ ����ֵ���ݹ����е��м���ʱ��������ֱ���� m_rgRp �滻
   //double Pm_X1;         //++ ����ֵ���ݹ����е��м���ʱ��������ֱ���� m_X1 �滻
   //double Pm_X2;         //++ ����ֵ���ݹ����е��м���ʱ��������ֱ���� m_X2 �滻
   
   //�� m_rgQp��m_rgRp��m_X1��m_X2 �ĺ����ǣ�
   //double  m_rgQp[21];     //++ FIFO ���������洢����������Ԫ����������
   //double m_rgRp[21];      //++ FIFO ���������洢����������Ԫ������ɺ�Ķ��η��������
   // double m_X1;          //++ ����ģ�͵�һ��ϵ��
   // double m_X2;          //++ ����ģ�͵ڶ���ϵ��

   //Define the largest variation of quantization parameters�õ�DDquant���м����pDuantQp
   PDuantQp=2;//���ƿ����ڶ�ͼ��Qpc=min{Qpp+2,max(Qpp-2,Qpc)}����clip����֤ͼ���Qpc��֮ǰͼ��Qpp+-2��Χ��
   
   /*basic unit layer rate control*/
   //������Ԫ�õ��ز��� ���������Ԫ��Ŀ����8 DDquant=2 С��ȡ1����Ϊ������ԪԽ�࣬ϣ����������Խ�� ����������ؽ���  
   //PAveHeaderBitsͷ���ı���������Ϊr-dģ����ȥ��ͷ����Ϣ��ģ��
    PAveHeaderBits1=0;//yynע��PAveHeaderBits1��ʾ��ǰ֡ͼ����ǰһ��basic unitΪֹ��ƽ��ͷ��Ϣ��������
    PAveHeaderBits3=0;      //yynע��PAveHeaderBits3��ʾǰһ֡ͼ���ƽ��ͷ��Ϣ��������
//yynע��PAveHeaderBits2��ʾ��ǰbasic unit��ƽ��ͷ��Ϣ�������Ĺ���ֵ��

    if(TotalNumberofBasicUnit>=9)
      DDquant=1;
    else
      DDquant=2;
	//���1�еĺ����
	MBPerRow=input->imagewidth/16;//Ϊ�˵�2 3���������1���ǻ����1��������Ԫȡǰһͼ�����л�����Ԫƽ��ֵ
	//��2 ����R<0ʱ�� ˵�������λ��û��ȫ������ ������������Ҫ����ǰ1������Ԫ���������Լ���λ��
		//������Ԫ�����ĿС��MBPerRow��m_Qc = MIN(m_Qc, PAveFrameQP+3);����m_Qc = MIN(m_Qc, PAveFrameQP+6);
	//Ϊ�������.���ͷ���ЧӦ.ƽ��ͼ��.��QP����
    
    /*adaptive field/frame coding*/
	img->FieldControl=0;//ȡ0����������Ч ��ֻ�Ў�����.qp��Χ��0-51��
    
    RC_MAX_QUANT = 51;  // clipping
    RC_MIN_QUANT = 0;//clipping
    
    /*compute thei initial QP*/
	//�ڵ�1��gopԤ�ȶ�����������Qp ��I�͵�1��P������GOPԤ�ȶ����QP
    bpp = 1.0*bit_rate /(frame_rate*img->width*img->height);
	if (img->width == 176) //�����Qcif��ʽ����L1 L2 L3ȡֵ����
    {
      L1 = 0.1;
      L2 = 0.3;
      L3 = 0.6;
	}else if (img->width == 352)//�����cif��ʽ����L1 L2 L3ȡֵ����
    {
      L1 = 0.2;
      L2 = 0.6;
      L3 = 1.2;
	}else //������ʽȡֵ
	{
      L1 = 0.6;
      L2 = 1.4;
      L3 = 2.4;
    }
	//����bpp���������Ӧ�ĳ�ʼQp
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

// Initiate one GOPĿ�ľ������GOP��QP��ʼֵ
void rc_init_GOP(int np, int nb)//����δ�����p��Ŀnp δ�����b��Ŀnb
{
  Boolean Overum=FALSE1;
  int OverBits;
  int OverDuantQp;
  int AllocatedBits;
  int GOPDquant;  //yynע��GOP��������

  /*check if the last GOP over uses its budget. If yes, the initial QP of the I frame in 
 the coming  GOP will be increased.*/

  if(R<0)//��ȻR<0˵������С���˱���λ��
    Overum=TRUE1;
  OverBits=-R;

  /*initialize the lower bound and the upper bound for the target bits of each frame, HRD consideration*/
  //���ǵ�hrd������ο����룩������R�п���λ�������޺�����
  LowerBound=(long)(R+bit_rate/frame_rate);
  UpperBound1=(long)(R+InitialDelayOffset);

  /*compute the total number of bits for the current GOP*/ //���㵱ǰgop�Ŀ���λ
  AllocatedBits = (int) floor((1 + np + nb) * bit_rate / frame_rate + 0.5);//np nbδ�����p b�� ����ǰ����������λ��
  R +=AllocatedBits;//R��ʾ��ǰgop����λ�� ����gop������ ������� 
  //����ǰR�����Ϊ��ǰR<0���ܴ��ڣ�˵����1��gop����δ���������� ���뽵�͵�ǰgopĿ��λ�� 
  Np  = np;
  Nb  = nb;

  OverDuantQp=(int)(8*OverBits/AllocatedBits+0.5);//ԭ�����R<0δ�������� �������OverDuantQp 8*Tr(i-1,Ngop)/Tr(i,0)�ϴ�ʣ���λ�����Ե�ǰ�����λ��
  GOPOverdue=FALSE1;//��������gop�Ƿ�ٵ������� ������������ԭ����Щgop������ �������յ�
  
/*field coding*/
  img->IFLAG=1;//��

/*Compute InitialQp for each GOP*/
  TotalPFrame=np;//P����Ŀ����δ�����np��Ŀ����Ϊ���ô˺��� gop������ͼ��δ����
  img->NumberofGOP++;//gop��Ŀ��1
  if(img->NumberofGOP==1)//��1��GOPȡ��bpp��ʼ�趨��QPֵ
  {
    MyInitialQp=input->SeinitialQP;
    PreviousQp2=MyInitialQp-1; //recent change -0;
    QPLastGOP=MyInitialQp;
  
  }
  else
  {//����GOP��Qp���㹫ʽ����QP=SumPqp/Np-8*Tr(i-1,Ngop)/Tr(i,0)-min{2,Ngop/15}
	//  =img->TotalQpforPPicture/img->NumberofPPicture-OverDuantQp-min{2,Ngop/15}

// 	  /*adaptive field/frame coding*///����Ӧ�峡����
    if((input->PicInterlace==ADAPTIVE_CODING)\
		||(input->MbInterlace))//�����ͼ������Ӧ�峡���룬���ߺ���ǳ�ģʽ���ߺ�����峡����Ӧģʽ
    {
		if (img->FieldFrame == 1)//����ǎ�����  TotalQpforPPicture�������FrameQPBuffer TotalQpforPPicture��ʾ����Pͼ���QP��
      {
        img->TotalQpforPPicture += FrameQPBuffer;
        QPLastPFrame = FrameQPBuffer;
      }
		else//������
      {
        img->TotalQpforPPicture += FieldQPBuffer;
        QPLastPFrame = FieldQPBuffer;
      }
      
    }
	/*compute the average QP of P frames in the previous GOP*///����֮ǰGOP����Pͼ���е�QPƽ��ֵ NumberofPPicture��ʾPͼ����Ŀ
    PAverageQp=(int)(1.0*img->TotalQpforPPicture/img->NumberofPPicture+0.5);

    GOPDquant=(int)(0.5+1.0*(np+nb+1)/15);
    if(GOPDquant>2)
        GOPDquant=2;

    PAverageQp-=GOPDquant;

	if (PAverageQp > (QPLastPFrame - 2))//Ϊ��ͼ������Կ��� ���ƽ����QP����֮ǰPͼ���QP-2������PAverageQp
      PAverageQp--;
	PAverageQp = MAX(QPLastGOP-2,  PAverageQp);//�޶���Χ���ϸ�GOP��QPֵ+-2��Χ��
    PAverageQp = MIN(QPLastGOP+2, PAverageQp);
	PAverageQp = MIN(RC_MAX_QUANT, PAverageQp);//�޶���Χ��0-51
    PAverageQp = MAX(RC_MIN_QUANT, PAverageQp);
  

    MyInitialQp=PAverageQp;
	QPLastGOP = MyInitialQp;//������һ��gop��QP
    Pm_Qp=PAverageQp;
    PAveFrameQP=PAverageQp;
    PreviousQp1=PreviousQp2;
    PreviousQp2=MyInitialQp-1;  
  }

  img->TotalQpforPPicture=0;//��TotalQpforPPictureΪ��ǰGOP����QP�����ۼ���
  img->NumberofPPicture=0;//��ǰGOP��û���� �����Pͼ����Ŀ NumberofPPictureΪ0
  NumberofBFrames=0; //��ǰGOP��û���� �����Bͼ����Ŀ NumberofPPictureΪ0
}


//��ʼ��ͼ�������fieldpic����������ͼ�񡡣��������ͼ��topfield �����Ƕ���
void rc_init_pict(int fieldpic,int topfield,int targetcomputation)//targetcomputation�Ƿ����Ŀ�껺��ȼ�Ϊ1��ʾ����
{
  int i;

/*compute the total number of basic units in a frame*///���㣱����BU��������Ԫ������Ŀ���ã������ܺ����Ŀ��BU�к����

 // ����˵��ͼ��img->BasicUnit�����������input->basicunit����Ͼ���
	//  �ڳ�ģʽ���峡����Ӧģʽ
	 // img->BasicUnit=input->basicunit*2;
  //ֻ����ģʽ
	 // img->BasicUnit=input->basicunit;
  //Ҳ����˵�ڳ�ģʽ���峡����Ӧģʽ��BU�к����Ŀ�ǣ���ֻ������ģʽ��BU�����Ŀ

  if(input->MbInterlace)//�ڳ�ģʽ���峡����Ӧģʽ
    TotalNumberofBasicUnit=img->Frame_Total_Number_MB/img->BasicUnit;
  img->NumberofCodedMacroBlocks=0;//��ǰͼ��δ���롡���Ա���ĺ����ĿΪ��

/*Normally, the bandwith for the VBR case is estimated by 
a congestion control algorithm. A bandwidth curve can be predefined if we only want to 
test the proposed algorithm*/
  //����Ϊ����VBR�������ʣ��߼����ԡ������������������á����ڣ������ѱ���Pͼ��ʱ�����ʡ�������1.5�������ڣ�����Pͼ�����ʻָ���ǰ
  if(input->channel_type==1)
  {
    if(img->NumberofCodedPFrame==58)
      bit_rate *=1.5;
    else if(img->NumberofCodedPFrame==59)
      PreviousBit_Rate=bit_rate;
  }

  /*predefine a target buffer level for each frame*///Ϊÿ��ͼ��.Ԥ������Ŀ�껺��ˮƽ

  if((fieldpic||topfield)&&targetcomputation)
  {
    switch (img->type)
    {
	case P_SLICE://Pͼ��

		 //   ��Ϊ������ʵʱ�仯�أ�����R�п���λ�����룱��ͼ�񣱸�ͼ��ظ���
		 //	  ����ԭ��ʽTr(i,j)=Tr(i,j-1)+ (u(i,j)-u(i,j-1))*(Ngop-j)/Fr-b(i,j-1)��
		 //   i�����i��gop��j�����j��ͼ��Tr��gop�пɷ����λԪ��u�Ǵ���,Fr�����ʣ�b���ѱ��������λ��
		 //   ���湫ʽ������⣬��ǰͼ����õ�λ����������һ��ͼ�����λ����ȥ����仯���ϸ�ͼ������ı���λ��
       /*Since the available bandwidth may vary at any time, the total number of 
       bits is updated picture by picture*/
        if(PreviousBit_Rate!=bit_rate)
			R +=(int) floor((bit_rate-PreviousBit_Rate)*(Np+Nb)/frame_rate+0.5);//�� (u(i,j)-u(i,j-1))*(Ngop-j)/Fr
              
       /* predefine the  target buffer level for each picture.
       frame layer rate control*///�趨������ʿ���Ϊÿ��ͼ���趨����ˮƽTargetBufferLevel
	   //	Tbl(i,2)=Bc(i,2)      Bc(i,2)�ǵ�i��gop��1��Pͼ�������Ϻ�ʵ�ʻ���ռ����
	   //   TBL(i,j+1)=TBL(i,j)-Tbl(i,2)/Np-1+AWp(i,j)*(L+1)*u(i,j)/Fr*AWp(i,j)*AWb(i,j)*L-u(i,j)/Fr AWp AWbΪ�˼��㸴�Ӷȵ�Ȩ

		if(img->BasicUnit==img->Frame_Total_Number_MB)//���BU�ں����Ŀ�������к����Ŀ
        {
			if(img->NumberofPPicture==1)//�Ѿ�����һ��Pͼ�� Ŀ�껺��ˮƽ���ڵ�ǰ����ռ���� Tbl(i,2)=Bc(i,2)
          {
			//    ��Ϊ��ʼCurrentBufferFullness=0 GOPTargetBufferLevel=0 ����  CurrentBufferFullness += nbits - bit_rate/frame_rate;����λ����ȥ����
			//	  ��ʵҲ����CurrentBufferFullness=BufferSize/8 GOPTargetBufferLevel=BufferSize/8 �ڴ����������������ԣ��
			//	  �������Ա�֤��������ͨ��ģ�ͣ�fluid flow traffic model��HRD�����ǲ������ϻ������������

            TargetBufferLevel=CurrentBufferFullness;
			DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/(TotalPFrame-1);//����DeltaP=-Tbl(i,2)/Np-1 GOPTargetBufferLevel=0��������
			TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
          }
          else if(img->NumberofPPicture>1)
			  TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
        }
        /*basic unit layer rate control*/
        else
        {
			if(img->NumberofCodedPFrame>0)//�ѱ���Pͼ�����Ŀ����0
          {
            /*adaptive frame/filed coding*/
            if(((input->PicInterlace==ADAPTIVE_CODING)||(input->MbInterlace))\
				&&(img->FieldControl==1))//��ͼ���峡����Ӧ �ͺ�鳡��������峡����Ӧ��FieldControl=1��������� һ�仰��鳡�������
			{//TotalNumberofBasicUnit��ͼ����BU��Ŀע�������ͼ�񼴿���ָ��ͼ��Ҳ����ָ��ͼ��
				for(i=0;i<TotalNumberofBasicUnit;i++)//FCBUPFMAD��֮ǰ�����Ʊ���mad,FCBUCFMAD�ǵ�ǰ�����Ʊ���mad
					FCBUPFMAD[i]=FCBUCFMAD[i];//����֮ǰ�����Ʊ���mad,��Ϊ��ͼ���Ѿ��������
            }
            else
            {
				for(i=0;i<TotalNumberofBasicUnit;i++)//�����ģʽ BUPFMAD��֮ǰmad,BUCFMAD�ǵ�ǰmad
                BUPFMAD[i]=BUCFMAD[i];
            }     
          }

			if(img->NumberofGOP==1)//��1��gop��2��pͼ��ʼ����tbl
          {
			  if(img->NumberofPPicture==1)//�Ѿ�����һ��Pͼ�� Ŀ�껺��ˮƽ���ڵ�ǰ����ռ���� Tbl(i,2)=Bc(i,2)
            {
			//	    ��Ϊ��ʼCurrentBufferFullness=0 GOPTargetBufferLevel=0 ����  CurrentBufferFullness += nbits - bit_rate/frame_rate;����λ����ȥ����
			//		��ʵҲ����CurrentBufferFullness=BufferSize/8 GOPTargetBufferLevel=BufferSize/8 �ڴ����������������ԣ��
			//		�������Ա�֤��������ͨ��ģ�ͣ�fluid flow traffic model��HRD�����ǲ������ϻ������������   
              TargetBufferLevel=CurrentBufferFullness;
			  DeltaP=(CurrentBufferFullness-GOPTargetBufferLevel)/(TotalPFrame-1);//����DeltaP=-Tbl(i,2)/Np-1 GOPTargetBufferLevel=0��������
			  TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
            }
			  else if(img->NumberofPPicture>1)//�ѵ�ǰGOP����Pͼ����Ŀ����1��
				TargetBufferLevel -=DeltaP;//TBL(i,j)-DeltaP
          }
			else if(img->NumberofGOP>1)//��2 3....��GOP��1��P��ʼ�ͼ���tbl������GOP��I������ĵ�ǰ����ռ����
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
		//AWp��AWb�ĸ��¹�ʽAWP=AWp/8+7*Wp/8  AWb=AWb/8+7*Wb/8  Wb=b(i,j)*QPb(i,j)/1.3636 Wp=b(i,j)*QPp(i,j) 
		if(img->NumberofCodedPFrame==1)//�����1��Pͼ���ѱ���AWp=Wp
          AWp=Wp;
        //����Ա���Pͼ��С��8�Ҵ���1 ���ѱ���ͼ��1/NumberofCodedPFrame����Ȩ �Ͼ��ѱ���P��ĿС��8,��1/8��Ȩƽ�������� 
        if((img->NumberofCodedPFrame<8)&&(img->NumberofCodedPFrame>1))
            AWp=Wp*(img->NumberofCodedPFrame-1)/img->NumberofCodedPFrame+\
              AWp/img->NumberofCodedPFrame;
		else if(img->NumberofCodedPFrame>1)//�ѱ���Pͼ�����8�Ļ����ѱ���ͼ��1/8����ȨAWP=AWp/8+7*Wp/8
            AWp=Wp/8+7*AWp/8;
          
        //compute the average complexity of B frames
		//����tbl+=AWp(i,j)*(L+1)*u(i,j)/Fr*AWp(i,j)*AWb(i,j)*L  L��ʾ����B����Ŀ ����2��Pͼ��֮��B����Ŀ�ͺ�
		//input->successive_Bframe����L
        if(input->successive_Bframe>0)
        {
          //compute the target buffer level
          TargetBufferLevel +=(AWp*(input->successive_Bframe+1)*bit_rate\
            /(frame_rate*(AWp+AWb*input->successive_Bframe))-bit_rate/frame_rate);
        }
        
        break;

         case B_SLICE:
         /* update the total number of bits if the bandwidth is changed*/
			 //����ԭ��ʽTr(i,j)=Tr(i,j-1)+ (u(i,j)-u(i,j-1))*(Ngop-j)/Fr-b(i,j-1)
			 if(PreviousBit_Rate!=bit_rate)//����䶯�����Gop�п���λ��R��������,����ֵu(i,j)-u(i,j-1))*(Ngop-j)/Fr
             R +=(int) floor((bit_rate-PreviousBit_Rate)*(Np+Nb)/frame_rate+0.5);
			 if((img->NumberofCodedPFrame==1)&&(img->NumberofCodedBFrame==1))//�����1��P�ͱ����1��Bͼ������³�ʼ�����ӶȲ���AWp

          {
            AWp=Wp;
            AWb=Wb;
          }
			 else if(img->NumberofCodedBFrame>1)//����ѱ����B��Ŀ����1
          {
            //compute the average weight
			  //AWb�ĸ��¹�ʽAWb=AWb/8+7*Wb/8 ���ѱ����B��Ŀ��С��8 ��Ȩϵ��ȡ1/NumberofCodedBFrame
            if(img->NumberofCodedBFrame<8)
              AWb=Wb*(img->NumberofCodedBFrame-1)/img->NumberofCodedBFrame+\
                AWb/img->NumberofCodedBFrame;
            else
              AWb=Wb/8+7*AWb/8;
          }

            break;
    }
     /*Compute the target bit for each frame*/
	// ����ÿ�����Ŀ��λ��,ΪʲôҪ��һ����,�������
    //	���ʵ�ʻ���ĵȼ���ͬ�����ǵ�Ŀ�껺��ȼ�,��Ȼ���Ա�֤ÿ��GOP����ȷ���䵽λ��,
	//	����R-Dģ�ͺ�madԤ��ģ���ǲ���ȷ��,����ʵ�ʻ���ˮƽ��Ŀ��λԪ�ȼ������ϴ�,
	//	��ô������Ҫ����ÿ��ͼ���Ŀ��λ����������ʵ�ʻ���ˮƽ��Ŀ��λԪ֮������,�����ϸ�ڿ���
	//	��ͬ���������浹ˮ,���Ҫ���õ�1��ˮ, ��Ȼ1��1�ε������!
	if(img->type==P_SLICE)//Pͼ��ģʽ
    {
      /*frame layer rate control*/
		if(img->BasicUnit==img->Frame_Total_Number_MB)//BU�к����Ŀ�������к����Ŀ
      {
		  if(img->NumberofCodedPFrame>0)//�ѱ����P��Ŀ����1
        {

		//	    ͨ��Gop��R����λ������T��ͼ�����λ��,��ʽ����
		//		f(i,j)=Wp(i,j-1)*Tr/(Wp(i,j-1)*Np+Wb(i,j-1)*Nb) 
		//		Ҳ����T = Wp*R/(Np*Wp+Nb*Wb)

			T = (long) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);//T��ͼ�����λ��,��Ȼ��Gop��R����λ��
		//	    ���濼�ǵ���i��Gop�ĵ�J��ͼ���Ŀ�����λ����  
		//		Ŀ�껺��ȼ�Tbl=TargetBufferLevel ʵ�ʻ�����ռ����Bc(i,j)=CurrentBufferFullness Ƶ��u(i,j)=bit_rate    
		//		GAMMAP����Bͼ��������GAMMAP=0.25;  ��Bͼ��������GAMMAP=0.5;
		//	    Ϊʲô��Bͼ�����Bͼ��GAMMAPҪС,���������,B�������λ����,��Ȼ��Ҫ��ȥ��ʵ�ʺ�Ŀ������ϵ��С
		//		��ʽf1=u(i,j)/Fr+gamma*(Tbl(i,j)-Bc(i,j))

                
          T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1=MAX(0,T1);
		//  Ŀ�����λ��f(i,j)��f(i,j)��f1��Ȩ�õ�
		//	  f(i,j)=beta*f(i,j)+(1-beta)*f1
		//	  beta��û��Bͼ�������beta=0.5          ��Bͼ����������0.9��˵����Bͼ�������Ҫ��f(i,j)�������
          T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
       }
		/*basic unit layer rate control*///������Ԫ�����
      else
      {
		  if((img->NumberofGOP==1)&&(img->NumberofCodedPFrame>0))//��1��GOP�����ѱ���P��Ŀ����0
        {
			//  ͨ��Gop��R����λ������T��ͼ�����λ��,��ʽ����
			//	:f(i,j)=Wp(i,j-1)*Tr/(Wp(i,j-1)*Np+Wb(i,j-1)*Nb) 

          T = (int) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (int) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1=MAX(0,T1);
          T = (int)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
		  else if(img->NumberofGOP>1)//��2  3  .....��GOP
        {
          T = (long) floor(Wp*R/(Np*Wp+Nb*Wb) + 0.5);
          T1 = (long) floor(bit_rate/frame_rate-GAMMAP*(CurrentBufferFullness-TargetBufferLevel)+0.5);
          T1 = MAX(0,T1);
          T = (long)(floor(BETAP*T+(1.0-BETAP)*T1+0.5));
        }
      }

		/*reserve some bits for smoothing*///����һЩλ��Ϊ��ƽ��
		//	�ܻ��ư�0.0*input->successive_Bframe����0.0�Ǹ��ɵ�ϵ�� ����T�п���λ��

      T=(long)((1.0-0.0*input->successive_Bframe)*T);
	  /*HRD consideration*///��֤ͼ�����λ��T��Զ��[LowerBound,UpperBound2]��Χ��
	  //�����LowerBound��UpperBound2��ͼ��������������rc_update_pict���µ�
	  //ͨ��������ȥ����λ��,���������������(ͼ�����λ��),ƿ�ӿյĸ����� ��ʾ�и�����õ�λ�� ��ô��ǰͼ�����λ��������Ӧ������
	  //���������������(ͼ�����λ��),��ʾ�п��õ�λ������ ƿ�ӿռ������ ��ô��ǰͼ�����λ��������Ӧ�ü���       
	  //LowerBound  +=(long)(bit_rate/frame_rate-nbits);
	  //UpperBound1 +=(long)(bit_rate/frame_rate-nbits); UpperBound2 = (long)(OMEGA*UpperBound1);         OMEGA=0.9��֤ԣ��
      T = MAX(T, (long) LowerBound);
        T = MIN(T, (long) UpperBound2);

      if((topfield)||(fieldpic&&((input->PicInterlace==ADAPTIVE_CODING)\
		  ||(input->MbInterlace))))//���������
        T_field=T;
    }
  }

  if(fieldpic||topfield)//�������߶�������
  {
	  /*frame layer rate control*///ͼ������ʿ��� Ҳ��������ʿ���
	  img->NumberofHeaderBits=0;//��δ���� �Ȱ�ͷ��Ϣλ��NumberofHeaderBits��������Ϣλ��NumberofTextureBits��0
    img->NumberofTextureBits=0;

	/*basic unit layer rate control*///ֻ����BU�к��С�����к����Ŀ�Ŵ���BU������Ԫ����
    if(img->BasicUnit<img->Frame_Total_Number_MB)
    {
		TotalFrameQP=0;//ͼ��������BU��QP��
		img->NumberofBasicUnitHeaderBits=0;//ͷ��Ϣλ��NumberofHeaderBits
		img->NumberofBasicUnitTextureBits=0;//������Ϣλ��NumberofTextureBits
		img->TotalMADBasicUnit=0;//ͼ��������BU��MAD��TotalMADBasicUnit
	  if(img->FieldControl==0)//�����
        NumberofBasicUnit=TotalNumberofBasicUnit;
      else//������
        NumberofBasicUnit=TotalNumberofBasicUnit/2;
    }
  }
    
  if((img->type==P_SLICE)&&(img->BasicUnit<img->Frame_Total_Number_MB)\
    &&(img->FieldControl==1))
  {
  /*top filed at basic unit layer rate control*///����������BU�����ʿ���
    if(topfield)
    {
      bits_topfield=0;
      T=(long)(T_field*0.6);//��������λ��//�����͵׳���һ���Ŀ���λ��T_field 0.6ϵ��˵������������Ҫ�� �������ο���
    }
  /*bottom filed at basic unit layer rate control*/
    else
    {
      T=T_field-bits_topfield;//�����õ�λ��T_field��ȥ�������������λ��bits_topfield�õ��׳����õ�λ��
      img->NumberofBasicUnitHeaderBits=0;//ͷ��Ϣλ��NumberofHeaderBits
      img->NumberofBasicUnitTextureBits=0;//������Ϣλ��NumberofTextureBits
      img->TotalMADBasicUnit=0;//ͼ��������BU��MAD��TotalMADBasicUnit
      NumberofBasicUnit=TotalNumberofBasicUnit/2;//������ÿ������BU��Ŀ��ȻΪ��BU��Ŀ��1��
    }
  }
}

//calculate MAD for the current macroblock ���㵱ǰ����MAD,diff������Yԭֵ��ȥԤ��ֵ
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

// update one picture after frame/field encoding���峡��������һ��ͼ�����
void rc_update_pict(int nbits)
{
  R-= nbits; /* remaining # of bits in GOP *///GOP��ʣ�����λ��R,�����ȥ�ϸ�ͼ��������λ��
  CurrentBufferFullness += nbits - bit_rate/frame_rate;//��ǰ���������������ȥ��������

  /*update the lower bound and the upper bound for the target bits of each frame, HRD consideration*/
  //ͨ��������ȥ����λ��,���������������(ͼ�����λ��),ƿ�ӿյĸ����� ��ʾ�и�����õ�λ�� ��ô��ǰͼ�����λ��������Ӧ������
  //���������������(ͼ�����λ��),��ʾ�п��õ�λ������ ƿ�ӿռ������ ��ô��ǰͼ�����λ��������Ӧ�ü���
  LowerBound  +=(long)(bit_rate/frame_rate-nbits);
  UpperBound1 +=(long)(bit_rate/frame_rate-nbits);
  UpperBound2 = (long)(OMEGA*UpperBound1);
  
  return;
}

// update after frame encoding�����������RC
void rc_update_pict_frame(int nbits)
{

/*update the
complexity weight of I, P, B frame*///����I P  B�ĸ��Ӷ�Ȩ
  int Avem_Qc;// Avem_Qc������BU��QPƽ��ֵ
  int X;
    
/*frame layer rate control*///ͼ������ʿ���
// I P  B�ĸ��Ӷ�ȨWb=b(i,j)*QPb(i,j)/1.3636 Wp=b(i,j)*QPp(i,j)  nbits=b(i,j)

  if(img->BasicUnit==img->Frame_Total_Number_MB)//BU�ĺ����Ŀ�������ں����Ŀ
    X = (int) floor(nbits*m_Qc+ 0.5);//ע����m_Qc nbits�������λ��
/*basic unit layer rate control*///������Ԫ�����ʿ���
  else
  {
    if(img->type==P_SLICE)//Pͼ��ģʽ
    {
      if(((img->IFLAG==0)&&(img->FieldControl==1))         //\�������Ҳ��ǵ�1��GOP�ĵ�1��I��ͼ��           IFLAG=1���ܴ���seq�ĵ�1����
        ||(img->FieldControl==0))//�����
      {
        Avem_Qc=TotalFrameQP/TotalNumberofBasicUnit;//���ƽ��������ԪQc
        X=(int)floor(nbits*Avem_Qc+0.5);
      }
    }
    else if(img->type==B_SLICE)//Bͼ��ģʽ
      X = (int) floor(nbits*m_Qc+ 0.5);
  }

//I P  B�ĸ��Ӷ�ȨWb=b(i,j)*QPb(i,j)/1.3636 Wp=b(i,j)*QPp(i,j)  nbits=b(i,j)
  switch (img->type)
  {
    case P_SLICE:
 /*filed coding*/
      if(((img->IFLAG==0)&&(img->FieldControl==1))      //\ �������Ҳ��ǵ�1��GOP�ĵ�1��I��ͼ��
        ||(img->FieldControl==0))//��ͼ��
      {
        Xp = X;
        Np--;//�ѱ���Pͼ��������δ����Pͼ���1
        Wp=Xp;
        Pm_Hp=img->NumberofHeaderBits;//ͷ��Ϣλ��
        img->NumberofCodedPFrame++;//�ѱ���Pͼ����ĿNumberofCodedPFrame
        img->NumberofPPicture++;//�ѱ���Pͼ����ĿNumberofPPicture

      }
      else if((img->IFLAG!=0)&&(img->FieldControl==1))//��1��GOP�ĵ�1��I��ͼ���Ѿ�����,��ô�˱�־��0
        img->IFLAG=0;
        break;
        case B_SLICE:
        Xb = X;
      Nb--;//�ѱ���Bͼ��������δ����Pͼ���1
        Wb=Xb/THETA; //����Wb=b(i,j)*QPb(i,j)/1.3636 THETA=1.3636
      img->NumberofCodedBFrame++;//�ѱ���Bͼ����ĿNumberofCodedBFrame
      NumberofBFrames++;//�ѱ���Pͼ����ĿNumberofCodedBFrame

        break;
    }
}

// coded bits for top field�������������λ��bits_topfield
void setbitscount(int nbits)
{
  bits_topfield = nbits;
}

//compute a  quantization parameter for each frame����ÿ����������

int updateQuantizationParameter(int topfield)
{
  double dtmp;//�ǵó���ѧ���Ķ��η���ʽ�� ��ʵ����delta=sqrt(b*b-4*a*c)

  int m_Bits;
  int BFrameNumber;
  int StepSize;
  int PAverageQP;
  int SumofBasicUnit;
  int i;
  
/*frame layer rate control*///������ʿ���

  if(img->BasicUnit==img->Frame_Total_Number_MB)//BU�ں����Ŀ�������ں����Ŀ ��Ȼ������BU����������

  {
/*fixed quantization parameter is used to coded I frame, the first P frame and the first B frame
  the quantization parameter is adjusted according the available channel bandwidth and 
  the type of vide*/  
/*top field*/
    if((topfield)||(img->FieldControl==0))//�Ƕ������߳������ź�Ϊ0

    {
      if(img->type==I_SLICE)//I������� m_qcȡ��ʼֵ

      {
        m_Qc=MyInitialQp;
        return m_Qc;
      }
      else if(img->type==B_SLICE)//B����

      {
        if(input->successive_Bframe==1)//P��P֮��ֻ��1��B���

        {//successive_Bframe=1��1�����κ���������0������ľͿ�Ц��,��������

          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
          if(BFrameNumber==0)//������������BFrameNumber=1˵��BFrameNumber��ʾ��2��Pͼ���еĵڼ�����ͼ�����������밢

            BFrameNumber=input->successive_Bframe;//��Ϊsuccessive_Bframe����˵����2��Pͼ����ֻ�У�����ͼ��BFrameNumberֻ����Ϊ��

          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {//��L=1��ʾ2��Pͼ��֮����1��Bͼ�� ��ʽ����:
        //  QB1=(QP1+QP2+2)/2  ��QP1!=QP2
         // QB1=QP1+2          ��QP1 =QP2

            if((input->PicInterlace==ADAPTIVE_CODING)\
              ||(input->MbInterlace))
            {
              if(img->FieldControl==0)//�����ƹر�

              {                   
                /*previous choice is frame coding*/
                if(img->FieldFrame==1)//֮ǰѡ���������

                {
                  PreviousQp1=PreviousQp2;//����P1B1P2B2P3��ǰ��B1��PreviousQp1����P1ͼ��,PreviousQp2=P2

                  PreviousQp2=FrameQPBuffer;
                }           
                /*previous choice is field coding*/
                else//֮ǰѡ���ǳ�����

                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }
//          m_Qc=QB1=(QP1+QP2+2)/2  ��QP1!=QP2
          //          m_Qc=QB1=QP1+2          ��QP1 =QP2

          if(PreviousQp1==PreviousQp2)
            m_Qc=PreviousQp1+2;
          else
            m_Qc=(PreviousQp1+PreviousQp2)/2+1;
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clippingǯλ0-51

          m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        else
        {
//L=successive_Bframe>1 ����������������� ˵��BFrameNumber��ʾ��2��Pͼ���еĵڼ�����ͼ��
    //      ����������L=2 PB1B2PB3B4P��ôNumberofBFrames=2�ѱ����B��Ŀ��2��ô��ǰ�ǵ�3��Bʵ��B3,BFrameNumber=1��ʵ��P��P��λ����1
       //    �ѱ����B��Ŀ��3��ô��ǰ�ǵ�3��Bʵ��B3,BFrameNumber=2
       //   ͨ���Ժ���B����Ŀ���� ��õ�ǰB��P ��P��λ��


          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
          if(BFrameNumber==0)//������NumberofBFrames=3ʱ�����ڵ�BFrameNumber������������Ϊ��

            BFrameNumber=input->successive_Bframe;
          /*adaptive field/frame coding*/
          else if(BFrameNumber==1)
          {//����֮ǰ����Pͼ���QPֵ��Ϊ�ǵڣ���B������Ҫ���£�����B��ڣ���B���õĲο�P��һ�µ�������

            if((input->PicInterlace==ADAPTIVE_CODING)\
              ||(input->MbInterlace))//ͼ������Ӧ�����֡����

            {
              if(img->FieldControl==0)//�������źŹرյ�ǰͼ�����ʹ�����Ҳ�����ǳ�����

              {
                /*previous choice is frame coding*///����֮ǰ����Pͼ���QPֵ

                if(img->FieldFrame==1)//ͼ��ѡ���������
                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FrameQPBuffer;
                }
                /*previous choice is field coding*/
                else//ͼ��ѡ���ǳ�����

                {
                  PreviousQp1=PreviousQp2;
                  PreviousQp2=FieldQPBuffer;
                }
              }
            }
          }
// ��L>1��ʾ����Pͼ��������Bͼ�� ��ʽ����:
//         QBi=QP1+alpha+max{min{(QP2-QP1)/(L-1),2*(i-1)},-2*(i-1)}
//         alphaȡֵ����:
//         alpha=-3 if QP2-QP1<=-2*L-3
//         alpha=-2 if QP2-QP1=-2*L-2
//         alpha=-1 if QP2-QP1=-2*L-1
//         alpha=0  if QP2-QP1=-2*L
//         alpha=1  if QP2-QP1=-2*L+1
//         alpha=2   �������
//         Ϊʲô��Ҫalpha=StepSize����QB,�����ǰ��գ���Pͼ���Qp��deltaQ��L�Ĺ�ϵ������ƽ����������ͼ������������

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
// ���水��max{min{(QP2-QP1)/(L-1),2*(i-1)},-2*(i-1)}���������Ʊ�֤����˳�򿿽��ڣ���P���������ӽӽ��ڣ���P��Qp,
//           �ӽ��ڣ���Pͼ����������ӽӽ��ڣ���P��Qp��Ҳ��Ϊ������ƽ����������ͼ������������

          m_Qc +=MIN(2*(BFrameNumber-1),MAX(-2*(BFrameNumber-1), \
            (BFrameNumber-1)*(PreviousQp2-PreviousQp1)/(input->successive_Bframe-1)));
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        return m_Qc;
      }
      else if((img->type==P_SLICE)&&(img->NumberofPPicture==0))//��ǰ��Pͼ�����ѱ����Pͼ����Ŀ�ǣ�

      {
        m_Qc=MyInitialQp;//MyInitialQp��GOP��QP

        
        if(img->FieldControl==0)//�����ƹر�

        {
          if(active_sps->frame_mbs_only_flag)//ֻ����ģʽ�������

          {
            img->TotalQpforPPicture +=m_Qc;//TotalQpforPPicture��ʾGOP���ѱ���Pͼ��Qp��

            PreviousQp1=PreviousQp2;//����2��Pͼ��Qc

            PreviousQp2=m_Qc;
            Pm_Qp=m_Qc;
          }
          /*adaptive field/frame coding*///����Ӧ�峡ģʽ

          else
            FrameQPBuffer=m_Qc;//FrameQPBuffer��QP����

        }
        
        return m_Qc;  
      }
      else
      {
        /*adaptive field/frame coding*/
        if(((input->PicInterlace==ADAPTIVE_CODING)//\ͼ������Ӧ�峡ģʽ

          ||(input->MbInterlace))//\��鳡ģʽ������Ӧ�峡ģʽ

          &&(img->FieldControl==0))//�����ƹر�

        {
          /*previous choice is frame coding*/
          if(img->FieldFrame==1)//��ģʽ
          {
            img->TotalQpforPPicture +=FrameQPBuffer;//TotalQpforPPicture��ʾGOP���ѱ���Pͼ��Qp��

            Pm_Qp=FrameQPBuffer;
          }
          /*previous choice is field coding*/
          else//��ģʽ

          {
            img->TotalQpforPPicture +=FieldQPBuffer;
            Pm_Qp=FieldQPBuffer;
          }
        }
// ������2��r-dģ��ϵ�� ����ģ���Ƶ� ��taylor�����Ƶ� �Ժ���ʱ������д��RDģ���Ƶ�
//         (R-H)/MAD=C1/Qstep+C2/(Qstep*Qstep) ����m_X1 m_X2����C1 C2
//         ������(T-H)/MAD=MADPictureC1/Qstep+MADPictureC2/(Qstep*Qstep) 

        m_X1=Pm_X1;
        m_X2=Pm_X2;
        m_Hp=PPreHeader;//֮ǰ��ͷ��Ϣλ��

        m_Qp=Pm_Qp;//m_QP��m_���ܴ����飬��ΪBU������������Ŀ������m_qp����Pͼ��Pm_Qp

        DuantQp=PDuantQp;
        MADPictureC1=PMADPictureC1;//PMADPictureC1��Ԥ��MAD

        MADPictureC2=PMADPictureC2;//PMADPictureC2��Ԥ��MAD

        PreviousPictureMAD=PPictureMAD[0];//PreviousPictureMAD֮ǰͼ��MAD

        
        /* predict the MAD of current picture*///Ԥ�⵱ǰͼ��MAD
// MADc=a1*MADp+a2 MADPictureC1=a1 MADPictureC2=a2 MADp��ǰ1��ͼ���Ӧλ�õ�MAD,Ԥ������������Ԫ��MAD
//         ��ΪBU������������Ŀ������ֱ����ͼ���MADֵPreviousPictureMAD

        CurrentFrameMAD=MADPictureC1*PreviousPictureMAD+MADPictureC2;
        
        /*compute the number of bits for the texture*/      
        
        if(T<0)//T��ͼ�����λ�������С��0����Ҫ����ǰͼ��QP������λ��

        {
          m_Qc=m_Qp+DuantQp;//����ǰͼ��QP

          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping��֤�����������������51

        }
        else
        {//��ʽ(T-H)/MAD=MADPictureC1/Qstep+MADPictureC2/(Qstep*Qstep)

          m_Bits =T-m_Hp;//ȥ��ͷ��Ϣ

//MINVALUE=4.0;     �ڻ�����Ԫ�� ���������Ԫ��λ��R R= MAX(R, (int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));
          //����TotalNumberofBasicUnit=1���Բ����ڻ�����Ԫ�� ֱ��������


          m_Bits = MAX(m_Bits, (int)(bit_rate/(MINVALUE*frame_rate)));
          dtmp = CurrentFrameMAD * m_X1 * CurrentFrameMAD * m_X1 //\dtmp��2�η���ʽdelta����ⶼ���� ��ȥ������

            + 4 * m_X2 * CurrentFrameMAD * m_Bits;
//������2�η���ʽ��� ������ѧ ����ϸ����

          if ((m_X2 == 0.0) || (dtmp < 0) || ((sqrt (dtmp) - m_X1 * CurrentFrameMAD) <= 0.0)) // fall back 1st order mode
            m_Qstep = (float) (m_X1 * CurrentFrameMAD / (double) m_Bits);
          else // 2nd order mode
            m_Qstep = (float) ((2 * m_X2 * CurrentFrameMAD) / (sqrt (dtmp) - m_X1 * CurrentFrameMAD));
//           ��������m_Qstep����������QP��ϵ����������+6��ô����������������������+1��ô������������12.5%
//           �����ǵ�uniform��һ����������� ����Ϊʲô�������ڸ����ܶȺ����������������ʷֲ���ʧ��Ĺ�ϵ��
//           ���0-51����������m_qc��

          m_Qc=Qstep2QP(m_Qstep);
          //�����Ǳ�֤+-2��1-51����ǯλ ��֤ͼ������ǰ���𲻴���������������������

          m_Qc = MIN(m_Qp+DuantQp,  m_Qc);  // control variation
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clipping
          m_Qc = MAX(m_Qp-DuantQp, m_Qc); // control variation
          m_Qc = MAX(RC_MIN_QUANT, m_Qc);
        }
        
        if(img->FieldControl==0)//���źſ��ƹر�

        {
          /*frame coding*/
          if(active_sps->frame_mbs_only_flag)//ֻ�к����ģʽ

          {
            img->TotalQpforPPicture +=m_Qc;//TotalQpforPPicture��ʾGOP���ѱ���Pͼ��Qp��

            PreviousQp1=PreviousQp2;//2��Pͼ���QP�м��������

            PreviousQp2=m_Qc;
            Pm_Qp=m_Qc;
          }
          /*adaptive field/frame coding*///����Ӧ�峡ģʽ

          else
            FrameQPBuffer=m_Qc;
        }
        
        return m_Qc;
      }
   }
   /*bottom field*///�׳���˵�� д������

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
  /*basic unit layer rate control*///������Ԫ���ʿ��� �����ǰ��������Ʋ�����˵�� 
  else
  {
    /*top filed of I frame*///����

    if(img->type==I_SLICE)
    {
      m_Qc=MyInitialQp;//ֱ����GOP��QP

      return m_Qc;
    }
    /*bottom field of I frame*///��P���������ƿ���//���д��޸�ע��

    else if((img->type==P_SLICE)&&(img->IFLAG==1)&&(img->FieldControl==1))
    {
      m_Qc=MyInitialQp;
      return m_Qc;
    }
    else if(img->type==B_SLICE)//B����

    {
      /*top filed of B frame*///����
      if((topfield)||(img->FieldControl==0))//�����ƹر�

      {
        if(input->successive_Bframe==1)
        {//BFrameNumberֻ��Ϊ1ǰ���Ѿ���������

          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
        if(BFrameNumber==0)
          BFrameNumber=input->successive_Bframe;
        /*adaptive field/frame coding*/
        else if(BFrameNumber==1)
        {
          if((input->PicInterlace==ADAPTIVE_CODING)//\ͼ������Ӧ�峡����

            ||(input->MbInterlace))//�������Ӧ�򳡱���

          {
            if(img->FieldControl==0)
            {       // ����֮ǰ����Pͼ���QPֵ
     
        /*previous choice is frame coding*/
              if(img->FieldFrame==1)//��ǰͼ���������

              {
                PreviousQp1=PreviousQp2;
                PreviousQp2=FrameQPBuffer;
              }
        /*previous choice is field coding*/
              else//��ǰͼ���ǳ�����

              {
                PreviousQp1=PreviousQp2;
                PreviousQp2=FieldQPBuffer;
              }
            }
          }
        }
//          m_Qc=QB1=(QP1+QP2+2)/2  ��QP1!=QP2
          //          m_Qc=QB1=QP1+2          ��QP1 =QP2

          if(PreviousQp1==PreviousQp2)
            m_Qc=PreviousQp1+2;
          else
            m_Qc=(PreviousQp1+PreviousQp2)/2+1;
          m_Qc = MIN(m_Qc, RC_MAX_QUANT); // clippingǯλ0-51

          m_Qc = MAX(RC_MIN_QUANT, m_Qc);//clipping
        }
        else//����B�����1�����

        {////BFrameNumberֻ��Ϊ1-successive_Bframeǰ���Ѿ���������,����2��P֮��ڼ���Bͼ��

          BFrameNumber=(NumberofBFrames+1)%input->successive_Bframe;
        if(BFrameNumber==0)
          BFrameNumber=input->successive_Bframe;
        /*adaptive field/frame coding*/
        else if(BFrameNumber==1)
        {
          if((input->PicInterlace==ADAPTIVE_CODING)\

            ||(input->MbInterlace))//ͼ������Ӧ�����֡����
          {//����֮ǰ����Pͼ���QPֵ

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
// �μ�ǰ��д��alphaȡֵ��ΪʲôҪ����ȡֵ

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
		//         //��֤����˳�򿿽��ڣ���P���������ӽӽ��ڣ���P��Qp,
//�ӽ��ڣ���Pͼ����������ӽӽ��ڣ���P��Qp��Ҳ��Ϊ������ƽ����������ͼ������������

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
    else if(img->type==P_SLICE)//P����

    {
      if((img->NumberofGOP==1)&&(img->NumberofPPicture==0))//��1��GOP�����Ѿ�����P��ĿΪ0,Ҳ����GOP������Pͼ�񶼻�û����ʱ��

      {
        if((img->FieldControl==0)||((img->FieldControl==1)\
          &&(img->IFLAG==0)))
        {
        /*top field of the first P frame*///��1��P�Ķ���,ȡGOP��QP��ʼֵ,��ΪPδ����,
 //       ���԰ѻ�����Ԫ��ͷ��ϢNumberofBasicUnitHeaderBits��������ϢNumberofBasicUnitTextureBits������

          m_Qc=MyInitialQp;
          img->NumberofBasicUnitHeaderBits=0;
          img->NumberofBasicUnitTextureBits=0;
          NumberofBasicUnit--;//BU��Ŀ��1

        /*bottom field of the first P frame*///��1��P�ĵ׳�

          if((!topfield)&&(NumberofBasicUnit==0))
          {
            /*frame coding or field coding*/
            if((active_sps->frame_mbs_only_flag)||(input->PicInterlace==FIELD_CODING))//ֻ�к������� ����ͼ�񳡱���

            {
              img->TotalQpforPPicture +=m_Qc;//��ǰP��Qc�����͸���

              PreviousQp1=PreviousQp2;
              PreviousQp2=m_Qc;
              PAveFrameQP=m_Qc;
// ������PAveHeaderBits2 PAveHeaderBits3 PAveHeaderBits1�����ڼ������е�Ԫ��ͷ��Ϣ��λ��
//               
//               Mhdrl��l��������Ԫʵ�ʲ�����ͷ��Ϣλ��
//               AveMhdrl=AveMhdrl*(1-1/l)+Mhdrl/l �ɼ�l=1ʱAveMhdrl=Mhdrl
//               PMhdrl��ǰһ��ͼ�����л�����ԪԤ��õ��� PMhdrl=AveMhdrl/Nunit+PMhdrl*(1-1/Nunit)��ʼ��PMhdrl=0
//               PAveHeaderBits2�Ǳ�ʾԤ��BU��ͷ��Ϣλ��Ҳ����ƽ��BU��ͷ��Ϣλ��


              PAveHeaderBits3=PAveHeaderBits2;
            }
            /*adaptive frame/field coding*/
            else if((input->PicInterlace==ADAPTIVE_CODING)\
              ||(input->MbInterlace))
            {
              if(img->FieldControl==0)
              {
                FrameQPBuffer=m_Qc;
                FrameAveHeaderBits=PAveHeaderBits2;//ƽ��BU��ͷ��Ϣλ��

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

        if(img->FieldControl==0)//�����ƹر� ������Ԫ��Ŀ����

          SumofBasicUnit=TotalNumberofBasicUnit;
        else//�����ƿ��� ������Ԫ��Ŀ����

          SumofBasicUnit=TotalNumberofBasicUnit/2;

        /*the average QP of the previous frame is used to coded the first basic unit of the current frame or field*/
        if(NumberofBasicUnit==SumofBasicUnit)//ͼ���ڻ�����Ԫ������ϼ���֮ǰͼ���ƽ��Qp���ڵ�1��BU��Ԫ����

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

          if(T<=0)//��С�ڣ���ͼ��������λ��С�ڣ�

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
           {//�����l��BU��Ŀ��λ��Bi=Tr*MADl*MADl/(sum(k=l-Nunit)MADk)����sum(k=l-Nunit)MADk���ѱ���BU��MAD��

             /*predict the MAD of current picture*/
             if(((input->PicInterlace==ADAPTIVE_CODING)||(input->MbInterlace))\
               &&(img->FieldControl==1))
             {//Ԥ�⵱ǰMADͨ��֮ǰFCBUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]��FCBUPFMAD�ȶԵ�TotalNumberofBasicUnit������1���洢��

               CurrentFrameMAD=MADPictureC1*FCBUPFMAD[TotalNumberofBasicUnit-NumberofBasicUnit]+MADPictureC2;
               TotalBUMAD=0;
			   ////TotalNumberofBasicUnit-NumberofBasicUni��ʾ�ѱ����BU��Ŀ

               for(i=TotalNumberofBasicUnit-1; i>=(TotalNumberofBasicUnit-NumberofBasicUnit);i--)
               {//����ǰ�湫ʽ���µ�ǰ������ԪMADc=a1*MADp+a2

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
             /*compute the number of texture bits*///��ȥͷ��Ϣλ���õ�������Ϣ����λ��
// Mhdrl��l��������Ԫʵ�ʲ�����ͷ��Ϣλ��
//               AveMhdrl=AveMhdrl*(1-1/l)+Mhdrl/l �ɼ�l=1ʱAveMhdrl=Mhdrl
//               PMhdrl��ǰһ��ͼ�����л�����ԪԤ��õ��� PMhdrl=AveMhdrl/Nunit+PMhdrl*(1-1/Nunit)��ʼ��PMhdrl=0
//               ���� PAveHeaderBits2=PMhdrlҲ���Ǳ�ʾԤ���ͷ��Ϣλ��    
             m_Bits -=PAveHeaderBits2;
             // �ڻ�����Ԫ�� ���������Ԫ��λ��R R= MAX(R, (int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));

             m_Bits=MAX(m_Bits,(int)(bit_rate/(MINVALUE*frame_rate*TotalNumberofBasicUnit)));
        //m_qc������˵�� ������

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
             if(input->basicunit>=MBPerRow)//ǰ��˵��������Ԫ�ں����Ŀ���ڣ����к����Ŀ������QP��Ϊ��ƽ��ͼ������

               m_Qc = MAX(PAveFrameQP-6, m_Qc);
             else
               m_Qc = MAX(PAveFrameQP-3, m_Qc);

             m_Qc = MAX(RC_MIN_QUANT, m_Qc);
             TotalFrameQP +=m_Qc;
             Pm_Qp=m_Qc;
             NumberofBasicUnit--;
             if((NumberofBasicUnit==0)&&(img->type==P_SLICE))//NumberofBasicUnit��ʾ���1��BU

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
