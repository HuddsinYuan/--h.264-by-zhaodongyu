#include <stdio.h>
#include <math.h>
/*#include <Windows.h>*/
#define b(i) ((i%2)==0?(i/2):((i-1)/2))
#define BOOL int
#define TRUE1 1
#define FALSE1 0

int width=352;
int height=288;

//�Ҷȸ�ʴ�ӳ���
void grayerosion(unsigned char *imputimage, int *cake, int cakeheight, int cakewidth, unsigned char *outputimage)   

{ 
	int px,py;  
	int gray;//��ʾ��ʴģ��
	int temp;
	int i,j;
	unsigned char *buf=(unsigned char*)malloc(sizeof(unsigned char)*(width+cakewidth)*(height+cakeheight));
    for(i=0;i<height;i++)
		for(j=0;j<width;j++)
			outputimage[i*width+j]=imputimage[i*width+j];
	for(i=0;i<height+cakeheight;i++)
		for(j=0;j<width+cakewidth;j++)
		{
			if(i<cakeheight/2||i>=height+cakeheight/2||j<cakewidth/2||j>=width+cakewidth/2)
				*(buf+i*(width+cakewidth)+j)=255;
			else
				*(buf+i*(width+cakewidth)+j)=*(imputimage+(i-cakeheight/2)*width+j-cakewidth/2);
		}					  
		
		//������㸯ʴ          
		for(i=cakeheight/2;i<height+cakeheight/2;i++)
			for(j=cakewidth/2;j<width;j++)
			{  
				gray=255;
				for(py=-cakeheight/2;py<=cakeheight/2;py++)
					for(px=-cakewidth/2;px<=cakewidth/2;px++)
					{  
						temp=*(buf+(i+py)*(width+cakewidth)+j+px)-*(cake+(cakeheight/2+py)*cakewidth+cakewidth/2+px);
						if(temp<gray)
							gray=temp;
						
					}
					
					if(gray<0)  *(outputimage+(i-cakeheight/2)*width+j-cakewidth/2)=0;
					else  *(outputimage+(i-cakeheight/2)*width+j-cakewidth/2)=gray;
					
			}
			
			
}

//�Ҷ������ӳ���
void grayinflation(unsigned char *imputimage, int *cake, int cakeheight, int cakewidth, unsigned char *outputimage)   

{   
	int px,py;
	
	//��ʾ����ģ��
	int gray;
	int temp;
	
	///�����������
	int i,j;
	unsigned char *pbuf=(unsigned char*)malloc(sizeof(unsigned char)*(width+cakewidth)*(height+cakeheight));
	for(i=0;i<height;i++)
		for(j=0;j<width;j++)
			outputimage[i*width+j]=imputimage[i*width+j];

	for(i=0;i<height+cakeheight;i++)
		for(j=0;j<width+cakewidth;j++)
		{
	
			if((i>cakeheight/2)&&(i<height+cakeheight/2)&&(j>cakewidth/2)&&(j<width+cakewidth/2))    
				*(pbuf+i*(width+cakewidth)+j)=*(imputimage+(i-cakeheight/2)*width+j-cakewidth/2);
			else 
				*(pbuf+i*(width+cakewidth)+j)=0;
		}					  
		
		
		for(i=cakeheight/2;i<height+cakeheight/2;i++)
			for(j=cakewidth/2;j<width;j++)
			{  
				gray=0;
				for(py=-cakeheight/2;py<=cakeheight/2;py++)
					for(px=-cakewidth/2;px<=cakewidth/2;px++)
					{  
						temp=*(pbuf+(i+py)*(width+cakewidth)+j+px)+*(cake+(cakeheight/2+py)*cakewidth+cakewidth/2+px);
						if(temp>gray)
							gray=temp;
						
					}
					
					if(gray>255) *(outputimage+(i-cakeheight/2)*width+j-cakewidth/2)=255;
					else  *(outputimage+(i-cakeheight/2)*width+j-cakewidth/2)=gray;
					
			}

}

//��ֵ�˲��ӳ���
void median(unsigned char *inputbmp,unsigned char *outputbmp)
{    int mask[3][3];
     int maskw=3;
     int maskh=3;
	 int i,j;
	 int px,py;
	 for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			mask[i][j]=1;
	 

	   for(i=maskh/2;i<height-maskh/2;i++)
		   for(j=maskw/2;j<width-maskw/2;j++)
		   {  			
				   int m=8;
				   int a[9];				
				   int t;
				   int p,q;
				   
				   //����������ֵ����ָ������
				   for(py=-maskh/2;py<=maskh/2;py++)
					   for(px=-maskw/2;px<=maskw/2;px++)
					   { 
						   if(m>=0)				
							   a[m]=*(inputbmp+(i+py)*width+j+px);
						   m--;
						   
					   }		   
					   
					   //�������е�������ð������           
					   for(p=0;p<8;p++)
						   for(q=0;q<8-p;q++)
							   if(a[q]>a[q+1])
							   {
								   t=a[q];
								   a[q]=a[q+1];
								   a[q+1]=t;
							   }					   
							   
					  *(outputbmp+i*width+j)=a[4];							   
							   
			   }			   
			   
}


void VideoSegment(FILE*fp,unsigned char*plane_buffer,int tol)
{  
	
 //*******��YUVW�ļ�*****************************************************   
    int ysize=width*height; 
  
    unsigned char *colorbuf;
    unsigned char *graybuf;

    unsigned char *differbuf;

    unsigned char *tempbuf;
    unsigned char *movebuf;
    unsigned char *movecopy;
    unsigned char *edgebuf;
    unsigned char *binarybuf;
    FILE *fpwrite;

    int currframe;
    int r=0;//�ƶ�ָ��ı���
    int i,j;
    int differ1,differ2;
    int r1,r2;
   
   
    int temp;
   
    int maskw=3;
    int maskh=3;
    //ģ���趨
    int block[500];//�ṹԪ�صĴ�С
    int blockh=3;
    int blockw=3;

    int k,r3;
    unsigned char middle[101376];
    unsigned char left[101376];
    unsigned char right[101376];
    int seg_area;
    int l;
	
    BOOL flag;
    BOOL havepoint;
    int beginpoint,endpoint;
    int GX,GY;
    int totalframe=tol; //������Ƶ֡����
    binarybuf=plane_buffer;
    colorbuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);
   
    graybuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);
 
    differbuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);
 
    tempbuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);

    movebuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);
  
    movecopy=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);
  
    edgebuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);
  
    /*binarybuf=(unsigned char*)malloc(sizeof(unsigned char)*ysize*3/2*totalframe);*/

    fseek(fp,0,0);

 
 
    fread(colorbuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
    fread(graybuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
    fread(differbuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
    fread(tempbuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
    fread(movebuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
    fread(edgebuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
    fread(colorbuf,ysize*3/2*totalframe,1,fp);
    fseek(fp,0,0);
 //*********��YUV�ļ����д���******************************************** 

   
   
   //////////////////////////////////////////////////////////
   //��Ϊ�Ҷ�ͼ��
    for(currframe=1;currframe<=120;currframe++)
    {
	 	r1=height*width*(currframe-1)*3/2+height*width;
	    r2=height*width*(currframe-1)*3/2+height*width*3/2;
        //printf("%d %d\n",r1,r2);
	    for(i=0;i<height*width*3/2;i++)
		{	  
	  	    if(r>=r1&&r<=r2)
			{	
			   graybuf[r]=128;
			   differbuf[r]=128;
		 	   tempbuf[r]=128;
			   movebuf[r]=128;
			   movecopy[r]=128;
			   edgebuf[r]=128;
			   binarybuf[r]=128;
		   	}		   
		    r++;
		} 
    }
   ///////////////////////////////////////////////////////////
   //����֡��ͼ����Ľ׾�������˶�����
    
	
	for(i=0;i<100;i++)
	   block[i]=10;
    
  
    for(currframe=1;currframe<=100;currframe++)
	{
		r1=height*width*(currframe-1)*3/2;
	    r2=height*width*(currframe+2)*3/2;
	    r3=height*width*(currframe+5)*3/2;       
	
		
		for(k=0;k<height*width;k++)
		{ 
			differ1=abs(tempbuf[r1+k]-tempbuf[r2+k]);
			differ2=abs(tempbuf[r1+k]-tempbuf[r3+k]);
			if((differ1>=6)||(differ2>=6))
			{ 
				tempbuf[r1+k]=255;
		   	    movebuf[r1+k]=255;
			}
			
			else
			{ 
				tempbuf[r1+k]=0;
			    movebuf[r1+k]=0; 
			}			
			
		}
	   
         median(movebuf+r1,movebuf+r1);
		/*
		r1=height*width*(currframe-1)*3/2;
		r2=height*width*(currframe+2)*3/2;
	

		for(i=0;i<width*height-1;i++)
		{
		   differbuf[r1+i]=abs(graybuf[r2+i]-graybuf[r1+i]);//��������֡��֡��
		}
		  
		temp=0;
		for(i=0;i<12;i++)
		  for(j=0;j<12;j++)
		  {
			  temp=temp+differbuf[r1+i*width+j];		  

		  }
		 average=(int)(temp/(12*12));
		 printf("average=%d\n",average);

		 temp=0;
         for(i=0;i<12;i++)
		   for(j=0;j<12;j++)		   
			   temp=temp+(differbuf[r1+i*width+j]-average)*(differbuf[r1+i*width+j]-average);		  
		 T=temp/(12*12);//��������	
		 printf("%d",T);          

        for(i=1;i<height-1;i++)
		  for(j=1;j<width-1;j++)
		  {   
			  temp=0;
			  for(p=-1;p<=1;p++)
				 for(q=-1;q<=1;q++)
				 {
                     temp=temp+differbuf[r1+(i+p)*width+j+q];					 
				 }
               average=(int)(temp/(maskh*maskw));//������֡���źŵľ�ֵ

			   temp=0; 
			   for(p=-1;p<=1;p++)
				 for(q=-1;q<=1;q++)
				   {
					  temp=temp+(differbuf[r1+(i+p)*width+j+q]-average)*(differbuf[r1+(i+p)*width+j+q]-average)
					          *(differbuf[r1+(i+p)*width+j+q]-average)*(differbuf[r1+(i+p)*width+j+q]-average);
				 }
				average=(int)(temp/(maskh*maskw));//֡��ͼ����Ľ׾�
				
				if(average>=40)
					differbuf[r1+i*width+j]=255;
				else
                    differbuf[r1+i*width+j]=0;
				

		  }
      
		  for(i=0;i<width*height;i++)
			  tempbuf[r1+i]=differbuf[r1+i];         

         //grayinflation(differbuf+r1,&block[0],3,3,tempbuf+r1);
		// grayerosion(tempbuf+r1, &block[0],5,5,differbuf+r1);
         grayinflation(tempbuf+r1,&block[0],3,3,differbuf+r1);
         median(differbuf+r1,differbuf+r1);*/  
	}

///***********************�����շָ������ж�ʱ�����򲹳�******************************
	
		
	for(seg_area=0;seg_area<8;seg_area++)
	{  
			
		for(i=0;i<width*height;i++)
		{ 
			middle[i]=0;//��ǰʱ��θ����س��ִ���ͳ������
			left[i]=0;//ǰһʱ��θ����س��ִ���ͳ������
		    right[i]=0;//��һʱ��θ����س��ִ���ͳ������
		}	
		for(i=0;i<height;i++)
			for(j=0;j<width;j++)			
				for(l=(seg_area*12+1);l<=(seg_area*12+12);l++)
				{ 		 
					if(movebuf[height*width*(l-1)*3/2+i*width+j]==255)  middle[i*width+j]++;				
				}						
				for(currframe=seg_area*12+1;currframe<=seg_area*12+12;currframe++)
				{
					for(i=0;i<height;i++)
						for(j=0;j<width;j++)
						{
							//if(middle[i*width+j]>=4||left[i*width+j]>=4||right[i*width+j]>=4)
							if(middle[i*width+j]>=1)
								movebuf[height*width*(currframe-1)*3/2+i*width+j]=255;
							else
								movebuf[height*width*(currframe-1)*3/2+i*width+j]=0;
						}					  
				}	
	}		

//************************����ĸ��*************************************
    for(i=0;i<height;i++)
        for(j=0;j<width;j++)
		{    
        	temp=0;
	        for(currframe=1;currframe<=100;currframe++)
			{ 
	        	r1=height*width*(currframe-1)*3/2;
	        	if(movebuf[r1+i*width+j]==255)
	        		temp++;			
			}
	
	    if(temp>0)
	    	tempbuf[i*width+j]=255;
		}


    for(currframe=1;currframe<=100;currframe++)
	{ 
    	r1=height*width*(currframe-1)*3/2;
	
    	for(i=0;i<height;i++)
    		for(j=0;j<width;j++)
			{
		    	if(tempbuf[i*width+j]==255)
		    		movebuf[r1+i*width+j]=255;
		    	else
		 	    	movebuf[r1+i*width+j]=0;
			}
        median(movebuf+r1,movebuf+r1);
		
	}
//*******************��Ե�޲�***********************************
for(currframe=1;currframe<=100;currframe++)
{
     r1=height*width*(currframe-1)*3/2;
     for(i=0;i<2;i++)//�޲��ϱ�Ե
        for(j=0;j<width;j++)
          movebuf[r1+i*width+j]=0;
     //for(i=height-1;i<height;i++)//�޲��±�Ե
       //for(j=0;j<width;j++)
         // movebuf[r1+i*width+j]=0;
    for(i=0;i<height;i++)//�޲����Ե
        for(j=0;j<20;j++)
          movebuf[r1+i*width+j]=0;
    for(i=0;i<height;i++)//�޲��ұ�Ե
        for(j=width-1;j<width;j++)
          movebuf[r1+i*width+j]=0;	
}
//*********************���ĸ��***************************************

for(currframe=1;currframe<=100;currframe++)
{   
	r1=height*width*(currframe-1)*3/2;
	for(i=0;i<width*height;i++)
        movecopy[r1+i]=movebuf[r1+i];
}

for(currframe=1;currframe<=100;currframe++)
{
   r1=height*width*(currframe-1)*3/2;   
   
  ////////////�Լ��ģ�����ˮƽ���///////////
  for(i=0;i<height;i++)
  {  
   //��������Ƿ���Ŀ���
   havepoint=FALSE1; 
	for(j=0;j<width;j++)
	{
		if(movebuf[r1+i*width+j]==255) 
		{
			havepoint=TRUE1;
			break;
		}
	}

	//������Ŀ���
    if(havepoint)   
	{ 	
	   //�ҵ�ˮƽ�׸�Ŀ���
       j=0;
	   flag=TRUE1;    
	   while(flag)
	   {
	   	 if(movebuf[r1+i*width+j]==255)
		 { 
			 beginpoint=j;
			 flag=FALSE1;
		 }
		 else
		 {
			 j++;
		 }
	   }
	   //�ҵ�ˮƽĩβĿ���
	   j=width-1;
	   flag=TRUE1;
	   while(flag)
	   { 
		   if(movebuf[r1+i*width+j]==255)
		   { 
		     endpoint=j;
		     flag=FALSE1;
		   }
	      else
		  {
		      j--;
		  }
          
	   }
       //����׸�Ŀ��㵽ĩβĿ���֮�������
       for(j=beginpoint;j<=endpoint;j++)	   
		   movebuf[r1+i*width+j]=255;
	}
  }

 ////////////�Լ��ģ����д�ֱ���/////////////////
 
  for(j=0;j<width;j++)
  {  
	  //��������Ƿ���Ŀ���
	  havepoint=FALSE1; 
	  for(i=0;i<height;i++)
	  {
		  if(movecopy[r1+i*width+j]==255) 
		  {
			  havepoint=TRUE1;
			  break;
		  }
	  }
	  
	  //������Ŀ���
	  if(havepoint)   
	  { 	
		  //�ҵ���ֱ�׸�Ŀ���
		  i=0;
		  flag=TRUE1;    
		  while(flag)
		  {
			  if(movecopy[r1+i*width+j]==255)
			  { 
				  beginpoint=i;
				  flag=FALSE1;
			  }
			  else
			  {
				  i++;
			  }
		  }
		  //ȷ����ֱ��ĩβ��
		  endpoint=height-1;
		 
		  //����׸�Ŀ��㵽ĩβĿ���֮�������
		  for(i=beginpoint;i<=endpoint;i++)	   
			  movecopy[r1+i*width+j]=255;
	  }
  
  }

  /////////////ȡ�����õ����ģ��////////////////////
  for(i=0;i<height;i++)
	  for(j=0;j<width;j++)
	  {
		  if((movebuf[r1+i*width+j]==255)&&(movecopy[r1+i*width+j]==255))
		  {	 		
			  movebuf[r1+i*width+j]=255;				    
		  }
		  else
		  {
			  movebuf[r1+i*width+j]=0;
		  }
	  }

	  //for(i=0;i<width*height;i++)
			  //tempbuf[r1+i]=movebuf[r1+i];     
	  
     // grayinflation(tempbuf+r1, &block[0],5,5,movebuf+r1);	  

}
///***********************sobel���ӱ�Ե���*****************************

for(currframe=1;currframe<=100;currframe++)
{
    r1=height*width*(currframe-1)*3/2;   

	for(i=1;i<height-1;i++)
       for(j=1;j<width-1;j++)
	   {
		   
		   //sobel����
		   GX=graybuf[r1+(i-1)*width+(j+1)]+graybuf[r1+i*width+(j+1)]*2+graybuf[r1+(i+1)*width+(j+1)]
			   -graybuf[r1+(i-1)*width+(j-1)]-graybuf[r1+i*width+(j-1)]*2-graybuf[r1+(i+1)*width+(j-1)];
			   
		   GY=graybuf[r1+(i-1)*width+(j-1)]+graybuf[r1+(i-1)*width+j]*2+graybuf[r1+(i+1)*width+(j+1)]
			   -graybuf[r1+(i+1)*width+(j-1)]-graybuf[r1+(i+1)*width+j]*2-graybuf[r1+(i+1)*width+(j+1)];	   
			   
			   
		   temp=sqrt(GX*GX+GY*GY)+0.5;
		   if(temp>255) temp=255;
		        edgebuf[r1+i*width+j]=temp;

			if(edgebuf[r1+i*width+j]>20)
		        edgebuf[r1+i*width+j]=255;
			else
			    edgebuf[r1+i*width+j]=0;

	   }

}

///********************ʱ���ں���ȡ�˶������Ե*************************

for(currframe=1;currframe<=100;currframe++)
{
    r1=height*width*(currframe-1)*3/2; 
	for(i=0;i<height;i++)
      for(j=0;j<width;j++)
	  {
         if(movebuf[r1+i*width+j]==255&&edgebuf[r1+i*width+j]==255)
			   edgebuf[r1+i*width+j]=255;
		 else
			   edgebuf[r1+i*width+j]=0;

	  }

}
//************************�����ȡ���*********************************

for(currframe=1;currframe<=100;currframe++)
{
	r1=height*width*(currframe-1)*3/2;
     //median(edgebuf+r1,edgebuf+r1);	
	for(i=0;i<width*height;i++)
	{
		tempbuf[r1+i]=edgebuf[r1+i];
		binarybuf[r1+i]=edgebuf[r1+i];
	}


 for(i=0;i<height;i++)
  {  
   //��������Ƿ���Ŀ���
   havepoint=FALSE1; 
	for(j=0;j<width;j++)
	{
		if(edgebuf[r1+i*width+j]==255) 
		{
			havepoint=TRUE1;
			break;
		}
	}

	//������Ŀ���
    if(havepoint)   
	{ 	
	   //�ҵ�ˮƽ�׸�Ŀ���
       j=0;
	   flag=TRUE1;    
	   while(flag)
	   {
	   	 if(edgebuf[r1+i*width+j]==255)
		 { 
			 beginpoint=j;
			 flag=FALSE1;
		 }
		 else
		 {
			 j++;
		 }
	   }
	   //�ҵ�ˮƽĩβĿ���
	   j=width-1;
	   flag=TRUE1;
	   while(flag)
	   { 
		   if(edgebuf[r1+i*width+j]==255)
		   { 
		     endpoint=j;
		     flag=FALSE1;
		   }
	      else
		  {
		      j--;
		  }
          
	   }
       //����׸�Ŀ��㵽ĩβĿ���֮�������
       for(j=beginpoint;j<=endpoint;j++)	   
		   tempbuf[r1+i*width+j]=255;
	}
  }

 ////////////�Լ��ģ����д�ֱ���/////////////////
 
  for(j=0;j<width;j++)
  {  
	  //��������Ƿ���Ŀ���
	  havepoint=FALSE1; 
	  for(i=0;i<height;i++)
	  {
		  if(edgebuf[r1+i*width+j]==255) 
		  {
			  havepoint=TRUE1;
			  break;
		  }
	  }
	  
	  //������Ŀ���
	  if(havepoint)   
	  { 	
		  //�ҵ���ֱ�׸�Ŀ���
		  i=0;
		  flag=TRUE1;    
		  while(flag)
		  {
			  if(edgebuf[r1+i*width+j]==255)
			  { 
				  beginpoint=i;
				  flag=FALSE1;
			  }
			  else
			  {
				  i++;
			  }
		  }
		  //ȷ����ֱ��ĩβ��
		  endpoint=height-1;
		 
		  //����׸�Ŀ��㵽ĩβĿ���֮�������
		  for(i=beginpoint;i<=endpoint;i++)	   
			  binarybuf[r1+i*width+j]=255;
	  }
  
  }

  /////////////ȡ�����õ����ģ��////////////////////
  for(i=0;i<height;i++)
	  for(j=0;j<width;j++)
	  {
		  if((tempbuf[r1+i*width+j]==255)&&(binarybuf[r1+i*width+j]==255))
		  {	 		
			  binarybuf[r1+i*width+j]=255;				    
		  }
		  else
		  {
			  binarybuf[r1+i*width+j]=0;
		  }
	  }

	  for(i=0;i<width*height;i++)
		  tempbuf[r1+i]=binarybuf[r1+i];     
	  
	  //grayerosion(binarybuf+r1, &block[0],7,7,tempbuf+r1); 
	  //grayinflation(tempbuf+r1,&block[10],7,7,binarybuf+r1);
	  median(binarybuf+r1,binarybuf+r1);

}

for(currframe=1;currframe<=100;currframe++)
{
    r1=height*width*3/2*(currframe-1);
	for(i=0;i<height/2;i++)
		for(j=0;j<width/2;j++)
			   {
		   			 binarybuf[r1+width*height+i*width/2+j]= binarybuf[r1+(i*2)*width+j*2]; 
		 	   		 binarybuf[r1+width*height*5/4+i*width/2+j]= binarybuf[r1+(i*2)*width+j*2]; 				  
				}
	   	   

}






//*********����ֵģ��ͶӰ����ɫͼ������********************************
	  
for(currframe=1;currframe<=100;currframe++)
{
		  
		  r1=height*width*3/2*(currframe-1);
		  for(i=0;i<height;i++)
			  for(j=0;j<width;j++)
			  {
				  if(binarybuf[r1+i*width+j]==0)
				  {
					  colorbuf[r1+i*width+j]=255;
					  colorbuf[r1+width*height+b(i)*width/2+b(j)]=128;//��Ϊһ�������ڵ��ĸ���U��V������һ����
					  colorbuf[r1+width*height+width*height/4+b(i)*width/2+b(j)]=128;
					  
				  }
				  
			  }			  
		  			  
 }

  
//*************��YUV�ļ�*************************************************
   fpwrite=fopen("graybuf.yuv","wb");
   fwrite(graybuf,ysize*3/2*90,1,fpwrite);
   fpwrite=fopen("differbuf.yuv","wb");
   fwrite(differbuf,ysize*3/2*90,1,fpwrite);
   fpwrite=fopen("movebuf.yuv","wb");
   fwrite(movebuf,ysize*3/2*90,1,fpwrite);
   fpwrite=fopen("edgebuf.yuv","wb");
   fwrite(edgebuf,ysize*3/2*90,1,fpwrite);
   fpwrite=fopen("binarybuf.yuv","wb");
   fwrite(binarybuf,ysize*3/2*90,1,fpwrite);
   fpwrite=fopen("colorbuf.yuv","wb");
   fwrite(colorbuf,ysize*3/2*90,1,fpwrite);

}