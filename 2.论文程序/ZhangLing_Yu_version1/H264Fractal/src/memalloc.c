
#include "memalloc.h"
#include "huffman.h"

/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////


extern ImageParameters *img;


/*!
************************************************************************
* \brief
*    Allocate 2D memory array -> unsigned char array2D[rows][columns]
*
* \par Output:
*    memory size in bytes
************************************************************************/
// Change 9-Aug-2001 P. List: dont allocate independant row arrays anymore
// but one complete array and move row-pointers to array. Now you can step
// to the next line with an offset of img->width
int get_mem2D(byte ***array2D, int rows, int columns)
{
	int i;
	
	if((*array2D      = (byte**)calloc(rows,        sizeof(byte*))) == NULL)
		no_mem_exit("get_mem2D: array2D");
	if(((*array2D)[0] = (byte* )calloc(columns*rows,sizeof(byte ))) == NULL)
		no_mem_exit("get_mem2D: array2D");
	
	for(i=1;i<rows;i++)
		(*array2D)[i] = (*array2D)[i-1] + columns ;
	
	return rows*columns;
}
/*!
 ************************************************************************
 * \brief
 *    Allocate 2D memory array -> int array2D[rows][columns]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
// same change as in get_mem2Dint
int get_mem2Dint(int ***array2D, int rows, int columns)
{
	int i;
	
	if((*array2D      = (int**)calloc(rows,        sizeof(int*))) == NULL)
		no_mem_exit("get_mem2Dint: array2D");
	if(((*array2D)[0] = (int* )calloc(rows*columns,sizeof(int ))) == NULL)
		no_mem_exit("get_mem2Dint: array2D");
	
	for(i=1 ; i<rows ; i++)
		(*array2D)[i] =  (*array2D)[i-1] + columns  ;
	
	return rows*columns*sizeof(int);
}
int get_mem2Ddouble(double ***array2D, int rows, int columns)
{
	int i;
	
	if((*array2D      = (double**)calloc(rows,        sizeof(double*))) == NULL)
		no_mem_exit("get_mem2D: array2D");
	if(((*array2D)[0] = (double* )calloc(columns*rows,sizeof(double ))) == NULL)
		no_mem_exit("get_mem2D: array2D");
	
	for(i=1;i<rows;i++)
		(*array2D)[i] = (*array2D)[i-1] + columns ;
	
	return rows*columns;
}
/*!
 ************************************************************************
 * \brief
 *    Allocate 2D memory array -> int64 array2D[rows][columns]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
// same change as in get_mem2Dint
int get_mem2Dint64(int64 ***array2D, int rows, int columns)
{
  int i;

  if((*array2D      = (int64**)calloc(rows,        sizeof(int64*))) == NULL)
    no_mem_exit("get_mem2Dint64: array2D");
  if(((*array2D)[0] = (int64* )calloc(rows*columns,sizeof(int64 ))) == NULL)
    no_mem_exit("get_mem2Dint64: array2D");

  for(i=1 ; i<rows ; i++)
    (*array2D)[i] =  (*array2D)[i-1] + columns  ;

  return rows*columns*sizeof(int64);
}


/*!
************************************************************************
* \brief
*    Allocate 3D memory array -> unsigned char array3D[frames][rows][columns]
*
* \par Output:
*    memory size in bytes
************************************************************************
*/
// same change as in get_mem2Dint
int get_mem3D(byte ****array3D, int frames, int rows, int columns)
{
	int  j;
	
	if(((*array3D) = (byte***)calloc(frames,sizeof(byte**))) == NULL)
		no_mem_exit("get_mem3D: array3D");
	
	for(j=0;j<frames;j++)
		get_mem2D( (*array3D)+j, rows, columns ) ;
	
	return frames*rows*columns;
}

/*!
************************************************************************
* \brief
*    Allocate 3D memory array -> int array3D[frames][rows][columns]
*
* \par Output:
*    memory size in bytes
************************************************************************
*/
// same change as in get_mem2Dint
int get_mem3Dint(int ****array3D, int frames, int rows, int columns)
{
	int  j;
	
	if(((*array3D) = (int***)calloc(frames,sizeof(int**))) == NULL)
		no_mem_exit("get_mem3Dint: array3D");
	
	for(j=0;j<frames;j++)
		get_mem2Dint( (*array3D)+j, rows, columns ) ;
	
	return frames*rows*columns*sizeof(int);
}
/*!
************************************************************************
* \brief
*    Allocate 3D memory array -> int64 array3D[frames][rows][columns]
*
* \par Output:
*    memory size in bytes
************************************************************************
*/
// same change as in get_mem2Dint
int get_mem3Dint64(int64 ****array3D, int frames, int rows, int columns)
{
	int  j;
	
	if(((*array3D) = (int64***)calloc(frames,sizeof(int64**))) == NULL)
		no_mem_exit("get_mem3Dint64: array3D");
	
	for(j=0;j<frames;j++)
		get_mem2Dint64( (*array3D)+j, rows, columns ) ;
	
	return frames*rows*columns*sizeof(int64);
}

/*!
************************************************************************
* \brief
*    free 4D memory array 
*    which was alocated with get_mem4Dint()
************************************************************************
*/
int get_mem4Dint(int *****array4D, int idx, int frames, int rows, int columns )
{
	int  j;
	
	if(((*array4D) = (int****)calloc(idx,sizeof(int***))) == NULL)
		no_mem_exit("get_mem4Dint: array4D");
	
	for(j=0;j<idx;j++)
		get_mem3Dint( (*array4D)+j, frames, rows, columns ) ;
	
	return idx*frames*rows*columns*sizeof(int);
}

void free_mem2D(byte **array2D)
{
  if (array2D)
  {
    if (array2D[0])
      free (array2D[0]);

    free (array2D);
  } 
}


/*!
************************************************************************
* \brief
*    free 2D memory array
*    which was alocated with get_mem2Dint()
************************************************************************
*/
void free_mem2Dint(int **array2D)
{
	if (array2D)
	{
		if (array2D[0]) 
			free (array2D[0]);
		else error ("free_mem2D: trying to free unused memory",100);
		
		free (array2D);
		
	} else
	{
		error ("free_mem2D: trying to free unused memory",100);
	}
}

void free_mem2Ddouble(double **array2D)
{
  if (array2D)
  {
    if (array2D[0])
      free (array2D[0]);

    free (array2D);
  } 
}

/*!
************************************************************************
* \brief
*    free 2D memory array
*    which was alocated with get_mem2Dint64()
************************************************************************
*/
void free_mem2Dint64(int64 **array2D)
{
	if (array2D)
	{
		if (array2D[0]) 
			free (array2D[0]);
		else error ("free_mem2Dint64: trying to free unused memory",100);
		
		free (array2D);
		
	} else
	{
		error ("free_mem2Dint64: trying to free unused memory",100);
	}
}

/*!
************************************************************************
* \brief
*    free 3D memory array
*    which was alocated with get_mem3D()
************************************************************************
*/
void free_mem3D(byte ***array3D, int frames)
{
	int i;
	
	if (array3D)
	{
		for (i=0;i<frames;i++)
		{ 
			free_mem2D(array3D[i]);
		}
		free (array3D);
	} else
	{
		error ("free_mem3D: trying to free unused memory",100);
	}
}

/*!
************************************************************************
* \brief
*    free 3D memory array 
*    which was alocated with get_mem3Dint()
************************************************************************
*/
void free_mem3Dint(int ***array3D, int frames)
{
	int i;
	
	if (array3D)
	{
		for (i=0;i<frames;i++)
		{ 
			free_mem2Dint(array3D[i]);
		}
		free (array3D);
	} else
	{
		error ("free_mem3D: trying to free unused memory",100);
	}
}


/*!
 ************************************************************************
 * \brief
 *    free 3D memory array 
 *    which was alocated with get_mem3Dint64()
 ************************************************************************
 */
void free_mem3Dint64(int64 ***array3D, int frames)
{
  int i;

  if (array3D)
  {
    for (i=0;i<frames;i++)
    { 
      free_mem2Dint64(array3D[i]);
    }
   free (array3D);
  } else
  {
    error ("free_mem3Dint64: trying to free unused memory",100);
  }
}

/*!
************************************************************************
* \brief
*    free 4D memory array 
*    which was alocated with get_mem4Dint()
************************************************************************
*/
void free_mem4Dint(int ****array4D, int idx, int frames )
{
	int  j;
	
	if (array4D)
	{
		for(j=0;j<idx;j++)
			free_mem3Dint( array4D[j], frames) ;
		free (array4D);
	} else
	{
		error ("free_mem4D: trying to free unused memory",100);
	}
}


// int get_mem5Dint(int ******array5D, int refs, int blocktype, int rows, int columns, int component)
// {
// 	int  j;
// 	
// 	if(((*array5D) = (int*****)calloc(refs,sizeof(int****))) == NULL)
// 		no_mem_exit("get_mem5Dint: array5D");
// 	
// 	;
// 	for(j=0;j<refs;j++)
// 		get_mem4Dint( (*array5D)+j, blocktype, rows, columns, component) ;
// 	
// 	return refs*blocktype*rows*columns*component*sizeof(int);
// }
// 
/////+++++++++++++++++++++++++Zhang Ling+++++++++++++++++++++++++++++/////



void no_mem_exit(char *where)
{
	printf("Could not allocate memory: %s",where);
	exit(0);
}
void Memoalloc()
{
	int width,height,i,j;
	int trans_count,next_count;
	width = imageWidth;
	height = imageHeight;

	/////三目统一使用的变量
	imgUV_org = (byte***)calloc(2,sizeof(byte**));
    imgUV_rec = (byte***)calloc(2,sizeof(byte**));
	imgUV_ref = (byte***)calloc(2,sizeof(byte**));
	imgUV_ref_c = (byte***)calloc(2,sizeof(byte**));//////////////////////////////////////////////////////////////////////////

	imgUV_ref_r = (byte***)calloc(2,sizeof(byte**));
	imgUV_ref_l = (byte***)calloc(2,sizeof(byte**));
	imgUV_ref_h = (byte***)calloc(2,sizeof(byte**));
	imgUV_ref_m = (byte***)calloc(2,sizeof(byte**));
	imgUV_ref_n = (byte***)calloc(2,sizeof(byte**));
	get_mem2D(&imgY_org,height,width);
	get_mem2D(&imgUV_org[0],height/2,width/2);
	get_mem2D(&imgUV_org[1],height/2,width/2);

	get_mem2D(&imgY_ref,height,width);//////////////////////////////////////////////////////////////////////////
	get_mem2D(&imgUV_ref[0],height/2,width/2);//////////////////////////////////////////////////////////////////////////
	get_mem2D(&imgUV_ref[1],height/2,width/2);

// 	get_mem2D(&imgY_ref_c,height,width);//////////////////////////////////////////////////////////////////////////
	get_mem2D(&imgUV_ref_c[0],height/2,width/2);//////////////////////////////////////////////////////////////////////////
	get_mem2D(&imgUV_ref_c[1],height/2,width/2);
	get_mem2D(&imgY_reff,height,width);

	get_mem2D(&imgY_ref_r,height,width);
	get_mem2D(&imgUV_ref_r[0],height/2,width/2);
	get_mem2D(&imgUV_ref_r[1],height/2,width/2);

	get_mem2D(&imgY_ref_l,height,width);
	get_mem2D(&imgUV_ref_l[0],height/2,width/2);
	get_mem2D(&imgUV_ref_l[1],height/2,width/2);
    
	get_mem2D(&imgY_ref_h,height,width);
 	get_mem2D(&imgUV_ref_h[0],height/2,width/2);
 	get_mem2D(&imgUV_ref_h[1],height/2,width/2);

	get_mem2D(&imgY_ref_m,height,width);
	get_mem2D(&imgUV_ref_m[0],height/2,width/2);
 	get_mem2D(&imgUV_ref_m[1],height/2,width/2);

	get_mem2D(&imgY_ref_n,height,width);
	get_mem2D(&imgUV_ref_n[0],height/2,width/2);
 	get_mem2D(&imgUV_ref_n[1],height/2,width/2);

	get_mem2D(&imgY_rec,height,width);
	get_mem2D(&imgUV_rec[0],height/2,width/2);
	get_mem2D(&imgUV_rec[1],height/2,width/2);

    sourceframe.sourceframe_Y = (byte*)calloc(width*height,sizeof(byte));
    sourceframe.sourceframe_U = (byte*)calloc(width*height/4,sizeof(byte));
	sourceframe.sourceframe_V = (byte*)calloc(width*height/4,sizeof(byte));

	if(input->num_regions>1)
	{
		sourceframe.plane_Y = (byte*)calloc(width*height,sizeof(byte));
		sourceframe.plane_U = (byte*)calloc(width*height/4,sizeof(byte));
		sourceframe.plane_V = (byte*)calloc(width*height/4,sizeof(byte));

		sourceframe.plane_Y_domain = (byte*)calloc(width*height,sizeof(byte));
		sourceframe.plane_U_domain = (byte*)calloc(width*height/4,sizeof(byte));
		sourceframe.plane_V_domain = (byte*)calloc(width*height/4,sizeof(byte));

		plane_UV = (byte***)calloc(2,sizeof(byte**));
		plane_UV_domain = (byte***)calloc(2,sizeof(byte**));
		get_mem2D(&plane_Y,height,width);
		get_mem2D(&plane_UV[0],height/2,width/2);
		get_mem2D(&plane_UV[1],height/2,width/2);
		
		get_mem2D(&plane_Y_domain,height,width);
		get_mem2D(&plane_UV_domain[0],height/2,width/2);
		get_mem2D(&plane_UV_domain[1],height/2,width/2);
		
		if (input->right==1)
		{
			plane_UV_domain_r = (byte***)calloc(2,sizeof(byte**));
			get_mem2D(&plane_Y_domain_r,height,width);
			get_mem2D(&plane_UV_domain_r[0],height/2,width/2);
			get_mem2D(&plane_UV_domain_r[1],height/2,width/2);
		}
		if (input->left==1)
		{
			plane_UV_domain_l = (byte***)calloc(2,sizeof(byte**));
			get_mem2D(&plane_Y_domain_l,height,width);
			get_mem2D(&plane_UV_domain_l[0],height/2,width/2);
			get_mem2D(&plane_UV_domain_l[1],height/2,width/2);
		}

		imgY_rec_region=(byte***)malloc(sizeof(byte**)*2);
		for(j=0;j<2;j++)
		{
			get_mem2D(&imgY_rec_region[j],height,width);
			//get_mem2D(&imgY_rec_region[j],height,width);
		}
		imgUV_rec_region=(byte****)malloc(sizeof(byte***)*2);  //imgUV_rec_region[con][num_region][y][x]
		for(i=0;i<2;i++)
		{
			imgUV_rec_region[i]=(byte***)calloc(2,sizeof(byte**));
			for(j=0;j<2;j++)
			{
				get_mem2D(&imgUV_rec_region[i][j],height/2,width/2);
//				get_mem2D(&imgUV_rec_region[i][j],height/2,width/2);
			}
		}
	}	
//range sum,16*16,8*8,4*4
	get_mem2Ddouble(&sum_16_org,img->frmHeightInMbs,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_8_org,img->frmHeightInMbs*2,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum_4_org,img->frmHeightInMbs*4,img->frmWidthInMbs*4);
	get_mem2Ddouble(&sum2_16_org,img->frmHeightInMbs,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_8_org,img->frmHeightInMbs*2,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum2_4_org,img->frmHeightInMbs*4,img->frmWidthInMbs*4);

	get_mem2Ddouble(&sum_16_U_org,img->frmHeightInMbs/2,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum_8_U_org,img->frmHeightInMbs,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_4_U_org,img->frmHeightInMbs*2,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum2_16_U_org,img->frmHeightInMbs/2,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum2_8_U_org,img->frmHeightInMbs,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_4_U_org,img->frmHeightInMbs*2,img->frmWidthInMbs*2);

	get_mem2Ddouble(&sum_16_V_org,img->frmHeightInMbs/2,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum_8_V_org,img->frmHeightInMbs,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_4_V_org,img->frmHeightInMbs*2,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum2_16_V_org,img->frmHeightInMbs/2,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum2_8_V_org,img->frmHeightInMbs,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_4_V_org,img->frmHeightInMbs*2,img->frmWidthInMbs*2);
//     16*8,8*16,8*4,4*8
	get_mem2Ddouble(&sum_16_8_org,img->frmHeightInMbs*2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_8_16_org,img->frmHeightInMbs,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum_8_4_org,img->frmHeightInMbs*4,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum_4_8_org,img->frmHeightInMbs*2,img->frmWidthInMbs*4);

	get_mem2Ddouble(&sum_16_8_U_org,img->frmHeightInMbs,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum_8_16_U_org,img->frmHeightInMbs/2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_8_4_U_org,img->frmHeightInMbs*2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_4_8_U_org,img->frmHeightInMbs,img->frmWidthInMbs*2);

	get_mem2Ddouble(&sum_16_8_V_org,img->frmHeightInMbs,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum_8_16_V_org,img->frmHeightInMbs/2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_8_4_V_org,img->frmHeightInMbs*2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum_4_8_V_org,img->frmHeightInMbs,img->frmWidthInMbs*2);


	get_mem2Ddouble(&sum2_16_8_org,img->frmHeightInMbs*2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_8_16_org,img->frmHeightInMbs,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum2_8_4_org,img->frmHeightInMbs*4,img->frmWidthInMbs*2);
	get_mem2Ddouble(&sum2_4_8_org,img->frmHeightInMbs*2,img->frmWidthInMbs*4);

	get_mem2Ddouble(&sum2_16_8_U_org,img->frmHeightInMbs,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum2_8_16_U_org,img->frmHeightInMbs/2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_8_4_U_org,img->frmHeightInMbs*2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_4_8_U_org,img->frmHeightInMbs,img->frmWidthInMbs*2);

	get_mem2Ddouble(&sum2_16_8_V_org,img->frmHeightInMbs,img->frmWidthInMbs/2);
	get_mem2Ddouble(&sum2_8_16_V_org,img->frmHeightInMbs/2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_8_4_V_org,img->frmHeightInMbs*2,img->frmWidthInMbs);
	get_mem2Ddouble(&sum2_4_8_V_org,img->frmHeightInMbs,img->frmWidthInMbs*2);

//domain sum,16*16,8*8,4*4
	
	
///////////////domain sum

	get_mem2Ddouble(&sum_16_ref,input->imageheight-15,input->imagewidth-15);

	get_mem2Ddouble(&sum_8_ref,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_ref,input->imageheight-3,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_ref,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_ref,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_ref,input->imageheight-3,input->imagewidth-3);

	get_mem2Ddouble(&sum_16_U_ref,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_U_ref,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_U_ref,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_U_ref,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_U_ref,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_U_ref,input->imageheight/2-3,input->imagewidth/2-3);

	get_mem2Ddouble(&sum_16_V_ref,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_V_ref,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_V_ref,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_V_ref,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_V_ref,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_V_ref,input->imageheight/2-3,input->imagewidth/2-3);

	get_mem2Ddouble(&sum_16_8_ref,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_16_ref,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum_8_4_ref,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_8_ref,input->imageheight-7,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_8_ref,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_16_ref,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum2_8_4_ref,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_8_ref,input->imageheight-7,input->imagewidth-3);

	get_mem2Ddouble(&sum_16_8_U_ref,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_U_ref,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_U_ref,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_U_ref,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_U_ref,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_U_ref,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_U_ref,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_U_ref,input->imageheight/2-7,input->imagewidth/2-3);

	get_mem2Ddouble(&sum_16_8_V_ref,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_V_ref,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_V_ref,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_V_ref,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_V_ref,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_V_ref,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_V_ref,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_V_ref,input->imageheight/2-7,input->imagewidth/2-3);

    //////////////////////////
	get_mem2Ddouble(&sum_16_ref_H,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_ref_H,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_ref_H,input->imageheight-3,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_ref_H,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_ref_H,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_ref_H,input->imageheight-3,input->imagewidth-3);
	
	get_mem2Ddouble(&sum_16_U_ref_H,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_U_ref_H,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_U_ref_H,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_U_ref_H,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_U_ref_H,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_U_ref_H,input->imageheight/2-3,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_V_ref_H,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_V_ref_H,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_V_ref_H,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_V_ref_H,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_V_ref_H,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_V_ref_H,input->imageheight/2-3,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_8_ref_H,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_16_ref_H,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum_8_4_ref_H,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_8_ref_H,input->imageheight-7,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_8_ref_H,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_16_ref_H,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum2_8_4_ref_H,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_8_ref_H,input->imageheight-7,input->imagewidth-3);
	
	get_mem2Ddouble(&sum_16_8_U_ref_H,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_U_ref_H,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_U_ref_H,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_U_ref_H,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_U_ref_H,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_U_ref_H,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_U_ref_H,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_U_ref_H,input->imageheight/2-7,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_8_V_ref_H,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_V_ref_H,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_V_ref_H,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_V_ref_H,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_V_ref_H,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_V_ref_H,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_V_ref_H,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_V_ref_H,input->imageheight/2-7,input->imagewidth/2-3);
     ////////////////////////////
	get_mem2Ddouble(&sum_16_ref_M,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_ref_M,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_ref_M,input->imageheight-3,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_ref_M,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_ref_M,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_ref_M,input->imageheight-3,input->imagewidth-3);
	
	get_mem2Ddouble(&sum_16_U_ref_M,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_U_ref_M,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_U_ref_M,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_U_ref_M,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_U_ref_M,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_U_ref_M,input->imageheight/2-3,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_V_ref_M,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_V_ref_M,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_V_ref_M,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_V_ref_M,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_V_ref_M,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_V_ref_M,input->imageheight/2-3,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_8_ref_M,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_16_ref_M,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum_8_4_ref_M,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_8_ref_M,input->imageheight-7,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_8_ref_M,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_16_ref_M,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum2_8_4_ref_M,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_8_ref_M,input->imageheight-7,input->imagewidth-3);
	
	get_mem2Ddouble(&sum_16_8_U_ref_M,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_U_ref_M,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_U_ref_M,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_U_ref_M,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_U_ref_M,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_U_ref_M,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_U_ref_M,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_U_ref_M,input->imageheight/2-7,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_8_V_ref_M,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_V_ref_M,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_V_ref_M,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_V_ref_M,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_V_ref_M,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_V_ref_M,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_V_ref_M,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_V_ref_M,input->imageheight/2-7,input->imagewidth/2-3);
	////////////////////////
	get_mem2Ddouble(&sum_16_ref_N,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_ref_N,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_ref_N,input->imageheight-3,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_ref_N,input->imageheight-15,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_ref_N,input->imageheight-7,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_ref_N,input->imageheight-3,input->imagewidth-3);
	
	get_mem2Ddouble(&sum_16_U_ref_N,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_U_ref_N,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_U_ref_N,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_U_ref_N,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_U_ref_N,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_U_ref_N,input->imageheight/2-3,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_V_ref_N,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_V_ref_N,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_V_ref_N,input->imageheight/2-3,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_V_ref_N,input->imageheight/2-15,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_V_ref_N,input->imageheight/2-7,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_V_ref_N,input->imageheight/2-3,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_8_ref_N,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum_8_16_ref_N,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum_8_4_ref_N,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum_4_8_ref_N,input->imageheight-7,input->imagewidth-3);
	get_mem2Ddouble(&sum2_16_8_ref_N,input->imageheight-7,input->imagewidth-15);
	get_mem2Ddouble(&sum2_8_16_ref_N,input->imageheight-15,input->imagewidth-7);
	get_mem2Ddouble(&sum2_8_4_ref_N,input->imageheight-3,input->imagewidth-7);
	get_mem2Ddouble(&sum2_4_8_ref_N,input->imageheight-7,input->imagewidth-3);
	
	get_mem2Ddouble(&sum_16_8_U_ref_N,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_U_ref_N,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_U_ref_N,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_U_ref_N,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_U_ref_N,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_U_ref_N,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_U_ref_N,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_U_ref_N,input->imageheight/2-7,input->imagewidth/2-3);
	
	get_mem2Ddouble(&sum_16_8_V_ref_N,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum_8_16_V_ref_N,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_8_4_V_ref_N,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum_4_8_V_ref_N,input->imageheight/2-7,input->imagewidth/2-3);
	get_mem2Ddouble(&sum2_16_8_V_ref_N,input->imageheight/2-7,input->imagewidth/2-15);
	get_mem2Ddouble(&sum2_8_16_V_ref_N,input->imageheight/2-15,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_8_4_V_ref_N,input->imageheight/2-3,input->imagewidth/2-7);
	get_mem2Ddouble(&sum2_4_8_V_ref_N,input->imageheight/2-7,input->imagewidth/2-3);
	if (input->right==1)
	{
		get_mem2Ddouble(&sum_16_ref_r,input->imageheight-15,input->imagewidth-15);
		get_mem2Ddouble(&sum_8_ref_r,input->imageheight-7,input->imagewidth-7);
		get_mem2Ddouble(&sum_4_ref_r,input->imageheight-3,input->imagewidth-3);
		get_mem2Ddouble(&sum2_16_ref_r,input->imageheight-15,input->imagewidth-15);
		get_mem2Ddouble(&sum2_8_ref_r,input->imageheight-7,input->imagewidth-7);
		get_mem2Ddouble(&sum2_4_ref_r,input->imageheight-3,input->imagewidth-3);
		
		get_mem2Ddouble(&sum_16_U_ref_r,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_U_ref_r,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_U_ref_r,input->imageheight/2-3,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_U_ref_r,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_U_ref_r,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_U_ref_r,input->imageheight/2-3,input->imagewidth/2-3);
		
		get_mem2Ddouble(&sum_16_V_ref_r,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_V_ref_r,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_V_ref_r,input->imageheight/2-3,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_V_ref_r,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_V_ref_r,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_V_ref_r,input->imageheight/2-3,input->imagewidth/2-3);
		
		get_mem2Ddouble(&sum_16_8_ref_r,input->imageheight-7,input->imagewidth-15);
		get_mem2Ddouble(&sum_8_16_ref_r,input->imageheight-15,input->imagewidth-7);
		get_mem2Ddouble(&sum_8_4_ref_r,input->imageheight-3,input->imagewidth-7);
		get_mem2Ddouble(&sum_4_8_ref_r,input->imageheight-7,input->imagewidth-3);
		get_mem2Ddouble(&sum2_16_8_ref_r,input->imageheight-7,input->imagewidth-15);
		get_mem2Ddouble(&sum2_8_16_ref_r,input->imageheight-15,input->imagewidth-7);
		get_mem2Ddouble(&sum2_8_4_ref_r,input->imageheight-3,input->imagewidth-7);
		get_mem2Ddouble(&sum2_4_8_ref_r,input->imageheight-7,input->imagewidth-3);
		
		get_mem2Ddouble(&sum_16_8_U_ref_r,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_16_U_ref_r,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_8_4_U_ref_r,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_8_U_ref_r,input->imageheight/2-7,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_8_U_ref_r,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_16_U_ref_r,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_8_4_U_ref_r,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_8_U_ref_r,input->imageheight/2-7,input->imagewidth/2-3);
		
		get_mem2Ddouble(&sum_16_8_V_ref_r,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_16_V_ref_r,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_8_4_V_ref_r,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_8_V_ref_r,input->imageheight/2-7,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_8_V_ref_r,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_16_V_ref_r,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_8_4_V_ref_r,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_8_V_ref_r,input->imageheight/2-7,input->imagewidth/2-3);
	}

	if (input->left==1)
	{
		get_mem2Ddouble(&sum_16_ref_l,input->imageheight-15,input->imagewidth-15);
		get_mem2Ddouble(&sum_8_ref_l,input->imageheight-7,input->imagewidth-7);
		get_mem2Ddouble(&sum_4_ref_l,input->imageheight-3,input->imagewidth-3);
		get_mem2Ddouble(&sum2_16_ref_l,input->imageheight-15,input->imagewidth-15);
		get_mem2Ddouble(&sum2_8_ref_l,input->imageheight-7,input->imagewidth-7);
		get_mem2Ddouble(&sum2_4_ref_l,input->imageheight-3,input->imagewidth-3);
		
		get_mem2Ddouble(&sum_16_U_ref_l,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_U_ref_l,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_U_ref_l,input->imageheight/2-3,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_U_ref_l,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_U_ref_l,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_U_ref_l,input->imageheight/2-3,input->imagewidth/2-3);
		
		get_mem2Ddouble(&sum_16_V_ref_l,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_V_ref_l,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_V_ref_l,input->imageheight/2-3,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_V_ref_l,input->imageheight/2-15,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_V_ref_l,input->imageheight/2-7,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_V_ref_l,input->imageheight/2-3,input->imagewidth/2-3);
		
		get_mem2Ddouble(&sum_16_8_ref_l,input->imageheight-7,input->imagewidth-15);
		get_mem2Ddouble(&sum_8_16_ref_l,input->imageheight-15,input->imagewidth-7);
		get_mem2Ddouble(&sum_8_4_ref_l,input->imageheight-3,input->imagewidth-7);
		get_mem2Ddouble(&sum_4_8_ref_l,input->imageheight-7,input->imagewidth-3);
		get_mem2Ddouble(&sum2_16_8_ref_l,input->imageheight-7,input->imagewidth-15);
		get_mem2Ddouble(&sum2_8_16_ref_l,input->imageheight-15,input->imagewidth-7);
		get_mem2Ddouble(&sum2_8_4_ref_l,input->imageheight-3,input->imagewidth-7);
		get_mem2Ddouble(&sum2_4_8_ref_l,input->imageheight-7,input->imagewidth-3);
		
		get_mem2Ddouble(&sum_16_8_U_ref_l,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_16_U_ref_l,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_8_4_U_ref_l,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_8_U_ref_l,input->imageheight/2-7,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_8_U_ref_l,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_16_U_ref_l,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_8_4_U_ref_l,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_8_U_ref_l,input->imageheight/2-7,input->imagewidth/2-3);
		
		get_mem2Ddouble(&sum_16_8_V_ref_l,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum_8_16_V_ref_l,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_8_4_V_ref_l,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum_4_8_V_ref_l,input->imageheight/2-7,input->imagewidth/2-3);
		get_mem2Ddouble(&sum2_16_8_V_ref_l,input->imageheight/2-7,input->imagewidth/2-15);
		get_mem2Ddouble(&sum2_8_16_V_ref_l,input->imageheight/2-15,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_8_4_V_ref_l,input->imageheight/2-3,input->imagewidth/2-7);
		get_mem2Ddouble(&sum2_4_8_V_ref_l,input->imageheight/2-7,input->imagewidth/2-3);
	}
	
	//C
	trans_Y=(struct transformation_node**)malloc(sizeof(struct transformation_node*)*input->num_regions);
	trans_U=(struct transformation_node**)malloc(sizeof(struct transformation_node*)*input->num_regions);
	trans_V=(struct transformation_node**)malloc(sizeof(struct transformation_node*)*input->num_regions);

	for(i=0;i<input->num_regions;i++)
	{
		trans_Y[i]=(struct transformation_node*)malloc(sizeof(struct transformation_node)*img->frmSizeInMbs);
		for (trans_count=0;trans_count<img->frmSizeInMbs;trans_count++)
		{
			trans_Y[i][trans_count].next = 
				(struct transformation_node*)malloc(sizeof(struct transformation_node)*4);//分配四个trans的内存
			for (next_count=0;next_count<4;next_count++)
			{
				trans_Y[i][trans_count].next[next_count].next=
					(struct transformation_node*)malloc(sizeof(struct transformation_node)*4);
			}
		}
		trans_U[i]=(struct transformation_node*)malloc(sizeof(struct transformation_node)*img->frmSizeInMbs/4);
		for (trans_count=0;trans_count<img->frmSizeInMbs/4;trans_count++)
		{
			trans_U[i][trans_count].next = 
				(struct transformation_node*)malloc(sizeof(struct transformation_node)*4);//分配四个trans的内存
			for (next_count=0;next_count<4;next_count++)
			{
				trans_U[i][trans_count].next[next_count].next=
					(struct transformation_node*)malloc(sizeof(struct transformation_node)*4);
			}
		}
		trans_V[i]=(struct transformation_node*)malloc(sizeof(struct transformation_node)*img->frmSizeInMbs/4);
		for (trans_count=0;trans_count<img->frmSizeInMbs/4;trans_count++)
		{
			trans_V[i][trans_count].next = 
				(struct transformation_node*)malloc(sizeof(struct transformation_node)*4);//分配四个trans的内存
			for (next_count=0;next_count<4;next_count++)
			{
				trans_V[i][trans_count].next[next_count].next=
					(struct transformation_node*)malloc(sizeof(struct transformation_node)*4);
			}
		}
	}
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 	for(j=0;j<10;j++)
// 		for (i=0;i<input->num_regions;i++)
// 		{
// 			code[j][i] = (unsigned char *)malloc(width*height*MAX_BITS/(8*4)*sizeof(unsigned char));
// 			code_ptr[j][i]= code[j][i];
// 		}
// 
// 	x_trans = (int**)malloc(input->num_regions*sizeof(int*));	
// 	y_trans = (int**)malloc(input->num_regions*sizeof(int*));
// 	s_trans = (int**)malloc(input->num_regions*sizeof(int*));
// 	o_trans = (int**)malloc(input->num_regions*sizeof(int*));
// 
// 	redual_transle = (int**)malloc(input->num_regions*sizeof(int*));	
// 	redual_transru = (int**)malloc(input->num_regions*sizeof(int*));	
// 	
// 	redual_transle_u = (int**)malloc(input->num_regions*sizeof(int*));	
// 	redual_transru_u = (int**)malloc(input->num_regions*sizeof(int*));	
// 	redual_transle_v = (int**)malloc(input->num_regions*sizeof(int*));	
// 	redual_transru_v = (int**)malloc(input->num_regions*sizeof(int*));	
// 	for (i=0;i<input->num_regions;i++)
// 	{
// 		x_trans[i]=(int*)malloc(width*height/16*sizeof(int));
// 		y_trans[i]=(int*)malloc(width*height/16*sizeof(int));
// 		s_trans[i]=(int*)malloc(width*height/16*sizeof(int));
// 		o_trans[i]=(int*)malloc(width*height/16*sizeof(int));
// 
// 		redual_transle[i]=(int*)malloc(width*height/16*sizeof(int));
// 		redual_transru[i]=(int*)malloc(width*height/16*sizeof(int));
// 	
// 		redual_transle_u[i]=(int*)malloc(width*height/32*sizeof(int));
// 		redual_transru_u[i]=(int*)malloc(width*height/32*sizeof(int));
// 		redual_transle_v[i]=(int*)malloc(width*height/32*sizeof(int));
// 		redual_transru_v[i]=(int*)malloc(width*height/32*sizeof(int));
// 
// 	}
// 
// 	Huff_ptr = (HCODE***)malloc(sizeof(HCODE**)*input->num_regions);
// 	for (i=0;i<input->num_regions;i++)
// 	{
// 		Huff_ptr[i] = (HCODE**)malloc(sizeof(HCODE*)*9);
// 		for (j=0;j<9;j++)
// 		{
// 			Huff_ptr[i][j]=(HCODE*)malloc(sizeof(HCODE)*(1<<MAX_BITS));
// 		}
// 	}
// 
// 	xy_p=(long**)malloc(sizeof(long*)*input->num_regions);
// 	s_p=(long**)malloc(sizeof(long*)*input->num_regions);
// 	o_p=(long**)malloc(sizeof(long*)*input->num_regions);
////////////////////////////
// 	redual_ple=(long**)malloc(sizeof(long*)*input->num_regions);
// 	redual_pru=(long**)malloc(sizeof(long*)*input->num_regions);
// 
// 	redual_ple_u=(long**)malloc(sizeof(long*)*input->num_regions);
// 	redual_pru_u=(long**)malloc(sizeof(long*)*input->num_regions);
// 	redual_ple_v=(long**)malloc(sizeof(long*)*input->num_regions);
// 	redual_pru_v=(long**)malloc(sizeof(long*)*input->num_regions);
// 
// 
// 
// ////////////////////////////
// 	for (i=0;i<input->num_regions;i++)
// 	{
// 		xy_p[i]=(long*)malloc(sizeof(long)*(1<<huff_search_range));
// 		s_p[i]=(long*)malloc(sizeof(long)*((int)((MAX_ALPHA-MIN_ALPHA)*100)/5+1));
// 		o_p[i]=(long*)malloc(sizeof(long)*((MAX_BETA-MIN_BETA)/5+1));	
// 
// 		redual_ple[i]=(long*)malloc(sizeof(long)*(1<<10));//N_ple
// 		redual_pru[i]=(long*)malloc(sizeof(long)*(1<<10));//N_pru
// 
// 		redual_ple_u[i]=(long*)malloc(sizeof(long)*(1<<10));//N_ple_u
// 		redual_pru_u[i]=(long*)malloc(sizeof(long)*(1<<10));//N_pru_u
// 		redual_ple_v[i]=(long*)malloc(sizeof(long)*(1<<10));//N_ple_v
// 		redual_pru_v[i]=(long*)malloc(sizeof(long)*(1<<10));//N_pru_v
// 	}
}

void Memofree()
{
	int i,j;

	int trans_count,next_count;

	free_mem2D(imgY_org);
	free_mem2D(imgUV_org[0]);
	free_mem2D(imgUV_org[1]);

	free_mem2D(imgY_rec);
	free_mem2D(imgUV_rec[0]);
	free_mem2D(imgUV_rec[1]);

    free_mem2D(imgY_reff);

	free_mem2D(imgY_ref);
	free_mem2D(imgUV_ref[0]);
	free_mem2D(imgUV_ref[1]);
	free(imgUV_org);
	free(imgUV_rec);
	free(imgUV_ref);

	free(sourceframe.sourceframe_Y);
	free(sourceframe.sourceframe_U);
	free(sourceframe.sourceframe_V);
	if(input->num_regions>1)
	{
		free(sourceframe.plane_Y);
		free(sourceframe.plane_U);
		free(sourceframe.plane_V);
		free(sourceframe.plane_Y_domain);
		free(sourceframe.plane_U_domain);
		free(sourceframe.plane_V_domain);

		free_mem2D(plane_Y);
		free_mem2D(plane_UV[0]);
		free_mem2D(plane_UV[1]);
		free(plane_UV);

		free_mem2D(plane_Y_domain);
		free_mem2D(plane_UV_domain[0]);
		free_mem2D(plane_UV_domain[1]);
		free(plane_UV_domain);

		if (input->right==1)
		{
			free_mem2D(plane_Y_domain_r);
			free_mem2D(plane_UV_domain_r[0]);
			free_mem2D(plane_UV_domain_r[1]);
			free(plane_UV_domain_r);
		}
		if (input->left==1)
		{
			free_mem2D(plane_Y_domain_l);
			free_mem2D(plane_UV_domain_l[0]);
			free_mem2D(plane_UV_domain_l[1]);
			free(plane_UV_domain_l);
		}

		for(j=0;j<2;j++)
		{
			free_mem2D(imgY_rec_region[j]);
		}
		free(imgY_rec_region);

		for(i=0;i<2;i++)
		{	
			for(j=0;j<2;j++)
			{
				free_mem2D(imgUV_rec_region[i][j]);
			}
			free(imgUV_rec_region[i]);
		}
		free(imgUV_rec_region);
	}	

	free_mem2Ddouble(sum_16_org);
	free_mem2Ddouble(sum_8_org);
	free_mem2Ddouble(sum_4_org);
	free_mem2Ddouble(sum2_16_org);
	free_mem2Ddouble(sum2_8_org);
	free_mem2Ddouble(sum2_4_org);

	free_mem2Ddouble(sum_16_U_org);
	free_mem2Ddouble(sum_8_U_org);
	free_mem2Ddouble(sum_4_U_org);
	free_mem2Ddouble(sum2_16_U_org);
	free_mem2Ddouble(sum2_8_U_org);
	free_mem2Ddouble(sum2_4_U_org);

	free_mem2Ddouble(sum_16_V_org);
	free_mem2Ddouble(sum_8_V_org);
	free_mem2Ddouble(sum_4_V_org);
	free_mem2Ddouble(sum2_16_V_org);
	free_mem2Ddouble(sum2_8_V_org);
	free_mem2Ddouble(sum2_4_V_org);
//     16*8,8*16,8*4,4*8
	free_mem2Ddouble(sum_16_8_org);
	free_mem2Ddouble(sum_8_16_org);
	free_mem2Ddouble(sum_8_4_org);
	free_mem2Ddouble(sum_4_8_org);

	free_mem2Ddouble(sum_16_8_U_org);
	free_mem2Ddouble(sum_8_16_U_org);
	free_mem2Ddouble(sum_8_4_U_org);
	free_mem2Ddouble(sum_4_8_U_org);

	free_mem2Ddouble(sum_16_8_V_org);
	free_mem2Ddouble(sum_8_16_V_org);
	free_mem2Ddouble(sum_8_4_V_org);
	free_mem2Ddouble(sum_4_8_V_org);


	free_mem2Ddouble(sum2_16_8_org);
	free_mem2Ddouble(sum2_8_16_org);
	free_mem2Ddouble(sum2_8_4_org);
	free_mem2Ddouble(sum2_4_8_org);

	free_mem2Ddouble(sum2_16_8_U_org);
	free_mem2Ddouble(sum2_8_16_U_org);
	free_mem2Ddouble(sum2_8_4_U_org);
	free_mem2Ddouble(sum2_4_8_U_org);

	free_mem2Ddouble(sum2_16_8_V_org);
	free_mem2Ddouble(sum2_8_16_V_org);
	free_mem2Ddouble(sum2_8_4_V_org);
	free_mem2Ddouble(sum2_4_8_V_org);

//domain sum,16*16,8*8,4*4
	
	
///////////////domain sum
	free_mem2Ddouble(sum_16_ref);
	free_mem2Ddouble(sum_8_ref);
	free_mem2Ddouble(sum_4_ref);
	free_mem2Ddouble(sum2_16_ref);
	free_mem2Ddouble(sum2_8_ref);
	free_mem2Ddouble(sum2_4_ref);

	free_mem2Ddouble(sum_16_U_ref);
	free_mem2Ddouble(sum_8_U_ref);
	free_mem2Ddouble(sum_4_U_ref);
	free_mem2Ddouble(sum2_16_U_ref);
	free_mem2Ddouble(sum2_8_U_ref);
	free_mem2Ddouble(sum2_4_U_ref);

	free_mem2Ddouble(sum_16_V_ref);
	free_mem2Ddouble(sum_8_V_ref);
	free_mem2Ddouble(sum_4_V_ref);
	free_mem2Ddouble(sum2_16_V_ref);
	free_mem2Ddouble(sum2_8_V_ref);
	free_mem2Ddouble(sum2_4_V_ref);

	free_mem2Ddouble(sum_16_8_ref);
	free_mem2Ddouble(sum_8_16_ref);
	free_mem2Ddouble(sum_8_4_ref);
	free_mem2Ddouble(sum_4_8_ref);
	free_mem2Ddouble(sum2_16_8_ref);
	free_mem2Ddouble(sum2_8_16_ref);
	free_mem2Ddouble(sum2_8_4_ref);
	free_mem2Ddouble(sum2_4_8_ref);

	free_mem2Ddouble(sum_16_8_U_ref);
	free_mem2Ddouble(sum_8_16_U_ref);
	free_mem2Ddouble(sum_8_4_U_ref);
	free_mem2Ddouble(sum_4_8_U_ref);
	free_mem2Ddouble(sum2_16_8_U_ref);
	free_mem2Ddouble(sum2_8_16_U_ref);
	free_mem2Ddouble(sum2_8_4_U_ref);
	free_mem2Ddouble(sum2_4_8_U_ref);

	free_mem2Ddouble(sum_16_8_V_ref);
	free_mem2Ddouble(sum_8_16_V_ref);
	free_mem2Ddouble(sum_8_4_V_ref);
	free_mem2Ddouble(sum_4_8_V_ref);
	free_mem2Ddouble(sum2_16_8_V_ref);
	free_mem2Ddouble(sum2_8_16_V_ref);
	free_mem2Ddouble(sum2_8_4_V_ref);
	free_mem2Ddouble(sum2_4_8_V_ref);

	for(i=0;i<num_regions;i++)
	{
		for (trans_count=0;trans_count<img->frmSizeInMbs;trans_count++)
		{
			for (next_count=0;next_count<4;next_count++)
			{
				free(trans_Y[i][trans_count].next[next_count].next);
			}
			free(trans_Y[i][trans_count].next);
		}
		free(trans_Y[i]);
		for (trans_count=0;trans_count<img->frmSizeInMbs/4;trans_count++)
		{
			for (next_count=0;next_count<4;next_count++)
			{
				free(trans_U[i][trans_count].next[next_count].next);
			}
			free(trans_U[i][trans_count].next);
		}
		free(trans_U[i]);
		for (trans_count=0;trans_count<img->frmSizeInMbs/4;trans_count++)
		{
			for (next_count=0;next_count<4;next_count++)
			{
				free(trans_V[i][trans_count].next[next_count].next);
			}
			free(trans_V[i][trans_count].next);
		}
		free(trans_V[i]);
	}

	free(trans_Y);
	free(trans_U);
	free(trans_V);
//huffman//////ZZZZZZZZZZZZZZZZZZZ	
// 	for(j=0;j<4;j++)
// 		for (i=0;i<input->num_regions;i++)
// 			free(code_ptr[j][i]);
// 		
// 	for (i=0;i<input->num_regions;i++)
// 	{
// 		free(x_trans[i]);
// 		free(y_trans[i]);
// 		free(s_trans[i]);
// 		free(o_trans[i]);	
// 	}
// 	free(x_trans);
// 	free(y_trans);
// 	free(s_trans);
// 	free(o_trans);
// 
// 	for (i=0;i<input->num_regions;i++)
// 	{
// 		for (j=0;j<3;j++)
// 			free(Huff_ptr[i][j]);
// 		free(Huff_ptr[i]);
// 	}
// 	free(Huff_ptr);
// 
/*	for (i=0;i<input->num_regions;i++)
	{
		free(o_p[i]);
		free(xy_p[i]);
		free(s_p[i]);				
	}*/
//huffman//////ZZZZZZZZZZZZZZZZZZZ
// 	free(xy_p);
// 	free(s_p);
// 	free(o_p);
// 
// 	free(redual_ple);
// 	free(redual_pru);
// 
// 	free(redual_ple_u);
// 	free(redual_pru_v);
// 	free(redual_ple_u);
// 	free(redual_pru_v);
}

