#ifndef  __DCT__H__
#define  __DCT__H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#define PI 3.1415926
#define  SAFE_DELETE(p) \
	     if (p)         \
		 {              \
		    free(p); \
	        p = NULL;   \
		 }
	
	// DCT
    
	int FDCT_2D(double *pdata, int bitrow, int bitcol);
	int FDCT_1D(double *pdata, int bitlen);
	
	// IDCT
	int FIDCT_2D(double *pdata, int bitrow, int bitcol);
	int FIDCT_1D(double *pdata, int bitlen);
		
	int Initial_2D_Param(int row, int col);
	int Initial_1D_Param(int len);
	int FBitReverse(double *pdata, int bitlen);
	int BitReverse(int bit, int bitlen);
	
	// DCT functions
	int FDCT_1D_No_Param(double *pdata, int bitlen);
	int Initial_DCT_Param(int bitlen);	
	int DCTForward(double *pdata,  int bitlen);
	int DCTBackward(double *pdata, int bitlen);
	
	// IDCT functions
	int FIDCT_1D_No_Param(double *pdata, int bitlen);
	int Initial_IDCT_Param(int bitlen);
	int IDCTForward(double *pdata,  int bitlen);
	int IDCTBackward(double *pdata, int bitlen);	
	
	double *m_pdata_2D_Temp;
	double *m_pCos_Temp;
#endif   /* __DCT__H__ */