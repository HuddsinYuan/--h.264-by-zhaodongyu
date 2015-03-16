#include "dct.h"

void swap(double *a,double *b)
{
	double temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

int Initial_1D_Param(int len)
{
	int hr = 0;	
	return hr;
}


int FDCT_1D(double *pdata, int bitlen)
{
	int i;
	int hr = 0;
	int Len = 1<<bitlen;
	double Coef = sqrt(2.0/Len);

	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}
	
		
	FDCT_1D_No_Param(pdata, bitlen);
	
	for (i = 0; i < Len; ++i)
	{
		pdata[i] = Coef * pdata[i];
	}
	
	return hr;
}

int FDCT_2D(double *pdata, int bitrow, int bitcol)
{
	int hr = 0;
	int row = 1 << bitrow;
	int col = 1 << bitcol;
	int i_row,i_col;
	double coef = 2.0f / sqrt( (double)(row*col) );

	if (pdata == NULL || bitrow <= 0 || bitcol <= 0)
	{
		return 1;
	}

	Initial_2D_Param(row, col);

	for (i_row = 0; i_row < row; ++i_row)
	{
		FDCT_1D_No_Param( pdata + i_row*col, bitcol);///////////
	}

	for (i_col = 0; i_col < col; ++i_col)
	{
		for (i_row = 0; i_row < row; ++i_row)
		{
			m_pdata_2D_Temp[ i_row ] = pdata[ i_row*col + i_col ];
		}

		FDCT_1D_No_Param(m_pdata_2D_Temp, bitrow);

		for (i_row = 0; i_row < row; ++i_row)
		{
			pdata[ i_row*col + i_col ] = m_pdata_2D_Temp[i_row] * coef;
		}
	}
	return hr;
}

int Initial_2D_Param(int row, int col)
{
	int hr = 0;
	if (row <= 0 || col <= 0)
	{
		return 1;
	}

	SAFE_DELETE(m_pdata_2D_Temp);

	m_pdata_2D_Temp = (double*)malloc(sizeof(double)*row);
	if (m_pdata_2D_Temp == NULL)
	{
		hr = 1;
	}
	
	return hr;
}

int FDCT_1D_No_Param(double *pdata, int bitlen)
{
	int hr = 0;
	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	Initial_DCT_Param(bitlen);

	DCTForward(pdata, bitlen);

	DCTBackward(pdata, bitlen);

	FBitReverse(pdata, bitlen);

	pdata[0] = 1 / sqrt(2.0) * pdata[0];

	return hr;
}

int Initial_DCT_Param(int bitlen)
{
	int hr = 0,i,j;
	
	int EndStart = 0;
	int Factor = 0;
	
	int Len = 1 << bitlen;
	if (bitlen <= 0)
	{
		return 1;
	}

	SAFE_DELETE(m_pCos_Temp);
	m_pCos_Temp = (double*)malloc(sizeof(double)*Len);
	for (i = 0; i < Len>>1; ++i)
	{
		m_pCos_Temp[Len-i-1] = (double)(2*i + 1);
	}

	
	for (i = 0; i < bitlen - 1; ++i)
	{
		EndStart = 1 << (bitlen - 1 - i);
		Factor = 1 << (i + 1);
		for (j = 0; j < EndStart>>1; ++j )
		{
			m_pCos_Temp[EndStart-j-1] = Factor * m_pCos_Temp[Len-j-1];
		}		
	}

	for (i = 0; i < Len; ++i)
	{
		m_pCos_Temp[i] = 2.0 * cos(m_pCos_Temp[i] * PI / (Len<<1) );
	}	

	return hr;
}

int DCTForward(double *pdata, int bitlen)
{
	int hr = 0;
	
	int Len = 1<<bitlen;
	double Temp1 = 0.0;
	double Temp2 = 0.0;
	int Wings, HalfWing, WingLen;
	int i_bitlen,i_Wing,i_HalfWing;

	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	for (i_bitlen = 0; i_bitlen < bitlen; ++i_bitlen)
	{
		Wings    = 1<<i_bitlen;
		WingLen  = Len>>i_bitlen;
		HalfWing = WingLen>>1;
		
		for (i_Wing = 0; i_Wing < Wings; ++i_Wing)
		{
			for (i_HalfWing = 0; i_HalfWing < HalfWing; ++i_HalfWing)
			{
				Temp1 = pdata[ i_Wing * WingLen + i_HalfWing ];
				Temp2 = pdata[ (i_Wing+1) * WingLen - 1 - i_HalfWing];
				if (i_Wing%2)
				{
					swap(&Temp1, &Temp2);
				}

				pdata[ i_Wing*WingLen + i_HalfWing ] = Temp1 + Temp2;
				pdata[ (i_Wing+1) * WingLen - 1 - i_HalfWing ] = (Temp1 - Temp2) * m_pCos_Temp[WingLen-1-i_HalfWing];
			}
		}
	}
	return hr;
}

int DCTBackward(double *pdata, int bitlen)
{
	int hr = 0;
	int i_bitlen,i_HalfWing;
	int i_Wing;
	int Len = 1<<bitlen;
	int Wings, HalfWing, WingLen, Temp1, Temp2;

	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	for ( i_bitlen= bitlen-1; i_bitlen >= 0; --i_bitlen)
	{
		Wings    = 1<<i_bitlen;
		WingLen  = 1<<(bitlen-i_bitlen);
		HalfWing = WingLen>>1;

		for ( i_Wing= 0; i_Wing < Wings; ++i_Wing)
		{
			for (i_HalfWing = 0; i_HalfWing < HalfWing; ++i_HalfWing)
			{
				if (i_HalfWing == 0)
				{
					pdata[HalfWing + i_Wing*WingLen] *= 0.5;
				}
				else
				{
					Temp1 = BitReverse(i_HalfWing, bitlen-i_bitlen-1);
					Temp2 = BitReverse(i_HalfWing-1, bitlen-i_bitlen-1);

					pdata[HalfWing + i_Wing*WingLen + Temp1] -= pdata[HalfWing + i_Wing*WingLen + Temp2];
				}
			}
		}
	}
	return hr;
}





int FIDCT_1D(double *pdata, int bitlen)
{
	int hr = 0,i;
	
	int Len = 1<<bitlen;
	double Coef = sqrt((double)Len/2.0);

	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	for (i = 0; i < Len; ++i)
	{
		pdata[i] = Coef * pdata[i];
	}

	FIDCT_1D_No_Param(pdata, bitlen);

	return hr;
}

int FIDCT_2D(double *pdata, int bitrow, int bitcol)
{
	int hr = 0,i_row,i_col;
	
	int row = 1 << bitrow;
	int col = 1 << bitcol;
	double coef = sqrt((double)(row*col)) / 2.0;

	if (pdata == NULL || bitrow <= 0 || bitcol <= 0)
	{
		return 1;
	}

	Initial_2D_Param(row, col);

	for (i_row = 0; i_row < row; ++i_row)
	{
		FIDCT_1D_No_Param( pdata + i_row*col, bitcol);
	}

	for (i_col = 0; i_col < col; ++i_col)
	{
		for (i_row = 0; i_row < row; ++i_row)
		{
			m_pdata_2D_Temp[ i_row ] = pdata[ i_row*col + i_col ];
		}

		FIDCT_1D_No_Param(m_pdata_2D_Temp, bitrow);

		for (i_row = 0; i_row < row; ++i_row)
		{
			pdata[ i_row*col + i_col ] = m_pdata_2D_Temp[i_row] * coef;
		}
	}

	return hr;
}

int FIDCT_1D_No_Param(double *pdata, int bitlen)
{
	int hr = 0;
	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	Initial_IDCT_Param(bitlen);

	pdata[0] = sqrt(2.0) * pdata[0];

	FBitReverse(pdata, bitlen);	

	IDCTForward(pdata, bitlen);

	IDCTBackward(pdata, bitlen);	

	return hr;
}

int Initial_IDCT_Param(int bitlen)
{
	int hr = 0,i,j;
	int Len = 1 << bitlen;
	int EndStart = 0;
	int Factor = 0;

	if (bitlen <= 0)
	{
		return 1;
	}

	SAFE_DELETE(m_pCos_Temp);
	m_pCos_Temp = (double*)malloc(sizeof(double)*Len);
	for (i = 0; i < Len>>1; ++i)
	{
		m_pCos_Temp[Len-i-1] = (double)(2*i + 1);
	}

	for (i = 0; i < bitlen - 1; ++i)
	{
		EndStart = 1 << (bitlen - 1 - i);
		Factor = 1 << (i + 1);
		for (j = 0; j < EndStart>>1; ++j )
		{
			m_pCos_Temp[EndStart-j-1] = Factor * m_pCos_Temp[Len-j-1];
		}		
	}

	for (i = 0; i < Len; ++i)
	{
		m_pCos_Temp[i] = 1.0 / (2.0 * cos(m_pCos_Temp[i] * PI / (Len<<1) ));
	}	

	return hr;
}

int IDCTBackward(double *pdata, int bitlen)
{
	int hr = 0;
	
	int Len = 1<<bitlen;
	double Temp1 = 0.0;
	double Temp2 = 0.0;
	int Wings, HalfWing, WingLen;
	int i_bitlen,i_Wing,i_HalfWing;

	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	for (i_bitlen = bitlen-1; i_bitlen >= 0; --i_bitlen)
	{
		Wings    = 1<<i_bitlen;
		WingLen  = Len>>i_bitlen;
		HalfWing = WingLen>>1;

		for (i_Wing = 0; i_Wing < Wings; ++i_Wing)
		{
			for (i_HalfWing = 0; i_HalfWing < HalfWing; ++i_HalfWing)
			{
				Temp1 = pdata[ i_Wing * WingLen + i_HalfWing ];
				Temp2 = pdata[ (i_Wing+1) * WingLen - 1 - i_HalfWing] * m_pCos_Temp[WingLen-1-i_HalfWing];
				if (i_Wing%2)
				{
					pdata[ i_Wing*WingLen + i_HalfWing ] = (Temp1 - Temp2) * 0.5;
					pdata[ (i_Wing+1) * WingLen - 1 - i_HalfWing ] = (Temp1 + Temp2) * 0.5;
				}
				else
				{
					pdata[ i_Wing*WingLen + i_HalfWing ] = (Temp1 + Temp2) * 0.5;
					pdata[ (i_Wing+1) * WingLen - 1 - i_HalfWing ] = (Temp1 - Temp2) * 0.5;
				}
				
			}
		}
	}

	return hr;
}

int IDCTForward(double *pdata, int bitlen)
{
	int hr = 0;
	
	int Len = 1<<bitlen;
	int Wings, HalfWing, WingLen, Temp1, Temp2,i_bitlen,i_Wing,i_HalfWing;

	if (pdata == NULL || bitlen <= 0)
	{
		return 1;
	}

	for (i_bitlen = 0; i_bitlen <bitlen; ++i_bitlen)
	{
		Wings    = 1<<i_bitlen;
		WingLen  = 1<<(bitlen-i_bitlen);
		HalfWing = WingLen>>1;

		for (i_Wing = 0; i_Wing < Wings; ++i_Wing)
		{
			for (i_HalfWing = HalfWing-1; i_HalfWing >= 0; --i_HalfWing)
			{
				if (i_HalfWing == 0)
				{
					pdata[HalfWing + i_Wing*WingLen] *= 2.0;
				}
				else
				{
					Temp1 = BitReverse(i_HalfWing, bitlen-i_bitlen-1);
					Temp2 = BitReverse(i_HalfWing-1, bitlen-i_bitlen-1);

					pdata[HalfWing + i_Wing*WingLen + Temp1] += pdata[HalfWing + i_Wing*WingLen + Temp2];
				}
			}
		}
	}

	return hr;
}

int FBitReverse(double *pdata, int bitlen)
{
	int hr = 0;
	
	int len = (1<<bitlen) - 1;
	int i = 1;
	
	if (pdata == NULL || bitlen <= 1)
	{
		return 1;
	}

	while (i < len)
	{
		int ii = BitReverse(i, bitlen);
		if (ii > i)
		{
			double Temp = pdata[ii];
			pdata[ii] = pdata[i];
			pdata[i]  = Temp;
		}

		i++;
	}
	return hr;
}

int BitReverse(int bit, int bitlen)
{
	int Temp1 = 0;
	int Temp2 = 1;
	int Half_Len = 1<<(bitlen-1);

	if (bitlen == 0)
	{
		return bit;
	}

	while (Half_Len)
	{
		if (Half_Len&bit)
		{
			Temp1 += Temp2;
		}

		Temp2 <<= 1;
		Half_Len >>= 1;
	}

	return Temp1;
}
