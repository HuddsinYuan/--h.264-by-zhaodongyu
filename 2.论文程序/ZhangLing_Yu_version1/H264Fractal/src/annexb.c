
/*!
 *************************************************************************************
 * \file annexb.c
 *
 * \brief
 *    Annex B Byte Stream format NAL Unit writing routines
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *************************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "global.h"
#include "nalucommon.h"

#include "i_global.h"


/////////////////////////////////////>

FILE *bits = NULL;                //!< the bit stream file
static int FindStartCode (unsigned char *Buf, int zeros_in_startcode);

int IsFirstByteStreamNALU=1;
int LastAccessUnitExists=0;
int NALUCount=0;

/////////////////////////////////////<



static FILE *f = NULL;    // the output file


/*!
 ********************************************************************************************
 * \brief 
 *    Writes a NALU to the Annex B Byte Stream
 *
 * \return
 *    number of bits written
 *
 ********************************************************************************************
*/
int WriteAnnexbNALU (NALU_t *n)
{
  int BitsWritten = 0;

  assert (n != NULL);
  assert (n->forbidden_bit == 0);
  assert (f != NULL);
  assert (n->startcodeprefix_len == 3 || n->startcodeprefix_len == 4);

//  printf ("WriteAnnexbNALU: writing %d bytes w/ startcode_len %d\n", n->len+1, n->startcodeprefix_len); 
  if (n->startcodeprefix_len > 3)
  {
    putc (0, f);
    BitsWritten =+ 8;
//     printf("WriteAnnexbNALU:BitsWritten(0)=%d\n",BitsWritten);
  }
  putc (0, f);
  putc (0, f);
  putc (1, f);
  BitsWritten += 24;
//   printf("WriteAnnexbNALU:BitsWritten(1)=%d\n",BitsWritten);
 
  n->buf[0] =
    n->forbidden_bit << 7      |
    n->nal_reference_idc << 5  |
    n->nal_unit_type;

//  printf ("WriteAnnexbNALU:First Byte %x, nal_ref_idc %x, nal_unit_type %d\n", n->buf[0], n->nal_reference_idc, n->nal_unit_type);

  if (n->len != fwrite (n->buf, 1, n->len, f))
  {
    printf ("Fatal: cannot write %d bytes to bitstream file, exit (-1)\n", n->len);
    exit (-1);
  }
  BitsWritten += n->len * 8;//2603
//   printf("WriteAnnexbNALU:BitsWritten(2)=%d\n",BitsWritten);

  fflush (f);
#if TRACE
  fprintf (p_trace, "\n\nAnnex B NALU w/ %s startcode, len %d, forbidden_bit %d, nal_reference_idc %d, nal_unit_type %d\n\n",
    n->startcodeprefix_len == 4?"long":"short", n->len, n->forbidden_bit, n->nal_reference_idc, n->nal_unit_type);
  fflush (p_trace);
#endif
  return BitsWritten;
}


/*!
 ********************************************************************************************
 * \brief 
 *    Opens the output file for the bytestream    
 *
 * \param Filename
 *    The filename of the file to be opened
 *
 * \return
 *    none.  Function terminates the program in case of an error
 *
 ********************************************************************************************
*/
void OpenAnnexbFile(char *Filename)
{
  if ((f=fopen (Filename, "wb"))==NULL)
  {
    printf ("Fatal: cannot open Annex B bytestream file '%s', exit (-1)\n", Filename);
    exit (-1);
  }
}


/*!
 ********************************************************************************************
 * \brief 
 *    Closes the output bit stream file
 *
 * \return
 *    none.  Funtion trerminates the program in case of an error
 ********************************************************************************************
*/
void CloseAnnexbFile() {
  if (fclose (f))
  {
    printf ("Fatal: cannot close Annex B bytestream file, exit (-1)\n");
    exit (-1);
  }
}

/////////////////////////////////////>


/*!
************************************************************************
* \brief
*    Returns the size of the NALU (bits between start codes in case of
*    Annex B.  nalu->buf and nalu->len are filled.  Other field in
*    nalu-> remain uninitialized (will be taken care of by NALUtoRBSP.
*
* \return
*     0 if there is nothing any more to read (EOF)
*    -1 in case of any error
*
*  \note Side-effect: Returns length of start-code in bytes. 
*
* \note
*   GetAnnexbNALU expects start codes at byte aligned positions in the file
*
************************************************************************
*/

int GetAnnexbNALU (NALU_t *nalu)//第三次进入的信息
{
	int info2, info3, pos = 0;
	int StartCodeFound, rewind;
	char *Buf;
	int LeadingZero8BitsCount=0, TrailingZero8Bits=0;
    
	if ((Buf = (char*)calloc (nalu->max_size , sizeof(char))) == NULL) no_mem_exit("GetAnnexbNALU: Buf");
	
	while(!feof(bits) && (Buf[pos++]=fgetc(bits))==0);
	
	if(feof(bits))
	{
		if(pos==0)
			return 0;
		else
		{
// 			printf( "GetAnnexbNALU can't read start code\n");
			free(Buf);
			return -1;
		}
	}
	
	if(Buf[pos-1]!=1)//1!=1
	{
// 		printf ("GetAnnexbNALU: no Start Code at the begin of the NALU, return -1\n");
		free(Buf);
		return -1;
	}
	
	if(pos<3)
	{
// 		printf ("GetAnnexbNALU: no Start Code at the begin of the NALU, return -1\n");
		free(Buf);
		return -1;
	}
	else if(pos==3)
	{
		nalu->startcodeprefix_len = 3;
		LeadingZero8BitsCount = 0;
	}
	else//=4
	{
		LeadingZero8BitsCount = pos-4;
		nalu->startcodeprefix_len = 4;
	}
	
	//the 1st byte stream NAL unit can has leading_zero_8bits, but subsequent ones are not
	//allowed to contain it since these zeros(if any) are considered trailing_zero_8bits
	//of the previous byte stream NAL unit.
	if(!IsFirstByteStreamNALU && LeadingZero8BitsCount>0)//!0&&0>0
	{
// 		printf ("GetAnnexbNALU: The leading_zero_8bits syntax can only be present in the first byte stream NAL unit, return -1\n");
		free(Buf);
		return -1;
	}
	IsFirstByteStreamNALU=0;
	
	StartCodeFound = 0;
	info2 = 0;
	info3 = 0;
	
	while (!StartCodeFound)
	{
		if (feof (bits))//0
		{
			//Count the trailing_zero_8bits
			while(Buf[pos-2-TrailingZero8Bits]==0)
				TrailingZero8Bits++;
			nalu->len = (pos-1)-nalu->startcodeprefix_len-LeadingZero8BitsCount-TrailingZero8Bits;
			memcpy (nalu->buf, &Buf[LeadingZero8BitsCount+nalu->startcodeprefix_len], nalu->len);     
			nalu->forbidden_bit = (nalu->buf[0]>>7) & 1;
			nalu->nal_reference_idc = (nalu->buf[0]>>5) & 3;
			nalu->nal_unit_type = (nalu->buf[0]) & 0x1f;
			
			// printf ("GetAnnexbNALU, eof case: pos %d nalu->len %d, nalu->reference_idc %d, nal_unit_type %d \n", pos, nalu->len, nalu->nal_reference_idc, nalu->nal_unit_type);
			
#if TRACE
			fprintf (p_trace, "\n\nLast NALU in File\n\n");
			fprintf (p_trace, "Annex B NALU w/ %s startcode, len %d, forbidden_bit %d, nal_reference_idc %d, nal_unit_type %d\n\n",
				nalu->startcodeprefix_len == 4?"long":"short", nalu->len, nalu->forbidden_bit, nalu->nal_reference_idc, nalu->nal_unit_type);
			fflush (p_trace);
#endif
			free(Buf);
			return pos-1;
		}
		Buf[pos++] = fgetc (bits);//return 101//136//132
		info3 = FindStartCode(&Buf[pos-4], 3);//0
		if(info3 != 1)
			info2 = FindStartCode(&Buf[pos-3], 2);//0
		StartCodeFound = (info2 == 1 || info3 == 1);//0
	}
	
	//Count the trailing_zero_8bits
	if(info3==1)	//if the detected start code is 00 00 01, trailing_zero_8bits is sure not to be present
	{
		while(Buf[pos-5-TrailingZero8Bits]==0)
			TrailingZero8Bits++;
	}
	// Here, we have found another start code (and read length of startcode bytes more than we should
	// have.  Hence, go back in the file
	rewind = 0;
	if(info3 == 1)
		rewind = -4;
	else if (info2 == 1)
		rewind = -3;
	else
		printf(" Panic: Error in next start code search \n");
	
	if (0 != fseek (bits, rewind, SEEK_CUR))
	{
		printf (errortext, ET_SIZE, "GetAnnexbNALU: Cannot fseek %d in the bit stream file", rewind);
		free(Buf);
		error(errortext, 600);
	}
	
	// Here the leading zeros(if any), Start code, the complete NALU, trailing zeros(if any)
	// and the next start code is in the Buf.
	// The size of Buf is pos, pos+rewind are the number of bytes excluding the next
	// start code, and (pos+rewind)-startcodeprefix_len-LeadingZero8BitsCount-TrailingZero8Bits
	// is the size of the NALU.
	
	nalu->len = (pos+rewind)-nalu->startcodeprefix_len-LeadingZero8BitsCount-TrailingZero8Bits;
	memcpy (nalu->buf, &Buf[LeadingZero8BitsCount+nalu->startcodeprefix_len], nalu->len);
	nalu->forbidden_bit = (nalu->buf[0]>>7) & 1;
	nalu->nal_reference_idc = (nalu->buf[0]>>5) & 3;
	nalu->nal_unit_type = (nalu->buf[0]) & 0x1f;
	
	
	//printf ("GetAnnexbNALU, regular case: pos %d nalu->len %d, nalu->reference_idc %d, nal_unit_type %d \n", pos, nalu->len, nalu->nal_reference_idc, nalu->nal_unit_type);
#if TRACE
	fprintf (p_trace, "\n\nAnnex B NALU w/ %s startcode, len %d, forbidden_bit %d, nal_reference_idc %d, nal_unit_type %d\n\n",
		nalu->startcodeprefix_len == 4?"long":"short", nalu->len, nalu->forbidden_bit, nalu->nal_reference_idc, nalu->nal_unit_type);
	fflush (p_trace);
#endif
	
	free(Buf);
	
	return (pos+rewind);
}




/*!
************************************************************************
* \brief
*    Opens the bit stream file named fn
* \return
*    none
************************************************************************
*/
void OpenBitstreamFile (char *fn)
{
	if (NULL == (bits=fopen(fn, "rb")))
	{
		printf (errortext, ET_SIZE, "Cannot open Annex B ByteStream file '%s'", input->infile);
		error(errortext,500);
	}
}


/*!
************************************************************************
* \brief
*    Closes the bit stream file
************************************************************************
*/
void CloseBitstreamFile()
{
	fclose (bits);
}


/*!
************************************************************************
* \brief
*    returns if new start code is found at byte aligned position buf.
*    new-startcode is of form N 0x00 bytes, followed by a 0x01 byte.
*
*  \return
*     1 if start-code is found or                      \n
*     0, indicating that there is no start code
*
*  \param Buf
*     pointer to byte-stream
*  \param zeros_in_startcode
*     indicates number of 0x00 bytes in start-code.
************************************************************************
*/
static int FindStartCode (unsigned char *Buf, int zeros_in_startcode)
{
	int info;
	int i;
	
	info = 1;
	for (i = 0; i < zeros_in_startcode; i++)
		if(Buf[i] != 0)
			info = 0;
		
		if(Buf[i] != 1)
			info = 0;
		return info;
}

void CheckZeroByteNonVCL(NALU_t *nalu, int * ret)
{
	int CheckZeroByte=0;
	
	//This function deals only with non-VCL NAL units
	if(nalu->nal_unit_type>=1&&nalu->nal_unit_type<=5)
		return;
	
	//for SPS and PPS, zero_byte shall exist
	if(nalu->nal_unit_type==NALU_TYPE_SPS || nalu->nal_unit_type==NALU_TYPE_PPS)
		CheckZeroByte=1;
	//check the possibility of the current NALU to be the start of a new access unit, according to 7.4.1.2.3
	if(nalu->nal_unit_type==NALU_TYPE_AUD  || nalu->nal_unit_type==NALU_TYPE_SPS ||
		nalu->nal_unit_type==NALU_TYPE_PPS || nalu->nal_unit_type==NALU_TYPE_SEI ||
		(nalu->nal_unit_type>=13 && nalu->nal_unit_type<=18))
	{
		if(LastAccessUnitExists)
		{
			LastAccessUnitExists=0;		//deliver the last access unit to decoder
			NALUCount=0;
		}
	}
	NALUCount++;
	//for the first NAL unit in an access unit, zero_byte shall exists
	if(NALUCount==1)
		CheckZeroByte=1;
	if(CheckZeroByte && nalu->startcodeprefix_len==3)
	{
		printf("warning: zero_byte shall exist\n");
		//because it is not a very serious problem, we may not indicate an error by setting ret to -1
		//*ret=-1;
	}
}

void CheckZeroByteVCL(NALU_t *nalu, int * ret)
{
	int CheckZeroByte=0;
	
	//This function deals only with VCL NAL units
	if(!(nalu->nal_unit_type>=1&&nalu->nal_unit_type<=5))
		return;
	
	if(LastAccessUnitExists)
	{
		NALUCount=0;
	}
	NALUCount++;
	//the first VCL NAL unit that is the first NAL unit after last VCL NAL unit indicates 
	//the start of a new access unit and hence the first NAL unit of the new access unit.						(sounds like a tongue twister :-)
	if(NALUCount==1)
		CheckZeroByte=1;
	LastAccessUnitExists=1;
	if(CheckZeroByte && nalu->startcodeprefix_len==3)
	{
		printf("warning: zero_byte shall exist\n");
		//because it is not a very serious problem, we may not indicate an error by setting ret to -1
		//*ret=-1;
	}
}


/////////////////////////////////////<