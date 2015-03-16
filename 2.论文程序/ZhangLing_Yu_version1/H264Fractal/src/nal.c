
/*!
 **************************************************************************************
 * \file
 *    nal.c
 * \brief
 *    Handles the operations on converting String of Data Bits (SODB)
 *    to Raw Byte Sequence Payload (RBSP), and then 
 *    onto Encapsulate Byte Sequence Payload (EBSP).
 *  \date 14 June 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *      - Shankar Regunathan                  <shanre@microsoft.de>
 *      - Stephan Wenger                      <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */


// #include "contributors.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "i_defines.h"
#include "global.h"

 /*!
 ************************************************************************
 * \brief
 *    Converts String Of Data Bits (SODB) to Raw Byte Sequence 
 *    Packet (RBSP)
 * \param currStream
 *        Bitstream which contains data bits.
 * \return None
 * \note currStream is byte-aligned at the end of this function
 *    
 ************************************************************************
*/

static byte *NAL_Payload_buffer;

void SODBtoRBSP(Bitstream *currStream)
{
  currStream->byte_buf <<= 1;//0
  currStream->byte_buf |= 1;//1
  currStream->bits_to_go--;//7
  currStream->byte_buf <<= currStream->bits_to_go;//128
  currStream->streamBuffer[currStream->byte_pos++] = currStream->byte_buf;//204
  currStream->bits_to_go = 8;
  currStream->byte_buf = 0;
//   printf("SODBtoRBSP\n");
}


/*!
************************************************************************
*  \brief
*     This function converts a RBSP payload to an EBSP payload
*     
*  \param streamBuffer
*       pointer to data bits
*  \param begin_bytepos
*            The byte position after start-code, after which stuffing to
*            prevent start-code emulation begins.
*  \param end_bytepos
*           Size of streamBuffer in bytes.
*  \param min_num_bytes
*           Minimum number of bytes in payload. Should be 0 for VLC entropy
*           coding mode. Determines number of stuffed words for CABAC mode.
*  \return 
*           Size of streamBuffer after stuffing.
*  \note
*      NAL_Payload_buffer is used as temporary buffer to store data.
*
*
************************************************************************
*/

int RBSPtoEBSP(byte *streamBuffer, int begin_bytepos, int end_bytepos, int min_num_bytes)
{

  int i, j, count;
//   printf("RBSPtoEBSP\n");

  for(i = begin_bytepos; i < end_bytepos; i++)
    NAL_Payload_buffer[i] = streamBuffer[i];

  count = 0;
  j = begin_bytepos;//0
  for(i = begin_bytepos; i < end_bytepos; i++) 
  {
    if(count == ZEROBYTES_SHORTSTARTCODE && !(NAL_Payload_buffer[i] & 0xFC)) //0
    {
      streamBuffer[j] = 0x03;
      j++;
      count = 0;   
    }
    streamBuffer[j] = NAL_Payload_buffer[i];
    if(NAL_Payload_buffer[i] == 0x00) //0     
      count++;
    else 
      count = 0;
    j++;
  }
  while (j < begin_bytepos+min_num_bytes) {//0
    streamBuffer[j] = 0x00; // cabac stuffing word
    streamBuffer[j+1] = 0x00;
    streamBuffer[j+2] = 0x03;
    j += 3;
    stat->bit_use_stuffingBits[img->type]+=16;
  }
  return j;
}

 /*!
 ************************************************************************
 * \brief
 *    Initializes NAL module (allocates NAL_Payload_buffer)
 ************************************************************************
*/

void AllocNalPayloadBuffer()
{

  const int buffer_size = (input->imagewidth * input->imageheight * 5);//4 // AH 190202: There can be data expansion with 
                                                    // does not everflow. 4 is probably safe multiplier.
  FreeNalPayloadBuffer();

  NAL_Payload_buffer = (byte *) calloc(buffer_size, sizeof(byte));
  assert (NAL_Payload_buffer != NULL);
//   printf("AllocNalPayloadBuffer\n");
}


 /*!
 ************************************************************************
 * \brief
 *   Finits NAL module (frees NAL_Payload_buffer)
 ************************************************************************
*/

void FreeNalPayloadBuffer()
{
  if(NAL_Payload_buffer)
  {
    free(NAL_Payload_buffer);
    NAL_Payload_buffer=NULL;
  }
//   printf("FreeNalPayloadBuffer\n");
}




/////////////////////////////////////>






 /*!
 ************************************************************************
 * \brief
 *    Converts RBSP to string of data bits
 * \param streamBuffer
 *          pointer to buffer containing data
 *  \param last_byte_pos
 *          position of the last byte containing data.
 * \return last_byte_pos
 *          position of the last byte pos. If the last-byte was entirely a stuffing byte,
 *          it is removed, and the last_byte_pos is updated.
 *  
************************************************************************/

int RBSPtoSODB(byte *streamBuffer, int last_byte_pos)
{
  int ctr_bit, bitoffset;
  
  bitoffset = 0; 
  //find trailing 1
  ctr_bit = (streamBuffer[last_byte_pos-1] & (0x01<<bitoffset));   // set up control bit
  
  while (ctr_bit==0)
  {                 // find trailing 1 bit
    bitoffset++;
    if(bitoffset == 8) 
    {
      if(last_byte_pos == 0)
        printf(" Panic: All zero data sequence in RBSP \n");
      assert(last_byte_pos != 0);
      last_byte_pos -= 1;
      bitoffset = 0;
    }
    ctr_bit= streamBuffer[last_byte_pos-1] & (0x01<<(bitoffset));
  }
  
  
  // We keep the stop bit for now
/*  if (remove_stop)
  {
    streamBuffer[last_byte_pos-1] -= (0x01<<(bitoffset));
    if(bitoffset == 7)
      return(last_byte_pos-1);
    else
      return(last_byte_pos);
  }
*/
  return(last_byte_pos);
  
}

/*!
************************************************************************
* \brief
*    Converts Encapsulated Byte Sequence Packets to RBSP
* \param streamBuffer
*    pointer to data stream
* \param end_bytepos
*    size of data stream
* \param begin_bytepos
*    Position after beginning 
************************************************************************/


int EBSPtoRBSP(byte *streamBuffer, int end_bytepos, int begin_bytepos)
{
	int i, j, count;
	count = 0;
	
	if(end_bytepos < begin_bytepos)
		return end_bytepos;
	
	j = begin_bytepos;
	
	for(i = begin_bytepos; i < end_bytepos; i++) 
	{ //starting from begin_bytepos to avoid header information
		if(count == ZEROBYTES_SHORTSTARTCODE && streamBuffer[i] == 0x03) 
		{
			i++;
			count = 0;
		}
		streamBuffer[j] = streamBuffer[i];
		if(streamBuffer[i] == 0x00)
			count++;
		else
			count = 0;
		j++;
	}
	
	return j;
}




/////////////////////////////////////<





