
/*!
 *****************************************************************************
 *
 * \file rtp_.c
 *
 * \brief
 *    Functions to handle RTP headers and packets per RFC1889 and RTP NAL spec
 *    Functions support little endian systems only (Intel, not Motorola/Sparc)
 *
 * \date
 *    30 September 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <memory.h>
#include <malloc.h>

#include "rtp_.h"
#include "elements.h"
#include "i_defines.h"
#include "header.h"

#include "global.h"
#include "fmo.h"
#include "parsetcommon.h"
#include "parset.h"
#include "nalucommon.h"


// A little trick to avoid those horrible #if TRACE all over the source code
#if TRACE
#define SYMTRACESTRING(s) strncpy(sym.tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // to nothing
#endif


int CurrentRTPTimestamp = 0;      //! The RTP timestamp of the current packet,
                                  //! incremented with all P and I frames
int CurrentRTPSequenceNumber = 0; //! The RTP sequence number of the current packet
                                  //! incremented by one for each sent packet

FILE *f;



/////////////////////////////////////>


FILE *bits;

int RTPReadPacket (RTPpacket_t *p, FILE *bits);


/////////////////////////////////////<



/*!
 *****************************************************************************
 *
 * \brief 
 *    ComposeRTPpacket composes the complete RTP packet using the various
 *    structure members of the RTPpacket_t structure
 *
 * \return
 *    0 in case of success
 *    negative error code in case of failure
 *
 * \par Parameters
 *    Caller is responsible to allocate enough memory for the generated packet
 *    in parameter->packet. Typically a malloc of 12+paylen bytes is sufficient
 *
 * \par Side effects
 *    none
 *
 * \note
 *    Function contains assert() tests for debug purposes (consistency checks
 *    for RTP header fields
 *
 * \date
 *    30 Spetember 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/


int ComposeRTPPacket (RTPpacket_t *p)

{
  // Consistency checks through assert, only used for debug purposes
  assert (p->v == 2);
  assert (p->p == 0);
  assert (p->x == 0);
  assert (p->cc == 0);    // mixer designers need to change this one
  assert (p->m == 0 || p->m == 1);
  assert (p->pt < 128);
  assert (p->seq < 65536);
  assert (p->payload != NULL);
  assert (p->paylen < 65536 - 40);  // 2**16 -40 for IP/UDP/RTP header
  assert (p->packet != NULL);

  // Compose RTP header, little endian

  p->packet[0] = (   (p->v)
                  |  (p->p << 2)
                  |  (p->x << 3)
                  |  (p->cc << 4) );
  p->packet[1] = (   (p->m)
                  |  (p->pt << 1) );
  p->packet[2] = p->seq & 0xff;
  p->packet[3] = (p->seq >> 8) & 0xff;

  memcpy (&p->packet[4], &p->timestamp, 4);  // change to shifts for unified byte sex
  memcpy (&p->packet[8], &p->ssrc, 4);// change to shifts for unified byte sex

  // Copy payload 

  memcpy (&p->packet[12], p->payload, p->paylen);
  p->packlen = p->paylen+12;
  return 0;
}



/*!
 *****************************************************************************
 *
 * \brief 
 *    WriteRTPPacket writes the supplied RTP packet to the output file
 *
 * \return
 *    0 in case of access
 *    <0 in case of write failure (typically fatal)
 *
 * \param p
 *    the RTP packet to be written (after ComposeRTPPacket() )
 * \param f
 *    output file
 *
 * \date
 *    October 23, 2001
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

int WriteRTPPacket (RTPpacket_t *p, FILE *f)

{
  int intime = -1;

  assert (f != NULL);
  assert (p != NULL);


  if (1 != fwrite (&p->packlen, 4, 1, f))
    return -1;
  if (1 != fwrite (&intime, 4, 1, f))
    return -1;
  if (1 != fwrite (p->packet, p->packlen, 1, f))
    return -1;
  return 0;
}





/*!
 *****************************************************************************
 *
 * \brief 
 *    int RTPWriteNALU write a NALU to the RTP file
 *
 * \return
 *    Number of bytes written to output file
 *
 * \par Side effects
 *    Packet written, RTPSequenceNumber and RTPTimestamp updated
 *   
 * \date
 *    December 13, 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/


int WriteRTPNALU (NALU_t *n)
{
  RTPpacket_t *p;

  assert (f != NULL);
  assert (n != NULL);
  assert (n->len < 65000);

  n->buf[0] =
    n->forbidden_bit << 7      |
    n->nal_reference_idc << 5  |
    n->nal_unit_type;

  // Set RTP structure elements and alloca() memory foor the buffers
  if ((p = (RTPpacket_t *) malloc (sizeof (RTPpacket_t))) == NULL)
    no_mem_exit ("RTPWriteNALU-1");
  if ((p->packet = malloc (MAXRTPPACKETSIZE)) == NULL)
    no_mem_exit ("RTPWriteNALU-2");
  if ((p->payload = malloc (MAXRTPPACKETSIZE)) == NULL)
    no_mem_exit ("RTPWriteNALU-3");

  p->v=2;
  p->p=0;
  p->x=0;
  p->cc=0;
  p->m=(n->startcodeprefix_len==4)&1;     // a long startcode of Annex B sets marker bit of RTP
                                          // Not exactly according to the RTP paylaod spec, but
                                          // good enough for now (hopefully).
                                          //! For error resilience work, we need the correct
                                          //! marker bit.  Introduce a nalu->marker and set it in
                                          //! terminate_slice()?
  p->pt=H26LPAYLOADTYPE;
  p->seq=CurrentRTPSequenceNumber++;
  p->timestamp=CurrentRTPTimestamp;
  p->ssrc=H26LSSRC;
  p->paylen = n->len;
  memcpy (p->payload, n->buf, n->len);



  // Generate complete RTP packet
  if (ComposeRTPPacket (p) < 0)
  {
    printf ("Cannot compose RTP packet, exit\n");
    exit (-1);
  }
  if (WriteRTPPacket (p, f) < 0)
  {
    printf ("Cannot write %d bytes of RTP packet to outfile, exit\n", p->packlen);
    exit (-1);
  }
  free (p->packet);
  free (p->payload);
  free (p);
  return (n->len * 8);
}


/*!
 ********************************************************************************************
 * \brief 
 *    RTPUpdateTimestamp: patches the RTP timestamp depending on the TR
 *
 * \param 
 *    tr: TRof the following NALUs
 *
 * \return
 *    none.  
 *
 ********************************************************************************************
*/


void RTPUpdateTimestamp (int tr)
{
  int delta;
  static int oldtr = -1;

  if (oldtr == -1)            // First invocation
  {
    CurrentRTPTimestamp = 0;  //! This is a violation of the security req. of
                              //! RTP (random timestamp), but easier to debug
    oldtr = 0;
    return;
  }

  /*! The following code assumes a wrap around of TR at 256, and
      needs to be changed as soon as this is no more true.
      
      The support for B frames is a bit tricky, because it is not easy to distinguish
      between a natural wrap-around of the tr, and the intentional going back of the
      tr because of a B frame.  It is solved here by a heuristic means: It is assumed that
      B frames are never "older" than 10 tr ticks.  Everything higher than 10 is considered
      a wrap around.
  */

  delta = tr - oldtr;

  if (delta < -10)        // wrap-around
    delta+=256;

  CurrentRTPTimestamp += delta * RTP_TR_TIMESTAMP_MULT;
  oldtr = tr;
}


/*!
 ********************************************************************************************
 * \brief 
 *    Opens the output file for the RTP packet stream
 *
 * \param Filename
 *    The filename of the file to be opened
 *
 * \return
 *    none.  Function terminates the program in case of an error
 *
 ********************************************************************************************
*/

void OpenRTPFile (char *Filename)
{
  if ((f = fopen (Filename, "wb")) == NULL)
  {
    printf ("Fatal: cannot open bitstream file '%s', exit (-1)\n", Filename);
    exit (-1);
  }
}


/*!
 ********************************************************************************************
 * \brief 
 *    Closes the output file for the RTP packet stream
 *
 * \return
 *    none.  Function terminates the program in case of an error
 *
 ********************************************************************************************
*/

void CloseRTPFile ()
{
  fclose(f);
}








#if 0
/*!
 *****************************************************************************
 *
 * \brief 
 *    int aggregationRTPWriteBits (int marker) write the Slice header for the RTP NAL      
 *
 * \return
 *    Number of bytes written to output file
 *
 * \param marker
 *    marker bit,
 *
 * \par Side effects
 *    Packet written, RTPSequenceNumber and RTPTimestamp updated
 *   
 * \date
 *    September 10, 2002
 *
 * \author
 *    Dong Tian   tian@cs.tut.fi
 *****************************************************************************/

int aggregationRTPWriteBits (int Marker, int PacketType, int subPacketType, void * bitstream, 
                    int BitStreamLenInByte, FILE *out)
{
  RTPpacket_t *p;
  int offset;

//  printf( "writing aggregation packet...\n");
  assert (out != NULL);
  assert (BitStreamLenInByte < 65000);
  assert (bitstream != NULL);
  assert ((PacketType&0xf) == 4);

  // Set RTP structure elements and alloca() memory foor the buffers
  p = (RTPpacket_t *) alloca (sizeof (RTPpacket_t));
  p->packet=alloca (MAXRTPPACKETSIZE);
  p->payload=alloca (MAXRTPPACKETSIZE);
  p->v=2;
  p->p=0;
  p->x=0;
  p->cc=0;
  p->m=Marker&1;
  p->pt=H26LPAYLOADTYPE;
  p->seq=CurrentRTPSequenceNumber++;
  p->timestamp=CurrentRTPTimestamp;
  p->ssrc=H26LSSRC;

  offset = 0;
  p->payload[offset++] = PacketType; // This is the first byte of the compound packet

  // FIRST, write the sei message to aggregation packet, if it is available
  if ( HaveAggregationSEI() )
  {
    p->payload[offset++] = sei_message[AGGREGATION_SEI].subPacketType; // this is the first byte of the first subpacket
    *(short*)&(p->payload[offset]) = sei_message[AGGREGATION_SEI].payloadSize;
    offset += 2;
    memcpy (&p->payload[offset], sei_message[AGGREGATION_SEI].data, sei_message[AGGREGATION_SEI].payloadSize);
    offset += sei_message[AGGREGATION_SEI].payloadSize;

    clear_sei_message(AGGREGATION_SEI);
  }

  // SECOND, write other payload to the aggregation packet
  // to do ...

  // LAST, write the slice data to the aggregation packet
  p->payload[offset++] = subPacketType;  // this is the first byte of the second subpacket
  *(short*)&(p->payload[offset]) = BitStreamLenInByte;
  offset += 2;
  memcpy (&p->payload[offset], bitstream, BitStreamLenInByte);
  offset += BitStreamLenInByte;

  p->paylen = offset;  // 1 +3 +seiPayload.payloadSize +3 +BitStreamLenInByte

  // Now the payload is ready, we can ...
  // Generate complete RTP packet
  if (ComposeRTPPacket (p) < 0)
  {
    printf ("Cannot compose RTP packet, exit\n");
    exit (-1);
  }
  if (WriteRTPPacket (p, out) < 0)
  {
    printf ("Cannot write %d bytes of RTP packet to outfile, exit\n", p->packlen);
    exit (-1);
  }
  return (p->packlen);

}


/*!
 *****************************************************************************
 * \isAggregationPacket
 * \brief 
 *    Determine if current packet is normal packet or compound packet (aggregation
 *    packet)
 *
 * \return
 *    return TRUE, if it is compound packet.
 *    return FALSE, otherwise.
 *   
 * \date
 *    September 10, 2002
 *
 * \author
 *    Dong Tian   tian@cs.tut.fi
 *****************************************************************************/
Boolean isAggregationPacket()
{
  if (HaveAggregationSEI())
  {
    return TRUE;
  }
  // Until Sept 2002, the JM will produce aggregation packet only for some SEI messages

  return FALSE;
}

/*!
 *****************************************************************************
 * \PrepareAggregationSEIMessage
 * \brief 
 *    Prepare the aggregation sei message.
 *    
 * \date
 *    September 10, 2002
 *
 * \author
 *    Dong Tian   tian@cs.tut.fi
 *****************************************************************************/
void PrepareAggregationSEIMessage()
{
  Boolean has_aggregation_sei_message = FALSE;
  // prepare the sei message here
  // write the spare picture sei payload to the aggregation sei message
  if (seiHasSparePicture && img->type != B_SLICE)
  {
    FinalizeSpareMBMap();
    assert(seiSparePicturePayload.data->byte_pos == seiSparePicturePayload.payloadSize);
    write_sei_message(AGGREGATION_SEI, seiSparePicturePayload.data->streamBuffer, seiSparePicturePayload.payloadSize, SEI_SPARE_PICTURE);
    has_aggregation_sei_message = TRUE;
  }
  // write the sub sequence information sei paylaod to the aggregation sei message
  if (seiHasSubseqInfo)
  {
    FinalizeSubseqInfo(img->layer);
    write_sei_message(AGGREGATION_SEI, seiSubseqInfo[img->layer].data->streamBuffer, seiSubseqInfo[img->layer].payloadSize, SEI_SUBSEQ_INFORMATION);
    ClearSubseqInfoPayload(img->layer);
    has_aggregation_sei_message = TRUE;
  }
  // write the sub sequence layer information sei paylaod to the aggregation sei message
  if (seiHasSubseqLayerInfo && img->current_frame == 0)
  {
    FinalizeSubseqLayerInfo();
    write_sei_message(AGGREGATION_SEI, seiSubseqLayerInfo.data, seiSubseqLayerInfo.payloadSize, SEI_SUBSEQ_LAYER_CHARACTERISTICS);
    seiHasSubseqLayerInfo = FALSE;
    has_aggregation_sei_message = TRUE;
  }
  // write the sub sequence characteristics payload to the aggregation sei message
  if (seiHasSubseqChar)
  {
    FinalizeSubseqChar();
    write_sei_message(AGGREGATION_SEI, seiSubseqChar.data->streamBuffer, seiSubseqChar.payloadSize, SEI_SUBSEQ_CHARACTERISTICS);
    ClearSubseqCharPayload();
    has_aggregation_sei_message = TRUE;
  }
  // write the pan scan rectangle info sei playload to the aggregation sei message
  if (seiHasPanScanRectInfo)
  {
    FinalizePanScanRectInfo();
    write_sei_message(AGGREGATION_SEI, seiPanScanRectInfo.data->streamBuffer, seiPanScanRectInfo.payloadSize, SEI_PANSCAN_RECT);
    ClearPanScanRectInfoPayload();
    has_aggregation_sei_message = TRUE;
  }
  // write the arbitrary (unregistered) info sei playload to the aggregation sei message
  if (seiHasUser_data_unregistered_info)
  {
    FinalizeUser_data_unregistered();
    write_sei_message(AGGREGATION_SEI, seiUser_data_unregistered.data->streamBuffer, seiUser_data_unregistered.payloadSize, SEI_USER_DATA_UNREGISTERED);
    ClearUser_data_unregistered();
    has_aggregation_sei_message = TRUE;
  }
  // write the arbitrary (unregistered) info sei playload to the aggregation sei message
  if (seiHasUser_data_registered_itu_t_t35_info)
  {
    FinalizeUser_data_registered_itu_t_t35();
    write_sei_message(AGGREGATION_SEI, seiUser_data_registered_itu_t_t35.data->streamBuffer, seiUser_data_registered_itu_t_t35.payloadSize, SEI_USER_DATA_REGISTERED_ITU_T_T35);
    ClearUser_data_registered_itu_t_t35();
    has_aggregation_sei_message = TRUE;
  }
  //write RandomAccess info sei payload to the aggregation sei message
  if (seiHasRandomAccess_info)
  {
    FinalizeRandomAccess();
    write_sei_message(AGGREGATION_SEI, seiRandomAccess.data->streamBuffer, seiRandomAccess.payloadSize, SEI_RANDOM_ACCESS_POINT);
    ClearRandomAccess();
    has_aggregation_sei_message = TRUE;
  }
  // more aggregation sei payload is written here...

  // JVT-D099 write the scene information SEI payload
  if (seiHasSceneInformation)
  {
    FinalizeSceneInformation();
    write_sei_message(AGGREGATION_SEI, seiSceneInformation.data->streamBuffer, seiSceneInformation.payloadSize, SEI_SCENE_INFORMATION);
    has_aggregation_sei_message = TRUE;
  }
  // End JVT-D099

  // after all the sei payload is written
  if (has_aggregation_sei_message)
    finalize_sei_message(AGGREGATION_SEI);
}

/*!
 *****************************************************************************
 * \begin_sub_sequence_rtp
 * \brief 
 *    do some initialization for sub-sequence under rtp
 *    
 * \date
 *    September 10, 2002
 *
 * \author
 *    Dong Tian   tian@cs.tut.fi
 *****************************************************************************/

void begin_sub_sequence_rtp()
{
  if ( input->of_mode != PAR_OF_RTP || input->NumFramesInELSubSeq == 0 ) 
    return;

  // begin to encode the base layer subseq
  if ( IMG_NUMBER == 0 )
  {
//    printf("begin to encode the base layer subseq\n");
    InitSubseqInfo(0);
    if (1)
      UpdateSubseqChar();
  }
  // begin to encode the enhanced layer subseq
  if ( IMG_NUMBER % (input->NumFramesInELSubSeq+1) == 1 )
  {
//    printf("begin to encode the enhanced layer subseq\n");
    InitSubseqInfo(1);  // init the sub-sequence in the enhanced layer
//    add_dependent_subseq(1);
    if (1)
      UpdateSubseqChar();
  }
}

/*!
 *****************************************************************************
 * \end_sub_sequence_rtp
 * \brief 
 *    do nothing
 *    
 * \date
 *    September 10, 2002
 *
 * \author
 *    Dong Tian   tian@cs.tut.fi
 *****************************************************************************/
void end_sub_sequence_rtp()
{
  // end of the base layer:
  if ( img->current_frame == input->no_frames_h264-1 )
  {
//    printf("end of encoding the base layer subseq\n");
    CloseSubseqInfo(0);
//    updateSubSequenceBox(0);
  }
  // end of the enhanced layer:
  if ( ((IMG_NUMBER%(input->NumFramesInELSubSeq+1)==0) && (input->successive_Bframe != 0) && (IMG_NUMBER>0)) || // there are B frames
    ((IMG_NUMBER%(input->NumFramesInELSubSeq+1)==input->NumFramesInELSubSeq) && (input->successive_Bframe==0))   // there are no B frames
    )
  {
//    printf("end of encoding the enhanced layer subseq\n");
    CloseSubseqInfo(1);
//    add_dependent_subseq(1);
//    updateSubSequenceBox(1);
  }
}

#endif



/////////////////////////////////////>

/*!
************************************************************************
* \brief
*    Fills nalu->buf and nalu->len with the payload of an RTP packet.  
*    Other fields in nalu-> remain uninitialized (will be taken care of 
*    by NALUtoRBSP.
*
* \return
*     4 in case of ok (for compatibility with GetAnnexbNALU)
*     0 if there is nothing any more to read (EOF)
*    -1 in case of any error
*
************************************************************************
*/

int GetRTPNALU (NALU_t *nalu)
{
	RTPpacket_t *p;
	int ret;
	
	if ((p=malloc (sizeof (RTPpacket_t)))== NULL)
		no_mem_exit ("GetRTPNALU-1");
	if ((p->packet=malloc (MAXRTPPACKETSIZE))== NULL)
		no_mem_exit ("GetRTPNALU-2");
	if ((p->payload=malloc (MAXRTPPACKETSIZE))== NULL)
		no_mem_exit ("GetRTPNALU-3");
	
	ret = RTPReadPacket (p, bits);
	nalu->forbidden_bit = 1;
	nalu->len = 0;
	
	if (ret < 0)
		return -1;
	if (ret == 0)
		return 0;
	
	assert (p->paylen < nalu->max_size);
	
	nalu->len = p->paylen;
	memcpy (nalu->buf, p->payload, p->paylen);
	nalu->forbidden_bit = (nalu->buf[0]>>7) & 1;
	nalu->nal_reference_idc = (nalu->buf[0]>>5) & 3;
	nalu->nal_unit_type = (nalu->buf[0]) & 0x1f;
	
	free (p->payload);
	free (p->packet);
	free (p); 
//     printf ("Got an RTP NALU, len %d, first byte %x\n", nalu->len, nalu->buf[0]);
	return nalu->len;
}



/*!
*****************************************************************************
*
* \brief 
*    DecomposeRTPpacket interprets the RTP packet and writes the various
*    structure members of the RTPpacket_t structure
*
* \return
*    0 in case of success
*    negative error code in case of failure
*
* \param p
*    Caller is responsible to allocate enough memory for the generated payload
*    in parameter->payload. Typically a malloc of paclen-12 bytes is sufficient
*
* \par Side effects
*    none
*
* \date
*    30 Spetember 2001
*
* \author
*    Stephan Wenger   stewe@cs.tu-berlin.de
*****************************************************************************/

int DecomposeRTPpacket (RTPpacket_t *p)

{
	// consistency check 
	assert (p->packlen < 65536 - 28);  // IP, UDP headers
	assert (p->packlen >= 12);         // at least a complete RTP header
	assert (p->payload != NULL);
	assert (p->packet != NULL);
	
	// Extract header information
	
	p->v  = p->packet[0] & 0x3;
	p->p  = (p->packet[0] & 0x4) >> 2;
	p->x  = (p->packet[0] & 0x8) >> 3;
	p->cc = (p->packet[0] & 0xf0) >> 4;
	
	p->m  = p->packet[1] & 0x1;
	p->pt = (p->packet[1] & 0xfe) >> 1;
	
	p->seq = p->packet[2] | (p->packet[3] << 8);
	
	memcpy (&p->timestamp, &p->packet[4], 4);// change to shifts for unified byte sex
	memcpy (&p->ssrc, &p->packet[8], 4);// change to shifts for unified byte sex
	
	// header consistency checks
	if (     (p->v != 2)
        || (p->p != 0)
        || (p->x != 0)
        || (p->cc != 0) )
	{
		printf ("DecomposeRTPpacket, RTP header consistency problem, header follows\n");
		DumpRTPHeader (p);
		return -1;
	}
	p->paylen = p->packlen-12;
	memcpy (p->payload, &p->packet[12], p->paylen);
	return 0;
}

/*!
*****************************************************************************
*
* \brief 
*    DumpRTPHeader is a debug tool that dumps a human-readable interpretation
*    of the RTP header
*
* \return
*    n.a.
* \param p
*    the RTP packet to be dumped, after DecompositeRTPpacket()
*
* \par Side effects
*    Debug output to stdout
*
* \date
*    30 Spetember 2001
*
* \author
*    Stephan Wenger   stewe@cs.tu-berlin.de
*****************************************************************************/

void DumpRTPHeader (RTPpacket_t *p)

{
	int i;
	for (i=0; i< 30; i++)
		printf ("%02x ", p->packet[i]);
	printf ("Version (V): %d\n", p->v);
	printf ("Padding (P): %d\n", p->p);
	printf ("Extension (X): %d\n", p->x);
	printf ("CSRC count (CC): %d\n", p->cc);
	printf ("Marker bit (M): %d\n", p->m);
	printf ("Payload Type (PT): %d\n", p->pt);
	printf ("Sequence Number: %d\n", p->seq);
	printf ("Timestamp: %d\n", p->timestamp);
	printf ("SSRC: %d\n", p->ssrc);
}


/*!
*****************************************************************************
*
* \brief 
*    RTPReadPacket reads one packet from file
*
* \return
*    0: EOF
*    negative: error
*    positive: size of RTP packet in bytes
*
* \param p
*    packet data structure, with memory for p->packet allocated
*
* \param bits
*    target file
*
* \par Side effects:
*   - File pointer in bits moved
*   - p->xxx filled by reading and Decomposepacket()
*
* \date
*    04 November, 2001
*
* \author
*    Stephan Wenger, stewe@cs.tu-berlin.de
*****************************************************************************/

int RTPReadPacket (RTPpacket_t *p, FILE *bits)
{
	int Filepos, intime;
	
	assert (p != NULL);
	assert (p->packet != NULL);
	assert (p->payload != NULL);
	
	Filepos = ftell (bits);
	if (4 != fread (&p->packlen,1, 4, bits))
    {
		return 0;
    }
    
	if (4 != fread (&intime, 1, 4, bits))
    {
		fseek (bits, Filepos, SEEK_SET);
		printf ("RTPReadPacket: File corruption, could not read Timestamp, exit\n");
		exit (-1);
    }
	
	assert (p->packlen < MAXRTPPACKETSIZE);
	
	if (p->packlen != fread (p->packet, 1, p->packlen, bits))
    {
		printf ("RTPReadPacket: File corruption, could not read %d bytes\n", p->packlen);
		exit (-1);    // EOF inidication
    }
	
	if (DecomposeRTPpacket (p) < 0)
    {
		// this should never happen, hence exit() is ok.  We probably do not want to attempt
		// to decode a packet that obviously wasn't generated by RTP
		printf ("Errors reported by DecomposePacket(), exit\n");
		exit (-700);
    }
    assert (p->pt == H26LPAYLOADTYPE);
    assert (p->ssrc == 0x12345678);
	return p->packlen;
}






/////////////////////////////////////<