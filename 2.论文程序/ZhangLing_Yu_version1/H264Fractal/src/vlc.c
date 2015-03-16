
/*!
 ***************************************************************************
 * \file vlc.c
 *
 * \brief
 *    (CA)VLC coding functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Lang鴜               <inge.lille-langoy@telenor.com>
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 ***************************************************************************
 */

// #include "contributors.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "elements.h"
#include "vlc.h"

#if TRACE
#define SYMTRACESTRING(s) strncpy(sym.tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // do nothing
#endif

/////////////////////////////////////>

extern void tracebits(const char *trace_str,  int len,  int info,int value1);


int UsedBits;      // for internal statistics, is adjusted by se_v, ue_v, u_1


// Note that all NA values are filled with 0

//! for the linfo_levrun_inter routine
const byte NTAB1[4][8][2] =
{
	{{1,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
	{{1,1},{1,2},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
	{{2,0},{1,3},{1,4},{1,5},{0,0},{0,0},{0,0},{0,0}},
	{{3,0},{2,1},{2,2},{1,6},{1,7},{1,8},{1,9},{4,0}},
};
const byte LEVRUN1[16]=
{
	4,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,
};


const byte NTAB2[4][8][2] =
{
	{{1,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
	{{1,1},{2,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
	{{1,2},{3,0},{4,0},{5,0},{0,0},{0,0},{0,0},{0,0}},
	{{1,3},{1,4},{2,1},{3,1},{6,0},{7,0},{8,0},{9,0}},
};

//! for the linfo_levrun__c2x2 routine
const byte LEVRUN3[4] =
{
	2,1,0,0
};
const byte NTAB3[2][2][2] =
{
	{{1,0},{0,0}},
	{{2,0},{1,1}},
};


/////////////////////////////////////<

/*! 
 *************************************************************************************
 * \brief
 *    ue_v, writes an ue(v) syntax element, returns the length in bits
 *
 * \param tracestring
 *    the string for the trace file
 * \param value
 *    the value to be coded
 *  \param part
 *    the Data Partition the value should be coded into
 *
 * \return
 *    Number of bits used by the coded syntax element
 *
 * \ note
 *    This function writes always the bit buffer for the progressive scan flag, and
 *    should not be used (or should be modified appropriately) for the interlace crap
 *    When used in the context of the Parameter Sets, this is obviously not a
 *    problem.
 *
 *************************************************************************************
 */
int ue_v (char *tracestring, int value, DataPartition *part)
{
  SyntaxElement symbol, *sym=&symbol;
  sym->type = SE_HEADER;
  sym->mapping = ue_linfo;               // Mapping rule: unsigned integer
  sym->value1 = value;
  sym->value2 = 0;
#if TRACE
  strncpy(sym->tracestring,tracestring,TRACESTRING_SIZE);
#endif
  assert (part->bitstream->streamBuffer != NULL);
  return writeSyntaxElement_UVLC (sym, part);
}


/*! 
 *************************************************************************************
 * \brief
 *    se_v, writes an se(v) syntax element, returns the length in bits
 *
 * \param tracestring
 *    the string for the trace file
 * \param value
 *    the value to be coded
 *  \param part
 *    the Data Partition the value should be coded into
 *
 * \return
 *    Number of bits used by the coded syntax element
 *
 * \ note
 *    This function writes always the bit buffer for the progressive scan flag, and
 *    should not be used (or should be modified appropriately) for the interlace crap
 *    When used in the context of the Parameter Sets, this is obviously not a
 *    problem.
 *
 *************************************************************************************
 */
int se_v (char *tracestring, int value, DataPartition *part)
{
  SyntaxElement symbol, *sym=&symbol;
  sym->type = SE_HEADER;
  sym->mapping = se_linfo;               // Mapping rule: signed integer
  sym->value1 = value;
  sym->value2 = 0;
#if TRACE
  strncpy(sym->tracestring,tracestring,TRACESTRING_SIZE);
#endif
  assert (part->bitstream->streamBuffer != NULL);
  return writeSyntaxElement_UVLC (sym, part);
}


/*! 
 *************************************************************************************
 * \brief
 *    u_1, writes a flag u(1) syntax element, returns the length in bits, 
 *    always 1
 *
 * \param tracestring
 *    the string for the trace file
 * \param value
 *    the value to be coded
 *  \param part
 *    the Data Partition the value should be coded into
 *
 * \return
 *    Number of bits used by the coded syntax element (always 1)
 *
 * \ note
 *    This function writes always the bit buffer for the progressive scan flag, and
 *    should not be used (or should be modified appropriately) for the interlace crap
 *    When used in the context of the Parameter Sets, this is obviously not a
 *    problem.
 *
 *************************************************************************************
 */
int u_1 (char *tracestring, int value, DataPartition *part)
{
  SyntaxElement symbol, *sym=&symbol;

  sym->bitpattern = value;
  sym->len = 1;
  sym->type = SE_HEADER;
  sym->value1 = value;
  sym->value2 = 0;
#if TRACE
  strncpy(sym->tracestring,tracestring,TRACESTRING_SIZE);
#endif
  assert (part->bitstream->streamBuffer != NULL);
  return writeSyntaxElement_fixed(sym, part);
}


/*! 
 *************************************************************************************
 * \brief
 *    u_v, writes a a n bit fixed length syntax element, returns the length in bits, 
 *
 * \param n
 *    length in bits
 * \param tracestring
 *    the string for the trace file
 * \param value
 *    the value to be coded
 *  \param part
 *    the Data Partition the value should be coded into
 *
 * \return
 *    Number of bits used by the coded syntax element 
 *
 * \ note
 *    This function writes always the bit buffer for the progressive scan flag, and
 *    should not be used (or should be modified appropriately) for the interlace crap
 *    When used in the context of the Parameter Sets, this is obviously not a
 *    problem.
 *
 *************************************************************************************
 */

int u_v (int n, char *tracestring, int value, DataPartition *part)
{
  SyntaxElement symbol, *sym=&symbol;

  sym->bitpattern = value;
  sym->len = n;
  sym->type = SE_HEADER;
  sym->value1 = value;
  sym->value2 = 0;
#if TRACE
  strncpy(sym->tracestring,tracestring,TRACESTRING_SIZE);
#endif
  assert (part->bitstream->streamBuffer != NULL);
  return writeSyntaxElement_fixed(sym, part);
}


/*!
 ************************************************************************
 * \brief
 *    mapping for ue(v) syntax elements
 * \param ue
 *    value to be mapped
 * \param dummy
 *    dummy parameter
 * \param info
 *    returns mapped value
 * \param len
 *    returns mapped value length
 ************************************************************************
 */
void ue_linfo(int ue, int dummy, int *len,int *info)
{
  int i,nn;

  nn=(ue+1)/2;

  for (i=0; i < 16 && nn != 0; i++)
  {
    nn /= 2;
  }
  *len= 2*i + 1;
  *info=ue+1-(int)pow(2,i);
}


/*!
 ************************************************************************
 * \brief
 *    mapping for se(v) syntax elements//函数说明：有符号0阶指数哥伦布编码
 * \param se
 *    value to be mapped//待映射的值
 * \param dummy
 *    dummy argument//虚拟参数
 * \param len
 *    returns mapped value length//返回映射值的长度
 * \param info
 *    returns mapped value//返回映射的值
 ************************************************************************
 */
void se_linfo(int se, int dummy, int *len,int *info)
{

  int i,n,sign,nn;

  sign=0;

  if (se <= 0)
  {
    sign=1;
  }
  n=abs(se) << 1;

  /*
  n+1 is the number in the code table.  Based on this we find length and info
  */

  nn=n/2;
  for (i=0; i < 16 && nn != 0; i++)
  {
    nn /= 2;
  }
  *len=i*2 + 1;
  *info=n - (int)pow(2,i) + sign;
}


/*!
 ************************************************************************
 * \par Input:
 *    Number in the code table
 * \par Output:
 *    length and info
 ************************************************************************
 */
void cbp_linfo_intra(int cbp, int dummy, int *len,int *info)
{
  extern const int NCBP[48][2];
  ue_linfo(NCBP[cbp][0], dummy, len, info);
}


/*!
 ************************************************************************
 * \par Input:
 *    Number in the code table
 * \par Output:
 *    length and info
 ************************************************************************
 */
void cbp_linfo_inter(int cbp, int dummy, int *len,int *info)
{
  extern const int NCBP[48][2];
  ue_linfo(NCBP[cbp][1], dummy, len, info);
}


/*!
 ************************************************************************
 * \brief
 *    2x2 transform of chroma DC
 * \par Input:
 *    level and run for coefficients
 * \par Output:
 *    length and info
 * \note
 *    see ITU document for bit assignment
 ************************************************************************
 */
void levrun_linfo_c2x2(int level,int run,int *len,int *info)
{
  const int NTAB[2][2]=
  {
    {1,5},
    {3,0}
  };
  const int LEVRUN[4]=
  {
    2,1,0,0
  };

  int levabs,i,n,sign,nn;

  if (level == 0) //  check if the coefficient sign EOB (level=0)
  {
    *len=1;
    return;
  }
  sign=0;
  if (level <= 0)
  {
    sign=1;
  }
  levabs=abs(level);
  if (levabs <= LEVRUN[run])
  {
    n=NTAB[levabs-1][run]+1;
  }
  else
  {
    n=(levabs-LEVRUN[run])*8 + run*2;
  }

  nn=n/2;

  for (i=0; i < 16 && nn != 0; i++)
  {
    nn /= 2;
  }
  *len= 2*i + 1;
  *info=n-(int)pow(2,i)+sign;
}


/*!
 ************************************************************************
 * \brief
 *    Single scan coefficients
 * \par Input:
 *    level and run for coefficiets
 * \par Output:
 *    lenght and info
 * \note
 *    see ITU document for bit assignment
 ************************************************************************
 */
void levrun_linfo_inter(int level,int run,int *len,int *info)
{
  const byte LEVRUN[16]=
  {
    4,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0
  };
  const byte NTAB[4][10]=
  {
    { 1, 3, 5, 9,11,13,21,23,25,27},
    { 7,17,19, 0, 0, 0, 0, 0, 0, 0},
    {15, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {29, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  };

  int levabs,i,n,sign,nn;

  if (level == 0)           //  check for EOB
  {
    *len=1;
    return;
  }

  if (level <= 0)
    sign=1;
  else
    sign=0;

  levabs=abs(level);
  if (levabs <= LEVRUN[run])
  {
    n=NTAB[levabs-1][run]+1;
  }
  else
  {
    n=(levabs-LEVRUN[run])*32 + run*2;
  }

  nn=n/2;

  for (i=0; i < 16 && nn != 0; i++)
  {
    nn /= 2;
  }
  *len= 2*i + 1;
  *info=n-(int)pow(2,i)+sign;

}


/*!
 ************************************************************************
 * \brief
 *    Double scan coefficients
 * \par Input:
 *    level and run for coefficiets
 * \par Output:
 *    lenght and info
 * \note
 *    see ITU document for bit assignment
 ************************************************************************
 */
void levrun_linfo_intra(int level,int run,int *len,int *info)
{
  const byte LEVRUN[8]=
  {
    9,3,1,1,1,0,0,0
  };

  const byte NTAB[9][5] =
  {
    { 1, 3, 7,15,17},
    { 5,19, 0, 0, 0},
    { 9,21, 0, 0, 0},
    {11, 0, 0, 0, 0},
    {13, 0, 0, 0, 0},
    {23, 0, 0, 0, 0},
    {25, 0, 0, 0, 0},
    {27, 0, 0, 0, 0},
    {29, 0, 0, 0, 0},
  };

  int levabs,i,n,sign,nn;

  if (level == 0)     //  check for EOB
  {
    *len=1;
    return;
  }
  if (level <= 0)
    sign=1;
  else
    sign=0;

  levabs=abs(level);
  if (levabs <= LEVRUN[run])
  {
    n=NTAB[levabs-1][run]+1;
  }
  else
  {
    n=(levabs-LEVRUN[run])*16 + 16 + run*2;
  }

  nn=n/2;

  for (i=0; i < 16 && nn != 0; i++)
  {
    nn /= 2;
  }
  *len= 2*i + 1;
  *info=n-(int)pow(2,i)+sign;
}


/*!
 ************************************************************************
 * \brief
 *    Makes code word and passes it back//函数说明：生成码字，并传回
 *    A code word has the following format: 0 0 0 ... 1 Xn ...X2 X1 X0.
 *
 * \par Input:
 *    Info   : Xn..X2 X1 X0                                             \n
 *    Length : Total number of bits in the codeword
 ************************************************************************
 */
 // NOTE this function is called with sym->inf > (1<<(sym->len/2)).  The upper bits of inf are junk
int symbol2uvlc(SyntaxElement *sym)
{
  int suffix_len=sym->len/2;  
  sym->bitpattern = (1<<suffix_len)|(sym->inf&((1<<suffix_len)-1));
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    generates UVLC code and passes the codeword to the buffer//生成UVLC码，将码字传输到比特流
 ************************************************************************
 */
int writeSyntaxElement_UVLC(SyntaxElement *se, DataPartition *this_dataPart)
{
  se->mapping(se->value1,se->value2,&(se->len),&(se->inf));//语法元素到UVLC的映射
  symbol2uvlc(se);//生成码字

  writeUVLC2buffer(se, this_dataPart->bitstream);//将ULVC码字写入比特流

  if(se->type != SE_HEADER)//如果不是头
    this_dataPart->bitstream->write_flag = 1;//写标志位=1

#if TRACE
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    passes the fixed codeword to the buffer
 ************************************************************************
 */
int writeSyntaxElement_fixed(SyntaxElement *se, DataPartition *this_dataPart)
{  
  writeUVLC2buffer(se, this_dataPart->bitstream);
  
  if(se->type != SE_HEADER)
    this_dataPart->bitstream->write_flag = 1;

#if TRACE
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    generates code and passes the codeword to the buffer
 ************************************************************************
 */
int writeSyntaxElement_Intra4x4PredictionMode(SyntaxElement *se, DataPartition *this_dataPart)
{

  if (se->value1 == -1)
  {
    se->len = 1;
    se->inf = 1;
  }
  else 
  {
    se->len = 4;  
    se->inf = se->value1;
  }

  se->bitpattern = se->inf;
  writeUVLC2buffer(se, this_dataPart->bitstream);

  if(se->type != SE_HEADER)
    this_dataPart->bitstream->write_flag = 1;

#if TRACE
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    generates UVLC code and passes the codeword to the buffer
 * \author
 *  Tian Dong
 ************************************************************************
 */
int writeSyntaxElement2Buf_UVLC(SyntaxElement *se, Bitstream* this_streamBuffer )
{

  se->mapping(se->value1,se->value2,&(se->len),&(se->inf));

  symbol2uvlc(se);

  writeUVLC2buffer(se, this_streamBuffer );

#if TRACE
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    writes UVLC code to the appropriate buffer够了一个字节就进行输出
放到byte_buf
 ************************************************************************
 */
void  writeUVLC2buffer(SyntaxElement *se, Bitstream *currStream)
{

  int i;
  unsigned int mask = 1 << (se->len-1);

  // Add the new bits to the bitstream.
  // Write out a byte if it is full
  for (i=0; i<se->len; i++)
  {
    currStream->byte_buf <<= 1;
    if (se->bitpattern & mask)
      currStream->byte_buf |= 1;
    currStream->bits_to_go--;
    mask >>= 1;
    if (currStream->bits_to_go==0)
    {
      currStream->bits_to_go = 8;
      currStream->streamBuffer[currStream->byte_pos++]=currStream->byte_buf;
      currStream->byte_buf = 0;
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    generates UVLC code and passes the codeword to the buffer
 * \author
 *  Tian Dong
 ************************************************************************
 */
int writeSyntaxElement2Buf_Fixed(SyntaxElement *se, Bitstream* this_streamBuffer )
{

  writeUVLC2buffer(se, this_streamBuffer );

#if TRACE
  if(se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    Makes code word and passes it back//函数说明：生成码字并传回
 *
 * \par Input:
 *    Info   : Xn..X2 X1 X0                                             \n
 *    Length : Total number of bits in the codeword
 ************************************************************************
 */

int symbol2vlc(SyntaxElement *sym)
{
  int info_len = sym->len;

  // Convert info into a bitpattern int
  sym->bitpattern = 0;

  // vlc coding
  while(--info_len >= 0)
  {
    sym->bitpattern <<= 1;
    sym->bitpattern |= (0x01 & (sym->inf >> info_len));
  }
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    generates VLC code and passes the codeword to the buffer//生成VLC码字，并传输到比特流
 ************************************************************************
 */
int writeSyntaxElement_VLC(SyntaxElement *se, DataPartition *this_dataPart)
{

  se->inf = se->value1;
  se->len = se->value2;
  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for NumCoeff and TrailingOnes
 ************************************************************************
 */

int writeSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *se, DataPartition *this_dataPart)
{
  int lentab[3][4][17] = 
  {
    {   // 0702
      { 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},//对应trailingones=0,totalcoeff=1-16的长度值
      { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
      { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
      { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16},//对应trailingones=3,totalcoeff=1-16的长度值
    },                                                 
    {                                                  
      { 2, 6, 6, 7, 8, 8, 9,11,11,12,12,12,13,13,13,14,14},
      { 0, 2, 5, 6, 6, 7, 8, 9,11,11,12,12,13,13,14,14,14},
      { 0, 0, 3, 6, 6, 7, 8, 9,11,11,12,12,13,13,13,14,14},
      { 0, 0, 0, 4, 4, 5, 6, 6, 7, 9,11,11,12,13,13,13,14},
    },                                                 
    {                                                  
      { 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9,10,10,10,10},
      { 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10},
      { 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,10},
      { 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9,10,10,10},
    },

  };

  int codtab[3][4][17] = 
  {
    {
      { 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7,4}, //对应trailingones=0,totalcoeff=1-16的长度值
      { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10,6}, 
      { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9,5}, 
      { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12,8},//对应trailingones=3,totalcoeff=1-16的长度值
    },
    {
      { 3,11, 7, 7, 7, 4, 7,15,11,15,11, 8,15,11, 7, 9,7}, 
      { 0, 2, 7,10, 6, 6, 6, 6,14,10,14,10,14,10,11, 8,6}, 
      { 0, 0, 3, 9, 5, 5, 5, 5,13, 9,13, 9,13, 9, 6,10,5}, 
      { 0, 0, 0, 5, 4, 6, 8, 4, 4, 4,12, 8,12,12, 8, 1,4},
    },
    {
      {15,15,11, 8,15,11, 9, 8,15,11,15,11, 8,13, 9, 5,1}, 
      { 0,14,15,12,10, 8,14,10,14,14,10,14,10, 7,12, 8,4},
      { 0, 0,13,14,11, 9,13, 9,13,10,13, 9,13, 9,11, 7,3},
      { 0, 0, 0,12,11,10, 9, 8,13,12,12,12, 8,12,10, 6,2},
    },
  };
  int vlcnum;

  vlcnum = se->len;

  // se->value1 : numcoeff
  // se->value2 : numtrailingones

  if (vlcnum == 3)
  {
    se->len = 6;  // 4 + 2 bit FLC
    if (se->value1 > 0)
    {
      se->inf = ((se->value1-1) << 2) | se->value2;
    }
    else
    {
      se->inf = 3;
    }
  }
  else
  {
    se->len = lentab[vlcnum][se->value2][se->value1];//16
    se->inf = codtab[vlcnum][se->value2][se->value1];//10
  }
  //se->inf = 0;

  if (se->len == 0)
  {
    printf("ERROR: (numcoeff,trailingones) not valid: vlc=%d (%d, %d)\n", 
      vlcnum, se->value1, se->value2);
    exit(-1);
  }

  symbol2vlc(se);//生成VLC码字

  writeUVLC2buffer(se, this_dataPart->bitstream);//将码字写入比特流
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for NumCoeff and TrailingOnes for Chroma DC
 ************************************************************************
 */
int writeSyntaxElement_NumCoeffTrailingOnesChromaDC(SyntaxElement *se, DataPartition *this_dataPart)
{
  int lentab[4][5] = 
  {
    { 2, 6, 6, 6, 6,},          
    { 0, 1, 6, 7, 8,}, 
    { 0, 0, 3, 7, 8,}, 
    { 0, 0, 0, 6, 7,},
  };

  int codtab[4][5] = 
  {
    {1,7,4,3,2},
    {0,1,6,3,3},
    {0,0,1,2,2},
    {0,0,0,5,0},
  };

  // se->value1 : numcoeff
  // se->value2 : numtrailingones
  se->len = lentab[se->value2][se->value1];
  se->inf = codtab[se->value2][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (numcoeff,trailingones) not valid: (%d, %d)\n", 
      se->value1, se->value2);
    exit(-1);
  }

  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for TotalZeros
 ************************************************************************
 */
int writeSyntaxElement_TotalZeros(SyntaxElement *se, DataPartition *this_dataPart)
{
  int lentab[TOTRUN_NUM][16] = 
  {
	  { 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9}, //对应totalcoeff=1,total_zero=0-15
    { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},  
    { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},  
    { 5,3,4,4,3,3,3,4,3,4,5,5,5},  
    { 4,4,4,3,3,3,3,3,4,5,4,5},  
    { 6,5,3,3,3,3,3,3,4,3,6},  
    { 6,5,3,3,3,2,3,4,3,6},  
    { 6,4,5,3,2,2,3,3,6},  
    { 6,6,4,2,2,3,2,5},  
    { 5,5,3,2,2,2,4},  
    { 4,4,3,3,1,3},  
    { 4,4,2,1,3},  
    { 3,3,1,2},  
    { 2,2,1},  
	  { 1,1},  //对应totalcoeff=15,total_zero=0-1
  };

  int codtab[TOTRUN_NUM][16] = 
  {
    {1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
    {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
    {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
    {3,7,5,4,6,5,4,3,3,2,2,1,0},
    {5,4,3,7,6,5,4,3,2,1,1,0},
    {1,1,7,6,5,4,3,2,1,1,0},
    {1,1,5,4,3,3,2,1,1,0},
    {1,1,1,3,3,2,2,1,0},
    {1,0,1,3,2,1,1,1,},
    {1,0,1,3,2,1,1,},
    {0,1,1,2,1,3},
    {0,1,1,1,1},
    {0,1,1,1},
    {0,1,1},
    {0,1},  
  };
  int vlcnum;

  vlcnum = se->len;

  // se->value1 : TotalZeros
  se->len = lentab[vlcnum][se->value1];
  se->inf = codtab[vlcnum][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (TotalZeros) not valid: (%d)\n",se->value1);
    exit(-1);
  }

  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for TotalZeros for Chroma DC
 ************************************************************************
 */
int writeSyntaxElement_TotalZerosChromaDC(SyntaxElement *se, DataPartition *this_dataPart)
{
  int lentab[3][4] = 
  {
    { 1, 2, 3, 3,},
    { 1, 2, 2, 0,},
    { 1, 1, 0, 0,}, 
  };

  int codtab[3][4] = 
  {
    { 1, 1, 1, 0,},
    { 1, 1, 0, 0,},
    { 1, 0, 0, 0,},
  };
  int vlcnum;

  vlcnum = se->len;

  // se->value1 : TotalZeros
  se->len = lentab[vlcnum][se->value1];
  se->inf = codtab[vlcnum][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (TotalZeros) not valid: (%d)\n",se->value1);
    exit(-1);
  }

  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for Run Before Next Coefficient, VLC0
 ************************************************************************
 */
int writeSyntaxElement_Run(SyntaxElement *se, DataPartition *this_dataPart)
{
  int lentab[TOTRUN_NUM][16] = 
  {
    {1,1},
    {1,2,2},
    {2,2,2,2},
    {2,2,2,3,3},
    {2,2,3,3,3,3},
    {2,3,3,3,3,3,3},
    {3,3,3,3,3,3,3,4,5,6,7,8,9,10,11},
  };

  int codtab[TOTRUN_NUM][16] = 
  {
    {1,0},
    {1,1,0},
    {3,2,1,0},
    {3,2,1,1,0},
    {3,2,3,2,1,0},
    {3,0,1,3,2,5,4},
    {7,6,5,4,3,2,1,1,1,1,1,1,1,1,1},
  };
  int vlcnum;

  vlcnum = se->len;

  // se->value1 : run
  se->len = lentab[vlcnum][se->value1];
  se->inf = codtab[vlcnum][se->value1];

  if (se->len == 0)
  {
    printf("ERROR: (run) not valid: (%d)\n",se->value1);
    exit(-1);
  }

  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for Coeff Level (VLC1)
 ************************************************************************
 */
int writeSyntaxElement_Level_VLC1(SyntaxElement *se, DataPartition *this_dataPart)
{
  int level, levabs, sign;

  level = se->value1;
  levabs = abs(level);
  sign = (level < 0 ? 1 : 0);

  
  if (levabs < 8)
  {
    se->len = levabs * 2 + sign - 1;
    se->inf = 1;
  }
  else if (levabs < 8+8)
  {
    // escape code1
    se->len = 14 + 1 + 4;
    se->inf = (1 << 4) | ((levabs - 8) << 1) | sign;
  }
  else
  {
    // escape code2
    se->len = 14 + 2 + 12;
    se->inf = (0x1 << 12) | ((levabs - 16)<< 1) | sign;
  }


  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    write VLC for Coeff Level
 ************************************************************************
 */
int writeSyntaxElement_Level_VLCN(SyntaxElement *se, int vlc, DataPartition *this_dataPart)
{
  int iCodeword;
  int iLength;

  int level = se->value1;

  int levabs = abs(level);
  int sign = (level < 0 ? 1 : 0);  

  int shift = vlc-1;
  int escape = (15<<shift)+1;

  int numPrefix = (levabs-1)>>shift;

  int sufmask = ~((0xffffffff)<<shift);
  int suffix = (levabs-1)&sufmask;

  if (levabs < escape)
  {
    iLength = numPrefix + vlc + 1;
    iCodeword = (1<<(shift+1))|(suffix<<1)|sign;
  }
  else
  {
    iLength = 28;
    iCodeword = (1<<12)|((levabs-escape)<<1)|sign;
  }
  se->len = iLength;
  se->inf = iCodeword;

  symbol2vlc(se);

  writeUVLC2buffer(se, this_dataPart->bitstream);
#if TRACE
  if (se->type <= 1)
    trace2out (se);
#endif

  return (se->len);
}


/*!
 ************************************************************************
 * \brief
 *    Write out a trace string on the trace file
 ************************************************************************
 */
#if TRACE
void trace2out(SyntaxElement *sym)
{
  static int bitcounter = 0;
  int i, chars;

  if (p_trace != NULL)
  {
    putc('@', p_trace);
    chars = fprintf(p_trace, "%i", bitcounter);
    while(chars++ < 6)
      putc(' ',p_trace);

    chars += fprintf(p_trace, "%s", sym->tracestring);
    while(chars++ < 55)
      putc(' ',p_trace);

  // Align bitpattern
    if(sym->len<15)
    {
      for(i=0 ; i<15-sym->len ; i++)
        fputc(' ', p_trace);
    }
    
    // Print bitpattern
    bitcounter += sym->len;
    for(i=1 ; i<=sym->len ; i++)
    {
      if((sym->bitpattern >> (sym->len-i)) & 0x1)
        fputc('1', p_trace);
      else
        fputc('0', p_trace);
    }
    fprintf(p_trace, " (%3d) \n",sym->value1);
  }
  fflush (p_trace);
}
#endif


/*!
 ************************************************************************
 * \brief
 *    puts the less than 8 bits in the byte buffer of the Bitstream into
 *    the streamBuffer.  
 *
 * \param
 *   currStream: the Bitstream the alignment should be established
 *
 ************************************************************************
 */
void writeVlcByteAlign(Bitstream* currStream)
{
  if (currStream->bits_to_go < 8)//1
  { // trailing bits to process
    currStream->byte_buf = (currStream->byte_buf <<currStream->bits_to_go) | (0xff >> (8 - currStream->bits_to_go));
    stat->bit_use_stuffingBits[img->type]+=currStream->bits_to_go;
    currStream->streamBuffer[currStream->byte_pos++]=currStream->byte_buf;
    currStream->bits_to_go = 8;
  }
}





/////////////////////////////////////>

/*! 
 *************************************************************************************
 * \brief
 *    ue_v, reads an ue(v) syntax element, the length in bits is stored in 
 *    the global UsedBits variable
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int ue_v_dec (char *tracestring, Bitstream *bitstream)
{
  SyntaxElement symbol, *sym=&symbol;

  assert (bitstream->streamBuffer != NULL);
  sym->type = SE_HEADER;
  sym->mapping_dec = linfo_ue_dec;   // Mapping rule
  SYMTRACESTRING(tracestring);
  readSyntaxElement_VLC_dec (sym, bitstream);// // // //sym没有赋值.bitstream is currstream
  UsedBits+=sym->len;//25=24+1//26=25+1//27=26+1
  return sym->value1;//0//0//0
}


/*! 
 *************************************************************************************
 * \brief
 *    ue_v, reads an se(v) syntax element, the length in bits is stored in 
 *    the global UsedBits variable
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int se_v_dec (char *tracestring, Bitstream *bitstream)
{
  SyntaxElement symbol, *sym=&symbol;

  assert (bitstream->streamBuffer != NULL);
  sym->type = SE_HEADER;
  sym->mapping_dec = linfo_se_dec;   // Mapping rule: signed integer
  SYMTRACESTRING(tracestring);
  readSyntaxElement_VLC_dec (sym, bitstream);
  UsedBits+=sym->len;//28=19+9
  return sym->value1;//-9还是这里不一样啊
}


/*! 
 *************************************************************************************
 * \brief
 *    ue_v, reads an u(v) syntax element, the length in bits is stored in 
 *    the global UsedBits variable
 *
 * \param LenInBits
 *    length of the syntax element
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int u_v_dec (int LenInBits, char*tracestring, Bitstream *bitstream)
{
  SyntaxElement symbol, *sym=&symbol;

  assert (bitstream->streamBuffer != NULL);
  sym->type = SE_HEADER;//0
  sym->mapping_dec = linfo_ue_dec;   // Mapping rule
  sym->len = LenInBits;//4
  SYMTRACESTRING(tracestring);
  readSyntaxElement_FLC_dec (sym, bitstream);
  UsedBits+=sym->len;//18=14+4
  return sym->inf;//0
};

                
/*! 
 *************************************************************************************
 * \brief
 *    ue_v, reads an u(1) syntax element, the length in bits is stored in 
 *    the global UsedBits variable
 *
 * \param tracestring
 *    the string for the trace file
 *
 * \param bitstream
 *    the stream to be read from
 *
 * \return
 *    the value of the coded syntax element
 *
 *************************************************************************************
 */
int u_1_dec (char *tracestring, Bitstream *bitstream)
{
  return u_v_dec (1, tracestring, bitstream);
}



/*!
 ************************************************************************
 * \brief
 *    mapping rule for ue(v) syntax elements
 * \par Input:
 *    lenght and info
 * \par Output:
 *    number in the code table
 ************************************************************************
 */
void linfo_ue_dec(int len, int info, int *value1, int *dummy)
{
  *value1 = (int)pow(2,(len/2))+info-1; // *value1 = (int)(2<<(len>>1))+info-1;
}

/*!
 ************************************************************************
 * \brief
 *    mapping rule for se(v) syntax elements
 * \par Input:
 *    lenght and info
 * \par Output:
 *    signed mvd
 ************************************************************************
 */
void linfo_se_dec(int len,  int info, int *value1, int *dummy)
{
  int n;
  n = (int)pow(2,(len/2))+info-1;
  *value1 = (n+1)/2;
  if((n & 0x01)==0)                           // lsb is signed bit
    *value1 = -*value1;
}


/*!
 ************************************************************************
 * \par Input:
 *    lenght and info
 * \par Output:
 *    cbp (intra)
 ************************************************************************
 */
void linfo_cbp_intra_dec(int len,int info,int *cbp, int *dummy)
{
  extern const byte NCBP_dec[48][2];
    int cbp_idx;
  linfo_ue_dec(len,info,&cbp_idx,dummy);
    *cbp=NCBP_dec[cbp_idx][0];
}

/*!
 ************************************************************************
 * \par Input:
 *    lenght and info
 * \par Output:
 *    cbp (inter)
 ************************************************************************
 */
void linfo_cbp_inter_dec(int len,int info,int *cbp, int *dummy)
{
  extern const byte NCBP_dec[48][2];
  int cbp_idx;
  linfo_ue_dec(len,info,&cbp_idx,dummy);
    *cbp=NCBP_dec[cbp_idx][1];
}

/*!
 ************************************************************************
 * \par Input:
 *    lenght and info
 * \par Output:
 *    level, run
 ************************************************************************
 */
void linfo_levrun_inter_dec(int len, int info, int *level, int *irun)
{
  int l2;
  int inf;
  if (len<=9)
  {
    l2=max(0,len/2-1);
    inf=info/2;
    *level=NTAB1[l2][inf][0];
    *irun=NTAB1[l2][inf][1];
    if ((info&0x01)==1)
      *level=-*level;                   // make sign
  }
  else                                  // if len > 9, skip using the array
  {
    *irun=(info&0x1e)>>1;
    *level = LEVRUN1[*irun] + info/32 + (int)pow(2,len/2 - 5);
    if ((info&0x01)==1)
      *level=-*level;
  }
    if (len == 1) // EOB
        *level = 0;
}


/*!
 ************************************************************************
 * \par Input:
 *    lenght and info
 * \par Output:
 *    level, run
 ************************************************************************
 */
void linfo_levrun_c2x2_dec(int len, int info, int *level, int *irun)
{
  int l2;
  int inf;

  if (len<=5)
  {
    l2=max(0,len/2-1);
    inf=info/2;
    *level=NTAB3[l2][inf][0];
    *irun=NTAB3[l2][inf][1];
    if ((info&0x01)==1)
      *level=-*level;                 // make sign
  }
  else                                  // if len > 5, skip using the array
  {
    *irun=(info&0x06)>>1;
    *level = LEVRUN3[*irun] + info/8 + (int)pow(2,len/2 - 3);
    if ((info&0x01)==1)
      *level=-*level;
  }
  if (len == 1) // EOB
    *level = 0;
}

/*!
 ************************************************************************
 * \brief
 *    read next UVLC codeword from UVLC-partition and
 *    map it to the corresponding syntax element
 ************************************************************************
 */
int readSyntaxElement_VLC_dec(SyntaxElement *sym, Bitstream *currStream)
{
  int frame_bitoffset = currStream->frame_bitoffset;//24//25//26//27//19
  byte *buf = currStream->streamBuffer;//77
  int BitstreamLengthInBytes = currStream->bitstream_length;//8//8//8//8//7324

  sym->len =  GetVLCSymbol_dec (buf, frame_bitoffset, &(sym->inf), BitstreamLengthInBytes);//1//1//1//3//9
  //sym->inf的值: 0，0,查buf的定义
  if (sym->len == -1)
    return -1;
  currStream->frame_bitoffset += sym->len;//25=24+1//26=25+1//27=26+1//28=19+9
  sym->mapping_dec(sym->len,sym->inf,&(sym->value1),&(sym->value2));

#if TRACE
  tracebits(sym->tracestring, sym->len, sym->inf, sym->value1);
#endif

  return 1;
}


/*!
 ************************************************************************
 * \brief
 *    read next UVLC codeword from UVLC-partition and
 *    map it to the corresponding syntax element
 ************************************************************************
 */
int readSyntaxElement_UVLC_dec(SyntaxElement *sym, struct img_par *img, struct inp_par_dec *inp, struct datapartition *dP)
{
  Bitstream   *currStream = dP->bitstream;

  return (readSyntaxElement_VLC_dec(sym, currStream));
}

/*!
 ************************************************************************
 * \brief
 *    read next VLC codeword for 4x4 Intra Prediction Mode and
 *    map it to the corresponding Intra Prediction Direction
 ************************************************************************
 */
int readSyntaxElement_Intra4x4PredictionMode_dec(SyntaxElement *sym, struct img_par *img, struct inp_par_dec *inp, struct datapartition *dP)
{
  Bitstream   *currStream            = dP->bitstream;
  int         frame_bitoffset        = currStream->frame_bitoffset;
  byte        *buf                   = currStream->streamBuffer;
  int         BitstreamLengthInBytes = currStream->bitstream_length;

  sym->len = GetVLCSymbol_IntraMode_dec (buf, frame_bitoffset, &(sym->inf), BitstreamLengthInBytes);

  if (sym->len == -1)
    return -1;

  currStream->frame_bitoffset += sym->len;
  sym->value1                  = sym->len == 1 ? -1 : sym->inf;

#if TRACE
  tracebits2(sym->tracestring, sym->len, sym->value1);
#endif

  return 1;
}

int GetVLCSymbol_IntraMode_dec (byte buffer[],int totbitoffset,int *info, int bytecount)
{

  register int inf;
  long byteoffset;      // byte from start of buffer
  int bitoffset;      // bit from start of byte
  int ctr_bit=0;      // control bit for current bit posision
  int bitcounter=1;
  int len;
  int info_bit;

  byteoffset = totbitoffset/8;
  bitoffset  = 7-(totbitoffset%8);
  ctr_bit    = (buffer[byteoffset] & (0x01<<bitoffset));   // set up control bit
  len        = 1;

  //First bit
  if (ctr_bit)
  {
    *info = 0;
    return bitcounter;
  }
  else
    len=4;

  // make infoword
  inf=0;                          // shortest possible code is 1, then info is always 0
  for(info_bit=0;(info_bit<(len-1)); info_bit++)
  {
    bitcounter++;
    bitoffset-=1;
    if (bitoffset<0)
    {                 // finished with current byte ?
      bitoffset=bitoffset+8;
      byteoffset++;
    }
    if (byteoffset > bytecount)
    {
      return -1;
    }
    inf=(inf<<1);
    if(buffer[byteoffset] & (0x01<<(bitoffset)))
      inf |=1;
  }

  *info = inf;
  return bitcounter;           // return absolute offset in bit from start of frame
}


/*!
 ************************************************************************
 * \brief
 *    test if bit buffer contains only stop bit
 *
 * \param buffer
 *    buffer containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param bytecount
 *    buffer length
 * \return
 *    true if more bits available
 ************************************************************************
 */
int more_rbsp_data_dec (byte buffer[],int totbitoffset,int bytecount)
{

  long byteoffset;      // byte from start of buffer
  int bitoffset;      // bit from start of byte
  int ctr_bit=0;      // control bit for current bit posision

  int cnt=0;

  
  byteoffset= totbitoffset/8;
  bitoffset= 7-(totbitoffset%8);

  assert (byteoffset<bytecount);

  // there is more until we're in the last byte
  if (byteoffset<(bytecount-1)) return TRUE1;

  // read one bit
  ctr_bit = (buffer[byteoffset] & (0x01<<bitoffset));
  
  // a stop bit has to be one
  if (ctr_bit==0) return TRUE1;

  bitoffset--;

  while (bitoffset>=0)
  {
    ctr_bit = (buffer[byteoffset] & (0x01<<bitoffset));   // set up control bit
    if (ctr_bit>0) cnt++;
    bitoffset--;
  }

  return (0!=cnt);

}


/*!
 ************************************************************************
 * \brief
 *    Check if there are symbols for the next MB
 ************************************************************************
 */
int uvlc_startcode_follows_dec(struct img_par *img, struct inp_par_dec *inp, int dummy)
{
  int dp_Nr = assignSE2partition_dec[img->currentSlice->dp_mode][SE_MBTYPE];
  DataPartition *dP = &(img->currentSlice->partArr[dp_Nr]);
  Bitstream   *currStream = dP->bitstream;
  byte *buf = currStream->streamBuffer;

  //KS: new function test for End of Buffer
  return (!(more_rbsp_data_dec(buf, currStream->frame_bitoffset,currStream->bitstream_length)));
}



/*!
************************************************************************
* \brief
*    test if bit buffer contains only stop bit
*
* \param buffer
*    buffer containing VLC-coded data bits
* \param totbitoffset
*    bit offset from start of partition
* \param bytecount
*    buffer length
* \return
*    true if more bits available
************************************************************************
*/
int more_rbsp_data (byte buffer[],int totbitoffset,int bytecount)
{
	
	long byteoffset;      // byte from start of buffer
	int bitoffset;      // bit from start of byte
	int ctr_bit=0;      // control bit for current bit posision
	
	int cnt=0;
	
	
	byteoffset= totbitoffset/8;
	bitoffset= 7-(totbitoffset%8);
	
	assert (byteoffset<bytecount);
	
	// there is more until we're in the last byte
	if (byteoffset<(bytecount-1)) return TRUE1;
	
	// read one bit
	ctr_bit = (buffer[byteoffset] & (0x01<<bitoffset));
	
	// a stop bit has to be one
	if (ctr_bit==0) return TRUE1;
	
	bitoffset--;
	
	while (bitoffset>=0)
	{
		ctr_bit = (buffer[byteoffset] & (0x01<<bitoffset));   // set up control bit
		if (ctr_bit>0) cnt++;
		bitoffset--;
	}
	
	return (0!=cnt);
	
}


/*!
************************************************************************
* \brief
*    Check if there are symbols for the next MB
************************************************************************
*/
int uvlc_startcode_follows(struct img_par *img, struct inp_par_dec *inp, int dummy)
{
	int dp_Nr = assignSE2partition_dec[img->currentSlice->dp_mode][SE_MBTYPE];
	DataPartition *dP = &(img->currentSlice->partArr[dp_Nr]);
	Bitstream   *currStream = dP->bitstream;
	byte *buf = currStream->streamBuffer;
	
	//KS: new function test for End of Buffer
	return (!(more_rbsp_data(buf, currStream->frame_bitoffset,currStream->bitstream_length)));
}



/*!
 ************************************************************************
 * \brief
 *  read one exp-golomb VLC symbol
 *
 * \param buffer
 *    containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param  info 
 *    returns the value of the symbol
 * \param bytecount
 *    buffer length
 * \return
 *    bits read
 ************************************************************************
 */
int GetVLCSymbol_dec (byte buffer[],int totbitoffset,int *info, int bytecount)//totbitoffset：//24//25//26//27
{

  register int inf;
  long byteoffset;      // byte from start of buffer
  int bitoffset;      // bit from start of byte
  int ctr_bit=0;      // control bit for current bit posision
  int bitcounter=1;
  int len;
  int info_bit;

  byteoffset= totbitoffset/8;//3//3//3//3
  bitoffset= 7-(totbitoffset%8);//7//6//5//4
  ctr_bit = (buffer[byteoffset] & (0x01<<bitoffset));//128=232&128//64=232&64//32=232&32 //0  // set up control bit

  len=1;
  while (ctr_bit==0)
  {                 // find leading 1 bit
    len++;
    bitoffset-=1;           
    bitcounter++;
    if (bitoffset<0)
    {                 // finish with current byte ?
      bitoffset=bitoffset+8;
      byteoffset++;
    }
    ctr_bit=buffer[byteoffset] & (0x01<<(bitoffset));//8
  }
    // make infoword
  inf=0;                          // shortest possible code is 1, then info is always 0
  for(info_bit=0;(info_bit<(len-1)); info_bit++)
  {
    bitcounter++;//3
    bitoffset-=1;//2
    if (bitoffset<0)
    {                 // finished with current byte ?
      bitoffset=bitoffset+8;
      byteoffset++;
    }
    if (byteoffset > bytecount)//3》8
    {
      return -1;
    }
    inf=(inf<<1);
    if(buffer[byteoffset] & (0x01<<(bitoffset)))
      inf |=1;
  }

  *info = inf;
  return bitcounter;//1//3           // return absolute offset in bit from start of frame
}

extern void tracebits2_dec(const char *trace_str,  int len,  int info) ;

/*!
 ************************************************************************
 * \brief
 *    code from bitstream (2d tables)
 ************************************************************************
 */

int code_from_bitstream_2d_dec(SyntaxElement *sym,  
                           DataPartition *dP,
                           int *lentab,
                           int *codtab,
                           int tabwidth,
                           int tabheight,
                           int *code)
{
  Bitstream   *currStream = dP->bitstream;
  int frame_bitoffset = currStream->frame_bitoffset;
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;

  int i,j;
  int len, cod;

  // this VLC decoding method is not optimized for speed
  for (j = 0; j < tabheight; j++) {
    for (i = 0; i < tabwidth; i++)
    {
      len = lentab[i];
      if (!len)
        continue;
      cod = codtab[i];

      if ((ShowBits_dec(buf, frame_bitoffset, BitstreamLengthInBytes, len) == cod))
      {
        sym->value1 = i;
        sym->value2 = j;
        currStream->frame_bitoffset += len; // move bitstream pointer
        sym->len = len;
        goto found_code;
      }
    }
    lentab += tabwidth;
    codtab += tabwidth;
  }
  
  return -1;  // failed to find code

found_code:

  *code = cod;

  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    read FLC codeword from UVLC-partition 
 ************************************************************************
 */
int readSyntaxElement_FLC_dec(SyntaxElement *sym, Bitstream *currStream)
{
  int frame_bitoffset = currStream->frame_bitoffset;//14
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;//7245

  if ((GetBits_dec(buf, frame_bitoffset, &(sym->inf), BitstreamLengthInBytes, sym->len)) < 0)
    return -1;

  currStream->frame_bitoffset += sym->len; //18 move bitstream pointer
  sym->value1 = sym->inf;//0

#if TRACE
  tracebits2_dec(sym->tracestring, sym->len, sym->inf);
#endif

  return 1;
}



/*!
 ************************************************************************
 * \brief
 *    read NumCoeff/TrailingOnes codeword from UVLC-partition 
 ************************************************************************
 */

int readSyntaxElement_NumCoeffTrailingOnes_dec(SyntaxElement *sym,  DataPartition *dP,
                                           char *type)
{
  Bitstream   *currStream = dP->bitstream;
  int frame_bitoffset = currStream->frame_bitoffset;
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;

  int vlcnum, retval;
  int code, *ct, *lt;

  int lentab[3][4][17] = 
  {
    {   // 0702
      { 1, 6, 8, 9,10,11,13,13,13,14,14,15,15,16,16,16,16},
      { 0, 2, 6, 8, 9,10,11,13,13,14,14,15,15,15,16,16,16},
      { 0, 0, 3, 7, 8, 9,10,11,13,13,14,14,15,15,16,16,16},
      { 0, 0, 0, 5, 6, 7, 8, 9,10,11,13,14,14,15,15,16,16},
    },                                                 
    {                                                  
      { 2, 6, 6, 7, 8, 8, 9,11,11,12,12,12,13,13,13,14,14},
      { 0, 2, 5, 6, 6, 7, 8, 9,11,11,12,12,13,13,14,14,14},
      { 0, 0, 3, 6, 6, 7, 8, 9,11,11,12,12,13,13,13,14,14},
      { 0, 0, 0, 4, 4, 5, 6, 6, 7, 9,11,11,12,13,13,13,14},
    },                                                 
    {                                                  
      { 4, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9,10,10,10,10},
      { 0, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9,10,10,10},
      { 0, 0, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,10},
      { 0, 0, 0, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 9,10,10,10},
    },

  };

  int codtab[3][4][17] = 
  {
    {
      { 1, 5, 7, 7, 7, 7,15,11, 8,15,11,15,11,15,11, 7,4}, 
      { 0, 1, 4, 6, 6, 6, 6,14,10,14,10,14,10, 1,14,10,6}, 
      { 0, 0, 1, 5, 5, 5, 5, 5,13, 9,13, 9,13, 9,13, 9,5}, 
      { 0, 0, 0, 3, 3, 4, 4, 4, 4, 4,12,12, 8,12, 8,12,8},
    },
    {
      { 3,11, 7, 7, 7, 4, 7,15,11,15,11, 8,15,11, 7, 9,7}, 
      { 0, 2, 7,10, 6, 6, 6, 6,14,10,14,10,14,10,11, 8,6}, 
      { 0, 0, 3, 9, 5, 5, 5, 5,13, 9,13, 9,13, 9, 6,10,5}, 
      { 0, 0, 0, 5, 4, 6, 8, 4, 4, 4,12, 8,12,12, 8, 1,4},
    },
    {
      {15,15,11, 8,15,11, 9, 8,15,11,15,11, 8,13, 9, 5,1}, 
      { 0,14,15,12,10, 8,14,10,14,14,10,14,10, 7,12, 8,4},
      { 0, 0,13,14,11, 9,13, 9,13,10,13, 9,13, 9,11, 7,3},
      { 0, 0, 0,12,11,10, 9, 8,13,12,12,12, 8,12,10, 6,2},
    },
  };

  vlcnum = sym->value1;
  // vlcnum is the index of Table used to code coeff_token
  // vlcnum==3 means (8<=nC) which uses 6bit FLC

  if (vlcnum == 3)
  {
    // read 6 bit FLC
    code = ShowBits_dec(buf, frame_bitoffset, BitstreamLengthInBytes, 6);
    currStream->frame_bitoffset += 6;
    sym->value2 = code & 3;
    sym->value1 = (code >> 2);

    if (!sym->value1 && sym->value2 == 3)
    {
      // #c = 0, #t1 = 3 =>  #c = 0
      sym->value2 = 0;
    }
    else
      sym->value1++;

    sym->len = 6;

    retval = 0;
  }
  else

  {
    lt = &lentab[vlcnum][0][0];
    ct = &codtab[vlcnum][0][0];
    retval = code_from_bitstream_2d_dec(sym, dP, lt, ct, 17, 4, &code);
  }

  if (retval)
  {
    printf("ERROR: failed to find NumCoeff/TrailingOnes\n");
    exit(-1);
  }

#if TRACE
  printf(sym->tracestring, 
    TRACESTRING_SIZE, "%s # c & tr.1s vlc=%d #c=%d #t1=%d",
           type, vlcnum, sym->value1, sym->value2);
  tracebits2_dec(sym->tracestring, sym->len, code);

#endif

  return retval;
}


/*!
 ************************************************************************
 * \brief
 *    read NumCoeff/TrailingOnes codeword from UVLC-partition ChromaDC
 ************************************************************************
 */
int readSyntaxElement_NumCoeffTrailingOnesChromaDC_dec(SyntaxElement *sym,  DataPartition *dP)
{
  int retval;
  int code, *ct, *lt;

  int lentab[4][5] = 
  {
    { 2, 6, 6, 6, 6,},          
    { 0, 1, 6, 7, 8,}, 
    { 0, 0, 3, 7, 8,}, 
    { 0, 0, 0, 6, 7,},
  };

  int codtab[4][5] = 
  {
    {1,7,4,3,2},
    {0,1,6,3,3},
    {0,0,1,2,2},
    {0,0,0,5,0},
  };



  lt = &lentab[0][0];
  ct = &codtab[0][0];

  retval = code_from_bitstream_2d_dec(sym, dP, lt, ct, 5, 4, &code);

  if (retval)
  {
    printf("ERROR: failed to find NumCoeff/TrailingOnes ChromaDC\n");
    exit(-1);
  }


#if TRACE
    printf(sym->tracestring, 
      TRACESTRING_SIZE, "ChrDC # c & tr.1s  #c=%d #t1=%d",
              sym->value1, sym->value2);
    tracebits2_dec(sym->tracestring, sym->len, code);

#endif

  return retval;
}




/*!
 ************************************************************************
 * \brief
 *    read Level VLC0 codeword from UVLC-partition 
 ************************************************************************
 */
int readSyntaxElement_Level_VLC0_dec(SyntaxElement *sym, struct datapartition *dP)
{
  Bitstream   *currStream = dP->bitstream;
  int frame_bitoffset = currStream->frame_bitoffset;
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;
  int len, sign, level, code;

  len = 0;
  while (!ShowBits_dec(buf, frame_bitoffset+len, BitstreamLengthInBytes, 1))
    len++;

  len++;
  code = 1;
  frame_bitoffset += len;

  if (len < 15)
  {
    sign = (len - 1) & 1;
    level = (len-1) / 2 + 1;
  }
  else if (len == 15)
  {
    // escape code
    code = (code << 4) | ShowBits_dec(buf, frame_bitoffset, BitstreamLengthInBytes, 4);
    len += 4;
    frame_bitoffset += 4;
    sign = (code & 1);
    level = ((code >> 1) & 0x7) + 8;
  }
  else if (len == 16)
  {
    // escape code
    code = (code << 12) | ShowBits_dec(buf, frame_bitoffset, BitstreamLengthInBytes, 12);
    len += 12;
    frame_bitoffset += 12;
    sign =  (code & 1);
    level = ((code >> 1) & 0x7ff) + 16;
  }
  else
  {
    printf("ERROR reading Level code\n");
    exit(-1);
  }

  if (sign)
    level = -level;

  sym->inf = level;
  sym->len = len;

#if TRACE
  tracebits2_dec(sym->tracestring, sym->len, code);
#endif
  currStream->frame_bitoffset = frame_bitoffset;
  return 0;

}

/*!
 ************************************************************************
 * \brief
 *    read Level VLC codeword from UVLC-partition 
 ************************************************************************
 */
int readSyntaxElement_Level_VLCN_dec(SyntaxElement *sym, int vlc, struct datapartition *dP)  
{
  
  Bitstream   *currStream = dP->bitstream;
  int frame_bitoffset = currStream->frame_bitoffset;
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;
  
  int levabs, sign;
  int len = 0;
  int code, sb;
  
  int numPrefix;
  int shift = vlc-1;
  int escape = (15<<shift)+1;
  
  // read pre zeros
  numPrefix = 0;
  while (!ShowBits_dec(buf, frame_bitoffset+numPrefix, BitstreamLengthInBytes, 1))
    numPrefix++;
  
  
  len = numPrefix+1;
  code = 1;
  
  if (numPrefix < 15)
  {
    levabs = (numPrefix<<shift) + 1;
    
    // read (vlc-1) bits -> suffix
    if (vlc-1)
    {
      sb =  ShowBits_dec(buf, frame_bitoffset+len, BitstreamLengthInBytes, vlc-1);
      code = (code << (vlc-1) )| sb;
      levabs += sb;
      len += (vlc-1);
    }
    
    // read 1 bit -> sign
    sign = ShowBits_dec(buf, frame_bitoffset+len, BitstreamLengthInBytes, 1);
    code = (code << 1)| sign;
    len ++;
  }
  else  // escape
  {
    // read 11 bits -> levabs
    // levabs += escape
    sb = ShowBits_dec(buf, frame_bitoffset+len, BitstreamLengthInBytes, 11);
    code = (code << 11 )| sb;
    
    levabs =  sb + escape;
    len+=11;
    
    // read 1 bit -> sign
    sign = ShowBits_dec(buf, frame_bitoffset+len, BitstreamLengthInBytes, 1);
    code = (code << 1)| sign;
    len++;
  }
  
  sym->inf = (sign)?-levabs:levabs;
  sym->len = len;
  
  currStream->frame_bitoffset = frame_bitoffset+len;
  
#if TRACE
  tracebits2_dec(sym->tracestring, sym->len, code);
#endif
  
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *    read Total Zeros codeword from UVLC-partition 
 ************************************************************************
 */
int readSyntaxElement_TotalZeros_dec(SyntaxElement *sym,  DataPartition *dP)
{
  int vlcnum, retval;
  int code, *ct, *lt;

  int lentab[TOTRUN_NUM][16] = 
  {
    
    { 1,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9},  
    { 3,3,3,3,3,4,4,4,4,5,5,6,6,6,6},  
    { 4,3,3,3,4,4,3,3,4,5,5,6,5,6},  
    { 5,3,4,4,3,3,3,4,3,4,5,5,5},  
    { 4,4,4,3,3,3,3,3,4,5,4,5},  
    { 6,5,3,3,3,3,3,3,4,3,6},  
    { 6,5,3,3,3,2,3,4,3,6},  
    { 6,4,5,3,2,2,3,3,6},  
    { 6,6,4,2,2,3,2,5},  
    { 5,5,3,2,2,2,4},  
    { 4,4,3,3,1,3},  
    { 4,4,2,1,3},  
    { 3,3,1,2},  
    { 2,2,1},  
    { 1,1},  
  };

  int codtab[TOTRUN_NUM][16] = 
  {
    {1,3,2,3,2,3,2,3,2,3,2,3,2,3,2,1},
    {7,6,5,4,3,5,4,3,2,3,2,3,2,1,0},
    {5,7,6,5,4,3,4,3,2,3,2,1,1,0},
    {3,7,5,4,6,5,4,3,3,2,2,1,0},
    {5,4,3,7,6,5,4,3,2,1,1,0},
    {1,1,7,6,5,4,3,2,1,1,0},
    {1,1,5,4,3,3,2,1,1,0},
    {1,1,1,3,3,2,2,1,0},
    {1,0,1,3,2,1,1,1,},
    {1,0,1,3,2,1,1,},
    {0,1,1,2,1,3},
    {0,1,1,1,1},
    {0,1,1,1},
    {0,1,1},
    {0,1},  
  };
  vlcnum = sym->value1;

  lt = &lentab[vlcnum][0];
  ct = &codtab[vlcnum][0];

  retval = code_from_bitstream_2d_dec(sym, dP, lt, ct, 16, 1, &code);

  if (retval)
  {
    printf("ERROR: failed to find Total Zeros\n");
    exit(-1);
  }


#if TRACE
    tracebits2_dec(sym->tracestring, sym->len, code);

#endif

  return retval;
}    

/*!
 ************************************************************************
 * \brief
 *    read Total Zeros Chroma DC codeword from UVLC-partition 
 ************************************************************************
 */
int readSyntaxElement_TotalZerosChromaDC_dec(SyntaxElement *sym,  DataPartition *dP)
{
  int vlcnum, retval;
  int code, *ct, *lt;

  int lentab[3][4] = 
  {
    { 1, 2, 3, 3,},
    { 1, 2, 2, 0,},
    { 1, 1, 0, 0,}, 
  };

  int codtab[3][4] = 
  {
    { 1, 1, 1, 0,},
    { 1, 1, 0, 0,},
    { 1, 0, 0, 0,},
  };

  vlcnum = sym->value1;

  lt = &lentab[vlcnum][0];
  ct = &codtab[vlcnum][0];

  retval = code_from_bitstream_2d_dec(sym, dP, lt, ct, 4, 1, &code);

  if (retval)
  {
    printf("ERROR: failed to find Total Zeros\n");
    exit(-1);
  }


#if TRACE
    tracebits2_dec(sym->tracestring, sym->len, code);

#endif

  return retval;
}    


/*!
 ************************************************************************
 * \brief
 *    read  Run codeword from UVLC-partition 
 ************************************************************************
 */
int readSyntaxElement_Run_dec(SyntaxElement *sym,  DataPartition *dP)
{
  int vlcnum, retval;
  int code, *ct, *lt;

  int lentab[TOTRUN_NUM][16] = 
  {
    {1,1},
    {1,2,2},
    {2,2,2,2},
    {2,2,2,3,3},
    {2,2,3,3,3,3},
    {2,3,3,3,3,3,3},
    {3,3,3,3,3,3,3,4,5,6,7,8,9,10,11},
  };

  int codtab[TOTRUN_NUM][16] = 
  {
    {1,0},
    {1,1,0},
    {3,2,1,0},
    {3,2,1,1,0},
    {3,2,3,2,1,0},
    {3,0,1,3,2,5,4},
    {7,6,5,4,3,2,1,1,1,1,1,1,1,1,1},
  };

  vlcnum = sym->value1;

  lt = &lentab[vlcnum][0];
  ct = &codtab[vlcnum][0];

  retval = code_from_bitstream_2d_dec(sym, dP, lt, ct, 16, 1, &code);

  if (retval)
  {
    printf("ERROR: failed to find Run\n");
    exit(-1);
  }


#if TRACE
    tracebits2_dec(sym->tracestring, sym->len, code);
#endif

  return retval;
}    


/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer
 *
 * \param buffer
 *    containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param info
 *    returns value of the read bits
 * \param bytecount
 *    total bytes in bitstream
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */
int GetBits_dec (byte buffer[],int totbitoffset,int *info, int bytecount, 
             int numbits)
{

  register int inf;
  long byteoffset;      // byte from start of buffer
  int bitoffset;      // bit from start of byte

  int bitcounter=numbits;

  byteoffset= totbitoffset/8;
  bitoffset= 7-(totbitoffset%8);

  inf=0;
  while (numbits)
  {
    inf <<=1;
    inf |= (buffer[byteoffset] & (0x01<<bitoffset))>>bitoffset;
    numbits--;
    bitoffset--;
    if (bitoffset < 0)
    {
      byteoffset++;
      bitoffset += 8;
      if (byteoffset > bytecount)
      {
        return -1;
      }
    }
  }

  *info = inf;
  return bitcounter;           // return absolute offset in bit from start of frame
}     

/*!
 ************************************************************************
 * \brief
 *  Reads bits from the bitstream buffer
 *
 * \param buffer
 *    buffer containing VLC-coded data bits
 * \param totbitoffset
 *    bit offset from start of partition
 * \param bytecount
 *    total bytes in bitstream
 * \param numbits
 *    number of bits to read
 *
 ************************************************************************
 */

int ShowBits_dec (byte buffer[],int totbitoffset,int bytecount, int numbits)
{

  register int inf;
  long byteoffset;      // byte from start of buffer
  int bitoffset;      // bit from start of byte

  byteoffset= totbitoffset/8;
  bitoffset= 7-(totbitoffset%8);

  inf=0;
  while (numbits)
  {
    inf <<=1;
    inf |= (buffer[byteoffset] & (0x01<<bitoffset))>>bitoffset;
    numbits--;
    bitoffset--;
    if (bitoffset < 0)
    {
      byteoffset++;
      bitoffset += 8;
      if (byteoffset > bytecount)
      {
        return -1;
      }
    }
  }

  return inf;           // return absolute offset in bit from start of frame
}     


/*!
 ************************************************************************
 * \brief
 *    peek at the next 2 UVLC codeword from UVLC-partition to determine
 *    if a skipped MB is field/frame
 ************************************************************************
 */
int peekSyntaxElement_UVLC_dec(SyntaxElement *sym, struct img_par *img, struct inp_par_dec *inp, struct datapartition *dP)
{
  Bitstream   *currStream = dP->bitstream;
  int frame_bitoffset = currStream->frame_bitoffset;
  byte *buf = currStream->streamBuffer;
  int BitstreamLengthInBytes = currStream->bitstream_length;


  sym->len =  GetVLCSymbol_dec (buf, frame_bitoffset, &(sym->inf), BitstreamLengthInBytes);
  if (sym->len == -1)
    return -1;
  frame_bitoffset += sym->len;
  sym->mapping_dec(sym->len,sym->inf,&(sym->value1),&(sym->value2));


#if TRACE
  tracebits(sym->tracestring, sym->len, sym->inf, sym->value1);
#endif

  return 1;
}


/////////////////////////////////////<
