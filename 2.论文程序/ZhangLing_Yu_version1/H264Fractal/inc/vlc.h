
/*!
 *************************************************************************************
 * \file vlc.h
 *
 * \brief
 *    Prototypes for VLC coding funtions
 * \author
 *     Karsten Suehring
 *************************************************************************************
 */

#ifndef _VLC_H_
#define _VLC_H_

int u_1 (char *tracestring, int value, DataPartition *part);
int se_v (char *tracestring, int value, DataPartition *part);
int ue_v (char *tracestring, int value, DataPartition *part);
int u_v (int n, char *tracestring, int value, DataPartition *part);


void levrun_linfo_c2x2(int level,int run,int *len,int *info);
void levrun_linfo_intra(int level,int run,int *len,int *info);
void levrun_linfo_inter(int level,int run,int *len,int *info);

int   writeSyntaxElement_UVLC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_fixed(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement2Buf_UVLC(SyntaxElement *se, Bitstream* this_streamBuffer );
void  writeUVLC2buffer(SyntaxElement *se, Bitstream *currStream);
int   writeSyntaxElement2Buf_Fixed(SyntaxElement *se, Bitstream* this_streamBuffer );
int   symbol2uvlc(SyntaxElement *se);
void  ue_linfo(int n, int dummy, int *len,int *info);
void  se_linfo(int mvd, int dummy, int *len,int *info);
void  cbp_linfo_intra(int cbp, int dummy, int *len,int *info);
void  cbp_linfo_inter(int cbp, int dummy, int *len,int *info);

// CAVLC
void  CAVLC_init();
int writeCoeff4x4_CAVLC (int reference, int b8, int b4, int param);

int   writeSyntaxElement_VLC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_TotalZeros(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_TotalZerosChromaDC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_Run(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_NumCoeffTrailingOnes(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_NumCoeffTrailingOnesChromaDC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_Level_VLC1(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_Level_VLCN(SyntaxElement *se, int vlc, DataPartition *this_dataPart);
int   writeSyntaxElement_Intra4x4PredictionMode(SyntaxElement *se, DataPartition *this_dataPart);




/////////////////////////////////////>


int se_v_dec (char *tracestring, Bitstream *bitstream);
int ue_v_dec (char *tracestring, Bitstream *bitstream);
int u_1_dec (char *tracestring, Bitstream *bitstream);
int u_v_dec (int LenInBits, char *tracestring, Bitstream *bitstream);

// UVLC mapping
void linfo_ue_dec(int len, int info, int *value1, int *dummy);
void linfo_se_dec(int len, int info, int *value1, int *dummy);

void linfo_cbp_intra_dec(int len,int info,int *cbp, int *dummy);
void linfo_cbp_inter_dec(int len,int info,int *cbp, int *dummy);
void linfo_levrun_inter_dec(int len,int info,int *level,int *irun);
void linfo_levrun_c2x2_dec(int len,int info,int *level,int *irun);

int  readSyntaxElement_VLC_dec (SyntaxElement *sym, Bitstream *currStream);
int  readSyntaxElement_UVLC_dec(SyntaxElement *sym, struct img_par *img, struct inp_par_dec *inp, struct datapartition *dp);
int  readSyntaxElement_Intra4x4PredictionMode_dec(SyntaxElement *sym, struct img_par *img, struct inp_par_dec *inp, struct datapartition *dp);

int  GetVLCSymbol_dec (byte buffer[],int totbitoffset,int *info, int bytecount);
int  GetVLCSymbol_IntraMode_dec (byte buffer[],int totbitoffset,int *info, int bytecount);

int readSyntaxElement_FLC_dec(SyntaxElement *sym, Bitstream *currStream);
int readSyntaxElement_NumCoeffTrailingOnes_dec(SyntaxElement *sym,  DataPartition *dP,
                                           char *type);
int readSyntaxElement_NumCoeffTrailingOnesChromaDC_dec(SyntaxElement *sym,  DataPartition *dP);
int readSyntaxElement_Level_VLC0_dec(SyntaxElement *sym, struct datapartition *dP);
int readSyntaxElement_Level_VLCN_dec(SyntaxElement *sym, int vlc, struct datapartition *dP);
int readSyntaxElement_TotalZeros_dec(SyntaxElement *sym,  DataPartition *dP);
int readSyntaxElement_TotalZerosChromaDC_dec(SyntaxElement *sym,  DataPartition *dP);
int readSyntaxElement_Run_dec(SyntaxElement *sym,  DataPartition *dP);
int GetBits_dec (byte buffer[],int totbitoffset,int *info, int bytecount, 
             int numbits);
int ShowBits_dec (byte buffer[],int totbitoffset,int bytecount, int numbits);





/////////////////////////////////////<


#endif

