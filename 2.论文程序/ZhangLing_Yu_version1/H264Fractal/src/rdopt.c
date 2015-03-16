
/*!
 ***************************************************************************
 * \file rdopt.c
 *
 * \brief
 *    Rate-Distortion optimized mode decision
 *
 * \author
 *    Heiko Schwarz <hschwarz@hhi.de>
 *
 * \date
 *    12. April 2001
 **************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include "rdopt_coding_state.h"
#include "elements.h"
#include "refbuf.h"
#include "intrarefresh.h"
#include "vlc.h"
#include "mbuffer.h"
#include "image.h"
#include "mb_access.h"
#include "fast_me.h"
#include "ratectl.h"            // head file for rate control
#include "cabac.h"            // head file for rate control

#include "i_global.h"
#include "i_defines.h"

//Rate control

int QP,QP2;
int DELTA_QP,DELTA_QP2;
int diffy[16][16];
static int pred[16][16];

extern       int  QP2QUANT  [40];

//==== MODULE PARAMETERS ====
int   best_mode;
int   rec_mbY[16][16], rec_mbU[8][8], rec_mbV[8][8], rec_mbY8x8[16][16];    // reconstruction values
int   mpr8x8[16][16];
int   ****cofAC=NULL, ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
int   ***cofDC=NULL;                       // [yuv][level/run][scan_pos]
int   **cofAC4x4=NULL, ****cofAC4x4intern=NULL; // [level/run][scan_pos]
int   cbp, cbp8x8, cnt_nonz_8x8;
int   cbp_blk, cbp_blk8x8;
int   frefframe[4][4], brefframe[4][4], b8mode[4], b8pdir[4];
int   best8x8mode [4];                // [block]
int   best8x8pdir [MAXMODE][4];       // [mode][block]
int   best8x8fwref  [MAXMODE][4];       // [mode][block]
int   b8_ipredmode[16], b8_intra_pred_modes[16];
CSptr cs_mb=NULL, cs_b8=NULL, cs_cm=NULL, cs_imb=NULL, cs_ib8=NULL, cs_ib4=NULL, cs_pc=NULL;
int   best_c_imode;
int   best_i16offset;

int   best8x8bwref     [MAXMODE][4];       // [mode][block]
int   abp_typeframe[4][4];

/*!
 ************************************************************************
 * \brief
 *    delete structure for RD-optimized mode decision
 ************************************************************************
 */
void clear_rdopt ()
{
  free_mem_DCcoeff (cofDC);
  free_mem_ACcoeff (cofAC);
  free_mem_ACcoeff (cofAC8x8);
  free_mem_ACcoeff (cofAC4x4intern);

  // structure for saving the coding state
  delete_coding_state (cs_mb);
  delete_coding_state (cs_b8);
  delete_coding_state (cs_cm);
  delete_coding_state (cs_imb);
  delete_coding_state (cs_ib8);
  delete_coding_state (cs_ib4);
  delete_coding_state (cs_pc);
}


/*!
 ************************************************************************
 * \brief
 *    create structure for RD-optimized mode decision
 ************************************************************************
 */
void init_rdopt ()
{
  get_mem_DCcoeff (&cofDC);
  get_mem_ACcoeff (&cofAC);
  get_mem_ACcoeff (&cofAC8x8);
  get_mem_ACcoeff (&cofAC4x4intern);
  cofAC4x4 = cofAC4x4intern[0][0];

  // structure for saving the coding state
  cs_mb  = create_coding_state ();
  cs_b8  = create_coding_state ();
  cs_cm  = create_coding_state ();
  cs_imb = create_coding_state ();
  cs_ib8 = create_coding_state ();
  cs_ib4 = create_coding_state ();
  cs_pc  = create_coding_state ();
}



/*! 
 *************************************************************************************
 * \brief
 *    Updates the pixel map that shows, which reference frames are reliable for
 *    each MB-area of the picture.
 *
 * \note
 *    The new values of the pixel_map are taken from the temporary buffer refresh_map
 *
 *************************************************************************************
 */
void UpdatePixelMap()
{
  int mx,my,y,x,i,j;
  if (img->type==I_SLICE)
  {
    for (y=0; y<img->height; y++)
    for (x=0; x<img->width; x++)
    {
      pixel_map[y][x]=1;
    }
  }
  else
  {
    for (my=0; my<img->height/8; my++)
    for (mx=0; mx<img->width/8;  mx++)
    {
      j = my*8 + 8;
      i = mx*8 + 8;
      if (refresh_map[my][mx])
      {
        for (y=my*8; y<j; y++)
        for (x=mx*8; x<i; x++)  pixel_map[y][x] = 1;
      }
      else
      {
        for (y=my*8; y<j; y++)
        for (x=mx*8; x<i; x++)  pixel_map[y][x] = min(pixel_map[y][x]+1, input->num_reference_frames+1);
      }
    }
  }
}

/*! 
 *************************************************************************************
 * \brief
 *    Checks if a given reference frame is reliable for the current 
 *    macroblock, given the motion vectors that the motion search has 
 *    returned.
 *
 * \return
 *    If the return value is 1, the reference frame is reliable. If it 
 *    is 0, then it is not reliable.
 *
 * \note
 *    A specific area in each reference frame is assumed to be unreliable
 *    if the same area has been intra-refreshed in a subsequent frame.
 *    The information about intra-refreshed areas is kept in the pixel_map.
 *
 *************************************************************************************
 */
int CheckReliabilityOfRef (int block, int list_idx, int ref, int mode)
{
  int y,x, block_y, block_x, dy, dx, y_pos, x_pos, yy, xx, pres_x, pres_y;
  int maxold_x  = img->width-1;
  int maxold_y  = img->height-1;
  int ref_frame = ref+1;

  int by0 = (mode>=4?2*(block/2):mode==2?2*block:0);
  int by1 = by0 + (mode>=4||mode==2?2:4);
  int bx0 = (mode>=4?2*(block%2):mode==3?2*block:0);
  int bx1 = bx0 + (mode>=4||mode==3?2:4);

  for (block_y=by0; block_y<by1; block_y++)
    for (block_x=bx0; block_x<bx1; block_x++)
    {
      y_pos  = img->all_mv[block_x][block_y][list_idx][ref][mode][1];
      y_pos += (img->block_y+block_y) * BLOCK_SIZE * 4;
      x_pos  = img->all_mv[block_x][block_y][list_idx][ref][mode][0];
      x_pos += (img->block_x+block_x) * BLOCK_SIZE * 4;
      
      /* Here we specify which pixels of the reference frame influence
         the reference values and check their reliability. This is
         based on the function Get_Reference_Pixel */
      
      dy = y_pos & 3;
      dx = x_pos & 3;

      y_pos = (y_pos-dy)/4;
      x_pos = (x_pos-dx)/4;

      if (dy==0 && dx==0) //full-pel
      {
        for (y=0 ; y < BLOCK_SIZE ; y++)
          for (x=0 ; x < BLOCK_SIZE ; x++)
            if (pixel_map[max(0,min(maxold_y,y_pos+y))][max(0,min(maxold_x,x_pos+x))] < ref_frame)
              return 0;
      }
      else  /* other positions */
      {
        if (dy == 0) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              pres_y = max(0,min(maxold_y,y_pos+y));
              for(xx=-2;xx<4;xx++) {
                pres_x = max(0,min(maxold_x,x_pos+x+xx));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }

        else if (dx == 0) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              pres_x = max(0,min(maxold_x,x_pos+x));
              for(yy=-2;yy<4;yy++) {
                pres_y = max(0,min(maxold_y,y_pos+yy+y));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }
        else if (dx == 2) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              for(yy=-2;yy<4;yy++) {
                pres_y = max(0,min(maxold_y,y_pos+yy+y));
                for(xx=-2;xx<4;xx++) {
                  pres_x = max(0,min(maxold_x,x_pos+xx+x));
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else if (dy == 2) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              for(xx=-2;xx<4;xx++) {
                pres_x = max(0,min(maxold_x,x_pos+xx+x));
                for(yy=-2;yy<4;yy++) {
                  pres_y = max(0,min(maxold_y,y_pos+yy+y));
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              pres_y = dy == 1 ? y_pos+y : y_pos+y+1;
              pres_y = max(0,min(maxold_y,pres_y));

              for(xx=-2;xx<4;xx++) {
                pres_x = max(0,min(maxold_x,x_pos+xx+x));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }

              pres_x = dx == 1 ? x_pos+x : x_pos+x+1;
              pres_x = max(0,min(maxold_x,pres_x));

              for(yy=-2;yy<4;yy++) {
                pres_y = max(0,min(maxold_y,y_pos+yy+y));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }

      }
    }

  return 1;
}



/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for an 4x4 Intra block
 *************************************************************************************
 */
double RDCost_for_4x4IntraBlocks (int*    nonzero,
                           int     b8,
                           int     b4,
                           int    ipmode,
                           double  lambda,
                           double  min_rdcost,
                           int mostProbableMode)
{
  double  rdcost;
  int     dummy, x, y, rate;
  int     distortion  = 0;
  int     block_x     = 8*(b8%2)+4*(b4%2);
  int     block_y     = 8*(b8/2)+4*(b4/2);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     pic_opix_y  = img->opix_y+block_y;
  byte    **imgY      = enc_picture->imgY;

  Slice          *currSlice    =  img->currentSlice;
  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];
  SyntaxElement  *currSE       = &img->MB_SyntaxElements[currMB->currSEnr];
  const int      *partMap      = assignSE2partition[input->partition_mode];
  DataPartition  *dataPart;

  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  dummy = 0;
  *nonzero = dct_luma (block_x, block_y, &dummy, 1);

  //===== get distortion (SSD) of 4x4 block =====计算误差平方和~~
  for (y=0; y<4; y++)
  for (x=pic_pix_x; x<pic_pix_x+4; x++)  
    distortion += img->quad [imgY_org[pic_opix_y+y][x] - imgY[pic_pix_y+y][x]];

  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
  currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

  //--- set position and type ---
  currSE->context = 4*b8 + b4;
  currSE->type    = SE_INTRAPREDMODE;

  //--- set function pointer ----
  if (input->symbol_mode != UVLC)    
  currSE->writing = writeIntraPredMode_CABAC;

  //--- choose data partition ---
  dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
  //--- encode and update rate ---
  if (input->symbol_mode == UVLC)    writeSyntaxElement_Intra4x4PredictionMode(currSE, dataPart);
  else                               dataPart->writeSyntaxElement (currSE, dataPart);

  rate = currSE->len;
  currSE++;
  currMB->currSEnr++;

  //===== RATE for LUMINANCE COEFFICIENTS =====
  if (input->symbol_mode == UVLC)
  {
    rate  += writeCoeff4x4_CAVLC (LUMA, b8, b4, 0);
  }
  else
  {
    rate  += writeLumaCoeff4x4_CABAC (b8, b4, 1);
  }
  rdcost = (double)distortion + lambda*(double)rate;

  return rdcost;
}

/*! 
 *************************************************************************************
 * \brief
 *    Mode Decision for an 4x4 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_4x4IntraBlocks (int  b8,  int  b4,  double  lambda,  int*  min_cost)
{
  int     ipmode, best_ipmode = 0, i, j, k, x, y, cost, dummy;
  int     c_nz, nonzero = 0, rec4x4[4][4], diff[16];
  double  rdcost;
  int     block_x     = 8*(b8%2)+4*(b4%2);
  int     block_y     = 8*(b8/2)+4*(b4/2);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     pic_opix_x   = img->opix_x+block_x;
  int     pic_opix_y   = img->opix_y+block_y;
  int     pic_block_x = pic_pix_x/4;
  int     pic_block_y = pic_pix_y/4;
  double  min_rdcost  = 1e30;

  int left_available, up_available, all_available;

  int     upMode;
  int     leftMode;
  int     mostProbableMode;

  PixelPos left_block;
  PixelPos top_block;

  getLuma4x4Neighbour(img->current_mb_nr, block_x/4, block_y/4, -1,  0, &left_block);
  getLuma4x4Neighbour(img->current_mb_nr, block_x/4, block_y/4,  0, -1, &top_block);


  // constrained intra pred
  if (input->UseConstrainedIntraPred)
  {
    left_block.available = left_block.available ? img->intra_block[ left_block.mb_addr ] : 0;//////img->intra_block
    top_block.available  = top_block.available  ? img->intra_block[ top_block.mb_addr ]  : 0;
  }
  
  upMode            = top_block.available ? img->ipredmode[top_block.pos_x ][top_block.pos_y ] : -1;
  leftMode          = left_block.available ? img->ipredmode[left_block.pos_x][left_block.pos_y] : -1;
  
  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;
  
  *min_cost = (1<<20);

  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  intrapred_luma (pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);

  //===== LOOP OVER ALL 4x4 INTRA PREDICTION MODES =====
  for (ipmode=0; ipmode<NO_INTRA_PMODE; ipmode++)
  {
    int available_mode =  (ipmode==DC_PRED) ||
        ((ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED) && up_available ) ||
        ((ipmode==HOR_PRED||ipmode==HOR_UP_PRED) && left_available ) ||(all_available);
    if( available_mode)
    {
      if (!input->rdopt)
      {
        for (k=j=0; j<4; j++)
          for (i=0; i<4; i++, k++)
          {
            diff[k] = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr[ipmode][j][i];
          }
        cost  = (ipmode == mostProbableMode) ? 0 : (int)floor(4 * lambda );
        cost += SATD (diff, input->hadamard);
        if (cost < *min_cost)
        {
          best_ipmode = ipmode;
          *min_cost   = cost;
        }
      }
      else
      {
        // get prediction and prediction error
        for (j=0; j<4; j++)
        for (i=0; i<4; i++)
        {
          img->mpr[block_x+i][block_y+j]  = img->mprr[ipmode][j][i];
          img->m7[i][j]                   = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr[ipmode][j][i];
        }

        //===== store the coding state =====
        store_coding_state (cs_cm);
        // get and check rate-distortion cost
        if ((rdcost = RDCost_for_4x4IntraBlocks (&c_nz, b8, b4, ipmode, lambda, min_rdcost, mostProbableMode)) < min_rdcost)
        {
          //--- set coefficients ---
          for (j=0; j<2; j++)
          for (i=0; i<18;i++)  cofAC4x4[j][i]=img->cofAC[b8][b4][j][i];

          //--- set reconstruction ---
          for (y=0; y<4; y++)
          for (x=0; x<4; x++)  rec4x4[y][x] = enc_picture->imgY[pic_pix_y+y][pic_pix_x+x];

          //--- flag if dct-coefficients must be coded ---
          nonzero = c_nz;

          //--- set best mode update minimum cost ---
          min_rdcost  = rdcost;
          best_ipmode = ipmode;
        }
        reset_coding_state (cs_cm);
      }
    }
  }

  //===== set intra mode prediction =====
  img->ipredmode[pic_block_x][pic_block_y] = best_ipmode;
  img->mb_data[img->current_mb_nr].intra_pred_modes[4*b8+b4] = mostProbableMode == best_ipmode ? -1 : best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1;

  if (!input->rdopt)
  {
    // get prediction and prediction error
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        img->mpr[block_x+i][block_y+j]  = img->mprr[best_ipmode][j][i];
        img->m7[i][j]                   = imgY_org[pic_opix_y+j][pic_opix_x+i] - img->mprr[best_ipmode][j][i];
      }
    nonzero = dct_luma (block_x, block_y, &dummy, 1);
  }
  else//1
  {
    //===== restore coefficients =====
    for (j=0; j<2; j++)
    for (i=0; i<18;i++)  img->cofAC[b8][b4][j][i]=cofAC4x4[j][i];
  
    //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
    for (y=0; y<4; y++)
    for (x=0; x<4; x++)
    {
      enc_picture->imgY[pic_pix_y+y][pic_pix_x+x] = rec4x4[y][x];
      img->mpr[block_x+x][block_y+y] = img->mprr[best_ipmode][y][x];
    }
  }

  return nonzero;
}


/*! 
 *************************************************************************************
 * \brief
 *    Mode Decision for an 8x8 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_8x8IntraBlocks(int b8,double lambda,int *cost)
{
  int  nonzero=0, b4;
  int  cost4x4;
  
  *cost = (int)floor(6.0 * lambda + 0.4999);

  for (b4=0; b4<4; b4++)
  {
    if (Mode_Decision_for_4x4IntraBlocks (b8, b4, lambda, &cost4x4))
    {
      nonzero        = 1;
    }
    *cost += cost4x4;
  }

  return nonzero;
}

/*! 
 *************************************************************************************
 * \brief
 *    4x4 Intra mode decision for an macroblock
 *************************************************************************************
 */
int Mode_Decision_for_Intra4x4Macroblock (double lambda,  int* cost)

{
  int  cbp=0, b8, cost8x8;

  for (*cost=0, b8=0; b8<4; b8++)
  {
    if (Mode_Decision_for_8x8IntraBlocks (b8, lambda, &cost8x8))
    {
      cbp |= (1<<b8);
    }
    *cost += cost8x8;
  }

  return cbp;
}


/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for an 8x8 Partition
 *************************************************************************************
 */
double RDCost_for_8x8blocks (int*    cnt_nonz,   // --> number of nonzero coefficients
                             int*    cbp_blk,    // --> cbp blk
                             double  lambda,     // <-- lagrange multiplier
                             int     block,      // <-- 8x8 block number
                             int     mode,       // <-- partitioning mode
                             int     pdir,       // <-- prediction direction
                             int     ref,        // <-- reference frame
                             int     bwd_ref)    // <-- abp type
{
  int  i, j, k;
  int  rate=0, distortion=0;
  int  dummy, mrate;
  int  fw_mode, bw_mode;
  int  cbp     = 0;
  int  pax     = 8*(block%2);
  int  pay     = 8*(block/2);
  int  i0      = pax/4;
  int  j0      = pay/4;
  int  bframe  = (img->type==B_SLICE);
  int  direct  = (bframe && mode==0);
  int  b8value = B8Mode2Value (mode, pdir);

  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice         *currSlice = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap   = assignSE2partition[input->partition_mode];

  EncodingEnvironmentPtr eep_dp;

  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  if (direct)
  {
    if (direct_pdir[img->block_x+i0][img->block_y+j0]<0) // mode not allowed
    {
      return (1e20);
    }
    else
    {
      *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, direct_pdir[img->block_x+i0][img->block_y+j0], 0, 0, max(0,direct_ref_idx[LIST_0][img->block_x+i0][img->block_y+j0]), direct_ref_idx[LIST_1][img->block_x+i0][img->block_y+j0]);
    }
  }
  else
  {
    fw_mode   = (pdir==0||pdir==2 ? mode : 0);
    bw_mode   = (pdir==1||pdir==2 ? mode : 0);
    *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, pdir, fw_mode, bw_mode, ref, bwd_ref);
  }

  //===== get residue =====
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_b8block (block, -1);
  }

  //=====
  //=====   GET DISTORTION
  //=====
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders ;k++)
    {
      decode_one_b8block (k, P8x8, block, mode, ref);
      for (j=img->opix_y+pay; j<img->opix_y+pay+8; j++)
      for (i=img->opix_x+pax; i<img->opix_x+pax+8; i++)
      {
        distortion += img->quad[imgY_org[j][i] - decs->decY[k][j][i]];
      }
    }
    distortion /= input->NoOfDecoders;
  }
  else
  {
    for (j=pay; j<pay+8; j++)
    for (i=img->pix_x+pax; i<img->pix_x+pax+8; i++)
    {
      distortion += img->quad [imgY_org[img->opix_y+j][i] - enc_picture->imgY[img->pix_y+j][i]];
    }
  }

  //=====
  //=====   GET RATE
  //=====
  //----- block 8x8 mode -----
  if (input->symbol_mode == UVLC)
  {
    ue_linfo (b8value, dummy, &mrate, &dummy);
    rate += mrate;
  }
  else
  {
    currSE->value1  = b8value;
    currSE->writing = writeB8_typeInfo_CABAC;
    currSE->type    = SE_MBTYPE;
    dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    dataPart->writeSyntaxElement (currSE, dataPart);
    rate += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }
////////////////////////////////先不要了 因为与参考帧有关
//   //----- motion information -----
  if (!direct)
  {
    if ((img->num_ref_idx_l0_active > 1 ) && (pdir==0 || pdir==2))
      rate  += writeReferenceFrame (mode, i0, j0, 1, ref);
    if(img->num_ref_idx_l1_active > 1 && img->type== B_SLICE)
    {
      if (pdir==1 || pdir==2)
      {
        rate  += writeReferenceFrame (mode, i0, j0, 0, bwd_ref);
      }
    }

    if (pdir==0 || pdir==2)
    {
      rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, ref, LIST_0, mode);
    }
    if (pdir==1 || pdir==2)
    {
      rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, bwd_ref, LIST_1, mode);
    }
  }

//   //----- coded block pattern (for CABAC only) -----
  if (input->symbol_mode == CABAC)
  {
    dataPart = &(currSlice->partArr[partMap[SE_CBP_INTER]]);
    eep_dp   = &(dataPart->ee_cabac);
    mrate    = arienco_bits_written (eep_dp);
    writeCBP_BIT_CABAC (block, ((*cnt_nonz>0)?1:0), cbp8x8, currMB, 1, eep_dp);
    mrate    = arienco_bits_written (eep_dp) - mrate;
    rate    += mrate;
  }

  //----- luminance coefficients -----
  if (*cnt_nonz)
  {
    rate += writeLumaCoeff8x8 (block, 0);
  }

  return (double)distortion + lambda * (double)rate;
}


/*! 
 *************************************************************************************
 * \brief
 *    Gets mode offset for intra16x16 mode
 *************************************************************************************
 */
int I16Offset (int cbp, int i16mode)
{
  return (cbp&15?13:1) + i16mode + ((cbp&0x30)>>2);
}


/*! 
 *************************************************************************************
 * \brief
 *    Sets modes and reference frames for an macroblock
 *************************************************************************************
 */
void SetModesAndRefframeForBlocks (int mode)
{
  int i,j,k,l;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int  bframe  = (img->type==B_SLICE);

  int list_offset   = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;

  //--- macroblock type ---
  currMB->mb_type = mode;
  
  //--- block 8x8 mode and prediction direction ---
  switch (mode)
  {
  case 0:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = 0;
      currMB->b8pdir[i] = (bframe?direct_pdir[img->block_x+(i%2)*2][img->block_y+(i/2)*2]:0);
    }
    break;
  case 1:
  case 2:
  case 3:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = mode;
      currMB->b8pdir[i] = best8x8pdir[mode][i];
    }
    break;
  case P8x8:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i]   = best8x8mode[i];
      currMB->b8pdir[i]   = best8x8pdir[mode][i];
    }
    break;
  case I4MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = IBLOCK;
      currMB->b8pdir[i] = -1;
    }
    break;
  case I16MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] =  0;
      currMB->b8pdir[i] = -1;
    }
    break;
  default:
    printf ("Unsupported mode in SetModesAndRefframeForBlocks!\n");
    exit (1);
  }
  
#define IS_FW ((best8x8pdir[mode][k]==0 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0 || !bframe))
#define IS_BW ((best8x8pdir[mode][k]==1 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0))
  //--- reference frame arrays ---//ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
  if (mode==0 || mode==I4MB || mode==I16MB)
  {
    if (bframe)
    { 
      for (j=0;j<4;j++)
        for (i=0;i<4;i++)
        {
          if(!mode)
          {     //direct mode
            enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = direct_ref_idx[LIST_0][img->block_x+i][img->block_y+j];
            enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = direct_ref_idx[LIST_1][img->block_x+i][img->block_y+j];
          }
          else
          {   //intra
            enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = -1;
            enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
            
          }
        }
    }
    else
    {
      for (j=0;j<4;j++)
        for (i=0;i<4;i++)
        {
          enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = (mode==0?0:-1);
        }
    }
  }
  else
  {
    if (bframe)
    {
      for (j=0;j<4;j++)
        for (i=0;i<4;i++)
        {
          k = 2*(j/2)+(i/2);
          l = 2*(j%2)+(i%2);
          
          if(mode == P8x8 && best8x8mode[k]==0)
          {           
            enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = direct_ref_idx[LIST_0][img->block_x+i][img->block_y+j];
            enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = direct_ref_idx[LIST_1][img->block_x+i][img->block_y+j];
          }
          else
          {
            enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = (IS_FW ? best8x8fwref[mode][k] : -1);
            enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = (IS_BW ? best8x8bwref[mode][k] : -1);
          }
        }
    }
    else
    {
      for (j=0;j<4;j++)
        for (i=0;i<4;i++)
        {
          k = 2*(j/2)+(i/2);
          l = 2*(j%2)+(i%2);
          enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = (IS_FW ? best8x8fwref[mode][k] : -1);
        }
    }
  }

  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
    {
      enc_picture->ref_pic_id [LIST_0][img->block_x+i][img->block_y+j] = 
        (enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]>=0 ? 
         enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]]:
         -1);
    }
  }
  if (bframe)
  {
    for (j=0;j<4;j++)
    {
      for (i=0;i<4;i++)
      {
        enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+j] = 
          (enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j]>=0 ? 
           enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j]]:
           -1);
      }
      
    }
  }
// 

#undef IS_FW
#undef IS_BW
}


/*! 
 *************************************************************************************
 * \brief
 *    Intra 16x16 mode decision
 *************************************************************************************
 */
void
Intra16x16_Mode_Decision (Macroblock* currMB, int* i16mode)
{
  intrapred_luma_16x16 ();   /* make intra pred for all 4 new modes */
  find_sad_16x16 (i16mode); //19408  /* get best new intra mode */
  currMB->cbp = dct_luma_16x16 (*i16mode);//15
}



/*! 
 *************************************************************************************
 * \brief
 *    Sets Coefficients and reconstruction for an 8x8 block
 *************************************************************************************
 */
void SetCoeffAndReconstruction8x8 (Macroblock* currMB)
{
  int block, k, j, i;
  int block_y = img->block_y;
  int **ipredmodes = img->ipredmode;

  //--- restore coefficients ---
  for (block=0; block<6; block++)
  for (k=0; k<4; k++)
  for (j=0; j<2; j++)
  for (i=0; i<65; i++)
    img->cofAC[block][k][j][i] = cofAC8x8[block][k][j][i];

  if (cnt_nonz_8x8<=5 && img->type!=SP_SLICE)
  {
    currMB->cbp     = 0;
    currMB->cbp_blk = 0;
    for (j=0; j<16; j++)
    for (i=0; i<16; i++)  enc_picture->imgY[img->pix_y+j][img->pix_x+i] = mpr8x8[j][i];
  }
  else
  {
    currMB->cbp     = cbp8x8;
    currMB->cbp_blk = cbp_blk8x8;
    for (j=0; j<16; j++)
    for (i=0; i<16; i++)  enc_picture->imgY[img->pix_y+j][img->pix_x+i] = rec_mbY8x8[j][i];
  }

  //===== restore intra prediction modes for 8x8+ macroblock mode =====
  for (k=0, j=block_y; j<block_y+4; j++)
  for (     i=img->block_x; i<img->block_x+4; i++, k++)
  {
    ipredmodes    [i][j] = b8_ipredmode       [k];
    currMB->intra_pred_modes[k] = b8_intra_pred_modes[k];
  }
}


/*! 
 *************************************************************************************
 * \brief
 *    Sets motion vectors for an macroblock
 *************************************************************************************
 */
void SetMotionVectorsMB (Macroblock* currMB, int bframe)
{
  int i, j, k, l, m, mode8, pdir8, ref, by, bx, bxr;
  int ******all_mv  = img->all_mv;
  int ******pred_mv = img->pred_mv;
  int  bw_ref;

  for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      mode8 = currMB->b8mode[k=2*(j/2)+(i/2)];
      l     = 2*(j%2)+(i%2);
      by    = img->block_y+j;
      bxr   = img->block_x+i;
      bx    = img->block_x+i+4;


      
        pdir8 = currMB->b8pdir[k];
        ref    = enc_picture->ref_idx[LIST_0][bxr][by];
        bw_ref = enc_picture->ref_idx[LIST_1][bxr][by];
      
      
      if (!bframe)//不是b帧，其实这个函数没啥用吧，我只要把向量直接保存到enc_picture就好了吧。
      {
        if (pdir8>=0) //(mode8!=IBLOCK)&&(mode8!=I16MB))  // && ref != -1)
        {
          enc_picture->mv[LIST_0][bxr][by][0] = all_mv [i][j][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][bxr][by][1] = all_mv [i][j][LIST_0][ ref][mode8][1];
        }
        else
        {
          enc_picture->mv[LIST_0][bxr][by][0] = 0;
          enc_picture->mv[LIST_0][bxr][by][1] = 0;
        }
      }
      else
      {
        if (pdir8==-1) // intra
        {
          enc_picture->mv[LIST_0][bxr][by][0] = 0;
          enc_picture->mv[LIST_0][bxr][by][1] = 0;
          enc_picture->mv[LIST_1][bxr][by][0] = 0;
          enc_picture->mv[LIST_1][bxr][by][1] = 0;
        }
        else if (pdir8==0) // forward
        {
          enc_picture->mv[LIST_0][bxr][by][0] = all_mv [i][j][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][bxr][by][1] = all_mv [i][j][LIST_0][ ref][mode8][1];
          enc_picture->mv[LIST_1][bxr][by][0] = 0;
          enc_picture->mv[LIST_1][bxr][by][1] = 0;
        }
        else if (pdir8==1) // backward
        {
          enc_picture->mv[LIST_0][bxr][by][0] = 0;
          enc_picture->mv[LIST_0][bxr][by][1] = 0;
          
          enc_picture->mv[LIST_1][bxr][by][0] = all_mv[i][j][LIST_1][bw_ref][mode8][0];
          enc_picture->mv[LIST_1][bxr][by][1] = all_mv[i][j][LIST_1][bw_ref][mode8][1];
        }
        else if (pdir8==2) // bidir
        {
          enc_picture->mv[LIST_0][bxr][by][0] = all_mv [i][j][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][bxr][by][1] = all_mv [i][j][LIST_0][ ref][mode8][1];
          

          enc_picture->mv[LIST_1][bxr][by][0] = all_mv[i][j][LIST_1][bw_ref][mode8][0];
          enc_picture->mv[LIST_1][bxr][by][1] = all_mv[i][j][LIST_1][bw_ref][mode8][1];
        }
        else 
        {
          error("invalid direction mode", 255);
        }
      }
  }
  
  // copy all the motion vectors into rdopt structure
  // Can simplify this by copying the MV's of the best mode (TBD)
  if(img->MbaffFrameFlag)
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        for (k=0;k<2;k++)
          for(l=0;l<img->max_num_references;l++)
            for(m=0;m<9;m++)
            {
              rdopt->all_mv [i][j][k][l][m][0]  = all_mv [i][j][k][l][m][0];
              rdopt->pred_mv[i][j][k][l][m][0]  = pred_mv[i][j][k][l][m][0];
              
              rdopt->all_mv [i][j][k][l][m][1]  = all_mv [i][j][k][l][m][1];
              rdopt->pred_mv[i][j][k][l][m][1]  = pred_mv[i][j][k][l][m][1];
            }
  }
}



/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for a macroblock
 *************************************************************************************
 */
int RDCost_for_macroblocks (double   lambda,      // <-- lagrange multiplier
                        int      mode,        // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                        double*  min_rdcost)  // <-> minimum rate-distortion cost
{
  int         i, j, k; //, k, ****ip4;
  int         i16mode, rate=0, distortion=0;
  double      rdcost;
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  Macroblock  *prevMB   = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1] : NULL;
  int         bframe    = (img->type==B_SLICE);//0
  int         tmp_cc;
  int         use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);//0
  int         cc_rate, dummy;
  
  //=====
  //=====  SET REFERENCE FRAMES AND BLOCK MODES
  //=====
  SetModesAndRefframeForBlocks (mode);

  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  if (bframe && mode==0)
  {
    int block_x=img->pix_x>>2;
    int block_y=img->pix_y>>2;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        if (direct_pdir[block_x+i][block_y+j]<0)
          return 0;
  }

  if (mode<P8x8)
  {
    LumaResidualCoding ();
  }
  else if (mode==P8x8)
  {
    SetCoeffAndReconstruction8x8 (currMB);
  }
  else if (mode==I4MB)
  {
    currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda, &dummy);
  }
  else if (mode==I16MB)
  {
    Intra16x16_Mode_Decision  (currMB, &i16mode);
  }

  if (input->rdopt==2 && img->type!=B_SLICE)//0
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_mb (mode==I16MB?i16mode:-1);
  }

  //Rate control
  if (mode == I16MB)
  {
	  for(i=0; i<16; i++)
	  for(j=0; j<16; j++)
	    pred[j][i] = img->mprr_2[i16mode][j][i];
  }else
  {
	  for(i=0; i<16; i++)
	  for(j=0; j<16; j++)
	    pred[j][i] = img->mpr[i][j];
  }

  img->i16offset = 0;
  dummy = 0;
  ChromaResidualCoding (&dummy);
  if (mode==I16MB)     img->i16offset = I16Offset  (currMB->cbp, i16mode);//23


  //=====
  //=====   GET DISTORTION
  //=====
  // LUMA
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders ;k++)
    {
      decode_one_mb (k, currMB);
      for (j=0; j<MB_BLOCK_SIZE; j++)
      for (i=img->opix_x; i<img->opix_x+MB_BLOCK_SIZE; i++)
      {
        distortion += img->quad[imgY_org[img->opix_y+j][i]-decs->decY[k][img->opix_y+j][i] ];
      }
    }
    distortion /= input->NoOfDecoders;
  }
  else
  {
    for (j=0; j<16; j++)
    for (i=img->opix_x; i<img->opix_x+16; i++)
    {
      distortion += img->quad[imgY_org[j+img->opix_y][i]-enc_picture->imgY[img->pix_y+j][i] ];
    }
  }

  // CHROMA
  for (j=0; j<8; j++)
  for (i=img->opix_c_x; i<img->opix_c_x+8; i++)
  {
    distortion += img->quad[imgUV_org[0][j+img->opix_c_y][i] - enc_picture->imgUV[0][img->pix_c_y+j][i] ];
    distortion += img->quad[imgUV_org[1][j+img->opix_c_y][i] - enc_picture->imgUV[1][img->pix_c_y+j][i] ];
  }


  //=====   S T O R E   C O D I N G   S T A T E   =====
  //---------------------------------------------------
  store_coding_state (cs_cm);
  

  //=====
  //=====   GET RATE
  //=====
  //----- macroblock header -----
  if (use_of_cc)//0
  {
    if (currMB->mb_type!=0 || (bframe && currMB->cbp!=0))
    {
      // cod counter and macroblock mode are written ==> do not consider code counter
      tmp_cc = img->cod_counter;
      rate   = writeMBHeader (1); 
      ue_linfo (tmp_cc, dummy, &cc_rate, &dummy);
      rate  -= cc_rate;
      img->cod_counter = tmp_cc;
    }
    else
    {
      // cod counter is just increased  ==> get additional rate
      ue_linfo (img->cod_counter+1, dummy, &rate,    &dummy);
      ue_linfo (img->cod_counter,   dummy, &cc_rate, &dummy);
      rate -= cc_rate;
    }
  }
  else
  {
    rate = writeMBHeader(1); //12//53
  }
  if (mode)
  {
    //----- motion information -----
    rate  += writeMotionInfo2NAL  ();//12//53
  }
  if (mode || (bframe && (currMB->cbp!=0 || input->symbol_mode==CABAC)))
  {
    rate  += writeCBPandLumaCoeff ();//946//714

    rate  += writeChromaCoeff     ();//1081//849
  }


  //=====   R E S T O R E   C O D I N G   S T A T E   =====
  //-------------------------------------------------------
  reset_coding_state (cs_cm);


  rdcost = (double)distortion + lambda * (double)rate;//

  if (rdcost >= *min_rdcost)//45975.710614268930//37879
  {
    return 0;
  }


  if ((img->MbaffFrameFlag) && (mode ? 0: ((img->type == B_SLICE) ? !currMB->cbp:1)))  // AFF and current is skip
  {
    if (img->current_mb_nr%2) //bottom
    {
      if (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1)) //top is skip
      {
        if (!(field_flag_inference() == currMB->mb_field)) //skip only allowed when correct inference
          return 0;
      }
    }
  }

 
  //=====   U P D A T E   M I N I M U M   C O S T   =====
  //-----------------------------------------------------
  *min_rdcost = rdcost;
  return 1;
}





/*! 
 *************************************************************************************
 * \brief
 *    Store macroblock parameters
 *************************************************************************************
 */
void store_macroblock_parameters (int mode)
{
  int  i, j, k, ****i4p, ***i3p;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  int        bframe   = (img->type==B_SLICE);

  //--- store best mode ---
  best_mode = mode;//10
  best_c_imode = currMB->c_ipred_mode;//0
  best_i16offset = img->i16offset;//23
  for (i=0; i<4; i++)
  {
    b8mode[i]   = currMB->b8mode[i];
    b8pdir[i]   = currMB->b8pdir[i];
  }

  //--- reconstructed blocks ----
  for (j=0; j<16; j++)
  for (i=0; i<16; i++)
  {
    rec_mbY[j][i] = enc_picture->imgY[img->pix_y+j][img->pix_x+i];
  }
  for (j=0; j<8; j++)
  for (i=0; i<8; i++)
  {
    rec_mbU[j][i] = enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x+i];
    rec_mbV[j][i] = enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x+i];
  }


  //--- store results of decoders ---
  if (input->rdopt==2 && img->type!=B_SLICE)//0
  {
    for (k=0; k<input->NoOfDecoders; k++)
    {
      for (j=img->pix_y; j<img->pix_y+16; j++)
      for (i=img->pix_x; i<img->pix_x+16; i++)
      {
        // Keep the decoded values of each MB for updating the ref frames
        decs->decY_best[k][j][i] = decs->decY[k][j][i];
      }
    }
  }

  //--- coeff, cbp, kac ---
  if (mode || bframe)//1
  {
    i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
    i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
    cbp     = currMB->cbp;
    cbp_blk = currMB->cbp_blk;
  }
  else
  {
    cbp = cbp_blk = 0;
  }
//////////////////////////////////////////////////////////////////////////ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {  
    frefframe[j][i]=enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j];

    if (bframe)
    {
      brefframe[j][i]=enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j];
    }
  }

}


/*! 
 *************************************************************************************
 * \brief
 *    Set stored macroblock parameters
 *************************************************************************************
 */
void set_stored_macroblock_parameters ()
{
  int  i, j, k, ****i4p, ***i3p,l,y,x;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode;
  int         bframe   = (img->type==B_SLICE);
  int         **ipredmodes = img->ipredmode;
  
  byte        **imgY  = enc_picture->imgY;
  byte       ***imgUV = enc_picture->imgUV;

  int        list_offset   = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;

  //===== reconstruction values =====
  for (j=0; j<16; j++)
  for (i=0; i<16; i++)
  {
    imgY[img->pix_y+j][img->pix_x+i] = rec_mbY[j][i];

    if(img->MbaffFrameFlag)//0
      rdopt->rec_mbY[j][i]       = rec_mbY[j][i]; 
  }
  for (j=0; j<8; j++)
  for (i=0; i<8; i++)
  {
    imgUV[0][img->pix_c_y+j][img->pix_c_x+i] = rec_mbU[j][i];
    imgUV[1][img->pix_c_y+j][img->pix_c_x+i] = rec_mbV[j][i];

    if(img->MbaffFrameFlag)//0
    {
      rdopt->rec_mbU[j][i]           = rec_mbU[j][i];   
      rdopt->rec_mbV[j][i]           = rec_mbV[j][i];   
    }
  }
//   for (y=0; y<1; y++)
//   {
// 	  for (x=0; x<16; x++)
// 	  {
// 		  if (x%2==0 && x%8!=0)
// 		  {
// 			  printf("%d(%d,%d)   ",imgY[y][x],imgUV[0][y/2][x/2],imgUV[1][y/2][x/2]);
// 		  } 
// 		  else if (x%2!=0 && x%8!=0)
// 		  {
// 			  printf("%d  ",imgY[y][x]);
// 		  }
// 		  else
// 		  {
// 			  printf("%d(%d,%d)\n",imgY[y][x],imgUV[0][y/2][x/2],imgUV[1][y/2][x/2]);
// 		  }
// 	  }
//   }
//   printf("\n");

  //===== coefficients and cbp =====
  i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
  i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
  currMB->cbp      = cbp;
  currMB->cbp_blk = cbp_blk;
  //==== macroblock type ====
  currMB->mb_type = mode;


  if(img->MbaffFrameFlag)//0
  {
    rdopt->mode = mode;
    rdopt->i16offset = img->i16offset;
    rdopt->cbp = cbp;
    rdopt->cbp_blk = cbp_blk;
    rdopt->mb_type  = mode;

    rdopt->prev_qp=currMB->prev_qp;
    rdopt->prev_delta_qp=currMB->prev_delta_qp;
    rdopt->qp=currMB->qp;

    for(i=0;i<6;i++)
      for(j=0;j<4;j++)
        for(k=0;k<2;k++)
          for(l=0;l<18;l++)
            rdopt->cofAC[i][j][k][l] =img->cofAC[i][j][k][l];

    for(i=0;i<3;i++)
        for(k=0;k<2;k++)
          for(l=0;l<18;l++)
            rdopt->cofDC[i][k][l]=img->cofDC[i][k][l];
  }


  for (i=0; i<4; i++)
  {
    currMB->b8mode[i]   = b8mode[i];
    currMB->b8pdir[i]   = b8pdir[i];
    if(img->MbaffFrameFlag)
    {
      rdopt->b8mode[i]  = b8mode[i];                  
      rdopt->b8pdir[i]  = b8pdir[i];                  
    }
  }

  if (input->rdopt==2 && img->type!=B_SLICE)//0
  {
    //! save the MB Mode of every macroblock
    decs->dec_mb_mode[img->mb_x][img->mb_y] = mode;
  }

  //==== reference frames =====
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    // backward prediction or intra
    if ((currMB->b8pdir[i/2+(j/2)*2] == 1) || IS_INTRA(currMB))
    {
      enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = -1;
      enc_picture->ref_pic_id [LIST_0][img->block_x+i][img->block_y+j] = -1;

      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0] =0;
      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1] = 0;
      if(img->MbaffFrameFlag)
        rdopt->refar[LIST_0][j][i] = -1;
    }
    else
    {
      enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = frefframe[j][i]; 
      enc_picture->ref_pic_id [LIST_0][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]];

      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0] = img->all_mv[i][j][LIST_0][frefframe[j][i]][currMB->b8mode[i/2+(j/2)*2]][0];
      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1] = img->all_mv[i][j][LIST_0][frefframe[j][i]][currMB->b8mode[i/2+(j/2)*2]][1];
      if(img->MbaffFrameFlag)
        rdopt->refar[LIST_0][j][i] = frefframe[j][i];
    }

    // forward prediction or intra
    if ((currMB->b8pdir[i/2+(j/2)*2] == 0) || IS_INTRA(currMB))
    {
      enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
      enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+j] = -1;  
      enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0] =0;
      enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1] = 0;
      if(img->MbaffFrameFlag)
        rdopt->refar[LIST_1][j][i] = -1;
    }
    
  }

  if (bframe)
  {
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        
        // forward
        if (IS_INTRA(currMB)||(currMB->b8pdir[i/2+(j/2)*2] == 0))
        {
          enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
          enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+j] = -1;
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0] = 0;
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1] = 0;
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_1][j][i] = -1;
        }
        else
        {
          enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = brefframe[j][i];
          enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j]];
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0] = img->all_mv[i][j][LIST_1][brefframe[j][i]][currMB->b8mode[i/2+(j/2)*2]][0];
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1] = img->all_mv[i][j][LIST_1][brefframe[j][i]][currMB->b8mode[i/2+(j/2)*2]][1];
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_1][j][i] = brefframe[j][i];
        }
      }
  }

  //==== intra prediction modes ====
  currMB->c_ipred_mode = best_c_imode;
  img->i16offset = best_i16offset;
  if (mode==P8x8)
  {
    for (k=0, j=img->block_y; j<img->block_y+4; j++)
      for (   i=img->block_x; i<img->block_x+4; i++, k++)
      {
        ipredmodes           [i][j] = b8_ipredmode       [k];
        currMB->intra_pred_modes[k] = b8_intra_pred_modes[k];
      }
  }
  else if (mode!=I4MB)
  {
    for (k=0, j=img->block_y; j<img->block_y+4; j++)
      for (   i=img->block_x; i<img->block_x+4; i++, k++)
      {
        ipredmodes           [i][j] = DC_PRED;
        currMB->intra_pred_modes[k] = DC_PRED;
      }
  }

  if(img->MbaffFrameFlag)
  {
    for (k=0, j=img->block_y; j<img->block_y+4; j++)
      for (   i=img->block_x; i<img->block_x+4; i++, k++)
      {
        rdopt->ipredmode[i][j]     = ipredmodes[i][j];
        rdopt->intra_pred_modes[k] = currMB->intra_pred_modes[k];
      }
    rdopt->c_ipred_mode = currMB->c_ipred_mode;
    rdopt->i16offset = img->i16offset;  // DH
  }

  //==== motion vectors =====
  SetMotionVectorsMB (currMB, bframe);
}



/*! 
 *************************************************************************************
 * \brief
 *    Set reference frames and motion vectors
 *************************************************************************************
 */
void SetRefAndMotionVectors (int block, int mode, int pdir, int fwref, int bwref)
{
  int     i, j=0;
  int     bslice  = (img->type==B_SLICE);
  int     pmode   = (mode==1||mode==2||mode==3?mode:4);
  int     j0      = ((block/2)<<1);
  int     i0      = ((block%2)<<1);
  int     j1      = j0 + (input->blc_size[pmode][1]>>2);
  int     i1      = i0 + (input->blc_size[pmode][0]>>2);

  int list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;


  if (pdir<0)
  {
    for (j=j0; j<j1; j++)
    for (i=i0; i<i1; i++)
    {
      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0] = 0;
      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1] = 0;
      enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0] = 0;
      enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1] = 0;
      enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = -1;
      enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
      enc_picture->ref_pic_id[LIST_0][img->block_x+i][img->block_y+j] = -1;
      enc_picture->ref_pic_id[LIST_1][img->block_x+i][img->block_y+j] = -1;
    }
    return;
  }

  if (!bslice)
  {
    for (j=j0; j<j1; j++)
    for (i=i0; i<i1; i++)
    {
      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0] = img->all_mv[i][j][LIST_0][fwref][mode][0];
      enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1] = img->all_mv[i][j][LIST_0][fwref][mode][1];
      enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = fwref;
      enc_picture->ref_pic_id[LIST_0][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_0+list_offset][enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]];
    }
  }
  else
  {
    for (j=j0; j<j1; j++)
      for (i=i0; i<i1; i++)
      {
        if (mode==0)
        {
          pdir  =direct_pdir[img->block_x+i][img->block_y+j];
          fwref =direct_ref_idx[LIST_0][img->block_x+i][img->block_y+j];
          bwref =direct_ref_idx[LIST_1][img->block_x+i][img->block_y+j];
        }
        
        if ((pdir==0 || pdir==2))
        {
          enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0] = img->all_mv[i][j][LIST_0][fwref][mode][0];
          enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1] = img->all_mv[i][j][LIST_0][fwref][mode][1];
          enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = fwref;
          enc_picture->ref_pic_id[LIST_0][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_0+list_offset][enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j]];
        }
        else
        {
          enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0] = 0;
          enc_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1] = 0;
          enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = -1;
          enc_picture->ref_pic_id[LIST_0][img->block_x+i][img->block_y+j] = -1;
        }

        if ((pdir==1 || pdir==2))
        {
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0] = img->all_mv[i][j][LIST_1][bwref][mode][0];
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1] = img->all_mv[i][j][LIST_1][bwref][mode][1];
          enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = bwref;
          enc_picture->ref_pic_id[LIST_1][img->block_x+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_1+list_offset][enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j]];
        }
        else
        {
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0] = 0;
          enc_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1] = 0;
          enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
          enc_picture->ref_pic_id[LIST_1][img->block_x+i][img->block_y+j] = -1;
        }
      }
  }
}

/*! 
 *************************************************************************************
 * \brief
 *    skip macroblock field inference
 * \return
 *    inferred field flag
 *************************************************************************************
 */
int field_flag_inference()
{
  int mb_field;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if (currMB->mbAvailA)
  {
    mb_field = img->mb_data[currMB->mbAddrA].mb_field;
  }
  else
  {
    // check top macroblock pair
    if (currMB->mbAvailB)
    {
      mb_field = img->mb_data[currMB->mbAddrB].mb_field;
    }
    else
      mb_field = 0;
  }

  return mb_field;
}

/*! 
 *************************************************************************************
 * \brief
 *    Mode Decision for a macroblock
 *    //+++该函数的作用是编码一个宏块（包括帧间、帧内、帧内预测的方式）。
 *    NOTE:从上面程序段中可以看出JM8.5中对7种宏块模式是采用全部遍历的方式，所以导致的计算复杂度很高。
 *************************************************************************************
 */
 void i_encode_one_macroblock ()
 {
   static const int  b8_mode_table[6]  = {0, 4, 5, 6, 7};         // DO NOT CHANGE ORDER !!! B片子宏块类型。0：Direct_8*8，4：8*8，5：8*4，6：4*8，7：4*4
   static const int  mb_mode_table[7]  = {0, 1, 2, 3, P8x8, I16MB, I4MB}; // DO NOT CHANGE ORDER !!!
   //+++//   0:16X16 Direct模式,在B帧中有效                                                                                                                                  
   //        1:Inter16X16,在帧间有效                            
   //        2:Inter16X8,在帧间有效                                
   //        3:Inter8X16,在帧间有效                                       
   //        亚宏块级(P8x8):帧间有效:包括Inter8x8,Inter8x4,Inter4x8,Inter4x4                                                 
   // I16MB:Intra16X16帧内有效                                                         
   // I4MB:Intra有效                                                     
   // I8MB:Intra有效                                                         
   // IPCM:Intra有效，不要预测，直接对RAW数据编码.   
   int         valid[MAXMODE];
   int         rerun,  index, mode,   i, j, k, ctr16x16, pdir,ref,i0, i1, j0, j1,block, dummy;// 
   double      qp, lambda_mode, lambda_motion, min_rdcost, rdcost = 0, max_rdcost=1e30;
   int         lambda_motion_factor;
   int          max_mcost=(1<<30),mcost, bw_mcost, bid_mcost,fw_mcost;//, 
   int          cnt_nonz = 0, best_cnt_nonz = 0, best_fw_ref = 0, best_pdir,curr_cbp_blk;//,
   int         cost=0;
   int         min_cost = max_mcost,   cost_direct=0, have_direct=0, min_cost8x8,cost8x8, i16mode;//
   //+++根据对下面代码的阅读,这儿有比较要说明一下上面声明的一些变量,这对于理解代码很重要.  
   //+++首先,要提到best_mode,我本来以为这个量是在本函数中声明的一个变量,然后仔细在上面寻找,都没有找到,后来才发现, 
   //+++wow,原来best_mode在这个文件(rdopt.c)的开头声明过,是本文件的一个全局变量,yes,终于明白这个了  
   //+++其次,对于上面的一些变量,介绍一下:  
   //+++min_rdcost是用来保存RDO模式下的最小代价;min_cost则是非RDO模式下的代价  
   //+++min_cost8x8是亚宏块级(P8x8)非RDO模式下最小模式代价;cost8x8是宏块(16x16)中4个8x8块的代价和(每一个8x8都是最佳代价min_cost8x8)  
   //+++fw_mcost,bw_mcost和bid_mcost是中间变量用来记录(前向,后向和双向)最佳参考帧对应的代价  
   //+++cost是一个中间变量,主要是在非RDO模式下,(帧间预测)记录宏块级或亚宏块级中某一个模式中所有块(16x16是1块,16x8是2块...)的所有代价之和(参考帧代价,运动估计)  
   //+++mcost也是一个中间变量,也是非RDO模式下,在进行最优参考帧时,用于累加代价之和(参考帧代价,运动估计)  

   int         intra1 = 0;
   
   int         intra       = (((img->type==P_SLICE||img->type==SP_SLICE) && img->mb_y==img->mb_y_upd && img->mb_y_upd!=img->mb_y_intra) || img->type==I_SLICE);
   int         spframe     = (img->type==SP_SLICE);
   int         siframe     = (img->type==SI_SLICE);
   int         bframe      = (img->type==B_SLICE);
   int         runs        = (input->RestrictRef==1 && input->rdopt==2 && (img->type==P_SLICE || img->type==SP_SLICE || (img->type==B_SLICE && img->nal_reference_idc>0)) ? 2 : 1);
   
   int         checkref    = (input->rdopt && input->RestrictRef && (img->type==P_SLICE || img->type==SP_SLICE));
   Macroblock* currMB      = &img->mb_data[img->current_mb_nr];//当前宏块
   Macroblock* prevMB      = img->current_mb_nr ? &img->mb_data[img->current_mb_nr-1]:NULL ;//前一个宏块
   
   int     **ipredmodes = img->ipredmode;//+++ ipredmode 帧内预测模式
   int     best_bw_ref  = -1;
   int     ******allmvs = img->all_mv;

   
   int  l,list_offset;

   int curr_mb_field = ((img->MbaffFrameFlag)&&(currMB->mb_field));

   // find out the correct list offsets
   if (curr_mb_field)//no
   {
     if(img->current_mb_nr%2)
       list_offset = 4; // bottom field mb
     else
       list_offset = 2; // top field mb
   }
   else//yes
   {
     list_offset = 0;  // no mb aff or frame mb
   }

   if(input->FMEnable)//no
     decide_intrabk_SAD();
   
   intra |= RandomIntra (img->current_mb_nr);    // Forced Pseudo-Random Intra

   //===== SET VALID MODES =====
   valid[I4MB]   = 1;//valid[9]++ 宏块可以采用intra4*4模式编码
   valid[I16MB]  = 1;//valid[10]++ 宏块可以采用intra16*16模式编码
   
   valid[0]      = (!intra );//++ 宏块是否采用inter模式编码
   valid[1]      = (!intra && input->InterSearch16x16);//++ 宏块是否采用inter16*16模式编码
   valid[2]      = (!intra && input->InterSearch16x8);
   valid[3]      = (!intra && input->InterSearch8x16);
   valid[4]      = (!intra && input->InterSearch8x8);
   valid[5]      = (!intra && input->InterSearch8x4);
   valid[6]      = (!intra && input->InterSearch4x8);
   valid[7]      = (!intra && input->InterSearch4x4);
   valid[P8x8]   = (valid[4] || valid[5] || valid[6] || valid[7]);	//++ 宏块是否采用亚宏块模式编码
   valid[12]     = (siframe);	//++ 宏块是否采用SI模式编码

   if (!img->MbaffFrameFlag)//+++//在非宏块级帧场自适应且是场模式下 进行色度矢量校正//yes
   {
     for (l=0+list_offset;l<(2+list_offset);l++)//遍历list0,list1
     {
       for(k = 0; k < listXsize[l]; k++)
       {
         listX[l][k]->chroma_vector_adjustment= 0;
         if(img->structure == TOP_FIELD && img->structure != listX[l][k]->structure)
           listX[l][k]->chroma_vector_adjustment = -2;
         if(img->structure == BOTTOM_FIELD && img->structure != listX[l][k]->structure)
           listX[l][k]->chroma_vector_adjustment = 2;
       }
     }
   }
   else
   {
     if (curr_mb_field)
     {
       for (l=0+list_offset;l<(2+list_offset);l++)
       {
         for(k = 0; k < listXsize[l]; k++)
         {
           listX[l][k]->chroma_vector_adjustment= 0;
           if(img->current_mb_nr % 2 == 0 && listX[l][k]->structure == BOTTOM_FIELD)
             listX[l][k]->chroma_vector_adjustment = -2;
           if(img->current_mb_nr % 2 == 1 && listX[l][k]->structure == TOP_FIELD)
             listX[l][k]->chroma_vector_adjustment = 2;
         }
       }
     }
     else//no
     {
       for (l=0+list_offset;l<(2+list_offset);l++)
       {
         for(k = 0; k < listXsize[l]; k++)
         {
           listX[l][k]->chroma_vector_adjustment= 0;
         }
       }
     }
     
   }
   //===== SET LAGRANGE PARAMETERS =====//+++ 设置RDO与非RDO下的lambda系数  
   if (input->rdopt)//yes
   {
     qp = (double)img->qp - SHIFT_QP;//16=28-12//输入的QP减去12

     if (input->successive_Bframe>0)//no
       lambda_mode   = 0.68 * pow (2, qp/3.0) * (img->type==B_SLICE? max(2.00,min(4.00,(qp / 6.0))):spframe?max(1.4,min(3.0,(qp / 12.0))):1.0);  
     //27.415882045712436
	 else
       lambda_mode   = 0.85 * pow (2, qp/3.0) * (img->type==B_SLICE? 4.0:spframe?max(1.4,min(3.0,(qp / 12.0))):1.0);  //34.27

     lambda_motion = sqrt (lambda_mode);//5.85
   }
   else//no
   {
     lambda_mode = lambda_motion = QP2QUANT[max(0,img->qp-SHIFT_QP)];
   }
   lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);//383651	//++ 使用该变量的原因主要是为了在运算过程中使用定点运算但又能提高精度
   
   
   for (rerun=0; rerun<runs; rerun++)//runs=1	//++ runs=2是针对loss rdo模式，其余的情况runs=1
   {
     if (runs==2)
     {
       if (rerun==0)   input->rdopt=1;
       else            input->rdopt=2;
     }
     
     // reset chroma intra predictor to default
     currMB->c_ipred_mode = DC_PRED_8;
	 //+++++++++++++++++++++++++++  
	 //+++帧间预测,计算其选择合适模式后的开销  
	 //+++(根据拉格朗日率失真准则)      
     if (!intra)//no
     {      
	   //+++//B帧的直接模式预测 
       //===== set direct motion vectors =====
       if (bframe)//+++P帧无direct模式 但有skip模式
       {
         Get_Direct_Motion_Vectors ();//+++不返回direct模式代价 
       }
       //+++//宏块级运动估计
       //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====
       for (min_cost=1<<20, best_mode=1, mode=1; mode<4; mode++)
       {
         if (valid[mode])//+++//对应于程序外部（即CFG文件中）的设置
         { //+++//对于16x16，MB只分一个块；对于16x8和8x16，MB被分成两个块  
           for (cost=0, block=0; block<(mode==1?1:2); block++)//++ 16*16 分割方式只需要计算一次，16*8 和 8*16 分割方式需要计算两次
           {  //+++/块匹配!!!lambda_motion用来求运动矢量消耗的码率
             PartitionMotionSearch (mode, block, lambda_motion);//++ 对 16*16、16*8、8*16 分割方式进行 ME

             //--- set 4x4 block indizes (for getting MV) ---
             j = (block==1 && mode==2 ? 2 : 0);//++ 如果现在处理的是 16*8 分割方式的第 2 个分割，则 j=2
             i = (block==1 && mode==3 ? 2 : 0);//++ 如果现在处理的是 8*16 分割方式的第 2 个分割，则 i=2
             //+++通过循环在list0中寻找前向预测的最佳参考帧(best_fw_ref),代价为fw_mcost,
             //--- get cost and reference frame for forward prediction ---
             for (fw_mcost=max_mcost, ref=0; ref<listXsize[LIST_0+list_offset]; ref++)
             {
               if (!checkref || ref==0 || CheckReliabilityOfRef (block, LIST_0, ref, mode))
               {//++ 参考帧索引号编码所需要使用的比特数作为编码代价的一部分
				   //+++这个地方就是利用拉格朗日进行计算开销的函数,事实上,这个地方要考虑很多东西的开销,可以发现,下面的代码中还要加上其他的  
				   //+++mcost->(1)REF_COST参考帧代价  
                 mcost  = (input->rdopt ? REF_COST (lambda_motion_factor, ref, LIST_0 + list_offset) : (int)(2*lambda_motion*min(ref,1)));
				 //+++motion_cost(运动搜索代价)是一个全局变量,其值在PartitionMotionSearch函数内进行过改变  
				 //+++so 此时motion_cost就是保存的相应的运动搜索代价  
				 //+++那么,mcost->(2)motion_cost运动搜索代价 
				 mcost += motion_cost[mode][LIST_0][ref][block];//+++可以发现motion_cost数组的参数mode,block和  
			                                            	 //+++PartitionMotionSearch(PMS)函数的参数是一致的,并且在PMS函数中  
			                                              //+++计算motion_cost时也是对ref(参考帧)进行循环的,所以将两者联系起来  
				                                      //+++我们就可以理解这儿motion_cost的数组中参数的含义了.  
				                                  //+++另外在[JM86核心函数研究]中有对motion_cost的详细介绍  
			                             	 //+++该数组存储了一个宏块的分块在不同帧间模式和参考帧下的运动估计的开销值。  
			                        	 //+++函数PartitionMotionSearch的主要任务就是更新该数组。同时,该数组也在决定最佳参考帧时被引用。
                 if (mcost < fw_mcost)//+++根据参考代价选出该模式下的最佳参考帧
                 {   //+++记录最佳参考帧及预测代价
                   fw_mcost    = mcost;
                   best_fw_ref = ref;
                 }
               }
             }
             //+++B帧还需记录向后及双向预测的各个参考帧与参考代价 
			 
			 //可以简化else是非b帧
             if (bframe)
             {
               //--- get cost for bidirectional prediction ---
               for (bw_mcost=max_mcost, ref=0; ref<listXsize[LIST_1 + list_offset]; ref++)
               {//+++B帧的计算方法与I,P帧不一样 
                 mcost  = (input->rdopt ? REF_COST (lambda_motion_factor, ref, LIST_1 + list_offset) : (int)(2*lambda_motion*min(ref,1)));
                 mcost += motion_cost[mode][LIST_1][ref][block];
                 if (mcost < bw_mcost)
                 {
                   bw_mcost    = mcost;
                   best_bw_ref = ref;
                 }
               }

               // search bidirectional between best forward and ref_idx=0 backward
               bid_mcost  = (input->rdopt ? (REF_COST (lambda_motion_factor, best_fw_ref,LIST_0+list_offset)+REF_COST (lambda_motion_factor, 0,LIST_1+list_offset)) : (int)(2*lambda_motion*min(best_fw_ref,1)));
               bid_mcost += BIDPartitionCost (mode, block, best_fw_ref, 0, lambda_motion_factor);

               //--- get prediction direction ----
               if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
               {
                 best_pdir = 0;
                 best_bw_ref = 0;
                 cost += fw_mcost;
               }
               else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
               {
                 best_pdir = 1;
                 cost += bw_mcost;
                 best_fw_ref = 0;
               }
               else
               {
                 best_pdir = 2;
                 cost += bid_mcost;
                 best_bw_ref = 0;
               }
             }
             else // if (bframe)
             {
				 best_pdir  = 0;//+++best_pdir是最佳预测方向,对于非B帧(P帧)预测方向就一个直接设置为0,而对于  
				                //+++B帧要从三种方式(前向,后向,双向)中选择一个,从上面if中的代码中我们就可以看出来 
               cost      += fw_mcost;
             }
			 //+++下面要记录参考帧和MV信息,说明一下:这些全局变量是作为BlockMotionSearch()中运动搜索时运动向量预测的依据  
			 //+++问:下一宏块的运动向量等信息也用当前宏块的运动向量预测吗？答:ME 算法是开放部分，你想怎么做就怎么做。JM 在 ME 时会用到相邻块的 MV；  
			 //+++对于MB是以4*4的块为单位来记录mv和ref的 
             if (mode==1)//+++16x16
             {
               if (best_pdir==1)//+++后向参考
               {
                 for (j=0; j<4; j++)
                 {
                   for (i=0; i<4; i++)
				   {   //+++保存参考信息至重建单元，B帧多做后项参考  
					   //+++enc_picture 中好像是最终数据，img 中存的是中间数据 
                     enc_picture->ref_idx[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = -1;
                     enc_picture->ref_pic_id [LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = -1; 
                     enc_picture->mv[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][0] = 0;
                     enc_picture->mv[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][1] = 0;
                   }
                 }
               }
               else
               {
                 for (j=0; j<4; j++)
                 {
                   for (i=0; i<4; i++)
                   {
                     enc_picture->ref_idx[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = best_fw_ref;
                     enc_picture->ref_pic_id [LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j]];  
                     enc_picture->mv[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][0] = img->all_mv[i][j][LIST_0][best_fw_ref][mode][0];
                     enc_picture->mv[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][1] = img->all_mv[i][j][LIST_0][best_fw_ref][mode][1];
                   }
                 }
               }

               if (bframe)
               {
                 if (best_pdir==0)
                 {
                   for (j=0; j<4; j++)
                   {
                     for (i=0; i<4; i++)
                     {
                       enc_picture->ref_idx[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = -1;
                       enc_picture->ref_pic_id [LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = -1;   
                       enc_picture->mv[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][0] = 0;
                       enc_picture->mv[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][1] = 0;
                     }
                   }
                 }
                 else
                 {
                   for (j=0; j<4; j++)
                   {
                     for (i=0; i<4; i++)
                     {
                       enc_picture->ref_idx[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = best_bw_ref;
                       enc_picture->ref_pic_id [LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j]];
                       if(best_bw_ref>=0)
                       {
                         enc_picture->mv[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][0] = img->all_mv[i][j][LIST_1][best_bw_ref][mode][0];
                         enc_picture->mv[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j][1] = img->all_mv[i][j][LIST_1][best_bw_ref][mode][1];
                       }
                     }
                   }
                 }
               }
             }
             else if (mode==2)
             {
               for (j=0; j<2; j++)
               {
                 for (i=0; i<4; i++)
                 {
                   if (best_pdir==1)
                   {
                     enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+block*2+j] = -1;
                     enc_picture->ref_pic_id [LIST_0][img->block_x+i][img->block_y+block*2+j] = -1;
                     enc_picture->mv[LIST_0][img->block_x+i][img->block_y+block*2+j][0] = 0;
                     enc_picture->mv[LIST_0][img->block_x+i][img->block_y+block*2+j][1] = 0;
                   }
                   else
                   {
                     enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+block*2+j] = best_fw_ref;
                     enc_picture->ref_pic_id [LIST_0][img->block_x+i][img->block_y+block*2+j] = enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+block*2+j]];
                     enc_picture->mv[LIST_0][img->block_x+i][img->block_y+block*2+j][0] = img->all_mv[i][j+block*2][LIST_0][best_fw_ref][mode][0];
                     enc_picture->mv[LIST_0][img->block_x+i][img->block_y+block*2+j][1] = img->all_mv[i][j+block*2][LIST_0][best_fw_ref][mode][1];
                   }

                   if (bframe)
                   {
                     if (best_pdir==0)
                     {
                       enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+block*2+j] = -1;
                       enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+block*2+j] = -1;
                       enc_picture->mv[LIST_1][img->block_x+i][img->block_y+block*2+j][0] = 0;
                       enc_picture->mv[LIST_1][img->block_x+i][img->block_y+block*2+j][1] = 0;
                     }
                     else
                     {
                       enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+block*2+j] = best_bw_ref;
                       enc_picture->ref_pic_id [LIST_1][img->block_x+i][img->block_y+block*2+j] = enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+block*2+j]];
                       if(best_bw_ref>=0)
                       {
                         enc_picture->mv[LIST_1][img->block_x+i][img->block_y+block*2+j][0] = img->all_mv[i][j][LIST_1][best_bw_ref][mode][0];
                         enc_picture->mv[LIST_1][img->block_x+i][img->block_y+block*2+j][1] = img->all_mv[i][j][LIST_1][best_bw_ref][mode][1];
                       }
                     }
                   }
                 }
               }
             }
             else
             {
               for (j=0; j<4; j++)
               {
                 for (i=0; i<2; i++)
                 {
                   if (best_pdir==1)
                   {
                     enc_picture->ref_idx[LIST_0][img->block_x+block*2+i][img->block_y+j] = -1;
                     enc_picture->ref_pic_id [LIST_0][img->block_x+block*2+i][img->block_y+j] = -1;
                     enc_picture->mv[LIST_0][img->block_x+block*2+i][img->block_y+j][0] = 0;
                     enc_picture->mv[LIST_0][img->block_x+block*2+i][img->block_y+j][1] = 0;
                   }
                   else
                   {
                     enc_picture->ref_idx[LIST_0][img->block_x+block*2+i][img->block_y+j] = best_fw_ref;
                     enc_picture->ref_pic_id [LIST_0][img->block_x+block*2+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+block*2+i][img->block_y+j]];
                     enc_picture->mv[LIST_0][img->block_x+block*2+i][img->block_y+j][0] = img->all_mv[i][j][LIST_0][best_fw_ref][mode][0];
                     enc_picture->mv[LIST_0][img->block_x+block*2+i][img->block_y+j][1] = img->all_mv[i][j][LIST_0][best_fw_ref][mode][1];
                   }

                   if (bframe)
                   {
                     if (best_pdir==0)
                     {
                       enc_picture->ref_idx[LIST_1][img->block_x+block*2+i][img->block_y+j] = -1;
                       enc_picture->ref_pic_id [LIST_1][img->block_x+block*2+i][img->block_y+j] = -1;
                       enc_picture->mv[LIST_1][img->block_x+block*2+i][img->block_y+j][0] = 0;
                       enc_picture->mv[LIST_1][img->block_x+block*2+i][img->block_y+j][1] = 0;
                     }
                     else
                     {
                       enc_picture->ref_idx[LIST_1][img->block_x+block*2+i][img->block_y+j] = best_bw_ref;
                       enc_picture->ref_pic_id [LIST_1][img->block_x+block*2+i][img->block_y+j] = enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+block*2+i][img->block_y+j]];
                       if(best_bw_ref>=0)
                       {
                         enc_picture->mv[LIST_1][img->block_x+block*2+i][img->block_y+j][0] = img->all_mv[i][j][LIST_1][best_bw_ref][mode][0];
                         enc_picture->mv[LIST_1][img->block_x+block*2+i][img->block_y+j][1] = img->all_mv[i][j][LIST_1][best_bw_ref][mode][1];
                       }
                     }
                   }
                 }
               }
             }
			 //+++以8x8块记录各个模式的预测信息 每个分割的8x8相同  
			 //+++此处解释一下best8x8fwref:： 该2维数组存储了一个宏块的4个8x8block在各个帧间预测模式下的最佳前向参考帧。  
			 //+++上面的代码完成了最佳参考帧的选择，所以这儿将最佳参考帧的信息存储在这个数组中。在计算预测值时,程序会引用该数组作为参考。  
             //----- set reference frame and direction parameters -----
             if (mode==3)//+++两个16x8:如果以8x8为单位,分别拥有(0,2)和(1,3) 
             {
               best8x8fwref   [3][block] = best8x8fwref   [3][block+2] = best_fw_ref;
               best8x8pdir    [3][block] = best8x8pdir    [3][block+2] = best_pdir;
               best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;
             }
             else if (mode==2)//+++两个8x16:如果以8x8为单位,分别拥有(0,1)和(2,3) 
             {
               best8x8fwref   [2][2*block] = best8x8fwref   [2][2*block+1] = best_fw_ref;
               best8x8pdir    [2][2*block] = best8x8pdir    [2][2*block+1] = best_pdir;
               best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;
             }
             else//+++一个16x16:如果以8x8为单位,(0,1,2,3)的设置是一样的 
             {
               best8x8fwref   [1][0] = best8x8fwref   [1][1] = best8x8fwref   [1][2] = best8x8fwref   [1][3] = best_fw_ref;
               best8x8pdir    [1][0] = best8x8pdir    [1][1] = best8x8pdir    [1][2] = best8x8pdir    [1][3] = best_pdir;
               best8x8bwref   [1][0] = best8x8bwref   [1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = best_bw_ref;
             }

             //--- set reference frames and motion vectors ---
			 if (mode>1 && block==0)  //+++mode>1,说明mode=2和3的时候执行,block==0,说明是对第一个分块.  
				 //+++这里说明一下:对于16x8，8x16来说，只有两个宏块分割,这里保存第一块分割的运动信息，目的是给下一个分割做运动矢量预测使用  
				 //+++由于前一个分割不需要下一个分割的运动矢量信息，因此只需要保存前一个分割即可.  
				 //+++另外，在这之前，事实上已经保存过，所以这里是重复保存了  
               SetRefAndMotionVectors (block, mode, best_pdir, best_fw_ref, best_bw_ref);

          } // for (block=0; block<(mode==1?1:2); block++)
		  //+++我在这个地方解释一下cost与min_cost:min_cost是为了寻找三种宏块模式(16x16,16x8,8x16)中代价最小的模式  
		  //+++由于宏块模式在计算的时候是进一步分成小块计算的,所以cost就是用来保存一种宏块模式下对应所有小块的代价之和  
		  //+++模式16x16下:cost对应一个16x16大小的块;模式16x8模式下:cost对应两个16x8大小的块......  
		  //+++所以下面这一个if语句就是在寻找代价最小的宏块模式(best_mode)  
          if (cost < min_cost)
          {
            best_mode = mode;
            min_cost  = cost;
          }
        } // if (valid[mode])
	  } // for (mode=1; mode<4; mode++)
	  //+++宏块级模式结束  
	  //+++现在宏块级运动估计结束,best_mode保存了最佳的宏块模式,min_cost记录其对应的代价  

	  //+++亚宏块级运动估计  
      if (valid[P8x8])
      {
        cost8x8 = 0;
        
        //===== store coding state of macroblock =====
        store_coding_state (cs_mb);
        //+++遍历16x16宏块的4个8x8块
        //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
        for (cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)
        {
          //--- set coordinates ---
          j0 = ((block/2)<<3);    j1 = (j0>>2);
          i0 = ((block%2)<<3);    i1 = (i0>>2);
		  //+++对宏块中的每个8x8块的模式循环  
		  //+++这个地方解释一下:min_cost8x8用于记录非RDO模式下的最小代价,而min_rdcost则是RDO模式下的最小代价  
          //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
          for (min_cost8x8=(1<<20), min_rdcost=1e30, index=(bframe?0:1); index<5; index++)
          {
            if (valid[mode=b8_mode_table[index]])
            {
              curr_cbp_blk = 0;//+++初始化当前cbp_blk
              
              if (mode==0)
              {
                //--- Direct Mode ---
                if (!input->rdopt)
                {
                  cost_direct += (cost = Get_Direct_Cost8x8 ( block, lambda_mode ));
                  if (cost==1<<30)
                    cost_direct = (1<<30);
                  have_direct ++;//+++记录是否做过direct模式的标记
                }
                best_fw_ref = direct_ref_idx[LIST_0][img->block_x+(block&1)*2][img->block_y+(block&2)];
                best_bw_ref = direct_ref_idx[LIST_1][img->block_x+(block&1)*2][img->block_y+(block&2)];
                best_pdir   = direct_pdir[img->block_x+(block&1)*2][img->block_y+(block&2)];
              } // if (mode==0)
              else//+++非direct模式下
              {
                //--- motion estimation for all reference frames ---
                PartitionMotionSearch (mode, block, lambda_motion);
                //+++这一部分与宏块级的类似 
                //--- get cost and reference frame for forward prediction ---
                for (fw_mcost=max_mcost, ref=0; ref<listXsize[LIST_0+list_offset]; ref++)
                {
                  if (!checkref || ref==0 || CheckReliabilityOfRef (block, LIST_0, ref, mode))
                  {
					  //+++(1)参考帧代价
                    mcost  = (input->rdopt ? REF_COST(lambda_motion_factor,ref,LIST_0+list_offset) : (int)(2*lambda_motion*min(ref,1)));
                     //+++(2)MV代价 
                    mcost += motion_cost[mode][LIST_0][ref][block];
                    if (mcost < fw_mcost)
                    {
                      fw_mcost    = mcost;
                      best_fw_ref = ref;
                    }
                  }
                }
                
                //store forward reference index for every block
                for (j=0; j<2; j++)
                  for (i=0; i<2; i++)
                  {
                    enc_picture->ref_idx[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = best_fw_ref;
                    enc_picture->ref_pic_id [LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = enc_picture->ref_pic_num[LIST_0 + list_offset][enc_picture->ref_idx[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j]];
                  }
                //+++B帧有所不同,还要考虑后向和双向模式,其实与上面的代码是相似的
                if (bframe)
                {
                  for (bw_mcost=max_mcost, ref=0; ref<listXsize[LIST_1+list_offset]; ref++)
                  {
                    mcost  = (input->rdopt ? REF_COST(lambda_motion_factor,ref,LIST_1+list_offset) : (int)(2*lambda_motion*min(ref,1)));
                    
                    mcost += motion_cost[mode][LIST_1][ref][block];
                    if (mcost < bw_mcost)
                    {
                      bw_mcost    = mcost;
                      best_bw_ref = ref;
                    }
                  }
                  
                  // bidirectional uses best forward and zero backward reference
                  bid_mcost  = (input->rdopt ? (REF_COST (lambda_motion_factor, best_fw_ref, LIST_0 + list_offset)+REF_COST (lambda_motion_factor, 0, LIST_1 + list_offset)) : (int)(2*lambda_motion*min(best_fw_ref,1)));
                  bid_mcost += BIDPartitionCost (mode, block, best_fw_ref, 0, lambda_motion_factor );
                  
                  //--- get prediction direction ----
                  if      (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
                  {
                    best_pdir = 0;
                    cost = fw_mcost;
                    best_bw_ref = -1;
                  }
                  else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
                  {
                    best_pdir = 1;
                    cost = bw_mcost;
                    best_fw_ref = -1;
                  }
                  else
                  {
                    best_pdir = 2;
                    cost = bid_mcost;
                    best_bw_ref = 0;
                  }
                    //store backward reference index for every block
                  for (j=0; j<2; j++)
                    for (i=0; i<2; i++)
                    {
                      enc_picture->ref_idx[LIST_0][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = best_fw_ref;
                      enc_picture->ref_idx[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = best_bw_ref;
                      //enc_picture->ref_pic_id [LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j] = enc_picture->ref_pic_num[LIST_1 + list_offset][enc_picture->ref_idx[LIST_1][img->block_x+(block&1)*2+i][img->block_y+(block&2)+j]];  
                    }
                } // if (bframe)
                else
                {
                  best_pdir = 0;
                  cost      = fw_mcost;
                }
              } // if (mode!=0)
              
              //--- store coding state before coding with current mode ---
              store_coding_state (cs_cm);

              //+++计算RDO模式代价rdcost或非RDO模式代价cost
              if (input->rdopt)
			  { //+++RDCost_for_8x8blocks是在4种亚宏块模式中选择最佳的模式(代价)  
				  //+++准确来说RDCost_for_8x8blocks主要是得到编号为block的8x8块,在mode模式下的代价
                //--- get and check rate-distortion cost ---
                rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,
                                               block, mode, best_pdir, best_fw_ref, best_bw_ref);
			                                    	//+++block:8x8块在16x16大小的宏块中的位置  
				                                    //+++mode:亚宏块模式  
              }
              else
              {
                cost += (REF_COST (lambda_motion_factor, B8Mode2Value (mode, best_pdir), list_offset + (best_pdir<1?0:1)) - 1);
              }
               //+++根据cost(非RDO模式)或rdcost(RDO模式)选择最佳模式 
              //--- set variables if best mode has changed ---
              if (( input->rdopt && rdcost < min_rdcost) ||
                  (!input->rdopt && cost   < min_cost8x8  )   )
              {
                min_cost8x8                  = cost;//+++非RDO代价
                min_rdcost                   = rdcost;//+++RDO代价
                best8x8mode          [block] = mode;//最佳的亚宏块模式
                best8x8pdir    [P8x8][block] = best_pdir;
                best8x8fwref   [P8x8][block] = best_fw_ref;
                best8x8bwref   [P8x8][block] = best_bw_ref;

                
                //--- store number of nonzero coefficients ---
                best_cnt_nonz  = cnt_nonz;
                //+++存储cbp、残差AC系数和重建帧中的预测值  
                if (input->rdopt)
                {
                  //--- store block cbp ---
                  cbp_blk8x8    &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
                  cbp_blk8x8    |= curr_cbp_blk;
                  
                  //--- store coefficients ---
                  for (k=0; k< 4; k++)
                    for (j=0; j< 2; j++)
                      for (i=0; i<65; i++)  cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT
                      
                      //--- store reconstruction and prediction ---
                  for (j=j0; j<j0+8; j++)
                    for (i=i0; i<i0+8; i++)
                    {
                      rec_mbY8x8[j][i] = enc_picture->imgY[img->pix_y+j][img->pix_x+i];
                      mpr8x8    [j][i] = img->mpr[i][j];
                    }
                }
                
                //--- store coding state ---
                store_coding_state (cs_b8);
              } // if (rdcost <= min_rdcost)
              
              //--- re-set coding state as it was before coding with current mode was performed ---
              reset_coding_state (cs_cm);
            } // if (valid[mode=b8_mode_table[index]])
		  } // for (min_rdcost=1e30, index=(bframe?0:1); index<6; index++)
		  //+++5种8x8模式循环结束  
		  //+++cost8x8:非RDO模式下的代价,其实cost8x8是一个宏块中的4个8x8块的代价和,  
		  //+++换句话说,其实cost8x8就是一个宏块(16x16)在亚宏块模式下的总的代价和.  
		  //+++这个量主要用来和上面计算的最佳宏块级模式(非RDO情况下)对应的最小代价min_cost来进行比较,  
		  //+++看宏块级模式或亚宏块级模式哪个更好,从下面的代码中我们也可以看出来. 
          
          cost8x8 += min_cost8x8;

//          if (!input->rdopt) cost8x8+= min_cost8x8;
//          else cost8x8 += min_rdcost;          
          
          if (!input->rdopt)
          {
            mode = best8x8mode[block];
            pdir = best8x8pdir[P8x8][block];
            //+++编码残差系数
            curr_cbp_blk  = 0;
            //+++cnt_nonz表示非零系数的编码代价；
            best_cnt_nonz = LumaResidualCoding8x8 (&dummy, &curr_cbp_blk, block, pdir,
                                                   (pdir==0||pdir==2?mode:0),
                                                   (pdir==1||pdir==2?mode:0),
                                                   (best8x8fwref[P8x8][block]),
                                                   (best8x8bwref[P8x8][block]));
            cbp_blk8x8   &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
            cbp_blk8x8   |= curr_cbp_blk;
            
            //--- store coefficients ---
            for (k=0; k< 4; k++)
              for (j=0; j< 2; j++)
                for (i=0; i<65; i++)  cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT
                
                //--- store reconstruction and prediction ---
            for (j=j0; j<j0+8; j++)
              for (i=i0; i<i0+8; i++)
              {
                rec_mbY8x8[j][i] = enc_picture->imgY[img->pix_y+j][img->pix_x+i];
                mpr8x8    [j][i] = img->mpr[i][j];
              }
          }
          
          //----- set cbp and count of nonzero coefficients ---
          if (best_cnt_nonz)
          {
            cbp8x8        |= (1<<block);
            cnt_nonz_8x8  += best_cnt_nonz;
          }
          
          mode=best8x8mode[block];//第2396行会用到mode,best8x8mode[block]应该是第block号8x8块所选取的最佳亚宏块级模式 
          //===== reset intra prediction modes (needed for prediction, must be stored after 8x8 mode dec.) =====
          j0 = img->block_y+2*(block/2);
          i0 = img->block_x+2*(block%2);
          for (j=j0; j<j0+2; j++)
            for (i=i0; i<i0+2; i++) 
              ipredmodes[i][j]         = DC_PRED;
          i0 = 4*block;
          for (i=i0; i<i0+4; i++)    currMB->intra_pred_modes[i]  = DC_PRED;
            
          if (block<3)
          {
            //===== re-set reconstructed block =====
            j0   = 8*(block/2);
            i0   = 8*(block%2);
            for (j=j0; j<j0+8; j++)
              for (i=i0; i<i0+8; i++)  
                enc_picture->imgY[img->pix_y+j][img->pix_x+i] = rec_mbY8x8[j][i];
          } // if (block<3)
           
		  //+++这儿与宏块级类似,要设置MV和参考帧 
          //===== set motion vectors and reference frames (prediction) =====
          SetRefAndMotionVectors (block, mode, best8x8pdir[P8x8][block], best8x8fwref[P8x8][block], best8x8bwref[P8x8][block]);
          //===== set the coding state after current block =====
          reset_coding_state (cs_b8);
        } // for (cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)
		//+++此刻完成了对一个宏块(16x16)在亚宏块模式下的计算,或者说对宏块中每个8x8块完成亚宏块模式的选择  
		//+++并且在亚宏块模式下的代价(4个8x8块代价之和)也计算完毕,在cost8x8中,下面会用到,和宏块级模式(min_cost)进行比较  
        //===== store intra prediction modes for 8x8+ macroblock mode =====
        for (k=0, j=img->block_y; j<img->block_y+4; j++)
        {
          for ( i=img->block_x; i<img->block_x+4; i++, k++)
          {
            b8_ipredmode       [k] = ipredmodes    [i][j];
            b8_intra_pred_modes[k] = currMB->intra_pred_modes[k];
          }
        }
        
        //--- re-set coding state (as it was before 8x8 block coding) ---
        reset_coding_state (cs_mb);
        for (i=0; i<16; i++)
          for(j=0; j<16; j++)
			//++ 原代码中这里 img->mpr 的下标写反了
//          diffy[j][i] = imgY_org[img->opix_y+j][img->opix_x+i]-img->mpr[j][i];
			diffy[j][i] = imgY_org[img->opix_y+j][img->opix_x+i]-img->mpr[i][j];

		//+++到这儿应该是一个结点了啊  
		//+++(在非RDO模式下)对一个宏块(当前宏块),此时宏块级(16x16,16x8,8x16)的最小开销已经计算出来,即min_cost  
		//+++子宏块级(8x8,8x4,4x8和4x4)的开销也计算出来(4个8x8块代价之和),即cost8x8  
		//+++这样可以比较这两种方式,选择开销更小的那种宏块分割方式  



        if(cost8x8 < min_cost)//+++// 帧间8x8模式与帧间16x16代价相比 
        {
           best_mode = P8x8;
           min_cost = cost8x8;
		 } 
		   //+++到这儿我们也可以看到,在非RDO模式下,best_mode保存的是一个宏块(当前宏块)帧间预测的最优模式,是(16x16,16x8,8x16和P8x8)之一.  
		   //+++(16x16,16x8,8x16)这三种宏块级的模式是选择最优模式后再与P8x8模式比较.min_cost现在保存的是一个宏块(当前宏块)帧间预测的最小代价  
		   //+++从下面帧内预测的代码我们也可以看到,min_cost作为帧间最小代价是要和帧内两种模式I4MB和I16MB进行比较的.  
        
      }
      else // if (valid[P8x8])
      {
        //+++这个地方是不用子宏块划分的方法时,将其开销设为很大
        cost8x8 = (1<<20);
      }
   
      // Find a motion vector for the Skip mode
      if((img->type == P_SLICE)||(img->type == SP_SLICE))
        FindSkipModeMotionVector ();//+++用于计算MV预测值 
    }
    else //帧内//yes
	{
       min_cost = (1<<20);
	}


	//===================================  
	//+++下面是对帧内预测进行求其开销  
    if (input->rdopt)//yes
    {//+++使用RDO准则的情况 
      int mb_available_up;
      int mb_available_left;
      int mb_available_up_left;
      
      min_rdcost = max_rdcost;//1.0000000000000000e+030
      
      // precompute all new chroma intra prediction modes
	  //++ 对色度进行帧内预测,并求出最优模式,  
	  //+++RDO模式下,需要色度预测,保存在currMB->c_ipred_mode中
      IntraChromaPrediction8x8(&mb_available_up, &mb_available_left, &mb_available_up_left);
	  //+++阅读完IntraChromaPrediction8x8函数的代码,并且结合下面的代码,我们可以看到:  
	  //+++这个函数不做任何运动估计和补偿，而仅仅是求出要进行帧内编码的色度宏块需要的参考块  
	  //+++本来在currMB->c_ipred_mode保存了最佳的模式,但是下面却被覆盖掉了 

	  //++ 分别在四种色度模式下进行 RDO 计算，如果是 inter 模式，因为色度预测模式与 SSD 计算
	  //++ 无关，因此只需要计算一次（利用 currMB->c_ipred_mode == DC_PRED_8 条件限制来实现）
      for (currMB->c_ipred_mode=DC_PRED_8; currMB->c_ipred_mode<=PLANE_8; currMB->c_ipred_mode++)
      {
        
        // bypass if c_ipred_mode is not allowed
          //+++若当前 c_ipred_mode 模式不被允许，则考虑下一个 c_ipred_mode 模式; 
        if ((currMB->c_ipred_mode==VERT_PRED_8 && !mb_available_up) ||
          (currMB->c_ipred_mode==HOR_PRED_8 && !mb_available_left) ||
          (currMB->c_ipred_mode==PLANE_8 && (!mb_available_left || !mb_available_up || !mb_available_up_left)))
          continue;
        
        
        //===== GET BEST MACROBLOCK MODE =====
        for (ctr16x16=0, index=0; index<7; index++)
        {
          mode = mb_mode_table[index];
          
          //--- for INTER16x16 check all prediction directions ---
          if (mode==1 && img->type==B_SLICE)//no
          {
            best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = ctr16x16;
            if (ctr16x16 < 2) index--;
            ctr16x16++;
          }
          
          img->NoResidueDirect = 0;
          
          if (valid[mode])
          {
            // bypass if c_ipred_mode not used
			//++ 设置当前宏块类型以及其中每个8*8块的分割方式和预测方向，每个4*4块的参考帧索引
			//++ 该函数在下面的 RDCost_for_macroblocks 函数内再次调用，进行了重复操作
            SetModesAndRefframeForBlocks (mode);
            if (currMB->c_ipred_mode == DC_PRED_8 ||	//++ 利用这个条件限制来实现 inter 模式时只计算一次 RDO
              (IS_INTRA(currMB) ))
            {
				//+++此刻我的理解是:RDCost_for_macroblocks是在模式0-3,P8x8,Intra4x4和Intra16x16之间选择一个最佳模式  
				//+++这是利用RDO的方法.上面再帧间预测时其实已经将在1-3中选择出了最佳模式,但那是在非RDO模式下的  
				//+++关于这个地方的RDCost_for_macroblocks我要说明一下:该函数是求出对应mode模式下的RDO代价,然后将这个代价与min_rdcost进行比较  
				//+++如果此时RDO的代价比较小,返回1,执行if,保存相关参数,并且修改min_rdcost,否则返回0,这样下面的if语句不执行,接着再比较下一个模式.  
				//+++所以我们可以看到(0-3,P8x8,Intra4x4和Intra16x16)这些模式的比较是随着外面的这个for循环,依次在RDCost_fo...函数中完成的,  
              if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))//++ 帧内模式时亮度存在重复计算情况：因为色度预测模式与亮度预测模式
                                   //++ 无关，所以在该色度模式循环中每次计算得到的 intra16*16 和 intra4*4 宏块类型的最佳亮度模式都是完全相同的
			  {
                //Rate control
                if(mode == P8x8)//0
                {
                  for (i=0; i<16; i++)
                    for(j=0; j<16; j++)
                      diffy[j][i] = imgY_org[img->opix_y+j][img->opix_x+i] - mpr8x8[j][i];
                }else
                {
                  for (i=0; i<16; i++)
                    for(j=0; j<16; j++)
                      diffy[j][i] = imgY_org[img->opix_y+j][img->opix_x+i] - pred[j][i];
                }
				//+++存储该模式下的一些参数:因为RDCost_for_macroblocks的计算结果是mode模式的预测开销小于min_rdcost,  
				//+++所以mode模式其实更好,所以现在要保存该模式的一些参数.
                store_macroblock_parameters (mode);
              }
            }
          }
          if (valid[0] && bframe && mode == 0 && currMB->cbp && (currMB->cbp&15) != 15) //g050
          {
            img->NoResidueDirect = 1;
            if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
            {
              //Rate control
              for (i=0; i<16; i++)
                for(j=0; j<16; j++)
                  diffy[j][i] = imgY_org[img->opix_y+j][img->opix_x+i] - pred[j][i];
                store_macroblock_parameters (mode);
            }
          }
        }
      }
    }
    else//+++不使用RDO准则的情况//no
    {
		//+++从下面的代码可以看到,在不同模式下计算出的cost都要和min_cost(帧间预测代价)进行比较  
		//+++然后决定是否更新min_cost和best_mode,所以我们也可以看到,最终的代价也是保存在min_cost中,  
		//+++best_mode对应保存最佳模式,从而用于最后的熵编码  

		//+++下面对Direct,I4MB和I16MB三种模式进行依次比较,代码是类似的.  
      if (valid[0] && bframe) // check DIRECT MODE
      {
        cost  = (have_direct?cost_direct:Get_Direct_CostMB (lambda_mode));
        cost -= (int)floor(16*lambda_motion+0.4999);
        if (cost <= min_cost)
        {
          //Rate control
          for (i=0; i<16; i++)
            for(j=0; j<16; j++)
              diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i]-img->mpr[i][j];
             //+++更新min_cost和best_mode
            min_cost  = cost;
            best_mode = 0;
        }
      }
      if (valid[I4MB]) // check INTRA4x4
      {
		  //+++从Intra4x4的9中模式中选择最佳的模式,函数Mode_Dec..中的cost就是当前宏块的代价  
		  //+++追踪各个函数,我们可以看到函数调用顺序如下:在非RDO模式下调用  
		  //+++Mode_Decision_for_Intra4x4Macroblock(得到一个宏块的代价)  
		  //+++->Mode_Decision_for_8x8IntraBlocks(得到一个8x8块代价)  
		  //+++->Mode_Decision_for_Intra4x4Macroblock(得到一个4x4块代价),  
		  //+++其实最后一个函数才是真正的Intra4x4模式选择的核心,在这个函数中包含9中模式的比较(RDO和非RDO) 
        currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda_mode, &cost);//++ 计算当前宏块采用intra4*4编码时的代价和CBP
        if (cost <= min_cost)
        {
          //Rate control
          for (i=0; i<16; i++)
            for(j=0; j<16; j++)
              diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i]-img->mpr[i][j];//++ 保存采用intra4*4编码时最佳模式的残差块
            //+++更新min_cost和best_mode            
            min_cost  = cost;
            best_mode = I4MB;
        }
      }
      if (valid[I16MB]) // check INTRA16x16
      {
		  //+++比较与RDO模式下代码可以发现,其实在RDO模式下进行Intra16x16预测时,调用的函数Intra16x16_Mode_Decision的代码  
		  //+++与这儿类似:也是先intrapred_luma_16x16,再find_sad_16x16,最后dct_luma_16x16  
		  //+++同时这也说明:intra16*16 在进行4种模式选择时候，无论是否在RDO情况下，其选择过程是相同的。

        intrapred_luma_16x16 ();	//++ 分别计算当前宏块在4种intra16*16帧内预测模式下的预测块
                                    //+++intrapred_luma_16x16执行完后,预测的结果保存在数组img->mprr_2中

        cost = find_sad_16x16 (&i16mode);	//++ 计算intra16*16类型的代价（以SATD作为判断标准，因为该段代码是不采用RDO的处理过程）
                                            //+++此时i16mode保存的是Intra16x16的最佳模式,cost是对应的代价

        if (cost < min_cost)	//++ 如果采用intra16*16类型编码的代价小于采用intra4*4类型编码的代价，则对当前宏块采用intra16*16类型编码
        {
          //Rate control
          for (i=0; i<16; i++)
            for(j=0; j<16; j++)
              diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i]-img->mprr_2[i16mode][j][i];
              //+++更新min_cost和best_mode          
            best_mode   = I16MB;
			//+++这里说明一下:find_sad_16x16函数中求sad时已经做过DCT变换，但在dct_luma_16x16函数中又要重新做一遍呢？  
			//+++答:前者是对每种预测模式都要做，后者是对最佳模式才做。每种模式做了之后并没有把结果都保存起来，所以要对最佳模式重新再做一次。 
            currMB->cbp = dct_luma_16x16 (i16mode);//++ 对当前宏块采用intra16*16类型进行编码，并计算CBP
        }
      }
    }
	//+++通过比较RDO模式与非RDO模式我们可以发现,  
	//+++在intra4*4与intra16*16之间做选择时，RDO情况下色度要参与计算，而非RDO情况下色度不参与计算  
	//+++并且我们可以看到在非RDO模式下,具体的帧内预测还是帧间预测的选择在这儿完成(可以在这儿看到具体的代码)  
	//+++但是对于RDO模式,具体的帧内还是帧间的选择却是在函数RDCost_for_macroblocks中完成的.  
	//+++帧内预测结束  
	//===================================  
    if (rerun==0)//yes
    {
      intra1 = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    }
  } // for (rerun=0; rerun<runs; rerun++)
  //+++到现在,帧内和帧间模式的选择终于结束  

  //+++下面的要进行存储编码参数,分RDO和非RDO  
  //+++印证了[h264率失真RDO研究]一文中的一句:  
  //+++"如果用RDO代价函数,因为计算代价时已经编码,择优后仅需对编码进行存储;  
  //+++如果用非RDO代价函数,则择优后还要进行编码过程,编码后存储参数."  
  //+++从下面的代码我们也可清楚的看到  
  if (input->rdopt)//yes
  {
    
    if ((cbp!=0 || best_mode==I16MB ))
      currMB->prev_cbp = 1;
    else if (cbp==0 && !input->RCEnable)
    {
      currMB->delta_qp = 0;
      currMB->qp = currMB->prev_qp;
      img->qp = currMB->qp;
      currMB->prev_cbp = 0;
    }
    //+++存储编码参数
    set_stored_macroblock_parameters ();
	//+++    (1)保存重建图像  
	//+++    rdopt->rec_mbY[j][i] = rec_mbY[j][i];  
	//+++    rdopt->rec_mbU[j][i] = rec_mbU[j][i];    
	//+++    rdopt->rec_mbV[j][i] = rec_mbV[j][i];   
	//+++    (2)变换系数和cbp  
	//+++    img->cofAC=i4p(cofAC);  
	//+++    img->cofDC=i3p(cofDC);  
	//+++    currMB->cbp      = cbp;  
	//+++    currMB->cbp_blk = cbp_blk;  
	//+++    (3)宏块模式  
	//+++    currMB->mb_type = mode(best_mode);  
	//+++    (4)参考帧保存  
	//+++    rdopt.c 1511行,对前向参考,后向参考帧分别进行了保存  
	//+++    (5)帧内预测模式  
	//+++    currMB->c_ipred_mode = best_c_imode;  
	//+++    (6)为当前块保存运动向量  
	//+++    SetMotionVectorsMB (currMB, bframe);  
	//+++    从上面其实可以看到"="右边的量在rdopt.c中都是全局变量,在完成i_encode_one_macroblock之前,对这些参数进行一下保存.
  }
  else//no//不使用RDO
  {
    //===== set parameters for chosen mode =====
    SetModesAndRefframeForBlocks (best_mode);//++ 设置当前宏块的参数，包括:宏块类型(mb_type)、
                                            //++ 4个8*8块的分割模式和预测方向(b8mode，b8pdir)、16个4*4块的参考帧索引（ref_idx，ref_pic_id）
	if (best_mode==P8x8)
    {
      SetCoeffAndReconstruction8x8 (currMB);
    }
    else
    {
      if (best_mode!=I4MB)
      {//+++    一进来先建立一个16*16的2维数组，不管当前编码的是16*8或者8*16，然后把当前要编码
	   //+++     的block复制到这个数组中，始终以数组的top/left为起点存放数据。
        for (k=0, j=img->block_y; j<img->block_y+4; j++)
          for (     i=img->block_x; i<img->block_x+4; i++, k++)
          {
            ipredmodes    [i][j] = DC_PRED;
            currMB->intra_pred_modes[k] = DC_PRED;
          }
          if (best_mode!=I16MB)
          {//+++    在RDO=0的情况下，前面只计算了预测值，没有计算DCT，所以这里做DCT
            LumaResidualCoding ();
            //Rate control
            for (i=0; i<16; i++)
              for(j=0; j<16; j++)
                diffy[j][i] = imgY_org[img->pix_y+j][img->pix_x+i]-img->mpr[i][j];
          }
      }
    }
    // precompute all chroma intra prediction modes
    IntraChromaPrediction8x8(NULL, NULL, NULL);//++ 对色度块进行处理，包括：分别计算4种色度帧内预测模式下的预测块、
                                               //++ 代价（该分支未采用RDO，因此以SATD作为判断标准）、最佳预测模式
    img->i16offset = 0;
    dummy = 0;
    ChromaResidualCoding (&dummy);//+++    如果是I帧，把上面得到的预测值复制到4个4*4的2维数组里，求残差，再做DCT；
	                              //+++    如果是P帧，则用由亮度为依据得到的MV得到色差的预测值，再做残差后做DCT。
    if (best_mode==I16MB)
    {
      img->i16offset = I16Offset  (currMB->cbp, i16mode);
    }
    SetMotionVectorsMB (currMB, bframe);
    
    //===== check for SKIP mode =====
    if ((img->type==P_SLICE || img->type==SP_SLICE) && best_mode==1 && currMB->cbp==0 &&
      enc_picture->ref_idx[LIST_0][img->block_x][img->block_y]==0 &&
      enc_picture->mv[LIST_0][img->block_x][img->block_y][0]==allmvs[0][0][0][0][0][0] &&
      enc_picture->mv[LIST_0][img->block_x][img->block_y][1]==allmvs[0][0][0][0][0][1]               )
    {
      currMB->mb_type=currMB->b8mode[0]=currMB->b8mode[1]=currMB->b8mode[2]=currMB->b8mode[3]=0;
    }
    
    if(img->MbaffFrameFlag)
      set_mbaff_parameters();
  }
  
//   // Rate control
  if(input->RCEnable)//no
  {   
    if(img->type==P_SLICE)
    {
      img->MADofMB[img->current_mb_nr] = calc_MAD();
      
      if(input->basicunit<img->Frame_Total_Number_MB)
      {
        img->TotalMADBasicUnit +=img->MADofMB[img->current_mb_nr];
        
        /* delta_qp is present only for non-skipped macroblocks*/
        if ((cbp!=0 || best_mode==I16MB))
          currMB->prev_cbp = 1;
        else
        {
          img->qp -= currMB->delta_qp;
          currMB->delta_qp = 0;
          currMB->qp = img->qp;
          currMB->prev_cbp = 0;
        }
        /* When MBAFF is used, delta_qp is only present for the first non-skipped macroblock of each 
        macroblock pair*/
        if (input->MbInterlace)
        {
          if(!currMB->mb_field)
          {
            DELTA_QP = currMB->delta_qp;
            QP      = currMB->qp;
          }
          else
          {
            DELTA_QP2 = currMB->delta_qp;
            QP2      = currMB->qp;
          }
        }       
      }
    }
  }
  
    //+++保存最佳代价在rdopt->min_rdcost中   
  if(input->rdopt)//yes
    rdopt->min_rdcost = min_rdcost;
  else
    rdopt->min_rdcost = min_cost; 

  if(img->MbaffFrameFlag)//no
  {
    if (img->current_mb_nr%2) //bottom
    {
      if ((currMB->mb_type ? 0:((img->type == B_SLICE) ? !currMB->cbp:1))  // bottom is skip
        &&(prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1))) // top is skip
      {
        if (!(field_flag_inference() == curr_mb_field))
        {
          rdopt->min_rdcost = 1e30;  // don't allow coding of an MB pair as skip if wrong inference
        }
      }
    }
  }
  
  //+++RestrictRefFrames应该是和容错有关系的。我的理解是，如果某个区域在后面帧里被intra方式update（从而可以避免错误传播),  
  //+++这个区域很可能是出现错误的区域，因此，不使用这个区域作参考以避免错误传播。请大家指正 
  //===== Decide if this MB will restrict the reference frames =====
  if (input->RestrictRef==1)//no
  {
    if (input->rdopt==1)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra ? 1 : 0);
    }
    else if (input->rdopt==2)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
    }
  }
  else if (input->RestrictRef==2)//no
  {
    refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
  }
  
  if(input->FMEnable)//no
    skip_intrabk_SAD(best_mode, listXsize[LIST_0+list_offset]);
}


void set_mbaff_parameters()
{
  int  i, j, k, l;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode;
  int         bframe   = (img->type==B_SLICE);
  int     **ipredmodes = img->ipredmode;

  if (!img->MbaffFrameFlag)
    return;

  //===== reconstruction values =====
  for (j=0; j<16; j++)
    for (i=0; i<16; i++)
      rdopt->rec_mbY[j][i]           = enc_picture->imgY[img->pix_y+j][img->pix_x+i]; 

  for (j=0; j<8; j++)
    for (i=0; i<8; i++)
    {
      rdopt->rec_mbU[j][i]           = enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x+i];   
      rdopt->rec_mbV[j][i]           = enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x+i];    
    }

  //===== coefficients and cbp =====
  rdopt->mode = mode;
  rdopt->i16offset = img->i16offset;  // For MBINTLC  -Rajeev
  rdopt->cbp = currMB->cbp;
  rdopt->cbp_blk = currMB->cbp_blk;
  rdopt->mb_type  =currMB->mb_type;

  if(rdopt->mb_type == 0 && mode != 0)
  {
    mode=0;
    rdopt->mode=0;
  }

  for(i=0;i<6;i++)
    for(j=0;j<4;j++)
      for(k=0;k<2;k++)
        for(l=0;l<18;l++)
          rdopt->cofAC[i][j][k][l] = img->cofAC[i][j][k][l];

  for(i=0;i<3;i++)
    for(k=0;k<2;k++)
      for(l=0;l<18;l++)
        rdopt->cofDC[i][k][l] = img->cofDC[i][k][l];


  for (i=0; i<4; i++)
  {
    rdopt->b8mode[i]  = currMB->b8mode[i];                  
    rdopt->b8pdir[i]  = currMB->b8pdir[i];                  
  }

  //==== reference frames =====
  for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      rdopt->refar[LIST_0][j][i]       = enc_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j];
    }

  if (bframe)
  {
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        rdopt->refar[LIST_1][j][i]     = enc_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j];
      }
  }


  for   (k=0, j=img->block_y; j<img->block_y+4; j++)
    for (     i=img->block_x; i<img->block_x+4; i++, k++)
    {
      rdopt->ipredmode[i][j]     = ipredmodes[i][j];
      rdopt->intra_pred_modes[k] = currMB->intra_pred_modes[k];
    }
}
