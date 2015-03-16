
/*!
 *************************************************************************************
 * \file mb_access.c
 *
 * \brief
 *    Functions for macroblock neighborhoods
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten String          <suehring@hhi.de>
 *************************************************************************************
 */

#include "global.h"

/////////////////////////////////////>


// extern StorablePicture *dec_picture;


/////////////////////////////////////<

/*!
 ************************************************************************
 * \brief
 *    returns 1 if the macroblock at the given address is available
 ************************************************************************
 */
int mb_is_available(int mbAddr, int currMbAddr)
{
  if ((mbAddr < 0) || (mbAddr > ((int)img->PicSizeInMbs - 1)))
    return 0;

  // the following line checks both: slice number and if the mb has been decoded
  if (!img->DeblockCall)
  {
    if (img->mb_data[mbAddr].slice_nr != img->mb_data[currMbAddr].slice_nr)
      return 0;
  }
  
  return 1;
}


/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 ************************************************************************
 */
// D   |B    |C
// -------------
// A   |cur  | 　　　检查临近宏块的可用性

void CheckAvailabilityOfNeighbors()
{
	const int mb_nr = img->current_mb_nr;//宏块编号
  Macroblock *currMB = &img->mb_data[mb_nr];//宏块指针指向当前宏块

  // mark all neighbors as unavailable
  currMB->mb_available_up   = NULL;//当前宏块上边宏块指针
  currMB->mb_available_left = NULL;//当前宏块左边宏块指针

  if (img->MbaffFrameFlag)//桢场自适应模式
  {
// 	  假设qcif 176*144那么img->PicWidthInMbs＝11
// 		 　0   2   4.... 　偶数代表顶场或者顶桢宏块
// 		   1   3   5....　　奇数代表底场或者底桢宏块
// 		   22  24  26...
// 		   23  25  27...
// 		   比如宏块4左边A就是宏块2,宏块2的A宏块就是宏块０注意下宏块３左边宏块A是宏块０，可能由于把32*16宏块对的缘故

    currMB->mbAddrA = 2 * (mb_nr/2 - 1);
// 	宏块24的上宏块B是宏块是宏块2,宏块26的上宏块B是宏块４
    currMB->mbAddrB = 2 * (mb_nr/2 - img->PicWidthInMbs);
    currMB->mbAddrC = 2 * (mb_nr/2 - img->PicWidthInMbs + 1);
    currMB->mbAddrD = 2 * (mb_nr/2 - img->PicWidthInMbs - 1);
    
	currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && (((mb_nr/2) % img->PicWidthInMbs)!=0);//保证当前宏块不是最左边的宏块
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
	currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr/2 +1) % img->PicWidthInMbs)!=0);//保证当前宏块不是最右边的宏块
	currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && (((mb_nr/2) % img->PicWidthInMbs)!=0);//保证当前宏块不是最左边的宏块
  }
  else
  {
// 	      D   |B    |C
// 		  -------------
// 		  A   |cur  | 
	  //当前宏块的左边宏块A
    currMB->mbAddrA = mb_nr - 1;
	//当前宏块的上边宏块B
	currMB->mbAddrB = mb_nr - img->PicWidthInMbs;//其中img->PicWidthInMbs指宏块为单位，1行中宏块数目
    currMB->mbAddrC = mb_nr - img->PicWidthInMbs + 1;
    currMB->mbAddrD = mb_nr - img->PicWidthInMbs - 1;
	//以上确定宏块A B C D位置下面代码确定A B C D的可用性
	currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);//mb_nr % img->PicWidthInMbs保证当前宏块不是最左边的宏块
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
	currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+1) % img->PicWidthInMbs)!=0);//保证当前宏块不是最右边的宏块
	currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);//保证当前宏块不是最左边的宏块
  }

  if (currMB->mbAvailA) currMB->mb_available_left = &(img->mb_data[currMB->mbAddrA]);//宏块左边宏块指针指向宏块A
  if (currMB->mbAvailB) currMB->mb_available_up   = &(img->mb_data[currMB->mbAddrB]);//宏块上边宏块指针指向宏块B
}


/*!
 ************************************************************************
 * \brief
 *    returns the x and y macroblock coordinates for a given MbAddress
 ************************************************************************
 */
void get_mb_block_pos (int mb_addr, int *x, int*y)
{

  if (img->MbaffFrameFlag)//no
  {
    *x = ((mb_addr/2)%img->PicWidthInMbs);
    *y = ( ((mb_addr/2)/img->PicWidthInMbs)*2+(mb_addr%2));
  }
  else
  {
    *x = (mb_addr%img->PicWidthInMbs);//宏块的光栅扫描坐标转换成二维x坐标（以宏块为单位）
    *y = (mb_addr/img->PicWidthInMbs);//宏块的光栅扫描坐标转换成二维y坐标（以宏块为单位）
  }
}


/*!
 ************************************************************************
 * \brief
 *    returns the x and y sample coordinates for a given MbAddress
 ************************************************************************
 */
void get_mb_pos (int mb_addr, int *x, int*y)
{
  get_mb_block_pos(mb_addr, x, y);
  
  (*x) *= MB_BLOCK_SIZE;
  (*y) *= MB_BLOCK_SIZE;
}


/*!
 ************************************************************************
 * \brief
 *    get neighbouring positions for non-aff coding求出邻接位置
 * \param curr_mb_nr//当前宏块编号
 *   current macroblock number (decoding order)
 * \param xN//相对当前宏块左上坐标
 *    input x position
 * \param yN//相对当前宏块左上坐标
 *    input y position
 * \param luma//亮度１　色度０
 *    1 if luma coding, 0 for chroma
 * \param pix//输出领接点位置信息
 *    returns position informations
 ************************************************************************
 *///在无自适应桢场情况下　求出邻接点在那个宏块和在图像中位置
void getNonAffNeighbour(unsigned int curr_mb_nr, int xN, int yN, int luma, PixelPos *pix)
{
	Macroblock *currMb = &img->mb_data[curr_mb_nr];//当前宏块指针
  int maxWH;

  if (luma)//亮度块最大宽度16
    maxWH = 16;
  else
    maxWH = 8;
  //xN yN领接点是和当前宏块左上角的相对位置
  if ((xN<0)&&(yN<0))//领接D宏块
  {
    pix->mb_addr   = currMb->mbAddrD;//23
    pix->available = currMb->mbAvailD;//0
  }
  else
	  if ((xN<0)&&((yN>=0)&&(yN<maxWH)))//领接点在邻接A宏块

  {
    pix->mb_addr  = currMb->mbAddrA;
    pix->available = currMb->mbAvailA;
  }
  else
	  if (((xN>=0)&&(xN<maxWH))&&(yN<0))//领接点在领接B宏块
  {
    pix->mb_addr  = currMb->mbAddrB;
    pix->available = currMb->mbAvailB;
  }
  else
	  if (((xN>=0)&&(xN<maxWH))&&((yN>=0)&&(yN<maxWH)))//领接点在当前宏块
  {
    pix->mb_addr  = curr_mb_nr;
    pix->available = 1;
  }
  else
	  if ((xN>=maxWH)&&(yN<0))//领接点在邻接C宏块
  {
    pix->mb_addr  = currMb->mbAddrC;
    pix->available = currMb->mbAvailC;
  }
  else 
  {
    pix->available = 0;
  }
	  if (pix->available || img->DeblockCall)//如果邻接点存在
  {
	  pix->x = (xN + maxWH) % maxWH;//求出邻接点在宏块内的位置
    pix->y = (yN + maxWH) % maxWH;
	get_mb_pos(pix->mb_addr, &(pix->pos_x), &(pix->pos_y));//求出邻接点所在宏块左上点在图像中坐标
	if (luma)//亮度，求出领接点在图像中的位置
    {
      pix->pos_x += pix->x;
      pix->pos_y += pix->y;
    }
	else//色度　pix->pos_x/2　pix->pos_y/2就说明了色度水平垂直各下采样
    {
      pix->pos_x = (pix->pos_x/2) + pix->x;
      pix->pos_y = (pix->pos_y/2) + pix->y;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    get neighbouring positions for aff coding
 * \param curr_mb_nr
 *   current macroblock number (decoding order)
 * \param xN
 *    input x position
 * \param yN
 *    input y position
 * \param luma
 *    1 if luma coding, 0 for chroma
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void getAffNeighbour(unsigned int curr_mb_nr, int xN, int yN, int luma, PixelPos *pix)
{
  Macroblock *currMb = &img->mb_data[curr_mb_nr];
  int maxWH;
  int yM = -1;

  if (luma)
    maxWH = 16;
  else
    maxWH = 8;

  // initialize to "not available"
  pix->available = 0;

  if(yN > (maxWH - 1))
  {
    return;
  }
  if ((xN > (maxWH -1)) && ((yN >= 0)&&(yN < (maxWH ))))
  {
    return;
  }

  if (xN < 0)
  {
    if (yN < 0)
    {
      if(!currMb->mb_field)
      {
        // frame
        if (curr_mb_nr%2 == 0)
        {
          // top
          pix->mb_addr  = currMb->mbAddrD  + 1;
          pix->available = currMb->mbAvailD;
           yM      = yN;
        }
        else
        {
          // bottom
          pix->mb_addr  = currMb->mbAddrA;
          pix->available = currMb->mbAvailA;
          if (currMb->mbAvailA)
          {
            if(!img->mb_data[currMb->mbAddrA].mb_field)
            {
               yM = yN;
            }
            else
            {
              (pix->mb_addr)++;
               yM = (yN + maxWH) >> 1;
            }
          }
        }
      }
      else
      {
        // field
        if(curr_mb_nr % 2 == 0)
        {
          // top
          pix->mb_addr  = currMb->mbAddrD;
          pix->available = currMb->mbAvailD;
          if (currMb->mbAvailD)
          {
            if(!img->mb_data[currMb->mbAddrD].mb_field)
            {
              (pix->mb_addr)++;
               yM = 2 * yN;
            }
            else
            {
               yM = yN;
            }
          }
        }
        else
        {
          // bottom
          pix->mb_addr  = currMb->mbAddrD+1;
          pix->available = currMb->mbAvailD;
           yM      = yN;
        }
      }
    }
    else
    { // xN < 0 && yN >= 0
      if ((yN >= 0) && (yN <maxWH))
      {
        if (!currMb->mb_field)
        {
          // frame
          if(curr_mb_nr % 2 == 0)
          {
            // top
            pix->mb_addr  = currMb->mbAddrA;
            pix->available = currMb->mbAvailA;
            if (currMb->mbAvailA)
            {
              if(!img->mb_data[currMb->mbAddrA].mb_field)
              {
                 yM = yN;
              }
              else
              {
                if (yN %2 == 0)
                {
                   yM = yN>> 1;
                }
                else
                {
                  (pix->mb_addr)++;
                   yM = yN>> 1;
                }
              }
            }
          }
          else
          {
            // bottom
            pix->mb_addr  = currMb->mbAddrA;
            pix->available = currMb->mbAvailA;
            if (currMb->mbAvailA)
            {
              if(!img->mb_data[currMb->mbAddrA].mb_field)
              {
                (pix->mb_addr)++;
                 yM = yN;
              }
              else
              {
                if (yN %2 == 0)
                {
                   yM = (yN + maxWH) >> 1;
                }
                else
                {
                  (pix->mb_addr)++;
                   yM = (yN + maxWH) >> 1;
                }
              }
            }
          }
        }
        else
        {
          // field
          if (curr_mb_nr % 2 == 0)
          {
            // top
            pix->mb_addr  = currMb->mbAddrA;
            pix->available = currMb->mbAvailA;
            if (currMb->mbAvailA)
            {
              if(!img->mb_data[currMb->mbAddrA].mb_field)
              {
                if (yN < (maxWH / 2))
                {
                   yM = yN << 1;
                }
                else
                {
                  (pix->mb_addr)++;
                   yM = (yN << 1 ) - maxWH;
                }
              }
              else
              {
                 yM = yN;
              }
            }
          }
          else
          {
            // bottom
            pix->mb_addr  = currMb->mbAddrA;
            pix->available = currMb->mbAvailA;
            if (currMb->mbAvailA)
            {
              if(!img->mb_data[currMb->mbAddrA].mb_field)
              {
                if (yN < (maxWH / 2))
                {
                   yM = (yN << 1) + 1;
                }
                else
                {
                  (pix->mb_addr)++;
                   yM = (yN << 1 ) + 1 - maxWH;
                }
              }
              else
              {
                 (pix->mb_addr)++;
                 yM = yN;
              }
            }
          }
        }
      }
    }
  }
  else
  {
    // xN >= 0
    if ((xN >= 0)&&(xN <maxWH))
    {
      if (yN<0)
      {
        if (!currMb->mb_field)
        {
          //frame
          if (curr_mb_nr % 2 == 0)
          {
            //top
            pix->mb_addr  = currMb->mbAddrB;
            // for the deblocker if the current MB is a frame and the one above is a field
            // then the neighbor is the top MB of the pair
            if (currMb->mbAvailB)
            {
              if (!(img->DeblockCall == 1 && (img->mb_data[currMb->mbAddrB]).mb_field))
                pix->mb_addr  += 1;
            }

            pix->available = currMb->mbAvailB;
             yM      = yN;
          }
          else
          {
            // bottom
            pix->mb_addr  = curr_mb_nr - 1;
            pix->available = 1;
             yM      = yN;
          }
        }
        else
        {
          // field
          if (curr_mb_nr % 2 == 0)
          {
            // top
            pix->mb_addr  = currMb->mbAddrB;
            pix->available = currMb->mbAvailB;
            if (currMb->mbAvailB)
            {
              if(!img->mb_data[currMb->mbAddrB].mb_field)
              {
                (pix->mb_addr)++;
                 yM = 2* yN;
              }
              else
              {
                 yM = yN;
              }
            }
          }
          else
          {
            // bottom
            pix->mb_addr  = currMb->mbAddrB + 1;
            pix->available = currMb->mbAvailB;
             yM      = yN;
          }
        }
      }
      else
      {
        // yN >=0
        // for the deblocker if this is the extra edge then do this special stuff
        if (yN == 0 && img->DeblockCall == 2)
        {
          pix->mb_addr  = currMb->mbAddrB + 1;
          pix->available = 1;
           yM      = yN - 1;
        }

        else if ((yN >= 0) && (yN <maxWH))
        {
          pix->mb_addr  = curr_mb_nr;
          pix->available = 1;
           yM      = yN;
        }
      }
    }
    else
    {
      // xN >= maxWH
      if(yN < 0)
      {
        if (!currMb->mb_field)
        {
          // frame
          if (curr_mb_nr % 2 == 0)
          {
            // top
            pix->mb_addr  = currMb->mbAddrC + 1;
            pix->available = currMb->mbAvailC;
             yM      = yN;
          }
          else
          {
            // bottom
            pix->available = 0;
          }
        }
        else
        {
          // field
          if (curr_mb_nr % 2 == 0)
          {
            // top
            pix->mb_addr  = currMb->mbAddrC;
            pix->available = currMb->mbAvailC;
            if (currMb->mbAvailC)
            {
              if(!img->mb_data[currMb->mbAddrC].mb_field)
              {
                (pix->mb_addr)++;
                 yM = 2* yN;
              }
              else
              {
                 yM = yN;
              }
            }
          }
          else
          {
            // bottom
            pix->mb_addr  = currMb->mbAddrC + 1;
            pix->available = currMb->mbAvailC;
             yM      = yN;
          }
        }
      }
    }
  }
  if (pix->available || img->DeblockCall)
  {
    pix->x = (xN + maxWH) % maxWH;
    pix->y = (yM + maxWH) % maxWH;
    get_mb_pos(pix->mb_addr, &(pix->pos_x), &(pix->pos_y));
    if (luma)
    {
      pix->pos_x += pix->x;
      pix->pos_y += pix->y;
    }
    else
    {
      pix->pos_x = (pix->pos_x/2) + pix->x;
      pix->pos_y = (pix->pos_y/2) + pix->y;
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    get neighbouring positions. MB AFF is automatically used from img structure
 * \param curr_mb_nr
 *   current macroblock number (decoding order)
 * \param xN
 *    input x position
 * \param yN
 *    input y position
 * \param luma
 *    1 if luma coding, 0 for chroma
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void getNeighbour(int curr_mb_nr, int xN, int yN, int luma, PixelPos *pix)
{
  if (curr_mb_nr<0)
     error ("getNeighbour: invalid macroblock number", 100);  //  原来的太复杂了 就改成如下，还少了结束编码把？
//      fprintf(stderr, "%s\n", text);                   不行 text没有定义  直接去掉算了

  if (img->MbaffFrameFlag)//0
    getAffNeighbour(curr_mb_nr, xN, yN, luma, pix);
  else
    getNonAffNeighbour(curr_mb_nr, xN, yN, luma, pix);
}


/*!
 ************************************************************************
 * \brief
 *    get neighbouring  get neighbouring 4x4 luma block
 * \param curr_mb_nr
 *   current macroblock number (decoding order)
 * \param block_x
 *    input x block position
 * \param block_y
 *    input y block position
 * \param rel_x
 *    relative x position of neighbor
 * \param rel_y
 *    relative y position of neighbor
 * \param pix
 *    returns position informations
 ************************************************************************
 */

// 求出4*4的邻接块，block_x block_y表示子块的大小，rel_x rel_y表示领接点与子块左上点的相对位置

void getLuma4x4Neighbour (int curr_mb_nr, int block_x, int block_y, int rel_x, int rel_y, PixelPos *pix)
{
	int x = 4* block_x + rel_x;//求出领接点在宏块内的位置，也就是和当前宏块左上角点的相对位置
  int y = 4* block_y + rel_y;

  getNeighbour(curr_mb_nr, x, y, 1, pix);//这里的１对应luma参数

  if (pix->available)
  {
    pix->x /= 4;
    pix->y /= 4;
    pix->pos_x /= 4;
    pix->pos_y /= 4;
  }
}


/*!
 ************************************************************************
 * \brief
 *    get neighbouring 4x4 chroma block
 * \param curr_mb_nr
 *   current macroblock number (decoding order)
 * \param block_x
 *    input x block position
 * \param block_y
 *    input y block position
 * \param rel_x
 *    relative x position of neighbor
 * \param rel_y
 *    relative y position of neighbor
 * \param pix
 *    returns position informations
 ************************************************************************
 */
void getChroma4x4Neighbour (int curr_mb_nr, int block_x, int block_y, int rel_x, int rel_y, PixelPos *pix)
{
  int x = 4* block_x + rel_x;
  int y = 4* block_y + rel_y;

  getNeighbour(curr_mb_nr, x, y, 0, pix);

  if (pix->available)
  {
    pix->x /= 4;
    pix->y /= 4;
    pix->pos_x /= 4;
    pix->pos_y /= 4;
  }
}
