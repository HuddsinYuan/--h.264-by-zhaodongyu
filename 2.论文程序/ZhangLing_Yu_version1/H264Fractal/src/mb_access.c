
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
// A   |cur  | ����������ٽ����Ŀ�����

void CheckAvailabilityOfNeighbors()
{
	const int mb_nr = img->current_mb_nr;//�����
  Macroblock *currMB = &img->mb_data[mb_nr];//���ָ��ָ��ǰ���

  // mark all neighbors as unavailable
  currMB->mb_available_up   = NULL;//��ǰ����ϱߺ��ָ��
  currMB->mb_available_left = NULL;//��ǰ�����ߺ��ָ��

  if (img->MbaffFrameFlag)//�峡����Ӧģʽ
  {
// 	  ����qcif 176*144��ôimg->PicWidthInMbs��11
// 		 ��0   2   4.... ��ż�����������߶�����
// 		   1   3   5....������������׳����ߵ�����
// 		   22  24  26...
// 		   23  25  27...
// 		   ������4���A���Ǻ��2,���2��A�����Ǻ�飰ע���º�飳��ߺ��A�Ǻ�飰���������ڰ�32*16���Ե�Ե��

    currMB->mbAddrA = 2 * (mb_nr/2 - 1);
// 	���24���Ϻ��B�Ǻ���Ǻ��2,���26���Ϻ��B�Ǻ�飴
    currMB->mbAddrB = 2 * (mb_nr/2 - img->PicWidthInMbs);
    currMB->mbAddrC = 2 * (mb_nr/2 - img->PicWidthInMbs + 1);
    currMB->mbAddrD = 2 * (mb_nr/2 - img->PicWidthInMbs - 1);
    
	currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && (((mb_nr/2) % img->PicWidthInMbs)!=0);//��֤��ǰ��鲻������ߵĺ��
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
	currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr/2 +1) % img->PicWidthInMbs)!=0);//��֤��ǰ��鲻�����ұߵĺ��
	currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && (((mb_nr/2) % img->PicWidthInMbs)!=0);//��֤��ǰ��鲻������ߵĺ��
  }
  else
  {
// 	      D   |B    |C
// 		  -------------
// 		  A   |cur  | 
	  //��ǰ������ߺ��A
    currMB->mbAddrA = mb_nr - 1;
	//��ǰ�����ϱߺ��B
	currMB->mbAddrB = mb_nr - img->PicWidthInMbs;//����img->PicWidthInMbsָ���Ϊ��λ��1���к����Ŀ
    currMB->mbAddrC = mb_nr - img->PicWidthInMbs + 1;
    currMB->mbAddrD = mb_nr - img->PicWidthInMbs - 1;
	//����ȷ�����A B C Dλ���������ȷ��A B C D�Ŀ�����
	currMB->mbAvailA = mb_is_available(currMB->mbAddrA, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);//mb_nr % img->PicWidthInMbs��֤��ǰ��鲻������ߵĺ��
    currMB->mbAvailB = mb_is_available(currMB->mbAddrB, mb_nr);
	currMB->mbAvailC = mb_is_available(currMB->mbAddrC, mb_nr) && (((mb_nr+1) % img->PicWidthInMbs)!=0);//��֤��ǰ��鲻�����ұߵĺ��
	currMB->mbAvailD = mb_is_available(currMB->mbAddrD, mb_nr) && ((mb_nr % img->PicWidthInMbs)!=0);//��֤��ǰ��鲻������ߵĺ��
  }

  if (currMB->mbAvailA) currMB->mb_available_left = &(img->mb_data[currMB->mbAddrA]);//�����ߺ��ָ��ָ����A
  if (currMB->mbAvailB) currMB->mb_available_up   = &(img->mb_data[currMB->mbAddrB]);//����ϱߺ��ָ��ָ����B
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
    *x = (mb_addr%img->PicWidthInMbs);//���Ĺ�դɨ������ת���ɶ�άx���꣨�Ժ��Ϊ��λ��
    *y = (mb_addr/img->PicWidthInMbs);//���Ĺ�դɨ������ת���ɶ�άy���꣨�Ժ��Ϊ��λ��
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
 *    get neighbouring positions for non-aff coding����ڽ�λ��
 * \param curr_mb_nr//��ǰ�����
 *   current macroblock number (decoding order)
 * \param xN//��Ե�ǰ�����������
 *    input x position
 * \param yN//��Ե�ǰ�����������
 *    input y position
 * \param luma//���ȣ���ɫ�ȣ�
 *    1 if luma coding, 0 for chroma
 * \param pix//�����ӵ�λ����Ϣ
 *    returns position informations
 ************************************************************************
 *///��������Ӧ�峡����¡�����ڽӵ����Ǹ�������ͼ����λ��
void getNonAffNeighbour(unsigned int curr_mb_nr, int xN, int yN, int luma, PixelPos *pix)
{
	Macroblock *currMb = &img->mb_data[curr_mb_nr];//��ǰ���ָ��
  int maxWH;

  if (luma)//���ȿ������16
    maxWH = 16;
  else
    maxWH = 8;
  //xN yN��ӵ��Ǻ͵�ǰ������Ͻǵ����λ��
  if ((xN<0)&&(yN<0))//���D���
  {
    pix->mb_addr   = currMb->mbAddrD;//23
    pix->available = currMb->mbAvailD;//0
  }
  else
	  if ((xN<0)&&((yN>=0)&&(yN<maxWH)))//��ӵ����ڽ�A���

  {
    pix->mb_addr  = currMb->mbAddrA;
    pix->available = currMb->mbAvailA;
  }
  else
	  if (((xN>=0)&&(xN<maxWH))&&(yN<0))//��ӵ������B���
  {
    pix->mb_addr  = currMb->mbAddrB;
    pix->available = currMb->mbAvailB;
  }
  else
	  if (((xN>=0)&&(xN<maxWH))&&((yN>=0)&&(yN<maxWH)))//��ӵ��ڵ�ǰ���
  {
    pix->mb_addr  = curr_mb_nr;
    pix->available = 1;
  }
  else
	  if ((xN>=maxWH)&&(yN<0))//��ӵ����ڽ�C���
  {
    pix->mb_addr  = currMb->mbAddrC;
    pix->available = currMb->mbAvailC;
  }
  else 
  {
    pix->available = 0;
  }
	  if (pix->available || img->DeblockCall)//����ڽӵ����
  {
	  pix->x = (xN + maxWH) % maxWH;//����ڽӵ��ں���ڵ�λ��
    pix->y = (yN + maxWH) % maxWH;
	get_mb_pos(pix->mb_addr, &(pix->pos_x), &(pix->pos_y));//����ڽӵ����ں�����ϵ���ͼ��������
	if (luma)//���ȣ������ӵ���ͼ���е�λ��
    {
      pix->pos_x += pix->x;
      pix->pos_y += pix->y;
    }
	else//ɫ�ȡ�pix->pos_x/2��pix->pos_y/2��˵����ɫ��ˮƽ��ֱ���²���
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
     error ("getNeighbour: invalid macroblock number", 100);  //  ԭ����̫������ �͸ĳ����£������˽�������ѣ�
//      fprintf(stderr, "%s\n", text);                   ���� textû�ж���  ֱ��ȥ������

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

// ���4*4���ڽӿ飬block_x block_y��ʾ�ӿ�Ĵ�С��rel_x rel_y��ʾ��ӵ����ӿ����ϵ�����λ��

void getLuma4x4Neighbour (int curr_mb_nr, int block_x, int block_y, int rel_x, int rel_y, PixelPos *pix)
{
	int x = 4* block_x + rel_x;//�����ӵ��ں���ڵ�λ�ã�Ҳ���Ǻ͵�ǰ������Ͻǵ�����λ��
  int y = 4* block_y + rel_y;

  getNeighbour(curr_mb_nr, x, y, 1, pix);//����ģ���Ӧluma����

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
