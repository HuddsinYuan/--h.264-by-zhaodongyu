
/*!
 **************************************************************************
 *  \file elements.h
 *  \brief  Header file for elements in H.264 streams
 *  \date 6.10.2000, 
 *  \version
 *      1.1
 *
 * \note
 *    Version 1.0 included three partition modes, no DP, 2 partitions per slice
 *      and 4 partitions per slice.  As per document VCEG-N72 this is changed
 *      in version 1.1 to only two partition modes, one without DP and one with 
 *      3 partition per slice
 *
 *  \author Sebastian Purreiter     <sebastian.purreiter@mch.siemens.de>
 *  \author Stephan Wenger          <stewe@cs.tu-berlin.de>
 *
 **************************************************************************
 */

#ifndef _ELEMENTS_H_
#define _ELEMENTS_H_

#include "global.h"

/*!
 *  definition of H.264 syntax elements
 *  order of elements follow dependencies for picture reconstruction
 */
/*!
 * \brief   Assignment of old TYPE or partition elements to new
 *          elements
 *
 *  old element     | new elements
 *  ----------------+-------------------------------------------------------------------
 *  TYPE_HEADER     | SE_HEADER, SE_PTYPE
 *  TYPE_MBHEADER   | SE_MBTYPE, SE_REFFRAME, SE_INTRAPREDMODE
 *  TYPE_MVD        | SE_MVD
 *  TYPE_CBP        | SE_CBP_INTRA, SE_CBP_INTER
 *  TYPE_COEFF_Y    | SE_LUM_DC_INTRA, SE_LUM_AC_INTRA, SE_LUM_DC_INTER, SE_LUM_AC_INTER
 *  TYPE_2x2DC      | SE_CHR_DC_INTRA, SE_CHR_DC_INTER
 *  TYPE_COEFF_C    | SE_CHR_AC_INTRA, SE_CHR_AC_INTER
 *  TYPE_EOS        | SE_EOS
*/




#define MAXPARTITIONMODES 2 //!< maximum possible partition modes as defined in assignSE2partition[][]

/*!
 *  \brief  lookup-table to assign different elements to partition
 *
 *  \note here we defined up to 6 different partitions similar to
 *      document Q15-k-18 described in the PROGFRAMEMODE.
 *      The Sliceheader contains the PSYNC information. \par
 *
 *      Elements inside a partition are not ordered. They are
 *      ordered by occurence in the stream.
 *      Assumption: Only partitionlosses are considered. \par
 *
 *      The texture elements luminance and chrominance are
 *      not ordered in the progressive form
 *      This may be changed in image.c \par
 *
 *  -IMPORTANT:
 *      Picture- or Sliceheaders must be assigned to partition 0. \par
 *      Furthermore partitions must follow syntax dependencies as
 *      outlined in document Q15-J-23.
 */


// A note on this table:
//
// While the assignment of values in enum data types is specified in C, it is not
// very ood style to have an "elementnumber", not even as a comment.
//
// Hence a copy of the relevant structure from global.h here
/*
typedef enum {
 0  SE_HEADER,
 1  SE_PTYPE,
 2  SE_MBTYPE,
 3  SE_REFFRAME,
 4  SE_INTRAPREDMODE,
 5  SE_MVD,
 6  SE_CBP_INTRA,
 7  SE_LUM_DC_INTRA,
 8  SE_CHR_DC_INTRA,
 9  SE_LUM_AC_INTRA,
10  SE_CHR_AC_INTRA,
11  SE_CBP_INTER,
12  SE_LUM_DC_INTER,
13  SE_CHR_DC_INTER,
14  SE_LUM_AC_INTER,
15  SE_CHR_AC_INTER,
16  SE_DELTA_QUANT_INTER,
17  SE_DELTA_QUANT_INTRA,
18  SE_BFRAME,
19  SE_EOS,
20  SE_MAX_ELEMENTS */ // number of maximum syntax elements
//} SE_type;


int * assignSE2partition[2];
int assignSE2partition_NoDP[SE_MAX_ELEMENTS];//没有数据分割，数据分割是给要传递的数据一个等级，等级越小越重要。
int assignSE2partition_DP[SE_MAX_ELEMENTS];
/////////////////////////////////////>

extern int assignSE2partition_dec[][SE_MAX_ELEMENTS];



/////////////////////////////////////<
#endif

