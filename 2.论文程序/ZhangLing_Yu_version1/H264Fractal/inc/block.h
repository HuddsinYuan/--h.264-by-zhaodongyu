
/*!
 ************************************************************************
 * \file block.h
 *
 * \author
 *  Inge Lille-Lang�y               <inge.lille-langoy@telenor.com>    \n
 *  Telenor Satellite Services                                         \n
 *  P.O.Box 6914 St.Olavs plass                                        \n
 *  N-0130 Oslo, Norway
 *
 ************************************************************************
 */

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "global.h"
/////////////////////////////////////>

#define DQ_BITS         6
#define DQ_ROUND        (1<<(DQ_BITS-1))

extern const byte QP_SCALE_CR[52] ;
extern const int  dequant_coef[6][4][4];




/////////////////////////////////////<

//! make chroma QP from quant
// const byte QP_SCALE_CR[52]=
// {
//     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
//    12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
//    28,29,29,30,31,32,32,33,34,34,35,35,36,36,37,37,
//    37,38,38,38,39,39,39,39
// };


//! single scan pattern
extern const byte SNGL_SCAN[16][2];
// const byte SNGL_SCAN[16][2] =
// {
//   {0,0},{1,0},{0,1},{0,2},
//   {1,1},{2,0},{3,0},{2,1},
//   {1,2},{0,3},{1,3},{2,2},
//   {3,1},{3,2},{2,3},{3,3}
// };
// 
//! field scan pattern

extern const byte FIELD_SCAN[16][2];

// const byte FIELD_SCAN[16][2] =
// {
//   {0,0},{0,1},{1,0},{0,2},
//   {0,3},{1,1},{1,2},{1,3},
//   {2,0},{2,1},{2,2},{2,3},
//   {3,0},{3,1},{3,2},{3,3}
// };


//! array used to find expencive coefficients
// extern const byte COEFF_COST[16];

// const byte COEFF_COST[16] =
// {
// 	3,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0
// };



//! bit cost for coefficients

extern const byte COEFF_BIT_COST[3][16][16];


// const byte COEFF_BIT_COST[3][16][16]=
// {
//   { // 2x2 scan (corrested per Gisle's Email 11/23/2000 by StW
//     { 3, 5, 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13},
//     { 5, 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13,13},
//     { 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13,13,15},
//     { 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13,13,15},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//     { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
//   },
//   {  // double scan
//     { 3, 5, 7, 7, 7, 9, 9, 9, 9,11,11,13,13,13,13,15},
//     { 5, 9, 9,11,11,13,13,13,13,15,15,15,15,15,15,15},
//     { 7,11,11,13,13,13,13,15,15,15,15,15,15,15,15,17},
//     { 9,11,11,13,13,13,13,15,15,15,15,15,15,15,15,17},
//     { 9,11,11,13,13,13,13,15,15,15,15,15,15,15,15,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//     {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
//   },
//   {    // single scan
//     { 3, 7, 9, 9,11,13,13,15,15,15,15,17,17,17,17,17},
//     { 5, 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17},
//     { 5, 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17},
//     { 7,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     { 7,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     { 7,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
//     {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
//     {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
//     {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
//     {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
//     {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
//     {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
//   },
// };
#endif
