/*Header file for dct-functions. Includes fdct-function for 2D forward   *
 *  dct transform, idct-function for 2D inverse dct transfrom and couple of *
 * constants which are needed in transformations.                           */
#ifndef _DCT_H_
#define _DCT_H_
#include <math.h>
#define W1 0.5*sqrt(2.0+sqrt(2.0)) /*cos(pi/8), sin(3pi/8)*/
#define W2 0.5*sqrt(2.0-sqrt(2.0)) /*cos(3pi/8), sin(pi/8)*/
#define W3 1/sqrt(2.0) /*cos(pi/4), sin(pi/4)*/
#define PI acos(-1.0)
#define W4 cos(PI/16.0) /*cos(pi/16), sin(7pi/16)*/
#define W5 sin(PI/16.0) /*cos(7pi/16), sin(pi/16)*/
#define W6 cos(3.0*PI/16.0) /*cos(3pi/16), sin(5pi/16)*/
#define W7 cos(5.0*PI/16.0) /*cos(5pi/16), sin(3pi/16)*/

/*extern void fdct(const double *inblock, double *outblock)                 *
 * *2D forward DCT transformation function for 8x8 blocks, uses  butterfly     *
 * *implementation.                                                            *
 * *Input parameters are:                                                      *
 * * const double *inblock: 8x8 block, contains pixel values to transform      *
 * * double *outblock: 8x8 block for transfrom coefficients                    *
 * *NOTE: inblock has to be scanned in horizontal order i.e. line by line      */ 
extern void fdct(const double *inblock, double *outblock);

/*extern void idct(const double *inblock, double *outblock)                 *
 * *2D inverse DCT transformation function for 8x8 blocks, uses  butterfly     *
 * *implementation.                                                            *
 * *Input parameters are:                                                      *
 * * const double *inblock: 8x8 block, contains coefficients for inverse       *
 * *                        transform                                          * 
 * * double *outblock: 8x8 block for pixel values                              *
 * *NOTE: transform coefficients have to scanned line by line                  *
 * *      e.g. 0 1  is 0123                                                    *
 * *           2 3                                                             */
extern void idct(const double *inblock, double *outblock);
#endif
