#include "atlas_misc.h"
#ifndef PFYDIST
   #define PFYDIST 4
#endif
#ifndef PFXDIST
   #define PFXDIST 0
#endif
#ifndef PFADIST
   #define PFADIST 0
#endif

#if !defined(ATL_3DNow) && !defined(ATL_SSE1) && \
    (PFADIST != 0 || PFXDIST != 0 || PFYDIST != 0)
   #include "atlas_prefetch.h"
#endif

#if defined(ATL_3DNow) || defined(ATL_SSE1)
   #ifndef PFIY
      #define PFIY prefetchnta
   #endif
   #ifndef PFIX
      #define PFIX prefetchnta
   #endif
   #ifndef PFIA
      #ifdef ATL_3DNow
         #define PFIA prefetchw
      #else
         #define PFIA prefetcht0
      #endif
   #endif
#endif
/*
 * X & A are prefetched in M loop PF[A,X]DIST (in bytes) ahead
 */
#if PFADIST == 0                /* flag for no prefetch */
   #define prefA(mem)
#else
   #if defined(ATL_3DNow) || defined(ATL_SSE1)
      #define prefA(mem) __asm__ __volatile__ \
         (Mstr(PFIA) " %0" : : "m" (*(((char *)(mem))+PFADIST)))
   #else
      #if PFLVL == 2
         #define prefA(mem) ATL_pfl2W(((char*)mem)+PFADIST)
      #else
         #define prefA(mem) ATL_pfl1W(((char*)mem)+PFADIST)
      #endif
   #endif
#endif
#if PFXDIST == 0                /* flag for no prefetch */
   #define prefX(mem)
#else
   #if defined(ATL_3DNow) || defined(ATL_SSE1)
      #define prefX(mem) __asm__ __volatile__ \
         (Mstr(PFIX) " %0" : : "m" (*(((char *)(mem))+PFXDIST)))
   #else
      #if PFLVL == 2
         #define prefX(mem) ATL_pfl2R(((char*)mem)+PFXDIST)
      #else
         #define prefX(mem) ATL_pfl1R(((char*)mem)+PFXDIST)
      #endif
   #endif
#endif
/*
 * Y is prefetched in N-loop, and always fetches next NU piece
 */
#if PFYDIST == 0                /* flag for no prefetch */
   #define prefY(mem)
#else
   #if defined(ATL_3DNow) || defined(ATL_SSE1)
      #define prefY(mem) \
         __asm__ __volatile__ (Mstr(PFIY) " %0" : : "m" (*(((char *)(mem)))))
   #else
      #if PFLVL == 2
         #define prefY(mem) ATL_pfl2R(mem)
      #else
         #define prefY(mem) ATL_pfl1R(mem)
      #endif
   #endif
#endif
#ifndef ATL_CINT
   #define ATL_CINT const int
#endif
#ifndef ATL_INT
   #define ATL_INT int
#endif
/* Need to handle BETA=0 case by assigning y to zero outside of loop */
void ATL_UGEMV(ATL_CINT M, ATL_CINT N, const TYPE *A, ATL_CINT lda,
               const TYPE *X, TYPE *Y)
/*
 * 8x4 unrolled mvn_c.
 * Extracted from ATLAS/tune/blas/gemv/atlas-l2g.base
 */
{
   ATL_CINT N4=(N/4)*4, M8=(M/8)*8, lda4=lda*4;
   ATL_INT j;

   #ifdef BETA0
      for (j=0; j < M; j++)
         Y[j] = ATL_rzero;
   #endif
   for (j=N4; j; j -= 4, A += lda4, X += 4)
   {
      const double *A0=A, *A1=A0+lda, *A2=A1+lda, *A3=A2+lda;
      const register double x0=X[0], x1=X[1], x2=X[2], x3=X[3];
      ATL_INT i;
      prefY(Y+4+4-1);
      for (i=0; i < M8; i += 8)
      {
         const double a0_0=A0[i+0], a1_0=A0[i+1], a2_0=A0[i+2], a3_0=A0[i+3],
                      a4_0=A0[i+4], a5_0=A0[i+5], a6_0=A0[i+6], a7_0=A0[i+7],
                      a0_1=A1[i+0], a1_1=A1[i+1], a2_1=A1[i+2], a3_1=A1[i+3],
                      a4_1=A1[i+4], a5_1=A1[i+5], a6_1=A1[i+6], a7_1=A1[i+7],
                      a0_2=A2[i+0], a1_2=A2[i+1], a2_2=A2[i+2], a3_2=A2[i+3],
                      a4_2=A2[i+4], a5_2=A2[i+5], a6_2=A2[i+6], a7_2=A2[i+7],
                      a0_3=A3[i+0], a1_3=A3[i+1], a2_3=A3[i+2], a3_3=A3[i+3],
                      a4_3=A3[i+4], a5_3=A3[i+5], a6_3=A3[i+6], a7_3=A3[i+7];
         register double y0=Y[i+0], y1=Y[i+1], y2=Y[i+2], y3=Y[i+3], y4=Y[i+4],
                         y5=Y[i+5], y6=Y[i+6], y7=Y[i+7];

         prefX(X);
         y0 += a0_0 * x0;
         prefA(A0);
         y1 += a1_0 * x0;
         y2 += a2_0 * x0;
         y3 += a3_0 * x0;
         y4 += a4_0 * x0;
         y5 += a5_0 * x0;
         y6 += a6_0 * x0;
         y7 += a7_0 * x0;
         y0 += a0_1 * x1;
         prefA(A1);
         y1 += a1_1 * x1;
         y2 += a2_1 * x1;
         y3 += a3_1 * x1;
         y4 += a4_1 * x1;
         y5 += a5_1 * x1;
         y6 += a6_1 * x1;
         y7 += a7_1 * x1;
         y0 += a0_2 * x2;
         prefA(A2);
         y1 += a1_2 * x2;
         y2 += a2_2 * x2;
         y3 += a3_2 * x2;
         y4 += a4_2 * x2;
         y5 += a5_2 * x2;
         y6 += a6_2 * x2;
         y7 += a7_2 * x2;
         y0 += a0_3 * x3;
         prefA(A3);
         y1 += a1_3 * x3;
         y2 += a2_3 * x3;
         y3 += a3_3 * x3;
         y4 += a4_3 * x3;
         y5 += a5_3 * x3;
         y6 += a6_3 * x3;
         y7 += a7_3 * x3;
         Y[i+0] = y0;
         Y[i+1] = y1;
         Y[i+2] = y2;
         Y[i+3] = y3;
         Y[i+4] = y4;
         Y[i+5] = y5;
         Y[i+6] = y6;
         Y[i+7] = y7;
      }
      for (i=M8; i < M; i++)
      {
         const double a0_0=A0[i], a0_1=A1[i], a0_2=A2[i], a0_3=A3[i];
         register double y0 = Y[i];

         y0 += a0_0 * x0;
         y0 += a0_1 * x1;
         y0 += a0_2 * x2;
         y0 += a0_3 * x3;
         Y[i] = y0;
      }
   }
/*
 * Do remaining columns with NU=1 cleanup
 */
   for (j=N-N4; j; j--, A += lda, X++)
   {
      const double *A0=A;
      const register double x0=X[0];
      ATL_INT i;
      prefY(Y+1+1-1);
      for (i=0; i < M8; i += 8)
      {
         const double a0_0=A0[i+0], a1_0=A0[i+1], a2_0=A0[i+2], a3_0=A0[i+3],
                      a4_0=A0[i+4], a5_0=A0[i+5], a6_0=A0[i+6], a7_0=A0[i+7];
         register double y0=Y[i+0], y1=Y[i+1], y2=Y[i+2], y3=Y[i+3], y4=Y[i+4],
                         y5=Y[i+5], y6=Y[i+6], y7=Y[i+7];

         prefX(X);
         y0 += a0_0 * x0;
         prefA(A0);
         y1 += a1_0 * x0;
         y2 += a2_0 * x0;
         y3 += a3_0 * x0;
         y4 += a4_0 * x0;
         y5 += a5_0 * x0;
         y6 += a6_0 * x0;
         y7 += a7_0 * x0;
         Y[i+0] = y0;
         Y[i+1] = y1;
         Y[i+2] = y2;
         Y[i+3] = y3;
         Y[i+4] = y4;
         Y[i+5] = y5;
         Y[i+6] = y6;
         Y[i+7] = y7;
      }
      for (i=M8; i < M; i++)
      {
         const double a0_0=A0[i];
         register double y0 = Y[i];

         y0 += a0_0 * x0;
         Y[i] = y0;
      }
   }
}
