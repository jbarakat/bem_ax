#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#include "atlas_prefetch.h"

void ATL_cJIK40x40x40TN0x0x0_aX_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=40, NB=40, KB=40, 
 * lda=0, ldb=0, ldc=0, mu=2, nu=5, ku=40, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const float *stM = A + (((lda) << 6)+((lda) << 4));
   const float *stN = B + (((ldb) << 6)+((ldb) << 4));
   const float *pfA = A + M;
   const float BetaAlpha = beta / alpha;
   #define incAk 80
   const int incAm = ((((((lda) << 1)) - 40) << 1)), incAn = -(((lda) << 6)+((lda) << 4));
   #define incBk 80
   const int incBm = -80, incBn = (((ldb) << 3)+((ldb) << 1));
   const int incAk0 = ((incAk) / 40), incBk0 = ((incBk) / 40);
   #define incCm 4
   const int incCn = (((ldc) << 3)+((ldc) << 1)) - 80;
   float *pC0=C, *pC1=pC0+(((ldc) << 1)), *pC2=pC1+(((ldc) << 1)), *pC3=pC2+(((ldc) << 1)), *pC4=pC3+(((ldc) << 1));
   const float *pA0=A, *pA1=pA0+(((lda) << 1));
   const float *pB0=B, *pB1=pB0+(((ldb) << 1)), *pB2=pB1+(((ldb) << 1)), *pB3=pB2+(((ldb) << 1)), *pB4=pB3+(((ldb) << 1));
   register int k;
   register float rA0, rA1;
   register float rB0, rB1, rB2, rB3, rB4;
   register float m0;
   register float rC0_0, rC1_0, rC0_1, rC1_1, rC0_2, rC1_2, rC0_3, rC1_3, rC0_4, rC1_4;
   do /* N-loop */
   {
      ATL_pfl1R(pfA+0);
      ATL_pfl1R(pfA+16);
      pfA += lda;
      do /* M-loop */
      {
         rA0 = BetaAlpha;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[2];
         rC1_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC1_1 = pC1[2];
         rC1_1 *= rA0;
         rC0_2 = *pC2;
         rC0_2 *= rA0;
         rC1_2 = pC2[2];
         rC1_2 *= rA0;
         rC0_3 = *pC3;
         rC0_3 *= rA0;
         rC1_3 = pC3[2];
         rC1_3 *= rA0;
         rC0_4 = *pC4;
         rC0_4 *= rA0;
         rC1_4 = pC4[2];
         rC1_4 *= rA0;
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = *pA1;
         m0 = rA0 * rB0;
         rB1 = *pB1;
         rB2 = *pB2;
         rB3 = *pB3;

/*
 *       Completely unrolled K-loop
 */
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rB4 = *pB4;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rA1 = pA1[2];
         rB1 = pB1[2];
         rB2 = pB2[2];
         rB3 = pB3[2];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[2];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rA1 = pA1[4];
         rB1 = pB1[4];
         rB2 = pB2[4];
         rB3 = pB3[4];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[4];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rA1 = pA1[6];
         rB1 = pB1[6];
         rB2 = pB2[6];
         rB3 = pB3[6];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[6];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rA1 = pA1[8];
         rB1 = pB1[8];
         rB2 = pB2[8];
         rB3 = pB3[8];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[8];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rA1 = pA1[10];
         rB1 = pB1[10];
         rB2 = pB2[10];
         rB3 = pB3[10];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[10];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rA1 = pA1[12];
         rB1 = pB1[12];
         rB2 = pB2[12];
         rB3 = pB3[12];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[12];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rA1 = pA1[14];
         rB1 = pB1[14];
         rB2 = pB2[14];
         rB3 = pB3[14];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[14];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rA1 = pA1[16];
         rB1 = pB1[16];
         rB2 = pB2[16];
         rB3 = pB3[16];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[16];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rA1 = pA1[18];
         rB1 = pB1[18];
         rB2 = pB2[18];
         rB3 = pB3[18];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[18];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rA1 = pA1[20];
         rB1 = pB1[20];
         rB2 = pB2[20];
         rB3 = pB3[20];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[20];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rA1 = pA1[22];
         rB1 = pB1[22];
         rB2 = pB2[22];
         rB3 = pB3[22];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[22];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rA1 = pA1[24];
         rB1 = pB1[24];
         rB2 = pB2[24];
         rB3 = pB3[24];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[24];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rA1 = pA1[26];
         rB1 = pB1[26];
         rB2 = pB2[26];
         rB3 = pB3[26];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[26];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rA1 = pA1[28];
         rB1 = pB1[28];
         rB2 = pB2[28];
         rB3 = pB3[28];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[28];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[30];
         rB0 = pB0[30];
         rA1 = pA1[30];
         rB1 = pB1[30];
         rB2 = pB2[30];
         rB3 = pB3[30];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[30];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[32];
         rB0 = pB0[32];
         rA1 = pA1[32];
         rB1 = pB1[32];
         rB2 = pB2[32];
         rB3 = pB3[32];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[32];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[34];
         rB0 = pB0[34];
         rA1 = pA1[34];
         rB1 = pB1[34];
         rB2 = pB2[34];
         rB3 = pB3[34];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[34];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[36];
         rB0 = pB0[36];
         rA1 = pA1[36];
         rB1 = pB1[36];
         rB2 = pB2[36];
         rB3 = pB3[36];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[36];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[38];
         rB0 = pB0[38];
         rA1 = pA1[38];
         rB1 = pB1[38];
         rB2 = pB2[38];
         rB3 = pB3[38];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[38];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[40];
         rB0 = pB0[40];
         rA1 = pA1[40];
         rB1 = pB1[40];
         rB2 = pB2[40];
         rB3 = pB3[40];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[40];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[42];
         rB0 = pB0[42];
         rA1 = pA1[42];
         rB1 = pB1[42];
         rB2 = pB2[42];
         rB3 = pB3[42];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[42];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[44];
         rB0 = pB0[44];
         rA1 = pA1[44];
         rB1 = pB1[44];
         rB2 = pB2[44];
         rB3 = pB3[44];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[44];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[46];
         rB0 = pB0[46];
         rA1 = pA1[46];
         rB1 = pB1[46];
         rB2 = pB2[46];
         rB3 = pB3[46];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[46];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[48];
         rB0 = pB0[48];
         rA1 = pA1[48];
         rB1 = pB1[48];
         rB2 = pB2[48];
         rB3 = pB3[48];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[48];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[50];
         rB0 = pB0[50];
         rA1 = pA1[50];
         rB1 = pB1[50];
         rB2 = pB2[50];
         rB3 = pB3[50];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[50];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[52];
         rB0 = pB0[52];
         rA1 = pA1[52];
         rB1 = pB1[52];
         rB2 = pB2[52];
         rB3 = pB3[52];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[52];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[54];
         rB0 = pB0[54];
         rA1 = pA1[54];
         rB1 = pB1[54];
         rB2 = pB2[54];
         rB3 = pB3[54];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[54];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[56];
         rB0 = pB0[56];
         rA1 = pA1[56];
         rB1 = pB1[56];
         rB2 = pB2[56];
         rB3 = pB3[56];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[56];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[58];
         rB0 = pB0[58];
         rA1 = pA1[58];
         rB1 = pB1[58];
         rB2 = pB2[58];
         rB3 = pB3[58];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[58];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[60];
         rB0 = pB0[60];
         rA1 = pA1[60];
         rB1 = pB1[60];
         rB2 = pB2[60];
         rB3 = pB3[60];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[60];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[62];
         rB0 = pB0[62];
         rA1 = pA1[62];
         rB1 = pB1[62];
         rB2 = pB2[62];
         rB3 = pB3[62];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[62];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[64];
         rB0 = pB0[64];
         rA1 = pA1[64];
         rB1 = pB1[64];
         rB2 = pB2[64];
         rB3 = pB3[64];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[64];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[66];
         rB0 = pB0[66];
         rA1 = pA1[66];
         rB1 = pB1[66];
         rB2 = pB2[66];
         rB3 = pB3[66];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[66];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[68];
         rB0 = pB0[68];
         rA1 = pA1[68];
         rB1 = pB1[68];
         rB2 = pB2[68];
         rB3 = pB3[68];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[68];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[70];
         rB0 = pB0[70];
         rA1 = pA1[70];
         rB1 = pB1[70];
         rB2 = pB2[70];
         rB3 = pB3[70];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[70];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[72];
         rB0 = pB0[72];
         rA1 = pA1[72];
         rB1 = pB1[72];
         rB2 = pB2[72];
         rB3 = pB3[72];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[72];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[74];
         rB0 = pB0[74];
         rA1 = pA1[74];
         rB1 = pB1[74];
         rB2 = pB2[74];
         rB3 = pB3[74];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[74];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[76];
         rB0 = pB0[76];
         rA1 = pA1[76];
         rB1 = pB1[76];
         rB2 = pB2[76];
         rB3 = pB3[76];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[76];
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rA0 = pA0[78];
         rB0 = pB0[78];
         rA1 = pA1[78];
         rB1 = pB1[78];
         rB2 = pB2[78];
         rB3 = pB3[78];
         rC1_4 += m0;
         m0 = rA0 * rB0;
         rB4 = pB4[78];
/*
 *       Drain pipe on last iteration of K-loop
 */
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rC1_1 += m0;
         m0 = rA0 * rB2;
         rC0_2 += m0;
         m0 = rA1 * rB2;
         rC1_2 += m0;
         m0 = rA0 * rB3;
         rC0_3 += m0;
         m0 = rA1 * rB3;
         rC1_3 += m0;
         m0 = rA0 * rB4;
         rC0_4 += m0;
         m0 = rA1 * rB4;
         rC1_4 += m0;
         pA0 += incAk;
         pA1 += incAk;
         pB0 += incBk;
         pB1 += incBk;
         pB2 += incBk;
         pB3 += incBk;
         pB4 += incBk;
         rB0 = alpha;
         rC0_0 *= rB0;
         rC0_1 *= rB0;
         rC0_2 *= rB0;
         rC0_3 *= rB0;
         rC0_4 *= rB0;
         rC1_0 *= rB0;
         rC1_1 *= rB0;
         rC1_2 *= rB0;
         rC1_3 *= rB0;
         rC1_4 *= rB0;
         *pC0 = rC0_0;
         pC0[2] = rC1_0;
         *pC1 = rC0_1;
         pC1[2] = rC1_1;
         *pC2 = rC0_2;
         pC2[2] = rC1_2;
         *pC3 = rC0_3;
         pC3[2] = rC1_3;
         *pC4 = rC0_4;
         pC4[2] = rC1_4;
         pC0 += incCm;
         pC1 += incCm;
         pC2 += incCm;
         pC3 += incCm;
         pC4 += incCm;
         pA0 += incAm;
         pA1 += incAm;
         pB0 += incBm;
         pB1 += incBm;
         pB2 += incBm;
         pB3 += incBm;
         pB4 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pC1 += incCn;
      pC2 += incCn;
      pC3 += incCn;
      pC4 += incCn;
      pA0 += incAn;
      pA1 += incAn;
      pB0 += incBn;
      pB1 += incBn;
      pB2 += incBn;
      pB3 += incBn;
      pB4 += incBn;
   }
   while(pB0 != stN);
}
#ifdef incAm
   #undef incAm
#endif
#ifdef incAn
   #undef incAn
#endif
#ifdef incAk
   #undef incAk
#endif
#ifdef incBm
   #undef incBm
#endif
#ifdef incBn
   #undef incBn
#endif
#ifdef incBk
   #undef incBk
#endif
#ifdef incCm
   #undef incCm
#endif
#ifdef incCn
   #undef incCn
#endif
#ifdef incCk
   #undef incCk
#endif
#ifdef Mb
   #undef Mb
#endif
#ifdef Nb
   #undef Nb
#endif
#ifdef Kb
   #undef Kb
#endif
