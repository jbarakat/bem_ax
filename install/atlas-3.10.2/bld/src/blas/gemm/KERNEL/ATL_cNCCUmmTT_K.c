#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#include "atlas_prefetch.h"

#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#include "atlas_prefetch.h"

static void ATL_cJIK0x0x24TT1x1x24_aX_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=T, MB=0, NB=0, KB=24, 
 * lda=0, ldb=0, ldc=0, mu=1, nu=1, ku=24, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   #define Nb N
   const float *stM = A + (((lda*Mb) << 1));
   const float *stN = B + (((Nb) << 1));
   const float *pfA = A + M;
   const float BetaAlpha = beta / alpha;
   #define incAk 48
   const int incAm = ((((lda) - 24) << 1)), incAn = -(((Mb*lda) << 1));
   const int incBk = (((ldb) << 1)), incBm = -(((ldb) << 5)+((ldb) << 4));
   #define incBn 2
   #define incCm 2
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   float *pC0=C;
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0;
   register float rB0;
   register float rC0_0;
   do /* N-loop */
   {
      ATL_pfl1R(pfA+0);
      pfA += lda;
      do /* M-loop */
      {
         rA0 = BetaAlpha;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[2];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[4];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[6];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[8];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[10];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[12];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[14];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[16];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[18];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[20];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[22];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[24];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[26];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[28];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[30];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[32];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[34];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[36];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[38];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[40];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[42];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[44];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         rA0 = pA0[46];
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         pB0 += incBk;
         pA0 += incAk;
         rB0 = alpha;
         rC0_0 *= rB0;
         *pC0 = rC0_0;
         pC0 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pA0 += incAn;
      pB0 += incBn;
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
void ATL_cJIK0x0x24TT0x0x0_aX_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=T, MB=0, NB=0, KB=24, 
 * lda=0, ldb=0, ldc=0, mu=12, nu=1, ku=24, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/12)*12;
   #define Nb N
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((lda*Mb) << 1));
   const float *stN = B + (((Nb) << 1));
   const float *pfA = A + M;
   const float BetaAlpha = beta / alpha;
   #define incAk 48
   const int incAm = ((((((lda) << 3)+((lda) << 2)) - 24) << 1)), incAn = -(((Mb*lda) << 1));
   const int incBk = (((ldb) << 1)), incBm = -(((ldb) << 5)+((ldb) << 4));
   #define incBn 2
   #define incCm 24
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   float *pC0=C;
   const float *pA0=A, *pA1=pA0+(((lda) << 1)), *pA2=pA1+(((lda) << 1)), *pA3=pA2+(((lda) << 1)), *pA4=pA3+(((lda) << 1)), *pA5=pA4+(((lda) << 1)), *pA6=pA5+(((lda) << 1)), *pA7=pA6+(((lda) << 1)), *pA8=pA7+(((lda) << 1)), *pA9=pA8+(((lda) << 1)), *pA10=pA9+(((lda) << 1)), *pA11=pA10+(((lda) << 1));
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9, rA10, rA11;
   register float rB0;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0, rC10_0, rC11_0;
   if (pA0 != stM)
   {
      do /* N-loop */
      {
         ATL_pfl1R(pfA+0);
         pfA += lda;
         do /* M-loop */
         {
            rA0 = BetaAlpha;
            rC0_0 = *pC0;
            rC0_0 *= rA0;
            rC1_0 = pC0[2];
            rC1_0 *= rA0;
            rC2_0 = pC0[4];
            rC2_0 *= rA0;
            rC3_0 = pC0[6];
            rC3_0 *= rA0;
            rC4_0 = pC0[8];
            rC4_0 *= rA0;
            rC5_0 = pC0[10];
            rC5_0 *= rA0;
            rC6_0 = pC0[12];
            rC6_0 *= rA0;
            rC7_0 = pC0[14];
            rC7_0 *= rA0;
            rC8_0 = pC0[16];
            rC8_0 *= rA0;
            rC9_0 = pC0[18];
            rC9_0 *= rA0;
            rC10_0 = pC0[20];
            rC10_0 *= rA0;
            rC11_0 = pC0[22];
            rC11_0 *= rA0;
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = *pA1;
            rA2 = *pA2;
            rA3 = *pA3;
            rA4 = *pA4;
            rA5 = *pA5;
            rA6 = *pA6;
            rA7 = *pA7;
            rA8 = *pA8;
            rA9 = *pA9;
            rA10 = *pA10;
            rA11 = *pA11;
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[2];
            rB0 = *pB0;
            rA1 = pA1[2];
            rA2 = pA2[2];
            rA3 = pA3[2];
            rA4 = pA4[2];
            rA5 = pA5[2];
            rA6 = pA6[2];
            rA7 = pA7[2];
            rA8 = pA8[2];
            rA9 = pA9[2];
            rA10 = pA10[2];
            rA11 = pA11[2];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[4];
            rB0 = *pB0;
            rA1 = pA1[4];
            rA2 = pA2[4];
            rA3 = pA3[4];
            rA4 = pA4[4];
            rA5 = pA5[4];
            rA6 = pA6[4];
            rA7 = pA7[4];
            rA8 = pA8[4];
            rA9 = pA9[4];
            rA10 = pA10[4];
            rA11 = pA11[4];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[6];
            rB0 = *pB0;
            rA1 = pA1[6];
            rA2 = pA2[6];
            rA3 = pA3[6];
            rA4 = pA4[6];
            rA5 = pA5[6];
            rA6 = pA6[6];
            rA7 = pA7[6];
            rA8 = pA8[6];
            rA9 = pA9[6];
            rA10 = pA10[6];
            rA11 = pA11[6];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[8];
            rB0 = *pB0;
            rA1 = pA1[8];
            rA2 = pA2[8];
            rA3 = pA3[8];
            rA4 = pA4[8];
            rA5 = pA5[8];
            rA6 = pA6[8];
            rA7 = pA7[8];
            rA8 = pA8[8];
            rA9 = pA9[8];
            rA10 = pA10[8];
            rA11 = pA11[8];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[10];
            rB0 = *pB0;
            rA1 = pA1[10];
            rA2 = pA2[10];
            rA3 = pA3[10];
            rA4 = pA4[10];
            rA5 = pA5[10];
            rA6 = pA6[10];
            rA7 = pA7[10];
            rA8 = pA8[10];
            rA9 = pA9[10];
            rA10 = pA10[10];
            rA11 = pA11[10];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[12];
            rB0 = *pB0;
            rA1 = pA1[12];
            rA2 = pA2[12];
            rA3 = pA3[12];
            rA4 = pA4[12];
            rA5 = pA5[12];
            rA6 = pA6[12];
            rA7 = pA7[12];
            rA8 = pA8[12];
            rA9 = pA9[12];
            rA10 = pA10[12];
            rA11 = pA11[12];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[14];
            rB0 = *pB0;
            rA1 = pA1[14];
            rA2 = pA2[14];
            rA3 = pA3[14];
            rA4 = pA4[14];
            rA5 = pA5[14];
            rA6 = pA6[14];
            rA7 = pA7[14];
            rA8 = pA8[14];
            rA9 = pA9[14];
            rA10 = pA10[14];
            rA11 = pA11[14];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[16];
            rB0 = *pB0;
            rA1 = pA1[16];
            rA2 = pA2[16];
            rA3 = pA3[16];
            rA4 = pA4[16];
            rA5 = pA5[16];
            rA6 = pA6[16];
            rA7 = pA7[16];
            rA8 = pA8[16];
            rA9 = pA9[16];
            rA10 = pA10[16];
            rA11 = pA11[16];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[18];
            rB0 = *pB0;
            rA1 = pA1[18];
            rA2 = pA2[18];
            rA3 = pA3[18];
            rA4 = pA4[18];
            rA5 = pA5[18];
            rA6 = pA6[18];
            rA7 = pA7[18];
            rA8 = pA8[18];
            rA9 = pA9[18];
            rA10 = pA10[18];
            rA11 = pA11[18];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[20];
            rB0 = *pB0;
            rA1 = pA1[20];
            rA2 = pA2[20];
            rA3 = pA3[20];
            rA4 = pA4[20];
            rA5 = pA5[20];
            rA6 = pA6[20];
            rA7 = pA7[20];
            rA8 = pA8[20];
            rA9 = pA9[20];
            rA10 = pA10[20];
            rA11 = pA11[20];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[22];
            rB0 = *pB0;
            rA1 = pA1[22];
            rA2 = pA2[22];
            rA3 = pA3[22];
            rA4 = pA4[22];
            rA5 = pA5[22];
            rA6 = pA6[22];
            rA7 = pA7[22];
            rA8 = pA8[22];
            rA9 = pA9[22];
            rA10 = pA10[22];
            rA11 = pA11[22];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[24];
            rB0 = *pB0;
            rA1 = pA1[24];
            rA2 = pA2[24];
            rA3 = pA3[24];
            rA4 = pA4[24];
            rA5 = pA5[24];
            rA6 = pA6[24];
            rA7 = pA7[24];
            rA8 = pA8[24];
            rA9 = pA9[24];
            rA10 = pA10[24];
            rA11 = pA11[24];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[26];
            rB0 = *pB0;
            rA1 = pA1[26];
            rA2 = pA2[26];
            rA3 = pA3[26];
            rA4 = pA4[26];
            rA5 = pA5[26];
            rA6 = pA6[26];
            rA7 = pA7[26];
            rA8 = pA8[26];
            rA9 = pA9[26];
            rA10 = pA10[26];
            rA11 = pA11[26];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[28];
            rB0 = *pB0;
            rA1 = pA1[28];
            rA2 = pA2[28];
            rA3 = pA3[28];
            rA4 = pA4[28];
            rA5 = pA5[28];
            rA6 = pA6[28];
            rA7 = pA7[28];
            rA8 = pA8[28];
            rA9 = pA9[28];
            rA10 = pA10[28];
            rA11 = pA11[28];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[30];
            rB0 = *pB0;
            rA1 = pA1[30];
            rA2 = pA2[30];
            rA3 = pA3[30];
            rA4 = pA4[30];
            rA5 = pA5[30];
            rA6 = pA6[30];
            rA7 = pA7[30];
            rA8 = pA8[30];
            rA9 = pA9[30];
            rA10 = pA10[30];
            rA11 = pA11[30];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[32];
            rB0 = *pB0;
            rA1 = pA1[32];
            rA2 = pA2[32];
            rA3 = pA3[32];
            rA4 = pA4[32];
            rA5 = pA5[32];
            rA6 = pA6[32];
            rA7 = pA7[32];
            rA8 = pA8[32];
            rA9 = pA9[32];
            rA10 = pA10[32];
            rA11 = pA11[32];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[34];
            rB0 = *pB0;
            rA1 = pA1[34];
            rA2 = pA2[34];
            rA3 = pA3[34];
            rA4 = pA4[34];
            rA5 = pA5[34];
            rA6 = pA6[34];
            rA7 = pA7[34];
            rA8 = pA8[34];
            rA9 = pA9[34];
            rA10 = pA10[34];
            rA11 = pA11[34];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[36];
            rB0 = *pB0;
            rA1 = pA1[36];
            rA2 = pA2[36];
            rA3 = pA3[36];
            rA4 = pA4[36];
            rA5 = pA5[36];
            rA6 = pA6[36];
            rA7 = pA7[36];
            rA8 = pA8[36];
            rA9 = pA9[36];
            rA10 = pA10[36];
            rA11 = pA11[36];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[38];
            rB0 = *pB0;
            rA1 = pA1[38];
            rA2 = pA2[38];
            rA3 = pA3[38];
            rA4 = pA4[38];
            rA5 = pA5[38];
            rA6 = pA6[38];
            rA7 = pA7[38];
            rA8 = pA8[38];
            rA9 = pA9[38];
            rA10 = pA10[38];
            rA11 = pA11[38];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[40];
            rB0 = *pB0;
            rA1 = pA1[40];
            rA2 = pA2[40];
            rA3 = pA3[40];
            rA4 = pA4[40];
            rA5 = pA5[40];
            rA6 = pA6[40];
            rA7 = pA7[40];
            rA8 = pA8[40];
            rA9 = pA9[40];
            rA10 = pA10[40];
            rA11 = pA11[40];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[42];
            rB0 = *pB0;
            rA1 = pA1[42];
            rA2 = pA2[42];
            rA3 = pA3[42];
            rA4 = pA4[42];
            rA5 = pA5[42];
            rA6 = pA6[42];
            rA7 = pA7[42];
            rA8 = pA8[42];
            rA9 = pA9[42];
            rA10 = pA10[42];
            rA11 = pA11[42];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[44];
            rB0 = *pB0;
            rA1 = pA1[44];
            rA2 = pA2[44];
            rA3 = pA3[44];
            rA4 = pA4[44];
            rA5 = pA5[44];
            rA6 = pA6[44];
            rA7 = pA7[44];
            rA8 = pA8[44];
            rA9 = pA9[44];
            rA10 = pA10[44];
            rA11 = pA11[44];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            rA0 = pA0[46];
            rB0 = *pB0;
            rA1 = pA1[46];
            rA2 = pA2[46];
            rA3 = pA3[46];
            rA4 = pA4[46];
            rA5 = pA5[46];
            rA6 = pA6[46];
            rA7 = pA7[46];
            rA8 = pA8[46];
            rA9 = pA9[46];
            rA10 = pA10[46];
            rA11 = pA11[46];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rC4_0 += rA4 * rB0;
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            pB0 += incBk;
            pA0 += incAk;
            pA1 += incAk;
            pA2 += incAk;
            pA3 += incAk;
            pA4 += incAk;
            pA5 += incAk;
            pA6 += incAk;
            pA7 += incAk;
            pA8 += incAk;
            pA9 += incAk;
            pA10 += incAk;
            pA11 += incAk;
            rB0 = alpha;
            rC0_0 *= rB0;
            rC1_0 *= rB0;
            rC2_0 *= rB0;
            rC3_0 *= rB0;
            rC4_0 *= rB0;
            rC5_0 *= rB0;
            rC6_0 *= rB0;
            rC7_0 *= rB0;
            rC8_0 *= rB0;
            rC9_0 *= rB0;
            rC10_0 *= rB0;
            rC11_0 *= rB0;
            *pC0 = rC0_0;
            pC0[2] = rC1_0;
            pC0[4] = rC2_0;
            pC0[6] = rC3_0;
            pC0[8] = rC4_0;
            pC0[10] = rC5_0;
            pC0[12] = rC6_0;
            pC0[14] = rC7_0;
            pC0[16] = rC8_0;
            pC0[18] = rC9_0;
            pC0[20] = rC10_0;
            pC0[22] = rC11_0;
            pC0 += incCm;
            pA0 += incAm;
            pA1 += incAm;
            pA2 += incAm;
            pA3 += incAm;
            pA4 += incAm;
            pA5 += incAm;
            pA6 += incAm;
            pA7 += incAm;
            pA8 += incAm;
            pA9 += incAm;
            pA10 += incAm;
            pA11 += incAm;
            pB0 += incBm;
         }
         while(pA0 != stM);
         pC0 += incCn;
         pA0 += incAn;
         pA1 += incAn;
         pA2 += incAn;
         pA3 += incAn;
         pA4 += incAn;
         pA5 += incAn;
         pA6 += incAn;
         pA7 += incAn;
         pA8 += incAn;
         pA9 += incAn;
         pA10 += incAn;
         pA11 += incAn;
         pB0 += incBn;
      }
      while(pB0 != stN);
   }
   if (k=M-Mb)
      ATL_cJIK0x0x24TT1x1x24_aX_bX(k, N, 24, alpha, ca + (((Mb*lda) << 1)), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
