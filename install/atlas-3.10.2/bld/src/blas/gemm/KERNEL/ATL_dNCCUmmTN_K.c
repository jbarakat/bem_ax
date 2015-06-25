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

static void ATL_dJIK0x0x24TN1x1x24_aX_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=24, 
 * lda=0, ldb=0, ldc=0, mu=1, nu=1, ku=24, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   #define Nb N
   const double *stM = A + (lda*Mb);
   const double *stN = B + (ldb*Nb);
   const double *pfA = A + M;
   const double BetaAlpha = beta / alpha;
   #define incAk 24
   const int incAm = ((lda) - 24), incAn = -(Mb*lda);
   #define incBk 24
   const int incBm = -24, incBn = (ldb);
   #define incCm 1
   const int incCn = (ldc) - (Mb);
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0;
   register double rC0_0;
   do /* N-loop */
   {
      ATL_pfl1R(pfA+0);
      ATL_pfl1R(pfA+8);
      ATL_pfl1R(pfA+16);
      pfA += lda;
      do /* M-loop */
      {
         rA0 = BetaAlpha;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rC0_0 += rA0 * rB0;
         pA0 += incAk;
         pB0 += incBk;
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
void ATL_dJIK0x0x24TN0x0x0_aX_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=24, 
 * lda=0, ldb=0, ldc=0, mu=12, nu=1, ku=24, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/12)*12;
   #define Nb N
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (lda*Mb);
   const double *stN = B + (ldb*Nb);
   const double *pfA = A + M;
   const double BetaAlpha = beta / alpha;
   #define incAk 24
   const int incAm = ((((lda) << 3)+((lda) << 2)) - 24), incAn = -(Mb*lda);
   #define incBk 24
   const int incBm = -24, incBn = (ldb);
   #define incCm 12
   const int incCn = (ldc) - (Mb);
   double *pC0=C;
   const double *pA0=A, *pA1=pA0+(lda), *pA2=pA1+(lda), *pA3=pA2+(lda), *pA4=pA3+(lda), *pA5=pA4+(lda), *pA6=pA5+(lda), *pA7=pA6+(lda), *pA8=pA7+(lda), *pA9=pA8+(lda), *pA10=pA9+(lda), *pA11=pA10+(lda);
   const double *pB0=B;
   register int k;
   register double rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9, rA10, rA11;
   register double rB0;
   register double rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0, rC10_0, rC11_0;
   if (pA0 != stM)
   {
      do /* N-loop */
      {
         ATL_pfl1R(pfA+0);
         ATL_pfl1R(pfA+8);
         ATL_pfl1R(pfA+16);
         pfA += lda;
         do /* M-loop */
         {
            rA0 = BetaAlpha;
            rC0_0 = *pC0;
            rC0_0 *= rA0;
            rC1_0 = pC0[1];
            rC1_0 *= rA0;
            rC2_0 = pC0[2];
            rC2_0 *= rA0;
            rC3_0 = pC0[3];
            rC3_0 *= rA0;
            rC4_0 = pC0[4];
            rC4_0 *= rA0;
            rC5_0 = pC0[5];
            rC5_0 *= rA0;
            rC6_0 = pC0[6];
            rC6_0 *= rA0;
            rC7_0 = pC0[7];
            rC7_0 *= rA0;
            rC8_0 = pC0[8];
            rC8_0 *= rA0;
            rC9_0 = pC0[9];
            rC9_0 *= rA0;
            rC10_0 = pC0[10];
            rC10_0 *= rA0;
            rC11_0 = pC0[11];
            rC11_0 *= rA0;
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = *pA1;
            rA2 = *pA2;
            rA3 = *pA3;
            rA4 = *pA4;
            rA5 = *pA5;
            rC0_0 += rA0 * rB0;
            rA6 = *pA6;
            rC1_0 += rA1 * rB0;
            rA7 = *pA7;
            rC2_0 += rA2 * rB0;
            rA8 = *pA8;
            rC3_0 += rA3 * rB0;
            rA9 = *pA9;
            rC4_0 += rA4 * rB0;
            rA10 = *pA10;
            rC5_0 += rA5 * rB0;
            rA11 = *pA11;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA1[1];
            rA2 = pA2[1];
            rA3 = pA3[1];
            rA4 = pA4[1];
            rA5 = pA5[1];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[1];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[1];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[1];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[1];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[1];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[1];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA1[2];
            rA2 = pA2[2];
            rA3 = pA3[2];
            rA4 = pA4[2];
            rA5 = pA5[2];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[2];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[2];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[2];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[2];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[2];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[2];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA1[3];
            rA2 = pA2[3];
            rA3 = pA3[3];
            rA4 = pA4[3];
            rA5 = pA5[3];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[3];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[3];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[3];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[3];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[3];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[3];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA1[4];
            rA2 = pA2[4];
            rA3 = pA3[4];
            rA4 = pA4[4];
            rA5 = pA5[4];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[4];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[4];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[4];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[4];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[4];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[4];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA1[5];
            rA2 = pA2[5];
            rA3 = pA3[5];
            rA4 = pA4[5];
            rA5 = pA5[5];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[5];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[5];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[5];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[5];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[5];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[5];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA1[6];
            rA2 = pA2[6];
            rA3 = pA3[6];
            rA4 = pA4[6];
            rA5 = pA5[6];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[6];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[6];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[6];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[6];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[6];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[6];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA1[7];
            rA2 = pA2[7];
            rA3 = pA3[7];
            rA4 = pA4[7];
            rA5 = pA5[7];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[7];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[7];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[7];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[7];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[7];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[7];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA1[8];
            rA2 = pA2[8];
            rA3 = pA3[8];
            rA4 = pA4[8];
            rA5 = pA5[8];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[8];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[8];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[8];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[8];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[8];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[8];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA1[9];
            rA2 = pA2[9];
            rA3 = pA3[9];
            rA4 = pA4[9];
            rA5 = pA5[9];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[9];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[9];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[9];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[9];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[9];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[9];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA1[10];
            rA2 = pA2[10];
            rA3 = pA3[10];
            rA4 = pA4[10];
            rA5 = pA5[10];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[10];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[10];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[10];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[10];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[10];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[10];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA1[11];
            rA2 = pA2[11];
            rA3 = pA3[11];
            rA4 = pA4[11];
            rA5 = pA5[11];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[11];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[11];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[11];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[11];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[11];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[11];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA1[12];
            rA2 = pA2[12];
            rA3 = pA3[12];
            rA4 = pA4[12];
            rA5 = pA5[12];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[12];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[12];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[12];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[12];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[12];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[12];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA1[13];
            rA2 = pA2[13];
            rA3 = pA3[13];
            rA4 = pA4[13];
            rA5 = pA5[13];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[13];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[13];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[13];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[13];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[13];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[13];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA1[14];
            rA2 = pA2[14];
            rA3 = pA3[14];
            rA4 = pA4[14];
            rA5 = pA5[14];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[14];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[14];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[14];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[14];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[14];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[14];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA1[15];
            rA2 = pA2[15];
            rA3 = pA3[15];
            rA4 = pA4[15];
            rA5 = pA5[15];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[15];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[15];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[15];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[15];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[15];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[15];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA1[16];
            rA2 = pA2[16];
            rA3 = pA3[16];
            rA4 = pA4[16];
            rA5 = pA5[16];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[16];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[16];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[16];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[16];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[16];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[16];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA1[17];
            rA2 = pA2[17];
            rA3 = pA3[17];
            rA4 = pA4[17];
            rA5 = pA5[17];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[17];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[17];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[17];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[17];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[17];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[17];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA1[18];
            rA2 = pA2[18];
            rA3 = pA3[18];
            rA4 = pA4[18];
            rA5 = pA5[18];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[18];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[18];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[18];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[18];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[18];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[18];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rA1 = pA1[19];
            rA2 = pA2[19];
            rA3 = pA3[19];
            rA4 = pA4[19];
            rA5 = pA5[19];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[19];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[19];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[19];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[19];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[19];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[19];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rA1 = pA1[20];
            rA2 = pA2[20];
            rA3 = pA3[20];
            rA4 = pA4[20];
            rA5 = pA5[20];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[20];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[20];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[20];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[20];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[20];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[20];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rA1 = pA1[21];
            rA2 = pA2[21];
            rA3 = pA3[21];
            rA4 = pA4[21];
            rA5 = pA5[21];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[21];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[21];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[21];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[21];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[21];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[21];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rA1 = pA1[22];
            rA2 = pA2[22];
            rA3 = pA3[22];
            rA4 = pA4[22];
            rA5 = pA5[22];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[22];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[22];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[22];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[22];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[22];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[22];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rA1 = pA1[23];
            rA2 = pA2[23];
            rA3 = pA3[23];
            rA4 = pA4[23];
            rA5 = pA5[23];
            rC0_0 += rA0 * rB0;
            rA6 = pA6[23];
            rC1_0 += rA1 * rB0;
            rA7 = pA7[23];
            rC2_0 += rA2 * rB0;
            rA8 = pA8[23];
            rC3_0 += rA3 * rB0;
            rA9 = pA9[23];
            rC4_0 += rA4 * rB0;
            rA10 = pA10[23];
            rC5_0 += rA5 * rB0;
            rA11 = pA11[23];
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rC10_0 += rA10 * rB0;
            rC11_0 += rA11 * rB0;
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
            pB0 += incBk;
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
            pC0[1] = rC1_0;
            pC0[2] = rC2_0;
            pC0[3] = rC3_0;
            pC0[4] = rC4_0;
            pC0[5] = rC5_0;
            pC0[6] = rC6_0;
            pC0[7] = rC7_0;
            pC0[8] = rC8_0;
            pC0[9] = rC9_0;
            pC0[10] = rC10_0;
            pC0[11] = rC11_0;
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
      ATL_dJIK0x0x24TN1x1x24_aX_bX(k, N, 24, alpha, ca + (Mb*lda), lda, cb, ldb, beta, cc + (Mb), ldc);
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
