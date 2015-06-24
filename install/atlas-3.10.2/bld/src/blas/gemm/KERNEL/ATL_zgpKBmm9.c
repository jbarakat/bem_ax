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

static void ATL_zJIK0x0x9TN1x1x9_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=9, 
 * lda=9, ldb=9, ldc=0, mu=1, nu=1, ku=9, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   #define Nb N
   const double *stM = A + (((Mb) << 3)+Mb);
   const double *stN = B + (((Nb) << 3)+Nb);
   const double *pfA = stM;
   const int incPFA0 = (((int)(stM - A))*1*1)/(M*N*sizeof(double));
   const int incPFA = (1 > incPFA0) ? 1 : incPFA0;
   #define incAk 9
   const int incAm = 0, incAn = -(((Mb) << 3)+Mb);
   #define incBk 9
   const int incBm = -9, incBn = 9;
   #define incCm 2
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0;
   register double rC0_0;
   do /* N-loop */
   {
      do /* M-loop */
      {
         ATL_pfl1R(pfA);
         pfA += incPFA;
         rA0 = beta;
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
         pA0 += incAk;
         pB0 += incBk;
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
void ATL_zJIK0x0x9TN9x9x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=9, 
 * lda=9, ldb=9, ldc=0, mu=10, nu=1, ku=9, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/10)*10;
   #define Nb N
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (((Mb) << 3)+Mb);
   const double *stN = B + (((Nb) << 3)+Nb);
   const double *pfA = stM;
   const int incPFA0 = (((int)(stM - A))*10*1)/(M*N*sizeof(double));
   const int incPFA = (1 > incPFA0) ? 1 : incPFA0;
   #define incAk 9
   const int incAm = 81, incAn = -(((Mb) << 3)+Mb);
   #define incBk 9
   const int incBm = -9, incBn = 9;
   #define incCm 20
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9;
   register double rB0;
   register double rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0;
   if (pA0 != stM)
   {
      do /* N-loop */
      {
         do /* M-loop */
         {
            ATL_pfl1R(pfA);
            pfA += incPFA;
            rA0 = beta;
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
            rA0 = *pA0;
            rB0 = *pB0;
            rC0_0 += rA0 * rB0;
            rA1 = pA0[9];
            rA2 = pA0[18];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[27];
            rA4 = pA0[36];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[45];
            rA6 = pA0[54];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[63];
            rA8 = pA0[72];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[81];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[10];
            rA2 = pA0[19];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[28];
            rA4 = pA0[37];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[46];
            rA6 = pA0[55];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[64];
            rA8 = pA0[73];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[82];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[11];
            rA2 = pA0[20];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[29];
            rA4 = pA0[38];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[47];
            rA6 = pA0[56];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[65];
            rA8 = pA0[74];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[83];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[12];
            rA2 = pA0[21];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[30];
            rA4 = pA0[39];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[48];
            rA6 = pA0[57];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[66];
            rA8 = pA0[75];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[84];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[13];
            rA2 = pA0[22];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[31];
            rA4 = pA0[40];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[49];
            rA6 = pA0[58];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[67];
            rA8 = pA0[76];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[85];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[14];
            rA2 = pA0[23];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[32];
            rA4 = pA0[41];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[50];
            rA6 = pA0[59];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[68];
            rA8 = pA0[77];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[86];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[15];
            rA2 = pA0[24];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[33];
            rA4 = pA0[42];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[51];
            rA6 = pA0[60];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[69];
            rA8 = pA0[78];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[87];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[16];
            rA2 = pA0[25];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[34];
            rA4 = pA0[43];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[52];
            rA6 = pA0[61];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[70];
            rA8 = pA0[79];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[88];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[17];
            rA2 = pA0[26];
            rC1_0 += rA1 * rB0;
            rA3 = pA0[35];
            rA4 = pA0[44];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[53];
            rA6 = pA0[62];
            rC3_0 += rA3 * rB0;
            rA7 = pA0[71];
            rA8 = pA0[80];
            rC4_0 += rA4 * rB0;
            rA9 = pA0[89];
            rC5_0 += rA5 * rB0;
            rC6_0 += rA6 * rB0;
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            pA0 += incAk;
            pB0 += incBk;
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
   if (k=M-Mb)
      ATL_zJIK0x0x9TN1x1x9_a1_bX(k, N, 9, alpha, ca + (((Mb) << 3)+Mb), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
