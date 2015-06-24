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

static void ATL_zJIK0x0x27TN1x1x27_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=27, 
 * lda=27, ldb=27, ldc=0, mu=1, nu=1, ku=27, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   #define Nb N
   const double *stM = A + (27*(Mb));
   const double *stN = B + (27*(Nb));
   const double *pfA = stM;
   const int incPFA0 = (((int)(stM - A))*1*1)/(M*N*sizeof(double));
   const int incPFA = (1 > incPFA0) ? 1 : incPFA0;
   #define incAk 27
   const int incAm = 0, incAn = -(27*(Mb));
   #define incBk 27
   const int incBm = -27, incBn = 27;
   const int incAk0 = ((incAk) / 27), incBk0 = ((incBk) / 27);
   #define incCm 2
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0;
   register double m0, m1, m2;
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
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         m0 = rA0 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         m1 = rA0 * rB0;
         rA0 = pA0[2];
         rB0 = pB0[2];
         m2 = rA0 * rB0;
         rA0 = pA0[3];
         rB0 = pB0[3];

/*
 *       Completely unrolled K-loop
 */
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rC0_0 += m1;
         m1 = rA0 * rB0;
         rA0 = pA0[26];
         rB0 = pB0[26];
/*
 *       Drain pipe on last iteration of K-loop
 */
         rC0_0 += m2;
         m2 = rA0 * rB0;
         rC0_0 += m0;
         rC0_0 += m1;
         rC0_0 += m2;
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
void ATL_zJIK0x0x27TN27x27x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=27, 
 * lda=27, ldb=27, ldc=0, mu=10, nu=1, ku=27, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/10)*10;
   #define Nb N
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (27*(Mb));
   const double *stN = B + (27*(Nb));
   const double *pfA = stM;
   const int incPFA0 = (((int)(stM - A))*10*1)/(M*N*sizeof(double));
   const int incPFA = (1 > incPFA0) ? 1 : incPFA0;
   #define incAk 27
   const int incAm = 243, incAn = -(27*(Mb));
   #define incBk 27
   const int incBm = -27, incBn = 27;
   const int incAk0 = ((incAk) / 27), incBk0 = ((incBk) / 27);
   #define incCm 20
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9;
   register double rB0;
   register double m0, m1, m2, m3, m4;
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
/*
 *          Start pipeline
 */
            rA0 = *pA0;
            rB0 = *pB0;
            m0 = rA0 * rB0;
            rA1 = pA0[27];
            rA2 = pA0[54];
            m1 = rA1 * rB0;
            rA3 = pA0[81];
            rA4 = pA0[108];
            m2 = rA2 * rB0;
            rA5 = pA0[135];
            rA6 = pA0[162];
            m3 = rA3 * rB0;
            rA7 = pA0[189];
            rA8 = pA0[216];
            m4 = rA4 * rB0;
            rA9 = pA0[243];

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[28];
            rA2 = pA0[55];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[82];
            rA4 = pA0[109];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[136];
            rA6 = pA0[163];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[190];
            rA8 = pA0[217];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[244];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[29];
            rA2 = pA0[56];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[83];
            rA4 = pA0[110];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[137];
            rA6 = pA0[164];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[191];
            rA8 = pA0[218];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[245];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[30];
            rA2 = pA0[57];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[84];
            rA4 = pA0[111];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[138];
            rA6 = pA0[165];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[192];
            rA8 = pA0[219];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[246];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[31];
            rA2 = pA0[58];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[85];
            rA4 = pA0[112];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[139];
            rA6 = pA0[166];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[193];
            rA8 = pA0[220];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[247];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[32];
            rA2 = pA0[59];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[86];
            rA4 = pA0[113];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[140];
            rA6 = pA0[167];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[194];
            rA8 = pA0[221];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[248];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[33];
            rA2 = pA0[60];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[87];
            rA4 = pA0[114];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[141];
            rA6 = pA0[168];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[195];
            rA8 = pA0[222];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[249];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[34];
            rA2 = pA0[61];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[88];
            rA4 = pA0[115];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[142];
            rA6 = pA0[169];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[196];
            rA8 = pA0[223];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[250];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA0[35];
            rA2 = pA0[62];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[89];
            rA4 = pA0[116];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[143];
            rA6 = pA0[170];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[197];
            rA8 = pA0[224];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[251];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA0[36];
            rA2 = pA0[63];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[90];
            rA4 = pA0[117];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[144];
            rA6 = pA0[171];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[198];
            rA8 = pA0[225];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[252];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA0[37];
            rA2 = pA0[64];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[91];
            rA4 = pA0[118];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[145];
            rA6 = pA0[172];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[199];
            rA8 = pA0[226];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[253];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA0[38];
            rA2 = pA0[65];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[92];
            rA4 = pA0[119];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[146];
            rA6 = pA0[173];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[200];
            rA8 = pA0[227];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[254];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA0[39];
            rA2 = pA0[66];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[93];
            rA4 = pA0[120];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[147];
            rA6 = pA0[174];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[201];
            rA8 = pA0[228];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[255];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA0[40];
            rA2 = pA0[67];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[94];
            rA4 = pA0[121];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[148];
            rA6 = pA0[175];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[202];
            rA8 = pA0[229];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[256];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA0[41];
            rA2 = pA0[68];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[95];
            rA4 = pA0[122];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[149];
            rA6 = pA0[176];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[203];
            rA8 = pA0[230];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[257];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA0[42];
            rA2 = pA0[69];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[96];
            rA4 = pA0[123];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[150];
            rA6 = pA0[177];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[204];
            rA8 = pA0[231];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[258];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA0[43];
            rA2 = pA0[70];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[97];
            rA4 = pA0[124];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[151];
            rA6 = pA0[178];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[205];
            rA8 = pA0[232];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[259];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA0[44];
            rA2 = pA0[71];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[98];
            rA4 = pA0[125];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[152];
            rA6 = pA0[179];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[206];
            rA8 = pA0[233];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[260];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA0[45];
            rA2 = pA0[72];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[99];
            rA4 = pA0[126];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[153];
            rA6 = pA0[180];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[207];
            rA8 = pA0[234];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[261];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rA1 = pA0[46];
            rA2 = pA0[73];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[100];
            rA4 = pA0[127];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[154];
            rA6 = pA0[181];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[208];
            rA8 = pA0[235];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[262];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rA1 = pA0[47];
            rA2 = pA0[74];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[101];
            rA4 = pA0[128];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[155];
            rA6 = pA0[182];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[209];
            rA8 = pA0[236];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[263];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rA1 = pA0[48];
            rA2 = pA0[75];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[102];
            rA4 = pA0[129];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[156];
            rA6 = pA0[183];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[210];
            rA8 = pA0[237];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[264];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rA1 = pA0[49];
            rA2 = pA0[76];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[103];
            rA4 = pA0[130];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[157];
            rA6 = pA0[184];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[211];
            rA8 = pA0[238];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[265];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rA1 = pA0[50];
            rA2 = pA0[77];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[104];
            rA4 = pA0[131];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[158];
            rA6 = pA0[185];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[212];
            rA8 = pA0[239];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[266];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rA1 = pA0[51];
            rA2 = pA0[78];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[105];
            rA4 = pA0[132];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[159];
            rA6 = pA0[186];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[213];
            rA8 = pA0[240];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[267];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rA1 = pA0[52];
            rA2 = pA0[79];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[106];
            rA4 = pA0[133];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[160];
            rA6 = pA0[187];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[214];
            rA8 = pA0[241];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[268];
            rC9_0 += m4;
            m4 = rA4 * rB0;
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rA1 = pA0[53];
            rA2 = pA0[80];
            rC5_0 += m0;
            m0 = rA0 * rB0;
            rA3 = pA0[107];
            rA4 = pA0[134];
            rC6_0 += m1;
            m1 = rA1 * rB0;
            rA5 = pA0[161];
            rA6 = pA0[188];
            rC7_0 += m2;
            m2 = rA2 * rB0;
            rA7 = pA0[215];
            rA8 = pA0[242];
            rC8_0 += m3;
            m3 = rA3 * rB0;
            rA9 = pA0[269];
            rC9_0 += m4;
            m4 = rA4 * rB0;
/*
 *          Drain pipe on last iteration of K-loop
 */
            rC0_0 += m0;
            m0 = rA5 * rB0;
            rC1_0 += m1;
            m1 = rA6 * rB0;
            rC2_0 += m2;
            m2 = rA7 * rB0;
            rC3_0 += m3;
            m3 = rA8 * rB0;
            rC4_0 += m4;
            m4 = rA9 * rB0;
            rC5_0 += m0;
            rC6_0 += m1;
            rC7_0 += m2;
            rC8_0 += m3;
            rC9_0 += m4;
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
      ATL_zJIK0x0x27TN1x1x27_a1_bX(k, N, 27, alpha, ca + (27*(Mb)), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
