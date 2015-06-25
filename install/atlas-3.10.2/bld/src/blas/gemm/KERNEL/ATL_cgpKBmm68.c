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

static void ATL_cJIK0x0x68TN1x1x8_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=68, 
 * lda=68, ldb=68, ldc=0, mu=1, nu=1, ku=8, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   #define Nb N
   const float *stM = A + (((Mb) << 6)+((Mb) << 2));
   const float *stN = B + (((Nb) << 6)+((Nb) << 2));
   const float *pfA = stM;
   const int incPFA0 = (((int)(stM - A))*1*1)/(M*N*sizeof(float));
   const int incPFA = (1 > incPFA0) ? 1 : incPFA0;
   #define incAk 8
   const int incAm = 4, incAn = -(((Mb) << 6)+((Mb) << 2));
   #define incBk 8
   const int incBm = -64, incBn = 68;
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
      do /* M-loop */
      {
         ATL_pfl1R(pfA);
         pfA += incPFA;
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         for (k=0; k < 8; k++) /* easy loop to unroll */
         {
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
            pA0 += incAk;
            pB0 += incBk;
         }
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
void ATL_cJIK0x0x68TN68x68x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=68, 
 * lda=68, ldb=68, ldc=0, mu=12, nu=1, ku=8, pf=1
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/12)*12;
   #define Nb N
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 6)+((Mb) << 2));
   const float *stN = B + (((Nb) << 6)+((Nb) << 2));
   const float *pfA = stM;
   const int incPFA0 = (((int)(stM - A))*12*1)/(M*N*sizeof(float));
   const int incPFA = (1 > incPFA0) ? 1 : incPFA0;
   #define incAk 8
   const int incAm = 752, incAn = -(((Mb) << 6)+((Mb) << 2));
   #define incBk 8
   const int incBm = -64, incBn = 68;
   #define incCm 24
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   float *pC0=C;
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9, rA10, rA11;
   register float rB0;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0, rC10_0, rC11_0;
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
            rC10_0 = pC0[20];
            rC10_0 *= rA0;
            rC11_0 = pC0[22];
            rC11_0 *= rA0;
            for (k=0; k < 8; k++) /* easy loop to unroll */
            {
               rA0 = *pA0;
               rB0 = *pB0;
               rA1 = pA0[68];
               rA2 = pA0[136];
               rA3 = pA0[204];
               rA4 = pA0[272];
               rA5 = pA0[340];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[408];
               rA7 = pA0[476];
               rA8 = pA0[544];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[612];
               rA10 = pA0[680];
               rA11 = pA0[748];
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
               rA0 = pA0[1];
               rB0 = pB0[1];
               rA1 = pA0[69];
               rA2 = pA0[137];
               rA3 = pA0[205];
               rA4 = pA0[273];
               rA5 = pA0[341];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[409];
               rA7 = pA0[477];
               rA8 = pA0[545];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[613];
               rA10 = pA0[681];
               rA11 = pA0[749];
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
               rA0 = pA0[2];
               rB0 = pB0[2];
               rA1 = pA0[70];
               rA2 = pA0[138];
               rA3 = pA0[206];
               rA4 = pA0[274];
               rA5 = pA0[342];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[410];
               rA7 = pA0[478];
               rA8 = pA0[546];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[614];
               rA10 = pA0[682];
               rA11 = pA0[750];
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
               rA0 = pA0[3];
               rB0 = pB0[3];
               rA1 = pA0[71];
               rA2 = pA0[139];
               rA3 = pA0[207];
               rA4 = pA0[275];
               rA5 = pA0[343];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[411];
               rA7 = pA0[479];
               rA8 = pA0[547];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[615];
               rA10 = pA0[683];
               rA11 = pA0[751];
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
               rA0 = pA0[4];
               rB0 = pB0[4];
               rA1 = pA0[72];
               rA2 = pA0[140];
               rA3 = pA0[208];
               rA4 = pA0[276];
               rA5 = pA0[344];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[412];
               rA7 = pA0[480];
               rA8 = pA0[548];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[616];
               rA10 = pA0[684];
               rA11 = pA0[752];
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
               rA0 = pA0[5];
               rB0 = pB0[5];
               rA1 = pA0[73];
               rA2 = pA0[141];
               rA3 = pA0[209];
               rA4 = pA0[277];
               rA5 = pA0[345];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[413];
               rA7 = pA0[481];
               rA8 = pA0[549];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[617];
               rA10 = pA0[685];
               rA11 = pA0[753];
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
               rA0 = pA0[6];
               rB0 = pB0[6];
               rA1 = pA0[74];
               rA2 = pA0[142];
               rA3 = pA0[210];
               rA4 = pA0[278];
               rA5 = pA0[346];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[414];
               rA7 = pA0[482];
               rA8 = pA0[550];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[618];
               rA10 = pA0[686];
               rA11 = pA0[754];
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
               rA0 = pA0[7];
               rB0 = pB0[7];
               rA1 = pA0[75];
               rA2 = pA0[143];
               rA3 = pA0[211];
               rA4 = pA0[279];
               rA5 = pA0[347];
               rC0_0 += rA0 * rB0;
               rA6 = pA0[415];
               rA7 = pA0[483];
               rA8 = pA0[551];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[619];
               rA10 = pA0[687];
               rA11 = pA0[755];
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
               pA0 += incAk;
               pB0 += incBk;
            }
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[68];
            rA2 = pA0[136];
            rA3 = pA0[204];
            rA4 = pA0[272];
            rA5 = pA0[340];
            rC0_0 += rA0 * rB0;
            rA6 = pA0[408];
            rA7 = pA0[476];
            rA8 = pA0[544];
            rC1_0 += rA1 * rB0;
            rA9 = pA0[612];
            rA10 = pA0[680];
            rA11 = pA0[748];
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
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[69];
            rA2 = pA0[137];
            rA3 = pA0[205];
            rA4 = pA0[273];
            rA5 = pA0[341];
            rC0_0 += rA0 * rB0;
            rA6 = pA0[409];
            rA7 = pA0[477];
            rA8 = pA0[545];
            rC1_0 += rA1 * rB0;
            rA9 = pA0[613];
            rA10 = pA0[681];
            rA11 = pA0[749];
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
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[70];
            rA2 = pA0[138];
            rA3 = pA0[206];
            rA4 = pA0[274];
            rA5 = pA0[342];
            rC0_0 += rA0 * rB0;
            rA6 = pA0[410];
            rA7 = pA0[478];
            rA8 = pA0[546];
            rC1_0 += rA1 * rB0;
            rA9 = pA0[614];
            rA10 = pA0[682];
            rA11 = pA0[750];
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
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[71];
            rA2 = pA0[139];
            rA3 = pA0[207];
            rA4 = pA0[275];
            rA5 = pA0[343];
            rC0_0 += rA0 * rB0;
            rA6 = pA0[411];
            rA7 = pA0[479];
            rA8 = pA0[547];
            rC1_0 += rA1 * rB0;
            rA9 = pA0[615];
            rA10 = pA0[683];
            rA11 = pA0[751];
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
      ATL_cJIK0x0x68TN1x1x8_a1_bX(k, N, 68, alpha, ca + (((Mb) << 6)+((Mb) << 2)), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
