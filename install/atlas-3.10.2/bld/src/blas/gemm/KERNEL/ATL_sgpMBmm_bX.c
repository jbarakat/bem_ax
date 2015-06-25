#ifdef ATL_UCLEANM
#define ATL_sJIK0x80x80TN80x80x0_a1_bX ATL_sgpMBmm_bX
#else
#define ATL_sJIK0x80x80TN80x80x0_a1_bX ATL_spMBmm_bX
#endif

#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
static void ATL_sJIK0x80x80TN1x1x6_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=80, KB=80, 
 * lda=80, ldb=80, ldc=0, mu=1, nu=1, ku=6, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   const float *stM = A + (((Mb) << 6)+((Mb) << 4));
   const float *stN = B + 6400;
   #define incAk 6
   const int incAm = 2, incAn = -(((Mb) << 6)+((Mb) << 4));
   #define incBk 6
   const int incBm = -78, incBn = 80;
   #define incCm 1
   const int incCn = (ldc) - (Mb);
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
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         for (k=0; k < 13; k++) /* easy loop to unroll */
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
            pA0 += incAk;
            pB0 += incBk;
         }
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
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
void ATL_sJIK0x80x80TN80x80x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=80, KB=80, 
 * lda=80, ldb=80, ldc=0, mu=10, nu=1, ku=6, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/10)*10;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 6)+((Mb) << 4));
   const float *stN = B + 6400;
   #define incAk 6
   const int incAm = 722, incAn = -(((Mb) << 6)+((Mb) << 4));
   #define incBk 6
   const int incBm = -78, incBn = 80;
   #define incCm 10
   const int incCn = (ldc) - (Mb);
   float *pC0=C;
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9;
   register float rB0;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0;
   if (pA0 != stM)
   {
      do /* N-loop */
      {
         do /* M-loop */
         {
            rA0 = beta;
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
            for (k=0; k < 13; k++) /* easy loop to unroll */
            {
               rA0 = *pA0;
               rB0 = *pB0;
               rA1 = pA0[80];
               rA2 = pA0[160];
               rC0_0 += rA0 * rB0;
               rA3 = pA0[240];
               rC1_0 += rA1 * rB0;
               rA4 = pA0[320];
               rC2_0 += rA2 * rB0;
               rA5 = pA0[400];
               rC3_0 += rA3 * rB0;
               rA6 = pA0[480];
               rC4_0 += rA4 * rB0;
               rA7 = pA0[560];
               rC5_0 += rA5 * rB0;
               rA8 = pA0[640];
               rC6_0 += rA6 * rB0;
               rA9 = pA0[720];
               rC7_0 += rA7 * rB0;
               rC8_0 += rA8 * rB0;
               rC9_0 += rA9 * rB0;
               rA0 = pA0[1];
               rB0 = pB0[1];
               rA1 = pA0[81];
               rA2 = pA0[161];
               rC0_0 += rA0 * rB0;
               rA3 = pA0[241];
               rC1_0 += rA1 * rB0;
               rA4 = pA0[321];
               rC2_0 += rA2 * rB0;
               rA5 = pA0[401];
               rC3_0 += rA3 * rB0;
               rA6 = pA0[481];
               rC4_0 += rA4 * rB0;
               rA7 = pA0[561];
               rC5_0 += rA5 * rB0;
               rA8 = pA0[641];
               rC6_0 += rA6 * rB0;
               rA9 = pA0[721];
               rC7_0 += rA7 * rB0;
               rC8_0 += rA8 * rB0;
               rC9_0 += rA9 * rB0;
               rA0 = pA0[2];
               rB0 = pB0[2];
               rA1 = pA0[82];
               rA2 = pA0[162];
               rC0_0 += rA0 * rB0;
               rA3 = pA0[242];
               rC1_0 += rA1 * rB0;
               rA4 = pA0[322];
               rC2_0 += rA2 * rB0;
               rA5 = pA0[402];
               rC3_0 += rA3 * rB0;
               rA6 = pA0[482];
               rC4_0 += rA4 * rB0;
               rA7 = pA0[562];
               rC5_0 += rA5 * rB0;
               rA8 = pA0[642];
               rC6_0 += rA6 * rB0;
               rA9 = pA0[722];
               rC7_0 += rA7 * rB0;
               rC8_0 += rA8 * rB0;
               rC9_0 += rA9 * rB0;
               rA0 = pA0[3];
               rB0 = pB0[3];
               rA1 = pA0[83];
               rA2 = pA0[163];
               rC0_0 += rA0 * rB0;
               rA3 = pA0[243];
               rC1_0 += rA1 * rB0;
               rA4 = pA0[323];
               rC2_0 += rA2 * rB0;
               rA5 = pA0[403];
               rC3_0 += rA3 * rB0;
               rA6 = pA0[483];
               rC4_0 += rA4 * rB0;
               rA7 = pA0[563];
               rC5_0 += rA5 * rB0;
               rA8 = pA0[643];
               rC6_0 += rA6 * rB0;
               rA9 = pA0[723];
               rC7_0 += rA7 * rB0;
               rC8_0 += rA8 * rB0;
               rC9_0 += rA9 * rB0;
               rA0 = pA0[4];
               rB0 = pB0[4];
               rA1 = pA0[84];
               rA2 = pA0[164];
               rC0_0 += rA0 * rB0;
               rA3 = pA0[244];
               rC1_0 += rA1 * rB0;
               rA4 = pA0[324];
               rC2_0 += rA2 * rB0;
               rA5 = pA0[404];
               rC3_0 += rA3 * rB0;
               rA6 = pA0[484];
               rC4_0 += rA4 * rB0;
               rA7 = pA0[564];
               rC5_0 += rA5 * rB0;
               rA8 = pA0[644];
               rC6_0 += rA6 * rB0;
               rA9 = pA0[724];
               rC7_0 += rA7 * rB0;
               rC8_0 += rA8 * rB0;
               rC9_0 += rA9 * rB0;
               rA0 = pA0[5];
               rB0 = pB0[5];
               rA1 = pA0[85];
               rA2 = pA0[165];
               rC0_0 += rA0 * rB0;
               rA3 = pA0[245];
               rC1_0 += rA1 * rB0;
               rA4 = pA0[325];
               rC2_0 += rA2 * rB0;
               rA5 = pA0[405];
               rC3_0 += rA3 * rB0;
               rA6 = pA0[485];
               rC4_0 += rA4 * rB0;
               rA7 = pA0[565];
               rC5_0 += rA5 * rB0;
               rA8 = pA0[645];
               rC6_0 += rA6 * rB0;
               rA9 = pA0[725];
               rC7_0 += rA7 * rB0;
               rC8_0 += rA8 * rB0;
               rC9_0 += rA9 * rB0;
               pA0 += incAk;
               pB0 += incBk;
            }
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[80];
            rA2 = pA0[160];
            rC0_0 += rA0 * rB0;
            rA3 = pA0[240];
            rC1_0 += rA1 * rB0;
            rA4 = pA0[320];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[400];
            rC3_0 += rA3 * rB0;
            rA6 = pA0[480];
            rC4_0 += rA4 * rB0;
            rA7 = pA0[560];
            rC5_0 += rA5 * rB0;
            rA8 = pA0[640];
            rC6_0 += rA6 * rB0;
            rA9 = pA0[720];
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[81];
            rA2 = pA0[161];
            rC0_0 += rA0 * rB0;
            rA3 = pA0[241];
            rC1_0 += rA1 * rB0;
            rA4 = pA0[321];
            rC2_0 += rA2 * rB0;
            rA5 = pA0[401];
            rC3_0 += rA3 * rB0;
            rA6 = pA0[481];
            rC4_0 += rA4 * rB0;
            rA7 = pA0[561];
            rC5_0 += rA5 * rB0;
            rA8 = pA0[641];
            rC6_0 += rA6 * rB0;
            rA9 = pA0[721];
            rC7_0 += rA7 * rB0;
            rC8_0 += rA8 * rB0;
            rC9_0 += rA9 * rB0;
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
      ATL_sJIK0x80x80TN1x1x6_a1_bX(k, 80, 80, alpha, ca + (((Mb) << 6)+((Mb) << 4)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
