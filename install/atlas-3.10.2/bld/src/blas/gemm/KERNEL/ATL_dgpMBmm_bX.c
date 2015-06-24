#ifdef ATL_UCLEANM
#define ATL_dJIK0x52x52TN52x52x0_a1_bX ATL_dgpMBmm_bX
#else
#define ATL_dJIK0x52x52TN52x52x0_a1_bX ATL_dpMBmm_bX
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
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
static void ATL_dJIK0x4x52TN1x4x1_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=4, KB=52, 
 * lda=52, ldb=52, ldc=0, mu=1, nu=4, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   const double *stM = A + (52*(Mb));
   const double *stN = B + 208;
   #define incAk 1
   const int incAm = 0, incAn = -(52*(Mb));
   #define incBk 1
   const int incBm = -52, incBn = 208;
   #define incCm 1
   const int incCn = (((ldc) << 2)) - (Mb);
   double *pC0=C, *pC1=pC0+(ldc), *pC2=pC1+(ldc), *pC3=pC2+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0, rB1, rB2, rB3;
   register double rC0_0, rC0_1, rC0_2, rC0_3;
   do /* M-loop */
   {
      rA0 = beta;
      rC0_0 = *pC0;
      rC0_0 *= rA0;
      rC0_1 = *pC1;
      rC0_1 *= rA0;
      rC0_2 = *pC2;
      rC0_2 *= rA0;
      rC0_3 = *pC3;
      rC0_3 *= rA0;
      for (k=0; k < 52; k++) /* easy loop to unroll */
      {
         rA0 = *pA0;
         rB0 = *pB0;
         rB1 = pB0[52];
         rB2 = pB0[104];
         rB3 = pB0[156];
         rC0_0 += rA0 * rB0;
         rC0_1 += rA0 * rB1;
         rC0_2 += rA0 * rB2;
         rC0_3 += rA0 * rB3;
         pA0 += incAk;
         pB0 += incBk;
      }
      *pC0 = rC0_0;
      *pC1 = rC0_1;
      *pC2 = rC0_2;
      *pC3 = rC0_3;
      pC0 += incCm;
      pC1 += incCm;
      pC2 += incCm;
      pC3 += incCm;
      pA0 += incAm;
      pB0 += incBm;
   }
   while(pA0 != stM);
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
static void ATL_dJIK0x4x52TN2x4x1_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=4, KB=52, 
 * lda=52, ldb=52, ldc=0, mu=2, nu=4, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (52*(Mb));
   const double *stN = B + 208;
   #define incAk 1
   const int incAm = 52, incAn = -(52*(Mb));
   #define incBk 1
   const int incBm = -52, incBn = 208;
   #define incCm 2
   const int incCn = (((ldc) << 2)) - (Mb);
   double *pC0=C, *pC1=pC0+(ldc), *pC2=pC1+(ldc), *pC3=pC2+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1;
   register double rB0, rB1, rB2, rB3;
   register double rC0_0, rC1_0, rC0_1, rC1_1, rC0_2, rC1_2, rC0_3, rC1_3;
   if (pA0 != stM)
   {
      do /* M-loop */
      {
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[1];
         rC1_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC1_1 = pC1[1];
         rC1_1 *= rA0;
         rC0_2 = *pC2;
         rC0_2 *= rA0;
         rC1_2 = pC2[1];
         rC1_2 *= rA0;
         rC0_3 = *pC3;
         rC0_3 *= rA0;
         rC1_3 = pC3[1];
         rC1_3 *= rA0;
         for (k=0; k < 52; k++) /* easy loop to unroll */
         {
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[52];
            rB1 = pB0[52];
            rB2 = pB0[104];
            rC0_0 += rA0 * rB0;
            rB3 = pB0[156];
            rC1_0 += rA1 * rB0;
            rC0_1 += rA0 * rB1;
            rC1_1 += rA1 * rB1;
            rC0_2 += rA0 * rB2;
            rC1_2 += rA1 * rB2;
            rC0_3 += rA0 * rB3;
            rC1_3 += rA1 * rB3;
            pA0 += incAk;
            pB0 += incBk;
         }
         *pC0 = rC0_0;
         pC0[1] = rC1_0;
         *pC1 = rC0_1;
         pC1[1] = rC1_1;
         *pC2 = rC0_2;
         pC2[1] = rC1_2;
         *pC3 = rC0_3;
         pC3[1] = rC1_3;
         pC0 += incCm;
         pC1 += incCm;
         pC2 += incCm;
         pC3 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
   }
   if (k=M-Mb)
      ATL_dJIK0x4x52TN1x4x1_a1_bX(k, 4, 52, alpha, ca + (52*(Mb)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
static void ATL_dJIK0x48x52TN1x6x1_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=48, KB=52, 
 * lda=52, ldb=52, ldc=0, mu=1, nu=6, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   const double *stM = A + (52*(Mb));
   const double *stN = B + 2496;
   #define incAk 1
   const int incAm = 0, incAn = -(52*(Mb));
   #define incBk 1
   const int incBm = -52, incBn = 312;
   #define incCm 1
   const int incCn = (((ldc) << 2)+((ldc) << 1)) - (Mb);
   double *pC0=C, *pC1=pC0+(ldc), *pC2=pC1+(ldc), *pC3=pC2+(ldc), *pC4=pC3+(ldc), *pC5=pC4+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0, rB1, rB2, rB3, rB4, rB5;
   register double rC0_0, rC0_1, rC0_2, rC0_3, rC0_4, rC0_5;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC0_2 = *pC2;
         rC0_2 *= rA0;
         rC0_3 = *pC3;
         rC0_3 *= rA0;
         rC0_4 = *pC4;
         rC0_4 *= rA0;
         rC0_5 = *pC5;
         rC0_5 *= rA0;
         for (k=0; k < 52; k++) /* easy loop to unroll */
         {
            rA0 = *pA0;
            rB0 = *pB0;
            rB1 = pB0[52];
            rB2 = pB0[104];
            rB3 = pB0[156];
            rC0_0 += rA0 * rB0;
            rB4 = pB0[208];
            rC0_1 += rA0 * rB1;
            rB5 = pB0[260];
            rC0_2 += rA0 * rB2;
            rC0_3 += rA0 * rB3;
            rC0_4 += rA0 * rB4;
            rC0_5 += rA0 * rB5;
            pA0 += incAk;
            pB0 += incBk;
         }
         *pC0 = rC0_0;
         *pC1 = rC0_1;
         *pC2 = rC0_2;
         *pC3 = rC0_3;
         *pC4 = rC0_4;
         *pC5 = rC0_5;
         pC0 += incCm;
         pC1 += incCm;
         pC2 += incCm;
         pC3 += incCm;
         pC4 += incCm;
         pC5 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pC1 += incCn;
      pC2 += incCn;
      pC3 += incCn;
      pC4 += incCn;
      pC5 += incCn;
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
void ATL_dJIK0x52x52TN52x52x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=52, KB=52, 
 * lda=52, ldb=52, ldc=0, mu=2, nu=6, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (52*(Mb));
   const double *stN = B + 2496;
   #define incAk 1
   const int incAm = 52, incAn = -(52*(Mb));
   #define incBk 1
   const int incBm = -52, incBn = 312;
   #define incCm 2
   const int incCn = (((ldc) << 2)+((ldc) << 1)) - (Mb);
   double *pC0=C, *pC1=pC0+(ldc), *pC2=pC1+(ldc), *pC3=pC2+(ldc), *pC4=pC3+(ldc), *pC5=pC4+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1;
   register double rB0, rB1, rB2, rB3, rB4, rB5;
   register double rC0_0, rC1_0, rC0_1, rC1_1, rC0_2, rC1_2, rC0_3, rC1_3, rC0_4, rC1_4, rC0_5, rC1_5;
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
            rC0_1 = *pC1;
            rC0_1 *= rA0;
            rC1_1 = pC1[1];
            rC1_1 *= rA0;
            rC0_2 = *pC2;
            rC0_2 *= rA0;
            rC1_2 = pC2[1];
            rC1_2 *= rA0;
            rC0_3 = *pC3;
            rC0_3 *= rA0;
            rC1_3 = pC3[1];
            rC1_3 *= rA0;
            rC0_4 = *pC4;
            rC0_4 *= rA0;
            rC1_4 = pC4[1];
            rC1_4 *= rA0;
            rC0_5 = *pC5;
            rC0_5 *= rA0;
            rC1_5 = pC5[1];
            rC1_5 *= rA0;
            for (k=0; k < 52; k++) /* easy loop to unroll */
            {
               rA0 = *pA0;
               rB0 = *pB0;
               rA1 = pA0[52];
               rB1 = pB0[52];
               rB2 = pB0[104];
               rC0_0 += rA0 * rB0;
               rB3 = pB0[156];
               rC1_0 += rA1 * rB0;
               rB4 = pB0[208];
               rC0_1 += rA0 * rB1;
               rB5 = pB0[260];
               rC1_1 += rA1 * rB1;
               rC0_2 += rA0 * rB2;
               rC1_2 += rA1 * rB2;
               rC0_3 += rA0 * rB3;
               rC1_3 += rA1 * rB3;
               rC0_4 += rA0 * rB4;
               rC1_4 += rA1 * rB4;
               rC0_5 += rA0 * rB5;
               rC1_5 += rA1 * rB5;
               pA0 += incAk;
               pB0 += incBk;
            }
            *pC0 = rC0_0;
            pC0[1] = rC1_0;
            *pC1 = rC0_1;
            pC1[1] = rC1_1;
            *pC2 = rC0_2;
            pC2[1] = rC1_2;
            *pC3 = rC0_3;
            pC3[1] = rC1_3;
            *pC4 = rC0_4;
            pC4[1] = rC1_4;
            *pC5 = rC0_5;
            pC5[1] = rC1_5;
            pC0 += incCm;
            pC1 += incCm;
            pC2 += incCm;
            pC3 += incCm;
            pC4 += incCm;
            pC5 += incCm;
            pA0 += incAm;
            pB0 += incBm;
         }
         while(pA0 != stM);
         pC0 += incCn;
         pC1 += incCn;
         pC2 += incCn;
         pC3 += incCn;
         pC4 += incCn;
         pC5 += incCn;
         pA0 += incAn;
         pB0 += incBn;
      }
      while(pB0 != stN);
   }
   ATL_dJIK0x4x52TN2x4x1_a1_bX(M, 4, 52, alpha, ca, lda, cb + 2496, ldb, beta, cc + (((ldc) << 5)+((ldc) << 4)), ldc);
   if (k=M-Mb)
      ATL_dJIK0x48x52TN1x6x1_a1_bX(k, 48, 52, alpha, ca + (52*(Mb)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
