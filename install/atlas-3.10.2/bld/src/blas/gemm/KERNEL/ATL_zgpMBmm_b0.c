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
static void ATL_zJIK0x56x56TN1x1x2_a1_b0
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=56, KB=56, 
 * lda=56, ldb=56, ldc=0, mu=1, nu=1, ku=2, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   #define Mb M
   const double *stM = A + (56*(Mb));
   const double *stN = B + 3136;
   #define incAk 2
   const int incAm = 0, incAn = -(56*(Mb));
   #define incBk 2
   const int incBm = -56, incBn = 56;
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
/*
 *       Feeble prefetch of C
 */
         rC0_0 = *pC0;
/*
 *       Peel 1st iter to assign C regs
 */
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 = rA0 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rC0_0 += rA0 * rB0;
         pA0 += incAk;
         pB0 += incBk;
/*
 *       Unpeeled K iterations
 */
         for (k=0; k < 27; k++) /* easy loop to unroll */
         {
            rA0 = *pA0;
            rB0 = *pB0;
            rC0_0 += rA0 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rC0_0 += rA0 * rB0;
            pA0 += incAk;
            pB0 += incBk;
         }
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
void ATL_zJIK0x56x56TN56x56x0_a1_b0
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=56, KB=56, 
 * lda=56, ldb=56, ldc=0, mu=12, nu=1, ku=2, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.10.2)
 */
{
   const int Mb = (M/12)*12;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (56*(Mb));
   const double *stN = B + 3136;
   #define incAk 2
   const int incAm = 616, incAn = -(56*(Mb));
   #define incBk 2
   const int incBm = -56, incBn = 56;
   #define incCm 24
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1, rA2, rA3, rA4, rA5, rA6, rA7, rA8, rA9, rA10, rA11;
   register double rB0;
   register double rC0_0, rC1_0, rC2_0, rC3_0, rC4_0, rC5_0, rC6_0, rC7_0, rC8_0, rC9_0, rC10_0, rC11_0;
   if (pA0 != stM)
   {
      do /* N-loop */
      {
         do /* M-loop */
         {
/*
 *          Feeble prefetch of C
 */
            rC0_0 = *pC0;
/*
 *          Peel 1st iter to assign C regs
 */
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[56];
            rA2 = pA0[112];
            rA3 = pA0[168];
            rA4 = pA0[224];
            rA5 = pA0[280];
            rA6 = pA0[336];
            rC0_0 = rA0 * rB0;
            rA7 = pA0[392];
            rA8 = pA0[448];
            rC1_0 = rA1 * rB0;
            rA9 = pA0[504];
            rA10 = pA0[560];
            rC2_0 = rA2 * rB0;
            rA11 = pA0[616];
            rC3_0 = rA3 * rB0;
            rC4_0 = rA4 * rB0;
            rC5_0 = rA5 * rB0;
            rC6_0 = rA6 * rB0;
            rC7_0 = rA7 * rB0;
            rC8_0 = rA8 * rB0;
            rC9_0 = rA9 * rB0;
            rC10_0 = rA10 * rB0;
            rC11_0 = rA11 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[57];
            rA2 = pA0[113];
            rA3 = pA0[169];
            rA4 = pA0[225];
            rA5 = pA0[281];
            rA6 = pA0[337];
            rC0_0 += rA0 * rB0;
            rA7 = pA0[393];
            rA8 = pA0[449];
            rC1_0 += rA1 * rB0;
            rA9 = pA0[505];
            rA10 = pA0[561];
            rC2_0 += rA2 * rB0;
            rA11 = pA0[617];
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
/*
 *          Unpeeled K iterations
 */
            for (k=0; k < 27; k++) /* easy loop to unroll */
            {
               rA0 = *pA0;
               rB0 = *pB0;
               rA1 = pA0[56];
               rA2 = pA0[112];
               rA3 = pA0[168];
               rA4 = pA0[224];
               rA5 = pA0[280];
               rA6 = pA0[336];
               rC0_0 += rA0 * rB0;
               rA7 = pA0[392];
               rA8 = pA0[448];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[504];
               rA10 = pA0[560];
               rC2_0 += rA2 * rB0;
               rA11 = pA0[616];
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
               rA1 = pA0[57];
               rA2 = pA0[113];
               rA3 = pA0[169];
               rA4 = pA0[225];
               rA5 = pA0[281];
               rA6 = pA0[337];
               rC0_0 += rA0 * rB0;
               rA7 = pA0[393];
               rA8 = pA0[449];
               rC1_0 += rA1 * rB0;
               rA9 = pA0[505];
               rA10 = pA0[561];
               rC2_0 += rA2 * rB0;
               rA11 = pA0[617];
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
      ATL_zJIK0x56x56TN1x1x2_a1_b0(k, 56, 56, alpha, ca + (56*(Mb)), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
#ifdef ATL_UCLEANM
#define ATL_zpMBmm_b0 ATL_zgpMBmm_b0
#endif

void ATL_zJIK0x56x56TN56x56x0_a1_b0(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x56x56TN56x56x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x56x56TN56x56x0_a1_b1(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);

void ATL_zpMBmm_b0(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
   ATL_zJIK0x56x56TN56x56x0_a1_b0(M, N, K, alpha, A, lda, B, ldb, -beta, C, ldc);
   ATL_zJIK0x56x56TN56x56x0_a1_b0(M, N, K, alpha, A, lda, B+N*ldb, ldb, beta, C+1, ldc);
   ATL_zJIK0x56x56TN56x56x0_a1_bX(M, N, K, alpha, A+M*lda, lda, B+N*ldb, ldb, -1.0, C, ldc);
   ATL_zJIK0x56x56TN56x56x0_a1_b1(M, N, K, alpha, A+M*lda, lda, B, ldb, 1.0, C+1, ldc);
}
