#define SREAL
#include "atlas_misc.h"
void ATL_supMBmm0_2_0_bX
   (const int M, const int N, const int K, const float alpha,
    const float *A, const int lda, const float *B, const int ldb,
    const float beta, float *C, const int ldc);
void ATL_supMBmm0_4_0_bX
   (const int M, const int N, const int K, const float alpha,
    const float *A, const int lda, const float *B, const int ldb,
    const float beta, float *C, const int ldc);
void ATL_sgpMBmm_bX
   (const int M, const int N, const int K, const float alpha,
    const float *A, const int lda, const float *B, const int ldb,
    const float beta, float *C, const int ldc);

void ATL_spMBmm_bX
   (const int M, const int N, const int K, const float alpha,
    const float *A, const int lda, const float *B, const int ldb,
    const float beta, float *C, const int ldc)
{

   if (M == (((((M) >> 2)) << 2)))
   {
      ATL_supMBmm0_4_0_bX(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   else if (M == (((((M) >> 1)) << 1)))
   {
      ATL_supMBmm0_2_0_bX(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   else ATL_sgpMBmm_bX(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
