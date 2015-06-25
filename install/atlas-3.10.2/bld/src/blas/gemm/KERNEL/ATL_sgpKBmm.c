#ifndef ATL_UCLEANK
   #define ATL_sgpKBmm ATL_spKBmm
#endif

void ATL_sJIK0x0x75TN75x75x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
void ATL_sJIK0x0x76TN76x76x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
void ATL_sJIK0x0x77TN77x77x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
void ATL_sJIK0x0x78TN78x78x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
void ATL_sJIK0x0x79TN79x79x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
void ATL_sJIK0x0x80TN80x80x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
void ATL_sJIK0x0x0TN0x0x0_a1_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);
typedef void (*MMfunc)(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc);

void ATL_sgpKBmm(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)
{
   static MMfunc mmfunc[  6] = {
                         ATL_sJIK0x0x75TN75x75x0_a1_bX,
                         ATL_sJIK0x0x76TN76x76x0_a1_bX,
                         ATL_sJIK0x0x77TN77x77x0_a1_bX,
                         ATL_sJIK0x0x78TN78x78x0_a1_bX,
                         ATL_sJIK0x0x79TN79x79x0_a1_bX,
                         ATL_sJIK0x0x80TN80x80x0_a1_bX,
                        };

   if (K <= 74) ATL_sJIK0x0x0TN0x0x0_a1_bX(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   else mmfunc[K-75](M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
#ifndef ATL_UCLEANK
void ATL_spKBmm_b0(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)
{
   ATL_sgpKBmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
void ATL_spKBmm_b1(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)
{
   ATL_sgpKBmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
void ATL_spKBmm_bX(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)
{
   ATL_sgpKBmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
#endif
