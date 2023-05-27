// Xiaozhi Zhu 521021910299
/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

// help function
void trans_1x8bits(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0 = A[i][j], a1 = A[i][j + 1], a2 = A[i][j + 2], a3 = A[i][j + 3], a4 = A[i][j + 4], a5 = A[i][j + 5], a6 = A[i][j + 6], a7 = A[i][j + 7];
    B[j][i] = a0, B[j + 1][i] = a1, B[j + 2][i] = a2, B[j + 3][i] = a3, B[j + 4][i] = a4, B[j + 5][i] = a5, B[j + 6][i] = a6, B[j + 7][i] = a7;
}

void trans_1x4bits(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0 = A[i][j], a1 = A[i][j + 1], a2 = A[i][j + 2], a3 = A[i][j + 3];
    B[j][i] = a0, B[j + 1][i] = a1, B[j + 2][i] = a2, B[j + 3][i] = a3;
}

void trans_2x4bits(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0 = A[i][j], a1 = A[i][j + 1], a2 = A[i][j + 2], a3 = A[i][j + 3];
    int b0 = A[i + 1][j], b1 = A[i + 1][j + 1], b2 = A[i + 1][j + 2], b3 = A[i + 1][j + 3];
    B[j][i] = a0, B[j + 1][i] = a1, B[j + 2][i] = a2, B[j + 3][i] = a3;
    B[j][i + 1] = b0, B[j + 1][i + 1] = b1, B[j + 2][i + 1] = b2, B[j + 3][i + 1] = b3;
}

void trans_1x5bits(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0 = A[i][j], a1 = A[i][j + 1], a2 = A[i][j + 2], a3 = A[i][j + 3], a4 = A[i][j + 4];
    B[j][i] = a0, B[j + 1][i] = a1, B[j + 2][i] = a2, B[j + 3][i] = a3, B[j + 4][i] = a4;
}

void trans_1x3bits(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0 = A[i][j], a1 = A[i][j + 1], a2 = A[i][j + 2];
    B[j][i] = a0, B[j + 1][i] = a1, B[j + 2][i] = a2;
}

//////////////////////////////////////// helper function for v3 ////////////////////////////////////////

void copy_8x8bits(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0, a1, a2, a3, a4, a5, a6, a7;
    for (int k = 0; k < 8; k++)
    {
        a0 = A[i + k][j], a1 = A[i + k][j + 1], a2 = A[i + k][j + 2], a3 = A[i + k][j + 3], a4 = A[i + k][j + 4], a5 = A[i + k][j + 5], a6 = A[i + k][j + 6], a7 = A[i + k][j + 7];
        B[j + k][i] = a0, B[j + k][i + 1] = a1, B[j + k][i + 2] = a2, B[j + k][i + 3] = a3, B[j + k][i + 4] = a4, B[j + k][i + 5] = a5, B[j + k][i + 6] = a6, B[j + k][i + 7] = a7;
    }
}

void self_trans_4x4bits_identity(int i, int j, int M, int N, int B[M][N]) // i,j与转置前保持一致
{
    int temp;
    temp = B[j][i + 1], B[j][i + 1] = B[j + 1][i], B[j + 1][i] = temp;
    temp = B[j][i + 2], B[j][i + 2] = B[j + 2][i], B[j + 2][i] = temp;
    temp = B[j][i + 3], B[j][i + 3] = B[j + 3][i], B[j + 3][i] = temp;
    temp = B[j + 1][i + 2], B[j + 1][i + 2] = B[j + 2][i + 1], B[j + 2][i + 1] = temp;
    temp = B[j + 1][i + 3], B[j + 1][i + 3] = B[j + 3][i + 1], B[j + 3][i + 1] = temp;
    temp = B[j + 2][i + 3], B[j + 2][i + 3] = B[j + 3][i + 2], B[j + 3][i + 2] = temp;
}

void self_trans_4x4bits(int i, int j, int M, int N, int B[M][N]) // i,j与转置前保持一致
{
    int a0, a1, a2, a3, b0, b1, b2, b3;
    for (int k = 0; k < 4; k++)
    {
        a0 = B[j + 4][i + k], a1 = B[j + 5][i + k], a2 = B[j + 6][i + k], a3 = B[j + 7][i + k];
        b0 = B[j + k][i + 4], b1 = B[j + k][i + 5], b2 = B[j + k][i + 6], b3 = B[j + k][i + 7];
        B[j + k][i + 4] = a0, B[j + k][i + 5] = a1, B[j + k][i + 6] = a2, B[j + k][i + 7] = a3;
        B[j + 4][i + k] = b0, B[j + 5][i + k] = b1, B[j + 6][i + k] = b2, B[j + 7][i + k] = b3;
    }
}

//////////////////////////////////////// helper function for v4 ////////////////////////////////////////
void copy(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0, a1, a2, a3, a4, a5, a6, a7;
    for (int k = 0; k < 4; k++)
    {
        a0 = A[i + k][j], a1 = A[i + k][j + 1], a2 = A[i + k][j + 2], a3 = A[i + k][j + 3], a4 = A[i + k][j + 4], a5 = A[i + k][j + 5], a6 = A[i + k][j + 6], a7 = A[i + k][j + 7];
        B[j][i + k] = a0, B[j + 1][i + k] = a1, B[j + 2][i + k] = a2, B[j + 3][i + k] = a3, B[j][i + k + 4] = a4, B[j + 1][i + k + 4] = a5, B[j + 2][i + k + 4] = a6, B[j + 3][i + k + 4] = a7;
    }
}

void transpose_with_trick(int i, int j, int M, int N, int A[N][M], int B[M][N])
{
    int a0, a1, a2, a3, a4, a5, a6, a7;
    int b0, b1, b2, b3;
    for (int k = 0; k < 4; k++)
    {
        a0 = A[i + 4][j + k], a1 = A[i + 5][j + k], a2 = A[i + 6][j + k], a3 = A[i + 7][j + k];
        a4 = A[i + 4][j + k + 4], a5 = A[i + 5][j + k + 4], a6 = A[i + 6][j + k + 4], a7 = A[i + 7][j + k + 4];
        b0 = B[j + k][i + 4], b1 = B[j + k][i + 5], b2 = B[j + k][i + 6], b3 = B[j + k][i + 7];
        B[j + k][i + 4] = a0, B[j + k][i + 5] = a1, B[j + k][i + 6] = a2, B[j + k][i + 7] = a3;
        B[j + 4 + k][i + 0] = b0, B[j + 4 + k][i + 1] = b1, B[j + 4 + k][i + 2] = b2, B[j + 4 + k][i + 3] = b3;
        B[j + 4 + k][i + 4] = a4, B[j + 4 + k][i + 5] = a5, B[j + 4 + k][i + 6] = a6, B[j + 4 + k][i + 7] = a7;
    }
}

//////////////////////////////////////// helper function for debug ////////////////////////////////////////
void traverse(int M, int N, int B[M][N])
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%d\t", B[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// different version of transpose
void transpose_32x32(int M, int N, int A[N][M], int B[N][M])
{
    const int STEP = 8;
    for (int i = 0; i * STEP < N; i++)
        for (int j = 0; j * STEP < M; j++)
            for (int k = i * STEP; k < i * STEP + STEP && k < N; k++)
                trans_1x8bits(k, j * STEP, M, N, A, B);
}

void transpose_64x64_v1(int M, int N, int A[N][M], int B[N][M])
{
    const int STEP = 4;
    for (int i = 0; i * STEP < N; i++)
        for (int j = 0; j * STEP < M; j++)
            for (int k = i * STEP; k < i * STEP + STEP && k < N; k++)
                trans_1x4bits(k, j * STEP, M, N, A, B);
}

void transpose_64x64_v2(int M, int N, int A[N][M], int B[M][N])
{
    const int STEP = 8;
    const int STEP_X = 4, STEP_Y = 2;
    for (int i = 0; i * STEP < N; i++)
        for (int j = 0; j * STEP < M; j++)
            for (int k = j * STEP; k < (j + 1) * STEP && k < M; k += STEP_X)
                for (int l = i * STEP; l < (i + 1) * STEP && l < N; l += STEP_Y)
                    trans_2x4bits(l, k, M, N, A, B);
}

void transpose_64x64_v3(int M, int N, int A[N][M], int B[M][N])
{
    const int STEP = 8;
    for (int i = 0; i * STEP < N; i++)
        for (int j = 0; j * STEP < M; j++)
        {
            copy_8x8bits(i * STEP, j * STEP, M, N, A, B);
            self_trans_4x4bits_identity(i * STEP + 4, j * STEP + 4, M, N, B);
            self_trans_4x4bits(i * STEP, j * STEP, M, N, B);
            self_trans_4x4bits_identity(i * STEP, j * STEP, M, N, B);
        }
}

void transpose_64x64_v4(int M, int N, int A[N][M], int B[M][N])
{
    const int STEP = 8;
    for (int i = 0; i < N; i += STEP)
        for (int j = 0; j < M; j += STEP)
        {
            copy(i, j, M, N, A, B);
            // traverse(M, N, B);
            transpose_with_trick(i, j, M, N, A, B);
            // traverse(M, N, B);
        }
}

void transpose_61x67(int M, int N, int A[N][M], int B[M][N])
{
    const int STEP = 8;
    for (int i = 0; i * STEP < N; i++)
        for (int j = 0; j * STEP < M; j++)
            for (int k = i * STEP; k < i * STEP + STEP && k < N; k++)
                if ((j + 1) * STEP <= M)
                    trans_1x8bits(k, j * STEP, M, N, A, B);
                else
                    trans_1x5bits(k, j * STEP, M, N, A, B);
}

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    if (M == 32 && N == 32)
        transpose_32x32(M, N, A, B);
    else if (M == 64 && N == 64)
        transpose_64x64_v4(M, N, A, B);
    else if (M == 61 && N == 67)
        transpose_61x67(M, N, A, B);
}

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++)
        {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    // registerTransFunction(trans, trans_desc);
}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; ++j)
        {
            if (A[i][j] != B[j][i])
            {
                return 0;
            }
        }
    }
    return 1;
}

// int main(int argc, char *argv[])
// {
//     int A[16][16], B[16][16];
//     // initialize A
//     for (int i = 0; i < 16; i++)
//         for (int j = 0; j < 16; j++)
//             A[i][j] = i * 16 + j;
//     // initialize B
//     for (int i = 0; i < 16; i++)
//         for (int j = 0; j < 16; j++)
//             B[i][j] = 0;
//     traverse(16, 16, A);
//     // printf("---------------------\n");
//     transpose_64x64_v4(16, 16, A, B);
//     traverse(16, 16, B);
// }
