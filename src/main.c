#include <cblas.h>
#include <stdio.h>

// Test OpenBLAS program
int main() {
    double A[] = {1.0, 0.0, 2.0, 3.0};
    double B[] = {0.0, 2.0, 3.0, 8.0};
    double C[] = {0.0, 0.0, 0.0, 0.0};

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1, A, 2, B, 2, 0, C, 2);
    printf("%.1f %.1f\n", C[0], C[2]);
    printf("%.1f %.1f",   C[1], C[3]);
    return 0;
}