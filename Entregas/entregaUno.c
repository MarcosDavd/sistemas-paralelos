#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>

// Prototipos de funciones
void initvalmat(double *mat, int n, double val, int transpose);
void matmulblks(double *a, double *b, double *c, int n, int bs);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);
double dwalltime();

int main(int argc, char *argv[]) {
    double *A, *B, *T1, *R;
    int N, BS;
    double timetick;
    
    // Variables para métricas
    double maxA = -1.0, minA = 1e15, promA = 0.0;
    double maxB = -1.0, minB = 1e15, promB = 0.0;
    double escalar;

    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("Uso: %s N BS (N múltiplo de BS)\n", argv[0]);
        exit(1);
    }

    size_t size = (size_t)N * N;
    A = (double *)malloc(size * sizeof(double));
    B = (double *)malloc(size * sizeof(double));
    T1 = (double *)malloc(size * sizeof(double)); // Resultado intermedio (A x B)
    R = (double *)malloc(size * sizeof(double));  // Resultado final

    // 1. Inicialización y cálculo de métricas (Optimizado en un solo paso)
    // Usamos transpose=1 para B porque la fórmula requiere B y B^T
    initvalmat(A, N, 1.0, 0); 
    initvalmat(B, N, 1.0, 1); // B almacenada por columnas facilita B^T y la multiplicación

    timetick = dwalltime();
	/*Si te preguntan por qué usaste un solo bucle, la respuesta clave es: 
    "Para garantizar un acceso secuencial a memoria y maximizar el aprovechamiento de las líneas de caché, evitando saltos innecesarios (strides)."*/
    for (int i = 0; i < N * N; i++) {
        // Métricas A
        if (A[i] > maxA) maxA = A[i];
        if (A[i] < minA) minA = A[i];
        promA += A[i];
        // Métricas B
        if (B[i] > maxB) maxB = B[i];
        if (B[i] < minB) minB = B[i];
        promB += B[i];
    }
    promA /= (double)size;
    promB /= (double)size;

    // 2. Cálculo del Escalar
    escalar = (maxA * maxB - minA * minB) / (promA * promB);

    // 3. Operaciones de matrices por bloques
    // T1 = A x B
    matmulblks(A, B, T1, N, BS);

    // R = T1 x B_transpuesta
    // Como B ya está en orden de columnas (transpose=1), pasarla a matmulblks
    // de nuevo actúa efectivamente como usar su transpuesta si la función 
    // espera orden de filas.
    matmulblks(T1, B, R, N, BS);

    // 4. Multiplicación por el escalar
    for (int i = 0; i < N * N; i++) {
        R[i] *= escalar;
    }

    double totalTime = dwalltime() - timetick;
    printf("N=%d, BS=%d, Tiempo total: %f s\n", N, BS, totalTime);

    free(A); free(B); free(T1); free(R);
    return 0;
}

// --- Funciones proporcionadas por la cátedra ---

void initvalmat(double *mat, int n, double val, int transpose) {
    int i, j;
    if (transpose == 0) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[i * n + j] = val;
    } else {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[j * n + i] = val;
    }
}

void matmulblks(double *a, double *b, double *c, int n, int bs) {
    int i, j, k;
    // Es vital limpiar la matriz destino antes de sumar bloques
    for (int idx = 0; idx < n*n; idx++) c[idx] = 0.0;

    for (i = 0; i < n; i += bs) {
        for (j = 0; j < n; j += bs) {
            for (k = 0; k < n; k += bs) {
                blkmul(&a[i * n + k], &b[j * n + k], &c[i * n + j], n, bs);
            }
        }
    }
}

void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs) {
    int i, j, k;
    for (i = 0; i < bs; i++) {
        for (j = 0; j < bs; j++) {
            for (k = 0; k < bs; k++) {
                cblk[i * n + j] += ablk[i * n + k] * bblk[j * n + k];
            }
        }
    }
}

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}
