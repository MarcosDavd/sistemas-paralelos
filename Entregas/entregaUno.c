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
    double *A, *B,*BT, *T1, *R;
    int N, BS;
    double timetick;
    
    // Variables para métricas
    double maxA = -1.0, minA = 99999.0, promA = 0.0;
    double maxB = -1.0, minB = 99999.0, promB = 0.0;
    double escalar;

    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
        exit(1);
    }


    A = (double *)malloc(N * N * sizeof(double));
    B = (double *)malloc(N * N * sizeof(double));
    BT = (double *)malloc(N * N * sizeof(double));
    T1 = (double *)malloc(N * N * sizeof(double)); // Resultado intermedio (A x B)
    R = (double *)malloc(N * N * sizeof(double));  // Resultado final

    // 1. Inicialización y cálculo de métricas (Optimizado en un solo paso)
    // Usamos transpose=1 para B porque la fórmula requiere B y B^T
    // Esto em ayuda a organizar las matrices en memoria fisica
    // Quedaria como un arreglo de una dimension
    initvalmat(A, N, 1.0, 0);
    initvalmat(B, N, 1.0, 0); 
    initvalmat(BT, N, 1.0, 1); // B^T

    timetick = dwalltime();
    //Para garantizar un acceso secuencial a memoria y maximizar el aprovechamiento de las líneas de caché, evitando saltos innecesarios (cahce miss)."*/
    //Puedo usar un solo for por la forma en que guarde las matrices 
    // con N*N cubro todo el recorrido
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
    promA =promA / (double)(N*N);
    promB = promB / (double)(N*N);
 
    // 2. Cálculo del Escalar
    escalar = (maxA * maxB - minA * minB) / (promA * promB);

    // T1 = A x B
    // Trabajo la ultiplicacion e las matrices A y B 
    // con bloques e tamaño BS aprovechando mejor el uso de 
    // Localidad espacial y remporal
    matmulblks(A, B, T1, N, BS);

    // R = T1 x B^t
    // B esta guardada por columnas en mem fisica que es lo mismo que 
    // tenerla transpuesta 
    matmulblks(T1, BT, R, N, BS);

    // escalar * R
    for (int i = 0; i < N * N; i++) {
        R[i] *= escalar;
    }

    double workTime = dwalltime() - timetick;

  printf("MMBLK-SEC;%d;%d;%lf;%lf\n",N,BS,workTime,((double)2*N*N*N)/(workTime*1000000000));

    free(A); 
    free(B); 
    free(T1); 
    free(R);
    return 0;
}

// cosas de la catedra

void initvalmat(double *mat, int n, double val, int transpose) {
    int i, j;
    if (transpose == 0) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[i * n + j] = val;// aca es como si leyera la memoria de corrido, no se prodcen saltos ni miss cache
    } else {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[j * n + i] = val;// esto invierto los inices para hacer la transpuesta de B
    }
}

void matmulblks(double *a, double *b, double *c, int n, int bs) {
    int i, j, k;
    // preparo la matriz acumuladora
    initvalmat(c,n,0.0,0);

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

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}