#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <pthread.h>

#define BS 64

double *A, *B, *C, *D, *R, *resultMatriz;
int N, T, blockSize;

void initvalmat(double *mat, int n, double val, int orden);

void matmulblks(double *a, double *b, double *c, int n, int bs, int inicio, int fin);

void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);

void *worker(void *arg)
{
    int id = *((int *)arg);
    int inicio = id * blockSize;
    int fin = inicio + blockSize;
    int i, j, offsetI, offsetJ;

    // Se realiza la multiplicacion de matrices A y B
    matmulblks(A, B, resultMatriz, N, BS, inicio, fin);

    // Se realiza la multiplicacion de matrices C y D
    matmulblks(C, D, R, N, BS, inicio, fin);


    // Se suman los resultados de las multiplicaciones de matrices y al mismo tiempo se multiplica por el escalar resultMatriz
    for (i = inicio; i < fin; i++)
    {
        offsetI = i * N;
        for (j = 0; j < N; j++)
        {
            R[offsetI + j] += resultMatriz[offsetI + j];
        }
    }

    pthread_exit(0);
}

double dwalltime()
{
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

int main(int argc, char *argv[])
{
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0))
    {
        printf("\nUsar: %s n t \n  n: Dimension de la matriz (nxn X nxn) \n  t: Cantidad de hilos \n", argv[0]);
        exit(1);
    }

    int i, j, offsetI, offsetJ;
    double timetick;
    int ids[T];
    pthread_attr_t attr;
    pthread_t threads[T];
    blockSize = N / T;
    int    size = N * N;

    A = (double *)malloc(size * sizeof(double));
    B = (double *)malloc(size * sizeof(double));
    C = (double *)malloc(size * sizeof(double));
    D = (double *)malloc(size * sizeof(double));
    R = (double *)malloc(size * sizeof(double));
    resultMatriz = (double *)malloc(size * sizeof(double));

    // Inicializa las matrices A, B, C y D en 1
    initvalmat(A, N, 1.0, 0); // Orden 0 para que se inicialice en orden de filas
    initvalmat(B, N, 1.0, 1); // Orden 1 para que se inicialice en orden de columnas
    initvalmat(C, N, 1.0, 0); // Orden 0 para que se inicialice en orden de filas
    initvalmat(D, N, 1.0, 1); // Orden 1 para que se inicialice en orden de columnas
    initvalmat(R, N, 0.0, 0);
    initvalmat(resultMatriz, N, 0.0, 0);


    // Crea los hilos
    for (i = 0; i < T; i++)
    {
        ids[i] = i;
        pthread_create(&threads[i], &attr, worker, &ids[i]);
    }

    timetick = dwalltime();

    for (i = 0; i < T; i++)
    {
        pthread_join(threads[i], NULL);
    }

    double workTime = dwalltime() - timetick;

    printf("Tiempo en segundos %f\n", workTime);


    free(A);
    free(B);
    free(C);
    free(D);
    free(R);
    free(resultMatriz);

    return (0);
}

// Inicializa una matriz de nxn con un valor val
void initvalmat(double *mat, int n, double val, int orden)
{
    int i, j, offsetI;

    if (orden == 0)
    {
        for (i = 0; i < n; i++)
        {
            offsetI = i * n;
            for (j = 0; j < n; j++)
            {
                mat[offsetI + j] = val;
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                mat[j * n + i] = val;
            }
        }
    }
}

// Realiza la multiplicacion de matrices en bloques de tamaño bs
void matmulblks(double *a, double *b, double *c, int n, int bs, int inicio, int fin)
{
    int i, j, k, offsetI, offsetJ;

    for (i = inicio; i < fin; i += bs)
    {
        offsetI = i * n;
        for (j = 0; j < n; j += bs)
        {
            offsetJ = j * n;
            for (k = 0; k < n; k += bs)
            {
                blkmul(&a[offsetI + k], &b[offsetJ + k], &c[offsetI + j], n, bs);
            }
        }
    }
}

void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs)
{
    int i, j, k, offsetI, offsetJ;
    double suma;

    for (i = 0; i < bs; i++)
    {
        offsetI = i * n;
        for (j = 0; j < bs; j++)
        {
            suma = 0;
            offsetJ = j * n;
            for (k = 0; k < bs; k++)
            {
                suma += ablk[offsetI + k] * bblk[offsetJ + k];
            }

            cblk[offsetI + j] += suma;
        }
    }
}