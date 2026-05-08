#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

void initvalmat(double *mat, int n, double val, int orden);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);

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
    double *A, *B, *BT, *T1, *R;
    // Matrices: A, B, BT, T1, R de tamaño NxN
    int N, T, BS;
    // argumentos: N (tamaño), T (threads), BS (block size)
    int i, j, k, offsetI, offsetJ;
    double timetick;
    double minA = INT_MAX, maxA = INT_MIN;
    double minB = INT_MAX, maxB = INT_MIN;
    double promA = 0.0, promB = 0.0;
    double escalar = 0.0;

    if ((argc != 4) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((T = atoi(argv[3])) <= 0) || ((N % BS) != 0))
    {
        printf("\nUsar: %s N T BS \n N: dimension, BS: tamaño de bloque, T: hilos\n", argv[0]);
        exit(1);
    }

    double size = (double)(N * N);
    omp_set_num_threads(T);

    A = (double *)malloc(size * sizeof(double));
    B = (double *)malloc(size * sizeof(double));
    BT = (double *)malloc(size * sizeof(double));
    T1 = (double *)malloc(size * sizeof(double));
    R = (double *)malloc(size * sizeof(double));

    initvalmat(A, N, 1.0, 0);
    initvalmat(B, N, 1.0, 1);
    initvalmat(T1, N, 0.0, 0);
    initvalmat(R, N, 0.0, 0);

    timetick = dwalltime();

#pragma omp parallel private(i, j, k, offsetI, offsetJ) // se abre region paralela, cada thread tiene copias privadas de i, j, k, offsetI, offsetJ
    {
// nowait no es necesario ya que cada resultado de los bloques for depende del bloque anterior, las barreras implicitas al final de cada for aseguran que todos los threads hayan terminado antes de avanzar
#pragma omp for reduction(min : minA, minB) reduction(max : maxA, maxB) reduction(+ : promA, promB) schedule(static)
        for (i = 0; i < N; i++)
        {
            offsetI = i * N;
            for (j = 0; j < N; j++)
            {
                double posA = A[offsetI + j];
                double posB = B[offsetI + j];

                if (posA < minA)
                    minA = posA;
                if (posA > maxA)
                    maxA = posA;
                promA += posA;

                if (posB < minB)
                    minB = posB;
                if (posB > maxB)
                    maxB = posB;
                promB += posB;
            }
        }

// cada thread calcula BT por separado para aprovechar el paralelismo
// con nowait porque no se necesita BT completa hasta  la segunda matmul
#pragma omp for nowait schedule(static)
        for (i = 0; i < N; i++)
        {
            offsetI = i * N;
            for (j = 0; j < N; j++)
            {
                BT[offsetI + j] = B[j * N + i];
            }
        }

// single porque no tiene sentido paralelizar el calculo del escalar
// igual que con los for, no se usa nowait porque se necesita el escalar para poder avanzar
#pragma omp single
        {
            promA = promA / size;
            promB = promB / size;
            escalar = (maxA * maxB - minA * minB) / (promA * promB);
        }

// calculo de T1 = A x B, sin nowait porque se necesita el resultado completo para poder hacer la siguiente multplicacion. Barrera implicita.
#pragma omp for schedule(static)
        for (i = 0; i < N; i += BS)
        {
            offsetI = i * N;
            for (j = 0; j < N; j += BS)
            {
                offsetJ = j * N;
                for (k = 0; k < N; k += BS)
                {
                    blkmul(&A[offsetI + k], &B[offsetJ + k], &T1[offsetI + j], N, BS);
                }
            }
        }

// igual que la primera multiplicacion, se necesita el resultado completo.
#pragma omp for schedule(static)
        for (i = 0; i < N; i += BS)
        {
            offsetI = i * N;
            for (j = 0; j < N; j += BS)
            {
                offsetJ = j * N;
                for (k = 0; k < N; k += BS)
                {
                    blkmul(&T1[offsetI + k], &BT[offsetJ + k], &R[offsetI + j], N, BS);
                }
            }
        }

// se aplica nowait porque no hay nada despues que dependa de R, no necesita barrera implicita
#pragma omp for nowait schedule(static)
        for (i = 0; i < N; i++)
        {
            offsetI = i * N;
            for (j = 0; j < N; j++)
            {
                R[offsetI + j] *= escalar;
            }
        }
    }

    double workTime = dwalltime() - timetick;
    printf("Matriz:%d, Bloque:%d, Hilos:%d, Tiempo: %f\n", N, BS, T, workTime);

    free(A);
    free(B);
    free(BT);
    free(T1);
    free(R);

    return 0;
}

void initvalmat(double *mat, int n, double val, int orden)
{
    int i, j, offsetI;
    if (orden == 0)
    {
        for (i = 0; i < n; i++)
        {
            offsetI = i * n;
            for (j = 0; j < n; j++)
                mat[offsetI + j] = val;
        }
    }
    else
    {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[j * n + i] = val;
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
