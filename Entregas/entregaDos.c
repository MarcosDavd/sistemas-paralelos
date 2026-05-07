#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>
#include <pthread.h>
#include <float.h>

//funciones
void initvalmat(double *mat, int n, double val, int transpose);
void matmulblks(double *a, double *b, double *c, int n, int bs, int inicio, int fin);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);
void transponer(double *b, double *bt, int n);
double dwalltime();

double *A, *B,*BT, *T1, *R;
int N;
int BS;
int T;
double maxA = -DBL_MAX, minA = DBL_MAX, promA = 0.0;
double maxB = -DBL_MAX, minB = DBL_MAX, promB = 0.0;
double escalar;
pthread_barrier_t barrera;
int blockSize;
pthread_mutex_t mut;
void * mi_funcion(void * arg){
    int id = *((int *)arg);
    int inicio = id * blockSize;
    int fin = inicio + blockSize;
    int i, j, offsetI, offsetJ;
    double valorA,valorB;
    // maximos del hilo
    double tmaxA=-DBL_MAX,tmaxB=-DBL_MAX;
    // minimos del hilo
    double tminA = DBL_MAX,tminB=DBL_MAX;
    // total acumulado del hilo para el promedio total de cada matriz
    double totalA=0, totalB=0;
    
    // Calculo el escalar recorriendo A y B
    int auxi;
    for(int i=inicio; i<fin;i++){
        auxi=i*N;
        for(int j=0; j< N;j++){
            valorA=A[auxi+j];
            valorB=B[auxi+j];
            //busco maximos
            if( valorA > tmaxA ){
                tmaxA=valorA;
            }
            if(valorB > tmaxB){
                tmaxB=valorB;
            }
            //Busco minimos
            if( valorA < tminA ){
                tminA=valorA;
            }
            if(valorB < tminB){
                tminB=valorB;
            }
            //acumulo 
            totalA=totalA + valorA;
            totalB=totalB + valorB;
        }
    }
    pthread_mutex_lock(&mut);
        if(minA > tminA){
            minA=tminA;
        }
        if(tmaxA > maxA){
            maxA=tmaxA;
        }
        if(minB > tminB){
            minB=tminB;
        }
        if(tmaxB > maxB){
            maxB=tmaxB;
        }
        promA+=totalA;
        promB+=totalB;
    pthread_mutex_unlock(&mut);
    // deberia usar un barrier ?
    pthread_barrier_wait(&barrera); // se espera para tenero todos los calores y poder calcular el escalar
    if(id==0){
        promA =promA / (double)(N*N);
        promB = promB / (double)(N*N);
        escalar = (maxA * maxB - minA * minB) / (promA * promB);
    }
    pthread_barrier_wait(&barrera);
    matmulblks(A, B, T1, N, BS, inicio,fin);
    pthread_barrier_wait(&barrera);
    matmulblks(T1, BT, R, N, BS, inicio,fin);
    // escalar * R
    for(int i = inicio; i< fin; i++){
        for (int j = 0; j < N; j++) {
            R[i*N+j] *= escalar;
        }
    }
    pthread_exit(0);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("\nUso: ./%s N BS T\n", argv[0]);
        exit(1);
    }

    N  = atoi(argv[1]);
    BS = atoi(argv[2]);
    T  = atoi(argv[3]);
    if (N <= 0 || BS <= 0 || T <= 0) {
        printf("\nError: N, BS y T deben ser mayores a 0\n");
        exit(1);
    }
    
    double timetick;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_barrier_init(&barrera, NULL, T);
    pthread_mutex_init(&mut,NULL); 
    blockSize= N/T;
    pthread_t hilos[T];
    int ids[T];


    A = (double *)malloc(N * N * sizeof(double));
    B = (double *)malloc(N * N * sizeof(double));
    BT = (double *)malloc(N * N * sizeof(double));
    T1 = (double *)malloc(N * N * sizeof(double)); 
    R = (double *)malloc(N * N * sizeof(double));  // Resultado final

    timetick = dwalltime();
    initvalmat(A, N, 1.0, 0);
    initvalmat(B, N, 1.0, 1); 
    initvalmat(T1,N,0.0,0);
    initvalmat(R,N,0.0,0);
    transponer(B, BT, N);
    for(int i=0; i<T; i++){
        ids[i]=i;
        pthread_create(&hilos[i], &attr,mi_funcion, &ids[i]);
    }
    for(int i=0; i<T;i++){
        pthread_join(hilos[i],NULL);
    }
    

    double workTime = dwalltime() - timetick;

printf("MMBLK-PTHREADS;"
       "N=%d;"
       "T=%d;"
       "BS=%d;"
       "TIME=%lf;"
       "GFLOPS=%lf\n",
       N,
       T,
       BS,
       workTime,
       ((double)4*N*N*N)/(workTime*1e9));
    pthread_mutex_destroy(&mut);
    pthread_barrier_destroy(&barrera);
    free(A); 
    free(B); 
    free(BT);
    free(T1);
    free(R);
    return 0;
}
void transponer(double *b, double *bt, int n) {
    int aux,i, j;
    for (i = 0; i < n; i++) {
        aux=i*n;
        for (j = 0; j < n; j++) {
            bt[j * n + i] = b[aux + j];
        }
    }
}
void initvalmat(double *mat, int n, double val, int transpose) {
    int i, j;
    if (transpose == 0) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[i * n + j] = val;
    } else {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                mat[j * n + i] = val;// con esto invierto los inices para hacer la transpuesta de B
    }
}

void matmulblks(double *a, double *b, double *c, int n, int bs, int inicio, int fin) {
    int i, j, k;

    for (i = inicio; i < fin; i += bs) {
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