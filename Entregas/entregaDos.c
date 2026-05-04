#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>
#include <pthread.h>
//funciones
void initvalmat(double *mat, int n, double val, int transpose);
void matmulblks(double *a, double *b, double *c, int n, int bs, int inicio, int fin);
void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs);
double dwalltime();

double *A, *B,*BT, *T1, *R;
int N, BS,T;
double maxA = -1.0, minA = 99999.0, promA = 0.0;
double maxB = -1.0, minB = 99999.0, promB = 0.0;
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
    double tmaxA=-1,tmaxB=-1;
    // minimos del hilo
    double tminA = 9999,tminB=9999;
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
    
    matmulblks(A, B, T1, N, BS, inicio,fin);
    pthread_barrier_wait(&barrera);// esto estaria bien ?
    matmulblks(T1, BT, R, N, BS, inicio,fin);
    pthread_barrier_wait(&barrera);
    // escalar * R
    for(int i = inicio; i< fin; i++){
        for (int j = 0; j < N * N; j++) {
            R[i*N+j] *= escalar;
        }
    }
    pthread_barrier_wait(&barrera);
    pthread_exit(0);
}

int main(int argc, char *argv[]) {
    double timetick;
    int ids[T];
    pthread_attr_t attr;
    pthread_t hilos[T];

    pthread_attr_init(&attr);
    pthread_mutex_init(&mut,NULL); 
    
    if ((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((BS = atoi(argv[2])) <= 0) || ((N % BS) != 0)) {
        printf("\nError en los parámetros. Usage: ./%s N BS (N debe ser multiplo de BS)\n", argv[0]);
        exit(1);
    }
     blockSize= N/T;


    A = (double *)malloc(N * N * sizeof(double));
    B = (double *)malloc(N * N * sizeof(double));
    BT = (double *)malloc(N * N * sizeof(double));
    T1 = (double *)malloc(N * N * sizeof(double)); // Resultado intermedio (B x BT)
    R = (double *)malloc(N * N * sizeof(double));  // Resultado final

    initvalmat(A, N, 1.0, 0);
    initvalmat(B, N, 1.0, 0); 
    initvalmat(T1,N,0.0,0);
    initvalmat(BT, N, 1.0, 1); // B^T
    initvalmat(R,N,0.0,0);
    timetick = dwalltime();
    //Puedo usar un solo for por la forma en que guarde las matrices 
    // con N*N cubro todo el recorrido
    /*
    for (int i = 0; i < N * N; i++) {
        if (A[i] > maxA) maxA = A[i];
        if (A[i] < minA) minA = A[i];
        promA += A[i];
        
        if (B[i] > maxB) maxB = B[i];
        if (B[i] < minB) minB = B[i];
        promB += B[i];
    }
    */
    pthread_barrier_init(&barrera, NULL, T);

    for(int i=0; i<T; i++){
        ids[i]=i;
        pthread_create(&hilos[i], &attr,mi_funcion, &ids[i]);
    }
    for(int i=0; i<T;i++){
        pthread_join(hilos[i],NULL);
    }
    

    double workTime = dwalltime() - timetick;

  printf("MMBLK-SEC;%d;%d;%lf;%lf\n",N,BS,workTime,((double)2*N*N*N)/(workTime*1000000000));
    
    pthread_mutex_destroy(&mut);
    pthread_barrier_destroy(&barrera);
    free(A); 
    free(B); 
    free(T1);
    free(BT);
    free(R);
    return 0;
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
    // preparo la matriz acumuladora

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