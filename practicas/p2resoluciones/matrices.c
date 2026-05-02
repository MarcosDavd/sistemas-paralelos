#include<stdio.h>
#include<stdlib.h>
#include <sys/time.h>
#include <pthread.h>
 int N,T, bloque;
  double *A,*B,*C;
void * funcion_matrices(void * ptr){
    int id = *(int *)ptr;
    int ini = bloque * id;  
    int fin = ini + bloque;
    int aux ;
    
    for(int i=ini; i < fin;i++){
        for(int j=0;j < N;j++){
            aux=0;
            for(int k=0;k < N;k++){
                aux+= A[i*N+k] * B[j*N+k];
            }
            C[i*N+j]= aux; 
        }
    }  
}


double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}

int main(int argc,char*argv[]){

 if((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0)){
        printf("\nUsar: %s n t \n n: Dimension del vector \n t: Cantidad de hilos \n", argv[0]);
        exit(1);
    }
 
  int ids[T];
  double timetick;
  bloque = N/T;
  pthread_attr_t attr;
  pthread_t threads[T];
  int i,j;

  A=(double*)malloc(sizeof(double)*N*N);
  B=(double*)malloc(sizeof(double)*N*N);
  C=(double*)malloc(sizeof(double)*N*N);

 for(i=0;i<N;i++){
   for(j=0;j<N;j++){
    A[i*N+j]=1;
    B[j*N+i]=1;
   }
  }   

  pthread_attr_init(&attr);// config de hilos
  for(int i = 0; i < T; i++){
    ids[i]=i;
    pthread_create(&threads[i],&attr,funcion_matrices,&ids[i]);
  }

  timetick = dwalltime();
  for(int i=0; i< T; i++){
    pthread_join(threads[i],NULL);
  }  
double time  =  dwalltime() -timetick;
printf("Multiplicacion de matrices de %dx%d. Tiempo en segundos %f\n",N,N,time);

 free(A);
 free(B);
 free(C);
 return(0);
}