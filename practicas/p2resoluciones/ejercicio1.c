#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

double *A,*B,*C;
int N,T,block_size;
void * nombre_funcion(void * ptr){
    // casteo el parametro ya que siempres se recibe como puntero
    int id = *(int*) ptr;
    int inicio = id * block_size; // trabaja donde le correponde, ej id 2, bs 2, inicio 4
    int fin = inicio + block_size;
    
    // para asignar el resto de trabajo al ultimo hilo
    if(id == T - 1){
        fin += N % T;
    }
    for(int i=inicio; i < fin; i++){
        C[i] = A[i]+B[i];
    }
    pthread_exit(0);
}

double dwalltime(){
    double sec;
    struct timeval tv;
    
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

void print_m(){
    for (int i= 0; i<N; i++){
        printf("Valor en la posicion %d = %f\n",i,C[i]);
    }
}
int main(int argc, char *argv[]){
    if((argc != 3) || ((N = atoi(argv[1])) <= 0) || ((T = atoi(argv[2])) <= 0)){
        printf("\nUsar: %s n t \n n: Dimension del vector \n t: Cantidad de hilos \n", argv[0]);
        exit(1);
    }
    int i;
    double timetick;
    int ids[T];
    block_size = N / T; // cuantas secciones le toca a cada hilo
    pthread_attr_t attr;
    pthread_t threads[T];
    
    A = (double *)malloc(sizeof(double)* N);
    B = (double *)malloc(sizeof(double)* N);
    C = (double *)malloc(sizeof(double)* N);
    
    for(int i=0; i< N; i++){
        A[i]=1;
        B[i]=1;
    }

    pthread_attr_init(&attr);

    for(int i=0; i < T;i++){
        ids[i]=i;
        pthread_create(&threads[i], &attr, nombre_funcion, &ids[i]);
    }
    timetick = dwalltime();// inicia
    for(int i=0; i<T; i++){
        pthread_join(threads[i], NULL);
    }
    double time = dwalltime() - timetick;// termina medicion de tiempo
    
    print_m();
    printf("Suma de vectores de dimension %d. Tiempo en segundos %f\n",N,time);
    free(A);
    free(B);
    free(C);
    return(0);
}