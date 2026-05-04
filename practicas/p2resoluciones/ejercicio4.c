#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

//enunciado
/*
Desarrolle un algoritmo paralelo que cuente la cantidad de veces que un elemento X aparece dentro de un
vector de N elementos enteros. Al finalizar, la cantidad de ocurrencias del elemento X debe quedar en una
variable llamada ocurrencias. Para la sincronización emplee mutex-locks. Pruebe con diversos tamaños de
N y T={2,4,8}. Analice el rendimiento
*/


// Vairiable global
int N = 100;
int miVector[100];
int ocurrencias = 0;
int x = 10;// El numero que se debe buscar 
int T=4; // cantidad de hilos
pthread_mutex_t mut;

void * analizar(void * ptr){
    int id = *(int *)ptr;
    int carga = N/T; // carga de trabajo que le toca a cada hilo 
    int inicio= id * carga;
    int fin = inicio + carga;
    for(int i=inicio; i<fin; i++){
        if(miVector[i]== x){
            pthread_mutex_lock(&mut);
                ocurrencias++;
            pthread_mutex_unlock(&mut);
        }
    }
    return NULL;
}

int main(){

    int ids[T];
    pthread_t hilo[T];
    pthread_mutex_init(&mut,NULL);
    
    srand(time(NULL));
    for(int i=0; i<N; i++){
        miVector[i]=rand() % N;// nuemro random del 0 a N 
    }
    
    for(int i=0; i<T; i++){
        pthread_create(&hilo[i],NULL,analizar,&ids[i]);
    }
    for(int i=0; i<T;i++){
        pthread_join(hilo[i],NULL);
    }
    printf("Cantidad de ocurrencias del numero 10 : %d\n",ocurrencias);
    pthread_mutex_destroy(&mut);
    return 0;
}
