# Apuntes de clase 3 
## Memory wall : Limitaciones en el rendimiento del sistema
Analogia de la manguera
* Latencia: es el tiempo que tarda en llegar el agua al pico desde que se abre la canilla
* Ancho de banda: es la cantidade de agua que puede salir dado el tamaño de la manguera
### Resumen de ideas para mejorar el rendimiento del sistema de memoria
* Explotar localidad espacil y temporal de datos es critico para amortizar la latencia e incrementar el ancho de banda efectivo
* La relacion que existe entre el numero de instrucciones y el numero de accesos a memoria es un buen indicador temprano del rendimiento efectivo del sistema. A mayor numero de instrucciones por numero de accesos a memoria  mejor ser al rendimiento de nuestro sistema.
* La organizacion de datos en la memoria y la forma en que se estrucutura el cosifo pueden impactar significativamente en el rendmiento final del sistema.
## Arreglos multidimiensionales y su organizacion en memoria
De las 4 opciones de arreglos la mejor es Arreglo dinamico como vector de elementos
* Ventajas
    + Favorece al aprovechamiento de la localidad de datos
    + Hace posible el uso de instrucciones SIMD, al permitir que todos los datos esten contiguos en memoria facilita el uso de instrucciones vectoriales
    + Facilita el intercambio de arreglos entre programas escritos en diferentes lenguajes
## Coherencia de cache en arquitecturas multiprocesador
