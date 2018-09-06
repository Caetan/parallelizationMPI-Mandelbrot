/****
*

Copyleft (C) 2018 Caetán Tojeiro Carpente

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/agpl-3.0.html

*
****/


/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mpi.h>

#define DEBUG 1
#define ANALISIS 1

#define         X_RESN  1024     /* x resolution */
#define         Y_RESN  1024       /* y resolution */
#define         X_MIN   -2.0	/*Boundaries of the mandelbrot set*/
#define         X_MAX    2.0
#define         Y_MIN   -2.0
#define         Y_MAX    2.0
#define		maxIterations	1000 /*Cuanto mayor, más detalle en la imagen y mayor coste computacional*/
			


typedef struct complextype
        {
        float real, imag;
        } Compl;



void analizar(int numprocs, int * flops, int * microseconds_comp, int * microseconds_comun){
    int i;
    int microseconds_comun_total = 0;
    int microseconds_comp_total = 0;
    int max_flops = 0;
    float bp = 0;
    int flops_total=0;
        
    FILE *f = fopen("resultados.txt", "w");
		if(f==NULL){
			printf("Error abriendo fichero\n");
			exit(1);
		}    
		
		
		fprintf(f, "TIEMPO DE COMUNICACIONES\n");
		for(i=0;i<numprocs;i++){
			fprintf(f, "El proceso %i ha tardado en comunicar %lf segundos\n", i, (double) microseconds_comun[i]/1E6);
			microseconds_comun_total += microseconds_comun[i];
		}
		
		
		fprintf(f, "\n\n");
		
		fprintf(f, "TIEMPO DE COMPUTACION\n");
		for(i=0;i<numprocs;i++){
			fprintf(f, "El tiempo de computacion del proceso %i es %lf segundos\n", i, (double) microseconds_comp[i]/1E6);
			microseconds_comp_total += microseconds_comp[i];
		}
		
		
		fprintf(f, "\n\n");
		
		fprintf(f, "NUMERO DE FLOPS\n");
		for(i=0;i<numprocs;i++){
			fprintf(f, "El proceso %i ha realizado %i flops\n", i, flops[i]);
			flops_total+=flops[i];
			if(flops[i]> max_flops){
				max_flops = flops[i];
			}
		}
		
		fprintf(f, "\n\n");
		
		fprintf(f, "El tiempo total de computacion total es %lf segundos\n", (double) microseconds_comp_total/1E6);
		fprintf(f, "El tiempo total de comunicaciones total es %lf segundos\n", (double) microseconds_comun_total/1E6);
		fprintf(f, "El total de flops realizados es de %i \n", flops_total);
		
		fprintf(f, "\n\n");
		
		bp = (double) flops_total/(max_flops*numprocs);
		fprintf(f, "siendo el numero máximo de flops que un proceso ha realizado %i\n", max_flops);
		fprintf(f, "Balanceo de carga = %lf \n", bp);

	fclose(f);

}


int main ( int argc, char *argv[] )
{

       /* Mandelbrot variables */
        int i, j, k;
        Compl   z, c;
        float   lengthsq, temp;
        int res[X_RESN][Y_RESN]; 
		struct timeval  ti, tf,tci, tcf;
		int microseconds,microseconds_c, flop = 0;
    	int *microseconds_total, *microseconds_c_total, *flop_total;
		
		//Inicializacion MPI
		int numprocs, miId;
		
		MPI_Status st;
		
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &miId);
		
		
		int proceso[numprocs], proceso_total[numprocs];
		int acc=0;
		int displs[numprocs];
		
		if (miId == 0) {
			if (ANALISIS){
				flop_total = malloc(numprocs * sizeof(double));
				microseconds_total = malloc(numprocs * sizeof(int));
				microseconds_c_total = malloc(numprocs * sizeof(int));
			}
		}
		
		
		//Calculamos cuantas filas tiene cada proceso
		if((X_RESN%numprocs)==0){  //en el caso de ser divisor exacto del numero de procesos, se asigna a cada uno el mismo numero de filas
			for(i=0;i<numprocs;i++){
				proceso[i] = X_RESN/numprocs;
				proceso_total[i] = proceso[i] * Y_RESN;  //paraa saber cuanto tiene que computar en total cada proceso, tambien la dimension Y
				displs[i] = i*proceso_total[i];  //desplazamiento dentro del rcvbuff (suma acumulada de los recvcount)
			}
		} else{				//en caso de no ser divisor exacto del numero de procesos
			for(i=0;i<(numprocs-1);i++){ //a todos los procesos menos al último se les asigna ceil(X_RESN/numprocs) filas
				proceso[i] = ((X_RESN/numprocs)+1);	 //No funciona ceil, hace floor
				proceso_total[i] = proceso[i] * Y_RESN;
				displs[i] = i*proceso_total[i];  //desplazamiento dentro del rcvbuff	
			}
			
			acc = proceso[0] * (numprocs-1);  //Para calcular cuantas filas llevamos ya asignadas a todos los procesos menos al ultimo
			proceso[numprocs-1] = X_RESN-acc;  //Para saber cuantas filas le corresponden al ultimo proceso
			proceso_total[numprocs-1] = proceso[numprocs-1] * Y_RESN;
			displs[numprocs-1] = displs[numprocs-2] + proceso_total[numprocs-2];  //el displs para el ultimo proceso sera la suma del displs anterior mas lo que ocupa en el "array" el proceso_total anterior
		}
		
		
		int res_aux[proceso[miId]] [Y_RESN];  //array con el numero de filas de cada proceso y la dimension de Y (el ultimo proceso puede tener un numero de filas diferente)
		


		/* Start measuring time */
		gettimeofday(&ti, NULL);
         
        /* Calculate and draw points */
        for(i=0; i < proceso[miId]; i++) {  //Cada proceso computa sus filas
			
			for(j=0; j < Y_RESN; j++) {  //y todas las columnas de sus filas
			  z.real = z.imag = 0.0;
			  c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
			  c.imag = Y_MAX - (i + (proceso[0] * miId)) * (Y_MAX - Y_MIN)/Y_RESN;  //Modificación para computar las filas que caen al proceso correspondiente, proceso[0] porque tiene el número de filas que caen a cada proceso (menos al último que puede variar)
			  
			  k = 0;

			  do  {    /* iterate for pixel color */

				temp = z.real*z.real - z.imag*z.imag + c.real;
				z.imag = 2.0*z.real*z.imag + c.imag;
				z.real = temp;
				lengthsq = z.real*z.real+z.imag*z.imag;
				k++;

			  } while (lengthsq < 4.0 && k < maxIterations);

			if (k >= maxIterations) res_aux[i][j] = 0;
			else res_aux[i][j] = k;

			}
			
			flop=flop+k*10;
			
		}
		
		/* End measuring time */
		gettimeofday(&tf, NULL);
		microseconds = (tf.tv_usec - ti.tv_usec)+ 1000000 * (tf.tv_sec - ti.tv_sec);
		/*printf ("(PERF) Time (seconds) = %lf\n", (double) microseconds/1E6);*/
		
		//Medicion del tiempo de comunicacion de cada proceso	
		gettimeofday(&tci, NULL);
		MPI_Gatherv(&res_aux, proceso_total[miId], MPI_INT, &res, proceso_total, displs, MPI_INT, 0, MPI_COMM_WORLD);  //El root siempre será el 0
		gettimeofday(&tcf, NULL);
		microseconds_c = (tcf.tv_usec - tci.tv_usec)+ 1000000 * (tcf.tv_sec - tci.tv_sec);
	
	if( DEBUG && miId == 0) {
		for(i=0;i<X_RESN;i++) {
			for(j=0;j<Y_RESN-1;j++)
		        	printf("%d\t", res[i][j]);
			printf( "%d\n", res[i][Y_RESN-1] );
		}
	}
	
	
        MPI_Gather(&flop, 1, MPI_INT, flop_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&microseconds, 1, MPI_INT, microseconds_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(&microseconds_c, 1, MPI_INT, microseconds_c_total, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (miId == 0){
            analizar(numprocs, flop_total, microseconds_total, microseconds_c_total);

            free(microseconds_c_total);
            free(microseconds_total);
            free(flop_total);
        }


		
/*	
	if(miId==0){
		FILE *f = fopen("resultados.txt", "w");
		if(f==NULL){
			printf("Error abriendo fichero\n");
			exit(1);
		}
		
		
		fprintf(f, "TIEMPO DE COMUNICACIONES\n");
		for(i=0;i<numprocs;i++){
			fprintf(f, "El proceso %i ha tardado en comunicar %lf microsegundos\n", i, (double) microseconds_comun[i]/1E6);
			microseconds_total += microseconds_comun[i];
		}
		
		fprintf(f, "\n\n");
		
		fprintf(f, "NUMERO DE FLOPS\n");
		for(i=0;i<numprocs;i++){
			fprintf(f, "El proceso %i ha realizado %lf miles de flops\n", i, (double) flops[i]/1E3);
			flops_total += flops[i];
		}
		
		fprintf(f, "\n\n");
		
		fprintf(f, "El tiempo total de comunicaciones es %lf microsegundos\n", (double) microseconds_total/1E6);
		fprintf(f, "El total de flops realizados es de %lf mil\n", (double) flops_total/1E3);
		
		fclose(f);
}	
*/	
	MPI_Finalize();

}
