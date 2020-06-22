#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define DEBUG 1            // comentar esta linha quando for medir tempo
#define ARRAY_SIZE 1000000  // trabalho final com o valores 10.000, 100.000, 1.000.000
#define PARTE 0.5
int vetor[ARRAY_SIZE];
/*
void bs(int n, int * vetor)
{
    int c=0, d, troca, trocou =1;

    while (c < (n-1) & trocou )
        {
        trocou = 0;
        for (d = 0 ; d < n - c - 1; d++)
            if (vetor[d] > vetor[d+1])
                {
                troca      = vetor[d];
                vetor[d]   = vetor[d+1];
                vetor[d+1] = troca;
                trocou = 1;
                }
        c++;
        }
}
*/
void intercala( int size, int* vet)
{
    int aux[size];
    int i = 0;
    int j = size/2;
    int total = 0;
    for (int i_aux = 0; i_aux < size; i_aux++) {
        if (((vet[i] <= vet[j]) && (i < (size / 2)))
            || (j == size))
            aux[i_aux] = vet[i++];
        else
            aux[i_aux] = vet[j++];
    }
    for(int i = 0; i < size; i++) { vet[i] = aux[i]; }
}
void bs(int size, int* vetor)
{
    for(int i=0;i<size-1;i++)
    {
        for(int j = i+1;j<size;j++)
        {
            if( vetor[i] > vetor[j])
            {
                int aux = vetor[i];
                vetor[i] = vetor[j];
                vetor[j] = aux;
            }
        }
    }
}


int main(int argc, char** argv)
{
    int my_rank, proc_n;
    bool pronto = false;
    MPI_Status status;
    int part;
    double t1,t2;
    MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_n);
    for(int i = 0; i < ARRAY_SIZE; i++) vetor[i] = ARRAY_SIZE-i-1;
    if(my_rank == 0 ) t1 = MPI_Wtime();  // inicia a contagem do tempo
    while(!pronto)
    {
        int estado[proc_n];
        for(int i = 0; i<proc_n;i++)
        {
            estado[i] = 0;
        }
        int tam = ARRAY_SIZE/proc_n;
        int resto = ARRAY_SIZE%proc_n;
        int maior_elem;
        int tam_aux = tam;
        int tam_aux_ultimo = tam_aux + resto;
        if(my_rank == proc_n-1)
        {
            tam_aux+=resto;
        }
        part = tam * PARTE;
        bs(tam_aux, &vetor[my_rank * tam]);
        if(my_rank != proc_n-1)
        {
            MPI_Send(&vetor[my_rank*tam + tam_aux - 1], 1, MPI_INT, my_rank + 1, 1, MPI_COMM_WORLD);
        }
        if(my_rank != 0)
        {
            MPI_Recv(&maior_elem, 1, MPI_INT, my_rank - 1, 1, MPI_COMM_WORLD, &status);
            if(maior_elem <= vetor[my_rank*tam])
            {
                estado[my_rank] = 1;
            }
        }
        pronto = true;
        for(int i = 1; i <proc_n; i++)
        {
            MPI_Bcast(&estado[i], 1, MPI_INT, i, MPI_COMM_WORLD);
            if(estado[i] == 0)
            {
                pronto = false;
                break;
            }
        }
        if(pronto) 
        {
            if(my_rank == 0)
            {

                t2 = MPI_Wtime();
                printf("\nTempo de execucao: %f\n\n", t2-t1);
                for(int i = 1; i <proc_n; i++)
                {
                    if(i == proc_n-1)
                    {
                        MPI_Recv(&vetor[i*tam], tam_aux_ultimo, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                    }else
                    {
                        MPI_Recv(&vetor[i*tam], tam_aux, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                    }                    
                    
                }                
                bool sorted = true;
                for (int i = 0; i < ARRAY_SIZE - 1; i++)
                {
                    if(vetor[i] >= vetor[i + 1])
                    {
                        sorted = false;
                        break;
                    }
                }                
                printf("VETOR SORTED: %s\n", sorted ? "true" : "false");               
            }
            else
            {
                MPI_Send(&vetor[my_rank*tam], tam_aux, MPI_INT, 0, 1, MPI_COMM_WORLD);
            }
            break;
        }
        if(my_rank != 0)
        {
            MPI_Send(&vetor[my_rank*tam], part, MPI_INT, my_rank - 1, 1, MPI_COMM_WORLD);
        }
        int vetor_aux[part * 2];
        if(my_rank != proc_n - 1)
        {
            for (int i = 0; i < part; i++)
            {
                vetor_aux[i] = vetor[my_rank*tam + tam_aux - part + i];
            }
            MPI_Recv(&vetor_aux[part], part, MPI_INT, my_rank + 1, 1, MPI_COMM_WORLD, &status);
            intercala(part*2, vetor_aux);
            MPI_Send(&vetor_aux[part], part, MPI_INT, my_rank + 1, 1, MPI_COMM_WORLD);
            for(int i = 0; i< part;i++)
            {
                vetor[my_rank*tam + tam_aux - part + i] = vetor_aux[i];
            }
        }
        if(my_rank != 0){
            MPI_Recv(&vetor[my_rank*tam], part, MPI_INT, my_rank - 1, 1, MPI_COMM_WORLD, &status);
        }
    }
    MPI_Finalize();
    return 0;
}