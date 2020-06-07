# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>
#include "mpi.h"
#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */
#define LIVE 0
#define KILL 1

#define VEC_SIZE 100

void sort(double* vetor,int size);
void intercala(double* vet, int size);

int main(int argc, char** argv)
{
    printf("INICIANDO");

	int my_rank;       // Identificador deste processo

	int proc_n; 

    int size = VEC_SIZE; 

    unsigned long seed = 1;

    double t1,t2;
	MPI_Status status;


    printf("Iniciando MPI INIT");
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_n);
    
    printf("Terminou MPI INIT");
    double vec[VEC_SIZE];
    printf("Inicializou vetor");

    if(my_rank != MASTER)
	{
        printf("Started Receiving size %d", my_rank);

        MPI_Recv(&size, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
        printf("Received size %d", my_rank);
        printf("Started Receiving vec %d", my_rank);
        MPI_Recv(&vec[0], size, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        printf("Received vec %d", my_rank);

        //MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
    }
    else
    {

		printf("Inicialização do vetor pelo mestre");
		t1 = MPI_Wtime();  // inicia a contagem do tempo

        double start = MPI_Wtime();

		/* Send matrix data to the worker tasks */
		
        //RandomNumberGenerator* rnd = new RandomNumberGenerator(seed);

        //CRIA VETOR RANDOM
        for(int i = 0; i < VEC_SIZE; i++)
        {
            vec[i] = (double)rand();//rnd->GetRandom(0.0, 1000000.0);
        }
        printf("Finalizou inicialização do vetor pelo mestre");
    }


    if(size <= VEC_SIZE/proc_n)
    {

        printf("Chegou na folha em %d",my_rank);
        sort(&vec[0], size);

        MPI_Send(&vec[0], size, MPI_DOUBLE, (my_rank - 1) / 2 , 1, MPI_COMM_WORLD);
    }
    else
    {
        printf("Processo %d iniciou", my_rank);
        int newSize = size/2;
        int newSize2 = newSize + size%2;
        
		MPI_Send(&newSize, 1, MPI_INT, my_rank * 2 + 1, 1, MPI_COMM_WORLD);
        MPI_Send(&vec[0], newSize, MPI_DOUBLE, my_rank * 2 + 1, 1, MPI_COMM_WORLD);
        printf("Processo %d mandou metade do vetor para processo %d", my_rank,my_rank * 2 + 1);
		MPI_Send(&newSize2, 1, MPI_INT, my_rank * 2 + 2, 1, MPI_COMM_WORLD);
        MPI_Send(&vec[newSize], newSize2, MPI_DOUBLE, my_rank * 2 + 2, 1, MPI_COMM_WORLD);
        printf("Processo %d mandou metade do vetor para processo %d", my_rank,my_rank * 2 + 2);

        MPI_Recv(&vec[0], newSize, MPI_DOUBLE, my_rank * 2 + 1, 1, MPI_COMM_WORLD, &status);
        printf("Processo %d RECEBEU metade do vetor para processo %d", my_rank,my_rank * 2 + 1);
        MPI_Recv(&vec[newSize], newSize2, MPI_DOUBLE, my_rank * 2 + 2, 1, MPI_COMM_WORLD, &status);
        printf("Processo %d RECEBEU metade do vetor para processo %d", my_rank,my_rank * 2 + 2);
        
        intercala(&vec[0], size);
        printf("Processo %d intercalou vetores", my_rank);

    }

    if(my_rank == MASTER)
    {
        t2 = MPI_Wtime(); // termina a contagem do tempo
        printf("\nTempo de execucao: %f\n\n", t2-t1);   


        for(int i = 0; i < size; i++) { std::cout << vec[i] << ", "; }


    }

    
}

void sort(double* vetor,int size)
{

    for(int i=0;i<size-1;i++)
    {
        for(int j = i+1;j<size;j++)
        {
            if( vetor[i] > vetor[j])
            {
                double aux = vetor[i];
                vetor[i] = vetor[j];
                vetor[j] = aux;
            }
        }
    }

}


void intercala(double* vet, int size)
{
    double aux[size];
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

    vet = aux;
}