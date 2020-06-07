# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>
#include "mpi.h"
#include "RandomNumberGenerator.h"
#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */
#define LIVE 0
#define KILL 1

#define VEC_SIZE 100


int main(int argc, char** argv)
{

	int my_rank;       // Identificador deste processo

	int proc_n; 

	MPI_Status status;


	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_n);
    
    double vec[] = new double[VEC_SIZE];
    

    if(my_rank != MASTER)
	{
        MPI_Recv(&size, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);

        MPI_Recv(&vec[0], size, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD, &status);
        //MPI_Recv(&offset, 1, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
    }
    else
    {
        timestamp();

		double t1,t2;
		t1 = MPI_Wtime();  // inicia a contagem do tempo

        double start = MPI_Wtime();

		/* Send matrix data to the worker tasks */
		
        RandomNumerGenerator rnd = new RandomNumerGenerator(1);

        //CRIA VETOR RANDOM
        for(int = 0; i < VEC_SIZE; i++)
        {
            vec[i] = rnd.GetRandom(0.0, 1000000.0);
        }

        



		int task_completed = 0;
		//printf("Started receiving results");

		
		
    }


    if(size <= VEC_SIZE/proc_n)
    {
        sort(&vec[0], size);

        MPI_Send(&vec[0], size, MPI_DOUBLE, (my_rank - 1) / 2 , 1, MPI_COMM_WORLD);
    }
    else
    {
        int newSize = size/2;
        int newSize2 = newSize + size%2;
        
		MPI_Send(&newSize, 1, MPI_INT, my_rank * 2 + 1, 1, MPI_COMM_WORLD);
        MPI_Send(&vec[0], newSize, MPI_DOUBLE, my_rank * 2 + 1, 1, MPI_COMM_WORLD);

		MPI_Send(&newSize2, 1, MPI_INT, my_rank * 2 + 2, 1, MPI_COMM_WORLD);
        MPI_Send(&vec[newSize], newSize2, MPI_DOUBLE, my_rank * 2 + 2, 1, MPI_COMM_WORLD);
     

        MPI_Recv(&vec[0], newSize, MPI_DOUBLE, my_rank * 2 + 1, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&vec[newSize], newSize2, MPI_DOUBLE, my_rank * 2 + 2, 1, MPI_COMM_WORLD, &status);

        intercala(&vec[0], size);

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
    double* aux;
    aux = new double[size];
    int i = 0;
    int j = size/2;
    int total = 0;

    for (i_aux = 0; i_aux < tam; i_aux++) {
        if (((vetor[i] <= vetor[j]) && (i < (size / 2)))
            || (j == tam))
            vetor_auxiliar[i_aux] = vetor[i++];
        else
            vetor_auxiliar[i_aux] = vetor[j++];
    }

    vet = aux;
}