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

#define VEC_SIZE 1000000

void sort(double* vetor,int size);
void intercala(double* vet, int size);
double vec[VEC_SIZE];

int main(int argc, char** argv)
{

	int my_rank;       // Identificador deste processo

	int proc_n;        // numero total de processos

    int size = VEC_SIZE; // tamanho atual do vetor


    double t1,t2;
	MPI_Status status;


	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_n);
    

    if(my_rank != MASTER)
	{

        //Se este processo nao eh o mestre, recebe de seu pai (my_rank - 1) / 2 o vetor a ser processado ou novamente dividido e o tamanho deste vetor

        MPI_Recv(&size, 1, MPI_INT, (my_rank - 1) / 2, 1, MPI_COMM_WORLD, &status);

        MPI_Recv(&vec[0], size, MPI_DOUBLE, (my_rank - 1) / 2, 1, MPI_COMM_WORLD, &status);

    }
    else
    {

		printf("Inicialização do vetor pelo mestre\n");
		t1 = MPI_Wtime();  // inicia a contagem do tempo

        double start = MPI_Wtime();

        //CRIA VETOR RANDOM
        for(int i = 0; i < VEC_SIZE; i++)
        {
            vec[i] = (double)rand();//rnd->GetRandom(0.0, 1000000.0);
        }

        printf("Finalizou inicialização do vetor pelo mestre\n");

        //Se este for o unico processo alocado, se faz a ordenacao total do vetor
        if(proc_n == 1)
        {
            printf("Mestre ordenando vetor sozinho\n");
            sort(&vec[0],size);
            t2 = MPI_Wtime(); // termina a contagem do tempo
            printf("\nTempo de execucao: %f\n\n", t2-t1);
            MPI_Finalize();
            return 0;   
        }

    }

    int delta =  VEC_SIZE/((proc_n+1)/2);
    if(size <= delta+1)
    {

        printf("Chegou na folha em %d\n",my_rank);
        sort(&vec[0], size);
        
        MPI_Send(&vec[0], size, MPI_DOUBLE, (my_rank - 1) / 2 , 1, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
       
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




    }else
    {
        MPI_Send(&vec[0], size, MPI_DOUBLE, (my_rank - 1) / 2 , 1, MPI_COMM_WORLD);
    }
    

    MPI_Finalize();
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
    


    for(int i = 0; i < size; i++) { vet[i] = aux[i]; }



}