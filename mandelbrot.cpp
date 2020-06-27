# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>
# include "mpi.h"
# include <omp.h>

#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */
#define TAREFAS 7 // Numero de tarefas no saco de trabalho para np = 8, processo 0 é o mestre
#define LIVE 0
#define KILL 1


using namespace std;
int main(int argc, char** argv);
int i4_min(int i1, int i2);
void i4pp_delete(int** a, int m, int n);
int** i4pp_new(int m, int n);
void timestamp();

//****************************************************************************80

int main(int argc, char** argv)


//****************************************************************************80
//
//  Purpose
//
//    MAIN is the main program for MANDELBROT_OPENMP.
//
//  Discussion:
//
//    MANDELBROT_OPENMP computes an image of the Mandelbrot set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    Local, int COUNT_MAX, the maximum number of iterations taken
//    for a particular pixel.
//
{

	int m = 1000; //tamanho da imagem em linhas
	int n = 1000; //tamanho da imagem em colunas

	int** b;
	int c;
	int** count;
	int count_max = 2000;
	int** g;
	int i;
	int j;
	int jhi;
	int jlo;
	int k;
	int kill_flag;
	string filename = "mandelbrot.ppm";
	ofstream output;
	int** r;
	double wtime;
	double x_max = 1.25;
	double x_min = -2.25;
	double x;
	double x1;
	double x2;
	double y_max = 1.75;
	double y_min = -1.75;
	double y;
	double y1;
	double y2;

	int offset, rows;  //Variaveis para controle de em qual linha comeca o processamento 
						//	e quantas linhas a partir desta serão processadas por essa thread
	int my_rank;       // Identificador deste processo

	int mtype, dest, extra, averow, source;
	int proc_n; 

	MPI_Status status;


	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_n);

	printf("MPI process rank %d started", my_rank);


	b = i4pp_new(m, n);
	count = i4pp_new(m, n);
	g = i4pp_new(m, n);
	r = i4pp_new(m, n);

	if(my_rank == MASTER)
	{
		timestamp();

		double t1,t2;
		t1 = MPI_Wtime();  // inicia a contagem do tempo


		std::cout << "\n";
		std::cout << "MANDELBROT_OPENMPI\n";
		std::cout << "  C++/OpenMP version\n";
		std::cout << "\n";
		std::cout << "  Create an ASCII PPM image of the Mandelbrot set.\n";
		std::cout << "\n";
		std::cout << "  For each point C = X + i*Y\n";
		std::cout << "  with X range [" << x_min << "," << x_max << "]\n";
		std::cout << "  and  Y range [" << y_min << "," << y_max << "]\n";
		std::cout << "  carry out " << count_max << " iterations of the map\n";
		std::cout << "  Z(n+1) = Z(n)^2 + C.\n";
		std::cout << "  If the iterates stay bounded (norm less than 2)\n";
		std::cout << "  then C is taken to be a member of the set.\n";
		std::cout << "\n";
		std::cout << "  An ASCII PPM image of the set is created using\n";
		std::cout << "    M = " << m << " pixels in the X direction and\n";
		std::cout << "    N = " << n << " pixels in the Y direction.\n";




		//SEND SEND SEND

		/* Measure start time */
		double start = MPI_Wtime();

		/* Send matrix data to the worker tasks */

		int tasks = (proc_n - 1); //numero de tasks, baseado na quantidade de workers (proc_n-1)

		kill_flag = LIVE; //flag enviada para matar os workers
		averow = m/tasks; //média de linhas por número de tasjs a seren criadas
		extra = m%tasks; //resto devido a divisao
		int last_sched_offset = 0;
		offset = 0;//posição atual no saco de trabalho
		
		mtype = FROM_MASTER;


		// Primeiramente enviamos estaticamente os pacotes de trabalho iniciais
		// Eh suficiente mandar a linha inicial, ou offset, e o numero linhas 
		// a serem processadas, ou rows, para cada thread worker para que aconteca o processamento
		for (dest=1; dest < proc_n; dest++)
		{
			if(extra > 0)
			{
				rows = averow + 1;
				extra--;
			}
			else
			{
				rows = averow;

			}
			MPI_Send(&kill_flag, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);

			MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);

			MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);


			offset = offset + rows;
			last_sched_offset = offset;
		}


		//O proximo passo eh esperar os primeiros workers terminarem o processamento, para
		// que possamos envia-los novos pacotes de trabalho dinamicamente. Isso se repete
		// ate a imagem inteira ter sido escalonada, ou seja, (last_sched_offset >= m).
		// E recebemos dos workers os pacotes processados ate que nao existam mais tasks, isto eh,
		// (task_completed >= tasks)

		mtype = FROM_WORKER;
		
		int task_completed = 0;
		printf("Started receiving results");

		while(task_completed < tasks)
		{
			
			mtype = FROM_WORKER;


			//Para que nao fiquemos presos esperadando um worker com id especifico, inicialmente
			// aceitamos mensagens de  MPI_ANY_SOURCE, ou seja, de qualquer id, desde que o tipo de
			// mensagem seja FROM_WORKER
			MPI_Recv(&offset, 1, MPI_INT, MPI_ANY_SOURCE, mtype, MPI_COMM_WORLD, &status);
			source = status.MPI_SOURCE;
			// Agora recebemos o resto das informacoes do mesmo id assinalado a source
			MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&r[offset][0], rows*n, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&g[offset][0], rows*n, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&b[offset][0], rows*n, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);

			printf("Master receibed offset: %d rows: %d from process %d", offset, rows, source);


			task_completed++;

			if(last_sched_offset < m){
				
				mtype = FROM_MASTER;

				printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

				//distribui linhas do resto se necessário
				if(extra > 0)
				{
					rows = averow + 1;
					extra--;
				}
				else
				{
					rows = averow;

				}

				//envia novas informações para o worker qual enviou os dados
				MPI_Send(&kill_flag, 1, MPI_INT, source, mtype, MPI_COMM_WORLD);

				MPI_Send(&last_sched_offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD);

				MPI_Send(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD);

				last_sched_offset = last_sched_offset + rows;
			}
			
		


		}

		kill_flag = KILL;
		//Por final informamos aos workers que nao ha mais trabalho a ser feito pelo sinal de kill
		for (dest=1; dest < proc_n; dest++)
		{
		
			MPI_Send(&kill_flag, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);

		}

		t2 = MPI_Wtime(); // termina a contagem do tempo
 		printf("\nTempo de execucao: %f\n\n", t2-t1);   


		///*
		//  Write data to an ASCII PPM file.
		//*/
		output.open(filename.c_str());

		output << "P3\n";
		output << n << "  " << m << "\n";
		output << 255 << "\n";
		for (i = 0; i < m; i++)
		{
			for (jlo = 0; jlo < n; jlo = jlo + 4)
			{
				jhi = i4_min(jlo + 4, n);
				for (j = jlo; j < jhi; j++)
				{
					output << "  " << r[i][j]
						<< "  " << g[i][j]
						<< "  " << b[i][j] << "\n";
				}
				output << "\n";
			}
		}

		output.close();
		std::cout << "\n";
		std::cout << "  Graphics data written to \"" << filename << "\".\n";

		/*
		Free memory.
		*/
		i4pp_delete(b, m, n);
		i4pp_delete(count, m, n);
		i4pp_delete(g, m, n);
		i4pp_delete(r, m, n);
		/*
		Terminate.
		*/
		std::cout << "\n";
		std::cout << "MANDELBROT_OPENMP\n";
		std::cout << "  Normal end of execution.\n";
		std::cout << "\n";
		timestamp();


	}
	else
	{

		while(true)
		{

		
			mtype = FROM_MASTER;
			//recebe a flag que informa se o worker deve finalizar ou continuar
			MPI_Recv(&kill_flag, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
			
			if(kill_flag == KILL)
			{
				break;
			}
			//recebe valores do mestre para o trabalho
			
			MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);

			omp_set_num_threads(16);
			
			printf("Process %d processing offset: %d rows: %d", my_rank, offset, rows);

			# pragma omp parallel \
			shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) \
			private ( i, j, k, x, x1, x2, y, y1, y2)
			{
			# pragma omp for schedule(dynamic)
				for (i = offset; i < offset + rows; i++)
				{


					for (j = 0; j < n; j++)
					{


						x = ((double)(j - 1) * x_max
						+ (double)(m - j) * x_min)
						/ (double)(m - 1);

						y = ((double)(i - 1) * y_max
							+ (double)(n - i) * y_min)
							/ (double)(n - 1);


						count[i][j] = 0;

						x1 = x;
						y1 = y;

						for (k = 1; k <= count_max; k++)
						{
							x2 = x1 * x1 - y1 * y1 + x;
							y2 = 2 * x1 * y1 + y;

							if (x2 < -2 || 2 < x2 || y2 < -1 || 1 < y2)
							{
								count[i][j] = k;
								break;
							}
							x1 = x2;
							y1 = y2;
						}

						if ((count[i][j] % 2) == 1)
						{
							r[i][j] = 255;
							g[i][j] = 255;
							b[i][j] = 255;
						}
						else
						{
							c = (int)(255.0 * sqrt(sqrt(sqrt(
								((double)(count[i][j]) / (double)(count_max))))));
							r[i][j] = 3 * c / 5;
							g[i][j] = 3 * c / 5;
							b[i][j] = c;
						}


					

					}
				}

			//SEND
			}

			printf("Process %d finished offset: %d rows: %d", my_rank, offset, rows);

			mtype = FROM_WORKER;

			//envia os dados para o mestre incluindo informações RGB calculadas e OFFSET e ROWS para a área da imagem calculada
			MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
			MPI_Send(&r[offset][0], rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
			MPI_Send(&g[offset][0], rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
			MPI_Send(&b[offset][0], rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);

			printf("Process %d sent results offset: %d rows: %d", my_rank, offset, rows);

			
		}


	}	

	
	MPI_Finalize();


	return 0;
}
//****************************************************************************80



//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
int i4_min(int i1, int i2)
{
	int value;

	if (i1 < i2)
	{
		value = i1;
	}
	else
	{
		value = i2;
	}
	return value;
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    I4PP_DELETE frees memory associated with an I4PP.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      int **a;
//      a = i4pp_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int **A, a pointer to the pointers to the array.
//
//    Input, int M, N, the number of rows and columns in the array.
//
void i4pp_delete(int** a, int m, int n)
{
	int i;

	for (i = 0; i < m; i++)
	{
		delete[] a[i];
	}

	delete[] a;

	return;
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    I4PP_NEW allocates a new I4PP.
//
//  Discussion:
//
//    A declaration of the form
//      int **a;
//    is necesary.  Then an assignment of the form:
//      a = i4pp_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, int **I4PP_NEW, a pointer to the pointers to the array.
//
int** i4pp_new(int m, int n)
{
	int** a;
	int i;

	a = new int* [m];

	if (a == NULL)
	{
		cerr << "\n";
		cerr << "I4PP_NEW - Fatal error!\n";
		cerr << "  Unable to allocate row pointer array.\n";
		exit(1);
	}

	for (i = 0; i < m; i++)
	{
		a[i] = new int[n];
		if (a[i] == NULL)
		{
			cerr << "\n";
			cerr << "I4PP_NEW - Fatal error!\n";
			cerr << "  Unable to allocate row array.\n";
			exit(1);
		}
	}

	return a;
}
//****************************************************************************80


//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
void timestamp()
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct std::tm* tm_ptr;
	std::time_t now;

	now = std::time(NULL);
	tm_ptr = std::localtime(&now);

	std::strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr);

	std::cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}
