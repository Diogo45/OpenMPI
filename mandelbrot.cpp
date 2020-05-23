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
#define TAREFAS 7; // Numero de tarefas no saco de trabalho para np = 8, processo 0 é o mestre



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

	// if (argc < 2) {
	// 	std::cout << "Number of threads missing!" << std::endl;
	// }

	// int xSize = atoi(argv[1]);
	// int ySize = atoi(argv[2]);

	// int debug = atoi(argv[3]);


	//std::cout << "Number of threads: " << nThreads << endl;
	int m = 1000;
	int n = 1000;

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

	int offset, rows; 
	int my_rank;       // Identificador deste processo

	int mtype, dest, extra, averow, source;
	int proc_n; 

	MPI_Status status;


	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&proc_n);

	b = i4pp_new(m, n);
	count = i4pp_new(m, n);
	g = i4pp_new(m, n);
	r = i4pp_new(m, n);

	if(my_rank == MASTER)
	{
		timestamp();
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
		averow = m/(proc_n-1);
		extra = m%(proc_n-1);
		offset = 0;
		mtype = FROM_MASTER;
		for (dest=1; dest < proc_n; dest++)
		{
			rows = (dest <= extra) ? averow+1 : averow;   

			MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);

			MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);

			offset = offset + rows;
		}


		//RECV RECV RECV

		mtype = FROM_WORKER;
		for (i=1; i < proc_n; i++)
		{
			source = i;
			printf("Receiving results from task %d\n",source);


			MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);

			
			MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);

			printf("Received offset and rows from task %d\n",source);
			MPI_Recv(&r[offset][0], rows*n, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			printf("Received r from task %d\n",source);
			MPI_Recv(&g[offset][0], rows*n, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			printf("Received g from task %d\n",source);
			MPI_Recv(&b[offset][0], rows*n, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			printf("Received b from task %d\n",source);
		}


		timestamp();

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

		//RECV
		mtype = FROM_MASTER;
		MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);

		printf("Task %d receiveing %d rows and offset=%d\n",my_rank,rows,offset);



		for (i = offset; i < offset + rows; i++)
		{

			//printf("Task %d starting processing line %d\n",my_rank,i);

			/*if (debug == 1)
			{
				int tid = omp_get_thread_num();
				printf("Hello %d\n", tid);
			}*/

			for (j = 0; j < n; j++)
			{


				// if(i == offset + rows - 1 && j >= 430 && j <= 438)
				// {
				// 	printf("1. TASK %d DOING SHIT IN J=%d\n", my_rank, j);
				// }

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


				if(i == offset + rows - 1 && j >= 430 && j <= 438)
				{
					printf("1. TASK %d r=%d g=%d b=%d \n", my_rank, r[i][j], g[i][j], b[i][j]);
				}
				

			}
		}

		//SEND
		
		mtype = FROM_WORKER;


		MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
		MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
		MPI_Send(&r[offset][0], rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
		MPI_Send(&g[offset][0], rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
		MPI_Send(&b[offset][0], rows*n, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);

		printf("Task %d Sending Back results\n",my_rank);

		
	}

	
	//	wtime = omp_get_wtime();
	//	/*
	//	  Carry out the iteration for each pixel, determining COUNT.
	//	*/
	//
	//	omp_set_num_threads(nThreads);
	//
	//# pragma omp parallel \
	//  shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) \
	//  private ( i, j, k, x, x1, x2, y, y1, y2 )
	//	{
	//# pragma omp for schedule(dynamic)

	//RECV(offset)
	//RECV(nRows)



	
//}

	//wtime = omp_get_wtime() - wtime;
	//std::cout << "\n";
	//std::cout << "  Time = " << wtime << " seconds.\n";
	///*
	//  Write data to an ASCII PPM file.
	//*/
	MPI_Finalize();
}
//****************************************************************************80

int i4_min(int i1, int i2)

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

void i4pp_delete(int** a, int m, int n)

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

int** i4pp_new(int m, int n)

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

void timestamp()

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
