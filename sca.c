#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SIZE 4

int main(int argc, char *argv[])
{
    int rank, size;     

    char rec_buf[10];          /* buffer where the received data should be stored*/

   /* the data to be distributed*/
    char data[8] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



   /* divide the data among processes as described by sendcounts and displs */
    MPI_Scatter(&data, 2, MPI_CHAR, &rec_buf, 2, MPI_CHAR, 0, MPI_COMM_WORLD);


    /* print what each process received */
    printf("rank = %d  file = %c : \n", rank, rec_buf[0]);
    printf("rank = %d  file = %c : \n", rank, rec_buf[1]);
    printf("rank = %d  file = %c : \n", rank, rec_buf[2]);

    printf("\n");



    MPI_Finalize();


    return 0;
}



