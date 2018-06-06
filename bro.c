

#include <stdio.h>
#include <mpi.h>

 int main(int argc, char *argv[]){

 int i,rank,size;
 int root,count; 
 int buffer[4];

 MPI_Status status;
 MPI_Request request;
 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD,&size);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 root=0;
 count=4;

 if(rank == root){

 for(i=0; i<count; i++){

 buffer[i]=i;
 }

 }

 MPI_Bcast(buffer,count,MPI_INT,root,MPI_COMM_WORLD);
 buffer[rank] = 8;

 printf("Rank is: %d, Value at buffer[%d] is: %d \n",
 rank, rank, buffer[rank]);

 printf("\n");
 MPI_Finalize();
 return 0;

 }



