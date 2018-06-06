

#include <stdio.h>
#include "mpi.h"

int main(int argc,char *argv[]){

  int i, sum, sumTotal, upToVal;
  int start, end, size, rank;
  double total_my_bcast_time;  

  upToVal = 10000;

  MPI_Init(&argc,&argv);

  total_my_bcast_time -= MPI_Wtime();
 
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  start = rank*(upToVal/size) + 1;

  if(rank==(size-1)){

    end = upToVal;

  }
  else{

    end = start + (upToVal/size)-1;

  }

  sum = 0;
  sumTotal=0;

  for(i=start; i<= end; i++){

    sum = sum +i;

  }

  printf("\n before Rank: %d, sum: %d, sumTotal: %d\n", rank, sum, sumTotal);

  MPI_Reduce (&sum, &sumTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

  printf("\nafter reduce Rank: %d, sum: %d, sumTotal: %d\n", rank, sum, sumTotal);

  total_my_bcast_time += MPI_Wtime();
  printf("\n total time: %d \n", total_my_bcast_time);

  MPI_Finalize();

  return 0;

}




