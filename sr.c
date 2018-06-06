

#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    int PING_PONG_LIMIT = 10;

    /*  Initialize the MPI environment  */
    MPI_Init(NULL, NULL);

    /* Get the number of processes */
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* Get the rank of the process */
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    /* Get the name of the processor  */
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    /* Print off a hello world message */
    /*printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);*/

    MPI_Status status; 


    int partner_rank = (world_rank + 1) % 2;

    int count = 10;


    while(count >0){
    
    if (world_rank == 0) {
         /*Increment the ping pong counit before you send it */
        MPI_Send(&count, 1, MPI_INT, 1, 9,MPI_COMM_WORLD);
        printf("%d sent and incremented ping_pong_count %d\n", world_rank, count);
        
        count --;
        MPI_Recv(&count, 1, MPI_INT, 1, 9,MPI_COMM_WORLD, &status);
        printf("%d received ping_pong_count %d\n",world_rank, count);

    } 

    else if (world_rank == 1) {
      
        MPI_Recv(&count, 1, MPI_INT, 0, 9,MPI_COMM_WORLD, &status);
        
        
        printf("%d received ping_pong_count %d\n",world_rank, count);
        count --;
        MPI_Send(&count, 1, MPI_INT, 0, 9,MPI_COMM_WORLD);
        printf("%d sent and incremented ping_pong_count " "%d\n", world_rank, count);
        
                                                                                                           
    }
    }
    


    /* Finalize the MPI environment. */
     MPI_Finalize();

     printf(" from  %d\n",world_rank);
 

}



