#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 
#include <mpi.h> 
#include <time.h>

int main(int argc, char **argv){
    //Initialize MPI
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Setup Cartesian Communicator without reordering
    int dims[2] = {0,0};
    MPI_Dims_create(size, 2, dims);
    int periods[2] = {1,1}; /*Periodicity in both dimensions*/
    int coords[2];
    MPI_Comm cart_comm_noReorder;
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm_noReorder);


    //Setup Cartesian Communicator with reordering
    MPI_Comm cart_comm_Reorder;
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm_Reorder);


    //Compare the two Cartesian Communicators created
    int result;
    MPI_Comm_compare(cart_comm_noReorder, cart_comm_Reorder, &result);

    if (result == 1){
        printf("%d The two Cartesian Communicators created are identical.\n", rank);
    }

}