#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 
#include <mpi.h> 
#include <time.h>

#define ALIVE '1'
#define DEAD '0'


int main(int argc, char **argv){
    //Initialize MPI
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Guarantee there exists the required input data
    int num_gen, num_rows, num_cols, num_rep;
    double probabilityAlive;
    if(argc != 6) {
        printf("Input data necessary: %s <number of generations> <rows of grid> <columns of grid> <probability of cell being alive [0, 1]> <number repetitions>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }
    num_gen = atoi(argv[1]);
    num_rows = atoi(argv[2]);   //nr
    num_cols = atoi(argv[3]);   //nc
    probabilityAlive = atof(argv[4]);
    num_rep = atoi(argv[5]);

    if (num_gen<2 || num_rows<1 || num_cols <1 || probabilityAlive<0 || probabilityAlive>1 || num_rep<1){
        printf("Invalid arguments...");
        MPI_Finalize();
        return -1;
    }


    //Setup Cartesian Communicator
    int dims[2] = {0,0};
    MPI_Dims_create(size, 2, dims);
    int periods[2] = {1,1}; /*Periodicity in both dimensions*/
    int coords[2];
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    const int num_rowsCC = dims[0]; 
    const int num_colsCC = dims[1]; 


    //Calculate rows,cols,displacement,extent for each process
    int *num_rows_proc; /* Number of rows for the i-th process [num_rows_local]*/
    int *num_cols_proc; /* Number of columns for the i-th process [num_cols_local]*/
    if (rank == 0){
		num_rows_proc = (int *) malloc(size * sizeof(int));
		num_cols_proc = (int *) malloc(size * sizeof(int));

        for (int i=0; i<size; i++) {
            num_rows_proc[i] = num_rows/num_rowsCC;
            num_cols_proc[i] = num_cols/num_colsCC;
        }
        for (int i=0; i<(num_rows%num_rowsCC); i++) {
            for (int j=0; j<num_colsCC; j++) {
                num_rows_proc[i*num_colsCC+j]++;
            }
        }
        for (int i=0; i<(num_cols%num_colsCC); i++) {
            for (int j=0; j<num_rowsCC; j++) {
                num_cols_proc[i+num_rowsCC*j]++;
            }
        }
	}


    //Scatter dimensions,displacement,extent of each process
    int num_rows_local;
    int num_cols_local;
	MPI_Scatter(num_rows_proc,1,MPI_INT,&num_rows_local,1,MPI_INT,0,cart_comm);
	MPI_Scatter(num_cols_proc,1,MPI_INT,&num_cols_local,1,MPI_INT,0,cart_comm);

	if(rank == 0){
        free(num_rows_proc);
        free(num_cols_proc);
    }


    //Allocate memory for the grids in each processs
    char** currentGrid = (char**)malloc((num_rows_local+2) * sizeof(char*));
    for (int i = 0; i < (num_rows_local+2); i++) {
        currentGrid[i] = (char*)malloc((num_cols_local+2) * sizeof(char));
    }

    char** nextGrid = (char**)malloc((num_rows_local+2) * sizeof(char*));
    for (int i = 0; i < (num_rows_local+2); i++) {
        nextGrid[i] = (char*)malloc((num_cols_local+2) * sizeof(char));
    }
	

    //Create grid to store timePast in each repetition
    double *timePast = (double*)malloc(num_rep * sizeof(double));


    //Start loop for the repetitions
    for (int r=1 ; r<=num_rep ; r++){
        //Initialize the sub-grid in each process
        srand((unsigned) time(NULL)+rank); //changes the random number generator, so that we can achieve different sequences in each run
            for (int i=1 ; i<=num_rows_local ; i++){     //first and last rows are our ghost rows, so we start in row 1 and end in row num_rows_localWithGhost-1
                for (int j=1 ; j<=num_cols_local ; j++){     //first and last columns are our ghost columns, so we start in column 1 and end in column num_columns_WithGhost-1
                    currentGrid[i][j] = ((rand() < probabilityAlive * RAND_MAX) ? ALIVE : DEAD);
                }
            }  

        int sources[8], destinations[8];
        int target[] = {0,1, 0,-1, -1,0, 1,0, -1,1, 1,1, 1,-1, -1,1};

        int vector[2];
        for (int i=0 ; i<8 ; i++){
            MPI_Cart_coords(cart_comm, rank, 2, vector);
            vector[0] += target[2*i];
            vector[1] += target[2*i+1];
            MPI_Cart_rank(cart_comm, vector, &destinations[i]);
            MPI_Cart_coords(cart_comm, rank, 2, vector);
            vector[0] -= target[2*i];
            vector[1] -= target[2*i+1];
            MPI_Cart_rank(cart_comm, vector, &sources[i]);
        }

        MPI_Comm new_cart_comm;
        MPI_Dist_graph_create_adjacent(cart_comm, 8, sources, MPI_UNWEIGHTED, 8, destinations, MPI_UNWEIGHTED, MPI_INFO_NULL, 1, &new_cart_comm);


        //Start counting time
        MPI_Barrier(new_cart_comm);
        double startTimelocal, endTimelocal, timePastlocal, timePastrep;
        startTimelocal = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000 

        
        //Create datatypes for the grids borders
        MPI_Datatype type_row, type_column, type_corner;
        int sendcount[8], recvcount[8];
        MPI_Aint senddisp[8], recvdisp[8];
        MPI_Datatype sendtype[8], recvtype[8];
        //Calculate values for the first column
        senddisp[0] = 1 * (num_cols_local+2) + num_rows_local;
        sendtype[0] = type_column;
        recvdisp[0] = 1 * (num_cols_local+2) + num_rows_local + 1;
        recvtype[0] = type_column;
        //Calculate values for the last column
        senddisp[1] = 1 * (num_cols_local+2) + num_rows_local + num_cols_local;
        sendtype[1] = type_column;
        recvdisp[1] = 1 * (num_cols_local+2) + num_rows_local + num_cols_local + 1;
        recvtype[1] = type_column;
        //Calculate values for the first row
        senddisp[2] = 1 * (num_cols_local+2) + 1;
        sendtype[2] = type_row;
        recvdisp[2] = 1;
        recvtype[2] = type_row;
        //Calculate values for the last row
        senddisp[3] = (num_rows_local+1) * (num_cols_local+2) + 1;
        sendtype[3] = type_row;
        recvdisp[3] = (num_rows_local+1) * (num_cols_local+2);
        recvtype[3] = type_row;
        //Calculate values for the top left corner
        senddisp[4] = num_rows_local * (num_cols_local+2) + 1;
        sendtype[4] = type_corner;
        recvdisp[4] = (num_rows_local+1) * (num_cols_local+2);
        recvtype[4] = type_corner;
        //Calculate values for the top right corner
        senddisp[5] = num_rows_local * (num_cols_local+2) + num_cols_local;
        sendtype[5] = type_corner;
        recvdisp[5] = (num_rows_local+1) * (num_cols_local+2) + num_cols_local + 1;
        recvtype[5] = type_corner;
        //Calculate values for the bottom left corner
        senddisp[6] = (num_rows_local+1) * (num_cols_local+2) + 1;
        sendtype[6] = type_corner;
        recvdisp[6] = num_rows_local * (num_cols_local+2);
        recvtype[6] = type_corner;
        //Calculate values for the bottom right corner
        senddisp[7] = (num_rows_local+1) * (num_cols_local+2) + num_cols_local;
        sendtype[7] = type_corner;
        recvdisp[7] = num_rows_local * (num_cols_local+2) + num_cols_local + 1;
        recvtype[7] = type_corner;

        for (int i=0 ; i<8 ; i++){
            sendcount[i] = 1;
            recvcount[i] = 1;
            senddisp[i] *= sizeof(char);
            recvdisp[i] *= sizeof(char);
        }
        

        //Start loop for the generations
        int numAliveNeighbors;
        for (int g=1 ; g<=num_gen ; g++){
            //Exchange the information prepared previously
            MPI_Neighbor_alltoallw(currentGrid, sendcount, senddisp, sendtype, currentGrid, recvcount, recvdisp, recvtype, new_cart_comm);


            //Calculate the inner grid
            for(int i=2 ; i<=(num_rows_local-1) ; i++){
                for(int j=2 ; j<=(num_cols_local-1) ; j++){
                    //Find the number of Alive neighbors
                    numAliveNeighbors = 0;
                    for(int i1=i-1 ; i1<=i+1 ; i1++){
                        for(int j1=j-1 ; j1<=j+1 ; j1++){
                            if (i1!=i || j1!=j){    
                                if (currentGrid[i1][j1] == ALIVE){  //if this neighbor is alive
                                    numAliveNeighbors = numAliveNeighbors + 1;
                                }
                            }
                        }
                    }
                    

                    //Decide if the cell will be ALive or Dead in the next generation
                    if (currentGrid[i][j]==ALIVE && (numAliveNeighbors==2 || numAliveNeighbors==3)){
                        nextGrid[i][j] = ALIVE;
                    }
                    else if (currentGrid[i][j]==DEAD && numAliveNeighbors==3){
                        nextGrid[i][j] = ALIVE;
                    }
                    else{
                        nextGrid[i][j] = DEAD;
                    }
                    
                }
            }

            //Calculate the outer grid
            for(int i=1 ; i<=num_rows_local ; i++){
                for(int j=1 ; j<=num_cols_local ; j++){
                    if(i==1 || i==num_rows_local || j==1 || j==num_cols_local){
                        //Find the number of Alive neighbors
                        numAliveNeighbors = 0;
                        for(int i1=i-1 ; i1<=i+1 ; i1++){
                            for(int j1=j-1 ; j1<=j+1 ; j1++){
                                if (i1!=i || j1!=j){    
                                    if (currentGrid[i1][j1] == ALIVE){  //if this neighbor is alive
                                        numAliveNeighbors = numAliveNeighbors + 1;
                                    }
                                }
                            }
                        }
                        

                        //Decide if the cell will be ALive or Dead in the next generation
                        if (currentGrid[i][j]==ALIVE && (numAliveNeighbors==2 || numAliveNeighbors==3)){
                            nextGrid[i][j] = ALIVE;
                        }
                        else if (currentGrid[i][j]==DEAD && numAliveNeighbors==3){
                            nextGrid[i][j] = ALIVE;
                        }
                        else{
                            nextGrid[i][j] = DEAD;
                        }
                        }
                }
            }
        
            //Swap current and next generation grids
            for (int i=1 ; i<=num_rows_local ; i++){
                for(int j=1 ; j<=num_cols_local ; j++){
                    currentGrid[i][j] = nextGrid[i][j];
                }
            }
        }

        //Stop counting the time and calculate time of the iteration
        endTimelocal = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000
        timePastlocal = endTimelocal - startTimelocal;
        MPI_Reduce(&timePastlocal,&timePastrep,1,MPI_DOUBLE,MPI_MAX,0,new_cart_comm);
        timePast[r-1] = timePastrep;


        //Count the number of alive and dead cells at the end of the last generation
        int numCellsAlivelocal = 0;
        int numCellsDeadlocal = 0;
        int numCellsAlive, numCellsDead;
        for(int i=1 ; i<=num_rows_local ; i++){
            for (int j=1 ; j<=num_cols_local ;j++){
                if (currentGrid[i][j] == ALIVE){
                    numCellsAlivelocal += 1;
                }
                else{
                    numCellsDeadlocal += 1;
                }
            }
        }
        MPI_Reduce(&numCellsAlivelocal,&numCellsAlive,1,MPI_INT,MPI_SUM,0,new_cart_comm);
        MPI_Reduce(&numCellsDeadlocal,&numCellsDead,1,MPI_INT,MPI_SUM,0,new_cart_comm);


        //Print the results
        if (rank == 0){
            printf("End of repetition: %d", r);
            printf("\nNumber of cells ALIVE= %d", numCellsAlive);
            printf("\nNumber of cells DEAD= %d", numCellsDead);
            printf("\nMicroseconds this repetition took = %.3f", timePast[r-1]);
            printf("\n\n");
        }
    }


    //Print out the mean of the timePast and the max and min timePast
    if (rank == 0){
        double sum = 0;
        double max = timePast[0];
        double min = timePast[0];
        for (int i=0 ; i<num_rep ; i++){
            sum += timePast[i];
            if(timePast[i]>max){
                max = timePast[i];
            }
            if(timePast[i]<min){
                min = timePast[i];
            }
        }
        double mean = sum/num_rep;
        printf("\nMean of the Time Past = %.3f microseconds\n", mean);
        printf("Maximum amount of Time Past in one iteration = %.3f microseconds\n", max);
        printf("Minimum amount of Time Past in one iteration = %.3f microseconds\n", min);
    }


    
    //Free allocated memory
    for (int i = 0; i < (num_rows_local+2); i++) {
        free(currentGrid[i]);
        free(nextGrid[i]);
    }
    free(currentGrid);
    free(nextGrid);
    free(timePast);
    

    //Finalize MPI and program
    MPI_Finalize();
    return 0;
}