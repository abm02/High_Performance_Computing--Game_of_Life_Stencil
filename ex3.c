#include "mpi.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define ALIVE '1'
#define DEAD '0'

int main(int argc, char* argv[]) {
    //Initialize MPI
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int mpiRoot = 0;


    //Guarantee there exists the required input data
    int num_gen, num_rows, num_cols, num_rep;
    double probabilityAlive;
    if (rank == mpiRoot){
        if(argc != 6) {
            printf("Input data necessary: %s <number of generations> <rows of grid> <columns of grid> <probability of cell being alive [0, 1]> <number repetitions>\n", argv[0]);
            return -1;
        }
        num_gen = atoi(argv[1]);
        num_rows = atoi(argv[2]);   //nr
        num_cols = atoi(argv[3]);   //nc
        probabilityAlive = atof(argv[4]);
        num_rep = atoi(argv[5]);

        if (num_gen<2 || num_rows<1 || num_cols <1 || probabilityAlive<0 || probabilityAlive>1 || num_rep<1){
            printf("Invalid arguments...");
            return -1;
        }
    }
    
    //Broadcast the data received from the command line to all the other processes
    MPI_Bcast(&num_gen, 1, MPI_INT, mpiRoot, MPI_COMM_WORLD);
    MPI_Bcast(&num_rows, 1, MPI_INT, mpiRoot, MPI_COMM_WORLD);
    MPI_Bcast(&num_cols, 1, MPI_INT, mpiRoot, MPI_COMM_WORLD);
    MPI_Bcast(&probabilityAlive, 1, MPI_DOUBLE, mpiRoot, MPI_COMM_WORLD);
    MPI_Bcast(&num_rep, 1, MPI_INT, mpiRoot, MPI_COMM_WORLD);


    //Start doing the partitioning
    int num_rows_local = num_rows/size;
    if (rank == size-1){    //we add the remaining rows to the last process
        num_rows_local += num_rows % size;
    }
    //Add the ghost rows and umns
    int num_rows_localWithGhost = num_rows_local + 2;
    int num_columns_WithGhost = num_cols + 2;


    //Create the grids
    char** currentGrid = (char**)malloc(num_rows_localWithGhost * sizeof(char*));
    for (int i = 0; i < num_rows_localWithGhost; i++) {
        currentGrid[i] = (char*)malloc(num_columns_WithGhost * sizeof(char));
    }

    char** nextGrid = (char**)malloc(num_rows_localWithGhost * sizeof(char*));
    for (int i = 0; i < num_rows_localWithGhost; i++) {
        nextGrid[i] = (char*)malloc(num_columns_WithGhost * sizeof(char));
    }


    //Create grid to store timePast in each repetition
    double *timePast = (double*)malloc(num_rep * sizeof(double));

    //Loop for the different repetitions
    for(int r=1 ; r<=num_rep ; r++){
        //Fill the initial grids
        srand((unsigned) time(NULL)+rank+r); //changes the random number generator, so that we can achieve different sequences in each run
        for (int i=1 ; i<=num_rows_local ; i++){     //first and last rows are our ghost rows, so we start in row 1 and end in row num_rows_localWithGhost-1
            for (int j=1 ; j<=num_cols ; j++){     //first and last columns are our ghost columns, so we start in column 1 and end in column num_columns_WithGhost-1
                currentGrid[i][j] = ((rand() < probabilityAlive * RAND_MAX) ? ALIVE : DEAD);
            }
        }  
        
        //Start counting time
        MPI_Barrier(MPI_COMM_WORLD);
        double startTimelocal, endTimelocal, timePastlocal, timePastrep;
        startTimelocal = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000 

        /*//----------------------------------------------------------------------------------------------------
        //DEBUG - print current grid 
        //All the ranks that are not the root rank are going to send its data to the root rank
        //This way, we make sure that the printing order is correct
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank != mpiRoot){
            for (int i=1 ; i<=num_rows_local ; i++){
                MPI_Send(&currentGrid[i][1], num_cols, MPI_CHAR, mpiRoot, 0, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD );

        if (rank == mpiRoot){
            printf("\n----------INITIAL GRID----------\n");
            for (int i=1 ; i<=num_rows_local ; i++){
                for (int j=1 ; j<=num_cols ; j++){
                    printf("%c", currentGrid[i][j]== ALIVE ? 'A' : 'D');
                }
                printf("\n");
            }
            for (int k=1 ; k<size ; k++){
                int aux = num_rows / size;
                if (k == size-1){   //if we are in the last process, we may have more rows
                    aux += num_rows % size;
                }
                int *buff = (int*)malloc(num_cols * sizeof(int));
                for (int i = 0; i < num_cols; i++) {
                    buff[i] = 0;
                }
                for (int i=0 ; i<aux ; i++){
                    MPI_Status statusPrint;
                    MPI_Recv(&buff[0], num_cols, MPI_CHAR, k, 0, MPI_COMM_WORLD, &statusPrint);
                    for (int j=0 ; j<num_cols ; j++){
                        printf("%c", buff[j]== ALIVE ? 'A' : 'D');
                    } 
                    printf("\n");
                }
                free(buff);
            }
        }
        printf("\n");
        //----------------------------------------------------------------------------------------------------*/
      
        //Loop for the different generations
        for(int g=1 ; g<=num_gen ; g++){
            //Find neighbors of the processes
            int upProcess,downProcess;   
            upProcess = (rank == 0) ? size-1 : rank-1;
            downProcess = (rank == size-1) ? 0 : rank+1;


            //Set up ghost rows - Send information to the neighbors
            MPI_Status status1, status2;
            //Send top row to the process above and Receive top row from below
            MPI_Sendrecv(&currentGrid[1][0], num_columns_WithGhost, MPI_CHAR, upProcess, 0, 
                        &currentGrid[num_rows_local+1][0], num_columns_WithGhost, MPI_CHAR, downProcess, 0, MPI_COMM_WORLD, &status1);
            //Send bottom row to the process below
            MPI_Sendrecv(&currentGrid[num_rows_local][0], num_columns_WithGhost, MPI_CHAR, downProcess, 0, 
                        &currentGrid[0][0], num_columns_WithGhost, MPI_CHAR, upProcess, 0, MPI_COMM_WORLD, &status2);

            //Set up ghost columns
            for (int i=0 ; i<num_rows_localWithGhost ; i++){
                currentGrid[i][0] = currentGrid[i][num_cols];
                currentGrid[i][num_cols+1] = currentGrid[i][1];
            }


            //Update the grid
            for (int i=1 ; i<=num_rows_local ; i++){ //we don't want to look at the ghost rows
                for (int j=1 ; j<=num_cols ; j++){ //we don't want to look at the ghost columns
                    int numAliveNeighbors = 0;
                    for (int i1=i-1 ;  i1<= i+1 ; i1++){
                        for (int j1=j-1 ; j1<=j+1 ; j1++){
                            if (i1!=i || j1!=j){    //We make sure that we are not examinating the central cell itself, just its neighbors
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


            //Swap current and next generation grids
            for (int i=1 ; i<=num_rows_local ; i++){
                for(int j=1 ; j<=num_cols ; j++){
                    currentGrid[i][j] = nextGrid[i][j];
                }
            }

            
            /*//----------------------------------------------------------------------------------------------------
            //DEBUG - print current grid 
            //All the ranks that are not the root rank are going to send its data to the root rank
            //This way, we make sure that the printing order is correct
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank != mpiRoot){
            for (int i=1 ; i<=num_rows_local ; i++){
                MPI_Send(&currentGrid[i][1], num_cols, MPI_CHAR, mpiRoot, 0, MPI_COMM_WORLD);
            }
            }
            MPI_Barrier(MPI_COMM_WORLD );

            if (rank == mpiRoot){
                printf("\n----------REPETITION: %d | GENERATION: %d----------\n", r, g);
                for (int i=1 ; i<=num_rows_local ; i++){
                    for (int j=1 ; j<=num_cols ; j++){
                        printf("%c", currentGrid[i][j]== ALIVE ? 'A' : 'D');
                    }
                    printf("\n");
                }
                for (int k=1 ; k<size ; k++){
                    int aux = num_rows / size;
                    if (k == size-1){   //if we are in the last process, we may have more rows
                        aux += num_rows % size;
                    }
                    int *buff = (int*)malloc(num_cols * sizeof(int));
                    for (int i = 0; i < num_cols; i++) {
                        buff[i] = 0;
                    }
                    for (int i=0 ; i<aux ; i++){
                        MPI_Status statusPrint;
                        MPI_Recv(&buff[0], num_cols, MPI_CHAR, k, 0, MPI_COMM_WORLD, &statusPrint);
                        for (int j=0 ; j<num_cols ; j++){
                            printf("%c", buff[i]== ALIVE ? 'A' : 'D');
                        } 
                        printf("\n");
                    }
                    free(buff);
                }
            }
            printf("\n");
            //----------------------------------------------------------------------------------------------------*/

        }
        //Stop counting the time and calculate time of the iteration
        endTimelocal = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000
        timePastlocal = endTimelocal - startTimelocal;
        MPI_Reduce(&timePastlocal,&timePastrep,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        timePast[r-1] = timePastrep;


        //Count the number of alive and dead cells at the end of the last generation
        int numCellsAlivelocal = 0;
        int numCellsDeadlocal = 0;
        int numCellsAlive, numCellsDead;
        for(int i=1 ; i<=num_rows_local ; i++){
            for (int j=1 ; j<=num_cols ;j++){
                if (currentGrid[i][j] == ALIVE){
                    numCellsAlivelocal += 1;
                }
                else{
                    numCellsDeadlocal += 1;
                }
            }
        }
        MPI_Reduce(&numCellsAlivelocal,&numCellsAlive,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&numCellsDeadlocal,&numCellsDead,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);


        //Print the results
        if (rank == mpiRoot){
            printf("End of repetition: %d", r);
            printf("\nNumber of cells ALIVE= %d", numCellsAlive);
            printf("\nNumber of cells DEAD= %d", numCellsDead);
            printf("\nMicroseconds this repetition took = %.3f", timePast[r-1]);
            printf("\n\n");
        }
    }

    //Print out the mean of the timePast and the max and min timePast
    if (rank == mpiRoot){
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

    MPI_Barrier(MPI_COMM_WORLD );
    //Free allocated memory
    for (int i = 0; i < num_rows_localWithGhost; i++) {
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

