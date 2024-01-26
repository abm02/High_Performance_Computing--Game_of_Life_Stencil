#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <getopt.h>
#include <assert.h>

#define ALIVE '1'
#define DEAD '0'

int main(int argc, char *argv[]) {
    //Guarantee there exists the required input data
    if(argc != 6) {
        printf("Input data necessary: %s <number of generations> <rows of grid> <collumns of grid> <probability of cell being alive [0, 1]> <number repetitions>\n", argv[0]);
        return -1;
    }
    int num_gen = atoi(argv[1]);
    int num_rows = atoi(argv[2]);
    int num_cols = atoi(argv[3]);
    double probabilityAlive = atof(argv[4]);
    int num_rep = atoi(argv[5]);

    if (num_gen<2 || num_rows<1 || num_cols <1 || probabilityAlive<0 || probabilityAlive>1 || num_rep<1){
        printf("Invalid arguments...");
        return -1;
    }

    //Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //Allocate memory for the grids
    //char* currentGrid = (char*)malloc(num_rows*num_cols*sizeof(char));
    char** currentGrid = (char**)malloc(num_rows * sizeof(char*));
    for (int i = 0; i < num_rows; i++) {
        currentGrid[i] = (char*)malloc(num_cols * sizeof(char));
    }

    char** nextGrid = (char**)malloc(num_rows * sizeof(char*));
    for (int i = 0; i < num_rows; i++) {
        nextGrid[i] = (char*)malloc(num_cols * sizeof(char));
    }


    //Create grid to store timePast in each repetition
    double *timePast = (double*)malloc(num_rep * sizeof(double));


    //Loop for the different repetitions
    for(int r=1 ; r<=num_rep ; r++){
        //Initialize the grid of the current generation, with random variables
        srand((unsigned) time(NULL)+r); //changes the random number generator, so that we can achieve different sequences in each run
        for(int i=0 ; i<num_rows ; i++){
            for(int j=0 ; j<num_cols ; j++){
                currentGrid[i][j] = ((rand() < probabilityAlive * RAND_MAX) ? ALIVE : DEAD);
            }
        }

        //Start counting the time 
        double startTime, endTime;
        startTime = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000

        //Loop for the different generations
        for(int g=1 ; g<=num_gen ; g++){
            for(int i=0 ; i<num_rows ; i++){
                for(int j=0 ; j<num_cols ; j++){
                    //Count number of alive neighbours
                    int numAliveNeighbors = 0;
                    for(int x1=-1 ; x1<=1 ; x1++){
                        for(int x2=-1 ; x2<=1 ; x2++){
                            //Make sure we are not considering the center(the cell itself) nor cells outside the bound
                            if((i+x1)>=0 && (i+x1)<num_rows && (j+x2)>=0 && (j+x2)<num_cols && (x1!=0 || x2!=0)){
                                if(currentGrid[i+x1][j+x2] == ALIVE){
                                    numAliveNeighbors += 1;
                                }
                            }
                        }
                    }

                    //Decide if in the next state, the cell will be Alive or Dead -NOT WORKING YET!!!!!!!!
                    if(currentGrid[i][j] == ALIVE){
                        if(numAliveNeighbors==2 || numAliveNeighbors==3){
                            nextGrid[i][j] = ALIVE;
                        }
                        else{
                            nextGrid[i][j] = DEAD;
                        }
                    }
                    else{
                        if(numAliveNeighbors==3){
                            nextGrid[i][j] = ALIVE;
                        }
                        else{
                            nextGrid[i][j] = DEAD;
                        }
                        
                    }              
                }
            }
            //Swap current and next generation grids
            for (int i=0 ; i<num_rows ; i++){
                for(int j=0 ; j<num_cols ; j++){
                    currentGrid[i][j] = nextGrid[i][j];
                }
            }

        /*//--------------------------------------------------
        //DEBUG - print grid at the end of the generation
        printf("Repetition: %d || Generation: %d\n", r, g);
        for(int i=0 ; i<num_rows ; i++){
            for(int j=0 ; j<num_cols ; j++){
                printf("%c", nextGrid[i][j]== ALIVE ? 'A' : 'D');
            }
            printf("\n");
        }
        printf("\n");
        //--------------------------------------------------*/
        }  

    //Stop counting the time and calculate time of the iteration
    endTime = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000
    timePast[r-1] = endTime - startTime; 
 
      //Count the number of alive and dead cells at the end of the last generation
      int numCellsAlive = 0;
      int numCellsDead = 0;
      for(int i=0 ; i<num_rows ; i++){
        for (int j=0 ; j<num_cols ;j++){
            if (currentGrid[i][j] == ALIVE){
                numCellsAlive += 1;
            }
            else{
                numCellsDead += 1;
            }
        }
        
      }
      printf("End of repetition: %d", r);
      printf("\nNumber of cells ALIVE= %d", numCellsAlive);
      printf("\nNumber of cells DEAD= %d", numCellsDead);
      printf("\nMicroseconds this repetition took = %.3f", timePast[r-1]);
      printf("\n\n");
    }


    //Print out the mean of the timePast and the max and min timePast
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

    

    //Free allocated memory
    for (int i = 0; i < num_rows; i++) {
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