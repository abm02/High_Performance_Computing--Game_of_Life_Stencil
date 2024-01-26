#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 
#include <mpi.h> 
#include <time.h>

#define ALIVE '1'
#define DEAD '0'
#define TAG 0

//Function to create a new datatype
void create_datatype(MPI_Datatype* type,int start1,int start2,int subsize1,int subsize2, int bigsize1, int bigsize2){
	const int array_of_subsizes[2] = {subsize1,subsize2};
	const int array_of_starts[2] = {start1,start2};
    const int array_of_bigsizes[2] = {bigsize1,bigsize2};

	MPI_Type_create_subarray(2,array_of_bigsizes,array_of_subsizes,array_of_starts,MPI_ORDER_C, MPI_CHAR,type);
	MPI_Type_commit(type);
}

//Function to find out the 8 neighbors 
void find_neighbours(MPI_Comm cart_comm,int rank,int num_rowsCC,int num_colsCC,int* left,int* right,int* top,int* bottom,int* topleft,int* topright,int* bottomleft,int* bottomright){
	int source,dest,disp=1;
	int my_coords[2];
	int corner_coords[2];
	int corner_rank;
	
	//Find top and bottom neighbors
	MPI_Cart_shift(cart_comm,0,disp,top,bottom);
	
	//Find left and right neighbors
	MPI_Cart_shift(cart_comm,1,disp,left,right);
	
    //Find top right corner
	MPI_Cart_coords(cart_comm,rank,2,my_coords);
	corner_coords[0] = my_coords[0] - 1; 
	corner_coords[1] = (my_coords[1]+1) % num_colsCC;
	if(corner_coords[0]<0)
		corner_coords[0] = num_rowsCC - 1;
	MPI_Cart_rank(cart_comm,corner_coords,topright);
	
	//Find top left corner
	MPI_Cart_coords(cart_comm,rank,2,my_coords);
	corner_coords[0] = my_coords[0] - 1; 
	corner_coords[1] = my_coords[1] - 1 ; 
	if(corner_coords[0]<0)
		corner_coords[0] = num_rowsCC - 1;
	if (corner_coords[1]<0)
		corner_coords[1] = num_colsCC - 1;
	MPI_Cart_rank(cart_comm,corner_coords,topleft);
	
	//Find bottom right corner
	MPI_Cart_coords(cart_comm,rank,2,my_coords);
	corner_coords[0] = (my_coords[0]+1) % num_rowsCC; 
	corner_coords[1] = (my_coords[1]+1) % num_colsCC; 
	MPI_Cart_rank(cart_comm,corner_coords,bottomright);
	
	//Find bottom left corner
	MPI_Cart_coords(cart_comm,rank,2,my_coords);
	corner_coords[0] = (my_coords[0]+1) % num_rowsCC; 
	corner_coords[1] = my_coords[1]-1 ;
	if (corner_coords[1]<0)
		corner_coords[1] = num_colsCC -1;
	MPI_Cart_rank(cart_comm,corner_coords,bottomleft);
}




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


        //Start counting time
        MPI_Barrier(cart_comm);
        double startTimelocal, endTimelocal, timePastlocal, timePastrep;
        startTimelocal = 1000000 * MPI_Wtime();       //MPI_Wtime() gives time in seconds, so to have microseconds, we multiply by 1000000 


        //Create 4 datatypes for sending
        MPI_Datatype firstcolumn_send,firstrow_send,lastcolumn_send,lastrow_send;
        create_datatype(&firstcolumn_send,1,1,num_rows_local,1, num_rows_local+2, num_cols_local+2);
        create_datatype(&firstrow_send,1,1,1,num_cols_local, num_rows_local+2, num_cols_local+2);
        create_datatype(&lastcolumn_send,1,num_cols_local,num_rows_local,1, num_rows_local+2, num_cols_local+2);
        create_datatype(&lastrow_send,num_rows_local,1,1,num_cols_local, num_rows_local+2, num_cols_local+2);

        //Create 4 datatypes for receiving
        MPI_Datatype firstcolumn_recv,firstrow_recv,lastcolumn_recv,lastrow_recv;
        create_datatype(&firstcolumn_recv,1,0,num_rows_local,1, num_rows_local+2, num_cols_local+2);
        create_datatype(&firstrow_recv,0,1,1,num_cols_local, num_rows_local+2, num_cols_local+2);
        create_datatype(&lastcolumn_recv,1,num_cols_local+1,num_rows_local,1, num_rows_local+2, num_cols_local+2);
        create_datatype(&lastrow_recv,num_rows_local+1,1,1,num_cols_local, num_rows_local+2, num_cols_local+2);


        //Find the ranks of the neighbors
        int left,right,bottom,top,topleft,topright,bottomleft,bottomright;
        find_neighbours(cart_comm,rank,num_rowsCC,num_colsCC,&left,&right,&top,&bottom,&topleft,&topright,&bottomleft,&bottomright);
        

        //Prepare exchange of necessary data between processes
        MPI_Request requests[16];
        MPI_Send_init(&(currentGrid[0][0]),1,firstcolumn_send,left,TAG,cart_comm,&requests[0]);
        MPI_Send_init(&(currentGrid[0][0]),1,firstrow_send,	top,TAG,cart_comm,&requests[1]);
        MPI_Send_init(&(currentGrid[0][0]),1,lastcolumn_send,right,TAG,cart_comm,&requests[2]);
        MPI_Send_init(&(currentGrid[0][0]),1,lastrow_send,	bottom,TAG,cart_comm,&requests[3]);
        MPI_Send_init(&(currentGrid[1][1]),1,MPI_CHAR,topleft,TAG,cart_comm,&requests[4]);
        MPI_Send_init(&(currentGrid[1][num_cols_local]),1,MPI_CHAR,topright,TAG,cart_comm,&requests[5]);
        MPI_Send_init(&(currentGrid[num_rows_local][num_cols_local]),1,MPI_CHAR,bottomright,TAG,cart_comm,&requests[6]);
        MPI_Send_init(&(currentGrid[num_rows_local][1]),1,MPI_CHAR,bottomleft,TAG,cart_comm,&requests[7]);
        
        MPI_Recv_init(&(currentGrid[0][0]),1,firstcolumn_recv,left,TAG,cart_comm,&requests[8]);
        MPI_Recv_init(&(currentGrid[0][0]),1,firstrow_recv,top,TAG,cart_comm,&requests[9]);
        MPI_Recv_init(&(currentGrid[0][0]),1,lastcolumn_recv,right,TAG,cart_comm,&requests[10]);
        MPI_Recv_init(&(currentGrid[0][0]),1,lastrow_recv,bottom,TAG,cart_comm,&requests[11]);
        MPI_Recv_init(&(currentGrid[0][0]),1,MPI_CHAR,topleft,TAG,cart_comm,&requests[12]);
        MPI_Recv_init(&(currentGrid[0][num_cols_local+1]),1,MPI_CHAR,topright,TAG,cart_comm,&requests[13]);
        MPI_Recv_init(&(currentGrid[num_rows_local+1][num_cols_local+1]),1,MPI_CHAR,bottomright,TAG,cart_comm,&requests[14]);
        MPI_Recv_init(&(currentGrid[num_rows_local+1][0]),1,MPI_CHAR,bottomleft,TAG,cart_comm,&requests[15]);


        
        //Start loop for the generations
        int numAliveNeighbors;
        for (int g=1 ; g<=num_gen ; g++){
            //Exchange the information prepared previously
            MPI_Startall(16,requests);


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

            //Make sure all requests have already been completed
            MPI_Status statuses[16];
            MPI_Waitall(16,requests,statuses);


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
        MPI_Reduce(&timePastlocal,&timePastrep,1,MPI_DOUBLE,MPI_MAX,0,cart_comm);
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
        MPI_Reduce(&numCellsAlivelocal,&numCellsAlive,1,MPI_INT,MPI_SUM,0,cart_comm);
        MPI_Reduce(&numCellsDeadlocal,&numCellsDead,1,MPI_INT,MPI_SUM,0,cart_comm);


        //Print the results
        if (rank == 0){
            printf("End of repetition: %d", r);
            printf("\nNumber of cells ALIVE= %d", numCellsAlive);
            printf("\nNumber of cells DEAD= %d", numCellsDead);
            printf("\nMicroseconds this repetition took = %.3f", timePast[r-1]);
            printf("\n\n");
        }

        //Free Allocated memory for the datatypes
        MPI_Type_free(&firstcolumn_send);
        MPI_Type_free(&firstrow_send);
        MPI_Type_free(&lastcolumn_send);
        MPI_Type_free(&lastrow_send);
        
        MPI_Type_free(&firstcolumn_recv);
        MPI_Type_free(&firstrow_recv);
        MPI_Type_free(&lastcolumn_recv);
        MPI_Type_free(&lastrow_recv);
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


    MPI_Barrier(cart_comm);
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