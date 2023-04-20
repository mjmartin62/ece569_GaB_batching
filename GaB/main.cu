/* ###########################################################################################################################
## Organization         : The University of Arizona
##                      :
## File name            : GaB.c
## Language             : C (ANSI)
## Short description    : Gallager-B Hard decision Bit-Flipping algorithm
##                      :
##                      :
##                      :
## History              : Modified 19/01/2016, Created by Burak UNAL
##                      :
## COPYRIGHT            : burak@email.arizona.edu
## ######################################################################################################################## */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include "GaB.h"
#include <iostream>
#include <time.h>
using namespace std;

#define arrondi(x) ((ceil(x)-x)<(x-floor(x))?(int)ceil(x):(int)floor(x))
#define min(x,y) ((x)<(y)?(x):(y))
#define signf(x) ((x)>=0?0:1)
#define	max(x,y) ((x)<(y)?(y):(x))
#define SQR(A) ((A)*(A))
#define BPSK(x) (1-2*(x))
#define PI 3.1415926536

// Verificaton Function ###############################################################################
int VerificationComputeSyndrome(int *Decide,int **Mat,int *RowDegree,int M)
{
	int Synd,k,l;
	for (k=0;k<M;k++)
	{
		Synd=0;
		for (l=0;l<RowDegree[k];l++) 
      Synd=Synd^Decide[Mat[k][l]];

		if (Synd==1) {
      printf("GPU error correction failed correction process \n");
      break;
    }
	}
	return(1-Synd);
}

// Codeword generator Function (Batch Processing) #####################################################
int CodewordBatchGenerator(int **Codeword, int **Receivedword, int **MatG, int *PermG, float alpha, int rank, int N, int *U, int Batch_size, int numWords) 
{
  int k, l, n;
  int bidx;
  
  // Outer loop is for each set of CWs packed into single 1D array
  for (bidx=0; bidx<Batch_size; bidx++){
    
    // Inner loop is for each CW
    // Initialize the relative position within the packed array to 0 then index by size of CW
    for (int packLoc = 0; packLoc < numWords * N; packLoc += N) {
      // Random generation and Encoding process
      for (k=0;k<rank;k++) 
        U[k] = 0;
	    for (k=rank;k<N;k++) 
        U[k] = floor(drand48()*2);
	    for (k=rank-1;k>=0;k--) { 
        for (l=k+1;l<N;l++) 
          U[k] = U[k]^(MatG[k][l]*U[l]); 
      }
	    for (k=0;k<N;k++) 
        Codeword[bidx][PermG[k]+ packLoc] = U[k];
    }


    // Add Noise across the packed CW 1D array and assign possibly corrupted set of Codewords to set of Receivedwords
    for (n=0; n<N*numWords; n++)  
      if (drand48()<alpha) 
        Receivedword[bidx][n]=1-Codeword[bidx][n]; 
      else 
        Receivedword[bidx][n]=Codeword[bidx][n];
  }

  return 0;
}

//#####################################################################################################
int GaussianElimination_MRB(int *Perm,int **MatOut,int **Mat,int M,int N)
{
	int k,n,m,m1,buf,ind,indColumn,nb,*Index,dep,Rank;
	Index=(int *)calloc(N,sizeof(int));
	// Triangularization
	indColumn=0;nb=0;dep=0;
	for (m=0;m<M;m++)
	{
		if (indColumn==N) { dep=M-m; break; }
		for (ind=m;ind<M;ind++) { if (Mat[ind][indColumn]!=0) break; }
		// If a "1" is found on the column, permutation of rows
		if (ind<M)
		{
			for (n=indColumn;n<N;n++) { buf=Mat[m][n]; Mat[m][n]=Mat[ind][n]; Mat[ind][n]=buf; }
		// bottom of the column ==> 0
			for (m1=m+1;m1<M;m1++)
			{
				if (Mat[m1][indColumn]==1) { for (n=indColumn;n<N;n++) Mat[m1][n]=Mat[m1][n]^Mat[m][n]; }
			}
			Perm[m]=indColumn;
		}
		// else we "mark" the column.
		else { Index[nb++]=indColumn; m--; }
		indColumn++;
	}
	Rank=M-dep;
	for (n=0;n<nb;n++) Perm[Rank+n]=Index[n];
	// Permutation of the matrix
	for (m=0;m<M;m++) { for (n=0;n<N;n++) MatOut[m][n]=Mat[m][Perm[n]]; }
	// Diagonalization
	for (m=0;m<(Rank-1);m++)
	{
		for (n=m+1;n<Rank;n++)
		{
			if (MatOut[m][n]==1) { for (k=n;k<N;k++) MatOut[m][k]=MatOut[n][k]^MatOut[m][k]; }
		}
	}
	free(Index);
	return(Rank);
}

//#####################################################################################################
// Function Inputs
// Argv 1 = matrix
// argv 2 = matrix
// argv 3 = codeword batch size
// argv 4 = reserved
// argv 5 = reserved
int main(int argc, char * argv[]) 
{
  // ----------------------------------------------------
  // Assign function inputs
  // ----------------------------------------------------
  int Batch_size = std::stoi(argv[3]);
  int Block_size = std::stoi(argv[4]);
  int numWords = std::stoi(argv[5]);
  int reserved_6 = std::stoi(argv[6]);

  
  // Variables Declaration
  FILE *f;
  int Graine,NbIter,nbtestedframes,NBframes;
  float alpha_max, alpha_min,alpha_step,alpha,NbMonteCarlo;
  clock_t start_time, end_time; 
  float total_time;
  long int total_words = 0;
  // ----------------------------------------------------
  // lecture des param de la ligne de commande
  // ----------------------------------------------------
  char *FileName,*FileMatrix,*FileResult,*name;
  FileName=(char *)malloc(200);
  FileMatrix=(char *)malloc(200);
  FileResult=(char *)malloc(200);
  name=(char *)malloc(200);

  // Copy input matrices
  strcpy(FileMatrix,argv[1]); 	// Matrix file
  strcpy(FileResult,argv[2]); 	// Results file
  
  // ----------------------------------------------------
  // Simulation input for GaB BF
  // ----------------------------------------------------
  NbMonteCarlo=100000;	    // Maximum nb of codewords sent
  NbIter=100; 	            // Maximum nb of iterations
  //alpha= 0.01;              // Channel probability of error
  NBframes=100;	            // Simulation stops when NBframes in error
  Graine=1;		            // Seed Initialization for Multiple Simulations

  // brkunl
  alpha_max= 0.0600;		    //Channel Crossover Probability Max and Min
  alpha_min= 0.0400;
  alpha_step=0.0100;
  
  // ----------------------------------------------------
  // Overrides for verification and testing runs
  // ----------------------------------------------------
  
  alpha_max = 0.03;
  alpha_min = 0.01;
  alpha_step = 0.01;
  NbMonteCarlo = 100;

  // ----------------------------------------------------
  // Load Matrix
  // ----------------------------------------------------
  int *ColumnDegree,*RowDegree,**Mat;
  int M,N,m,n,k,i,j;
  
  // Read size of an M x N H matrix from a file
  strcpy(FileName,FileMatrix);
  strcat(FileName,"_size");
  f=fopen(FileName,"r");
  fscanf(f,"%d",&M);
  fscanf(f,"%d",&N);

  // Initialize 1D array to track number of "1s" by Row and Column in H matrix
  ColumnDegree=(int *)calloc(N,sizeof(int));
  RowDegree=(int *)calloc(M,sizeof(int));

  fclose(f);
  strcpy(FileName,FileMatrix);
  strcat(FileName,"_RowDegree");
  f=fopen(FileName,"r");

  for (m=0;m<M;m++) 
   fscanf(f,"%d",&RowDegree[m]);
  fclose(f);

  Mat=(int **)calloc(M,sizeof(int *));for (m=0;m<M;m++) Mat[m]=(int *)calloc(RowDegree[m],sizeof(int));
  strcpy(FileName,FileMatrix);
  f=fopen(FileName,"r");for (m=0;m<M;m++) { for (k=0;k<RowDegree[m];k++) fscanf(f,"%d",&Mat[m][k]); }fclose(f);
  for (m=0;m<M;m++) { for (k=0;k<RowDegree[m];k++) ColumnDegree[Mat[m][k]]++; }

   int RowDegMax = RowDegree[0];
  // ----------------------------------------------------
  // Build Graph
  // ----------------------------------------------------
  int NbBranch,**NtoB,*Interleaver,*ind,numColumn,numBranch;
  
  // Count number of edges in H matrix (graph)
  NbBranch=0; 
  for (m=0;m<M;m++) 
    NbBranch=NbBranch+RowDegree[m];

  // Allocate memory for intermediary H matrix array and then populate
  NtoB=(int **)calloc(N,sizeof(int *)); 
  for (n=0;n<N;n++) 
    NtoB[n]=(int *)calloc(ColumnDegree[n],sizeof(int));
  
  ind=(int *)calloc(N,sizeof(int));
  numBranch=0;
  // Sweep down the rows
  for (m=0;m<M;m++) { 
    for (k=0;k<RowDegree[m];k++) { 
      // Read the column location from the H matrix and store in intermediary
      numColumn=Mat[m][k]; 
      NtoB[numColumn][ind[numColumn]++]=numBranch++; 
      } 
  }
  free(ind);
  
  // Build the interleaver array
  Interleaver=(int *)calloc(NbBranch,sizeof(int));
  numBranch=0;
  for (n=0;n<N;n++) { 
    for (k=0;k<ColumnDegree[n];k++) 
      Interleaver[numBranch++]=NtoB[n][k]; 
  }

  // ----------------------------------------------------
  // Decoder variables and Host memory allocation
  // ----------------------------------------------------
  int **CtoV,**VtoC,**Codeword,**Receivedword,**Decide,*U,l,kk;
  int numB;
  int RowDegreeConst, ColumnDegreeConst;
  int **IsCodeword;
  int iter[Batch_size][numWords];

  RowDegreeConst = 8;
  ColumnDegreeConst = 4;
  
  CtoV=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
    CtoV[k]=(int *)calloc(numWords * NbBranch,sizeof(int));
  
  VtoC=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
    VtoC[k]=(int *)calloc(numWords *NbBranch,sizeof(int));

  Codeword=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
    Codeword[k]=(int *)calloc(numWords * N,sizeof(int));

  Receivedword=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
    Receivedword[k]=(int *)calloc(numWords * N,sizeof(int));

  Decide=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
    Decide[k]=(int *)calloc(numWords * N,sizeof(int));

  IsCodeword=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++)
    IsCodeword[k]=(int *)calloc(numWords,sizeof(int));

 
  U=(int *)calloc(N,sizeof(int));
  srand48(time(0)+Graine*31+113);


  // ----------------------------------------------------
  // Gaussian Elimination for the Encoding Matrix (Full Representation)
  // ----------------------------------------------------
  int **MatFull,**MatG,*PermG;
  int rank;
  MatG=(int **)calloc(M,sizeof(int *));for (m=0;m<M;m++) MatG[m]=(int *)calloc(N,sizeof(int));
  MatFull=(int **)calloc(M,sizeof(int *));for (m=0;m<M;m++) MatFull[m]=(int *)calloc(N,sizeof(int));
  PermG=(int *)calloc(N,sizeof(int)); for (n=0;n<N;n++) PermG[n]=n;
  for (m=0;m<M;m++) { for (k=0;k<RowDegree[m];k++) { MatFull[m][Mat[m][k]]=1; } }
  rank=GaussianElimination_MRB(PermG,MatG,MatFull,M,N);
  //for (m=0;m<N;m++) printf("%d\t",PermG[m]);printf("\n");

  // Variables for Statistics
  int nb;
  int NiterMoy,NiterMax;
  int Dmin;
  int NbTotalErrors,NbBitError;
  int NbUnDetectedErrors;
  int NbError[Batch_size][numWords];

  strcpy(FileName,FileResult);
  f=fopen(FileName,"w");
  
  fprintf(f,"-------------------------Gallager B--------------------------------------------------\n");
  fprintf(f,"alpha\t\tNbEr(BER)\t\tNbFer(FER)\t\tNbtested\t\tIterAver(Itermax)\t\tNbUndec(Dmin)\n");

  printf("-------------------------Gallager B  Parallel code--------------------------------------------------\n");
  printf("alpha\t\t\tNbEr(BER)\t\t\t\tNbFer(FER)\t\t\tNbtested\t\tIterAver(Itermax)\t\tNbUndec(Dmin)\n");
  

  // ----------------------------------------------------
  // Constant Device Memory allocations
  // ----------------------------------------------------
  int *Dev_ColumnDegree, *Dev_RowDegree, *Dev_Interleaver, *Dev_Mat;

  if (cudaMalloc((void **) &Dev_Interleaver, NbBranch * sizeof(int)) != cudaSuccess) {
    printf("malloc error for *Dev_Interleaver \n");
    return 0;
  }
  if (cudaMalloc((void **) &Dev_ColumnDegree, sizeof(int)) != cudaSuccess) {
    printf("malloc error for *Col_Degree \n");
    return 0;
  }
  if (cudaMalloc((void **) &Dev_RowDegree, sizeof(int)) != cudaSuccess) {
    printf("malloc error for *Dev_RowDegree \n");
    return 0;
  }
  if(cudaMalloc((void **) &Dev_Mat, M * RowDegMax * sizeof(int)) != cudaSuccess) {
    printf("malloc error for *Dev_Mat \n");
    return 0;
  }
  
  // ----------------------------------------------------
  // Constant Device Memory Transfers
  // ----------------------------------------------------
  //Copy interleaver to Device 
  if (cudaMemcpy(Dev_Interleaver, Interleaver, NbBranch * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
    printf("data transfer error from host to device on Dev_Interleaver\n");
    return 0;
  }

  // Copy H-Matrix to Global memory 
  // Flatten Matrix
  int *Mat_flattened;
  Mat_flattened = (int *)calloc(RowDegreeConst*M,sizeof(int));
  for (int m=0; m<M; m++){
    for (int r=0; r<RowDegreeConst; r++){
      Mat_flattened[m*RowDegreeConst + r] = Mat[m][r];
    }
  }
  cudaMemcpy(Dev_Mat, Mat_flattened, RowDegreeConst * M *sizeof(int), cudaMemcpyHostToDevice);

  // ----------------------------------------------------
  // Stream based Device Memory allocations
  // ----------------------------------------------------
  // Declare stream sizes, state array and CUDA streams
  int stream_count = Batch_size;
  int stream_state[stream_count];
  cudaStream_t stream[stream_count];
  
  // Declare stream based device memory pointers
  // Definitions:
  // Dev_Receivedword is corrupted word
  // Dev_Decide is decoded word after decoder iteration(s)
  // Dev_Syndrome is the result to check if H*CW is valid
  // Dev_VtoC is the VN to CN message array
  // DEV_CtoV is the CN to VN message array
  int *Dev_Receivedword[stream_count], *Dev_Decide[stream_count], *Dev_Syndrome[stream_count];
  int *Dev_VtoC[stream_count], *Dev_CtoV[stream_count];

  
  // Allocate memory on the GPU for each stream
  for (int m=0; m<stream_count; m++) {
    cudaMalloc((void **) &Dev_Receivedword[m], numWords * N * sizeof(int));
    cudaMalloc((void **) &Dev_Decide[m], numWords * N * sizeof(int));
    cudaMalloc((void **) &Dev_Syndrome[m], numWords * sizeof(int));
    cudaMalloc((void **) &Dev_VtoC[m], numWords * NbBranch * sizeof(int));
    cudaMalloc((void **) &Dev_CtoV[m], numWords * NbBranch * sizeof(int));
  }

  // Assign host side pinned memory where each word and initialized message arrays will be assigned to a stream
  for (int m=0; m<stream_count; m++) {
    cudaHostAlloc((void**) &Receivedword[m], numWords * N * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &Decide[m], numWords * N * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &IsCodeword[m], numWords * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &VtoC[m], numWords * NbBranch * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &CtoV[m], numWords * NbBranch * sizeof(int), cudaHostAllocDefault);
  }

  // invoke CUDA streams
  for (int m=0; m<stream_count; m++) {
    cudaStreamCreate(&stream[m]);
  }

  // Loop for different channel bit error rates
  for(alpha=alpha_max;alpha>=alpha_min;alpha-=alpha_step) {
    NiterMoy=0;NiterMax=0;
    Dmin=1e5;
    NbTotalErrors=0;NbBitError=0;
    NbUnDetectedErrors=0;

    // Main loop for max number of codeword simulations per bit error rate
    for (nb=0, nbtestedframes=0; nb<NbMonteCarlo; nb += Batch_size*numWords) {
      
      // Fill codeword pseudo buffer
      CodewordBatchGenerator(Codeword, Receivedword, MatG, PermG, alpha, rank, N, U, Batch_size, numWords);
      total_words = total_words + Batch_size*numWords;


      // Initialize stream state array (1 = not complete, 0 = complete)
      // and batch state detector
      for (int m=0; m<Batch_size; m++){
        stream_state[m] = 1;
      }
      int batch_state = 1;
      
      // Initialize iteration trackers
      int iter_batch=0;
      for (int i=0; i<Batch_size; i++) {
        for (int j=0; j<numWords; j++) {
          iter[i][j] = 0;
        }
      }

      // ##############################################################
      //                         TIMER STARTS 
      start_time = clock();
      // ##############################################################

      // Stream corrupted set of codewords from pinned host memory to device
      for (int m=0; m<stream_count; m++) {
        if(cudaMemcpyAsync(Dev_Receivedword[m], Receivedword[m], numWords * N * sizeof(int), cudaMemcpyHostToDevice, stream[m]) != cudaSuccess){
          printf("Mem tx issue \n");
        }
      }


      // Outer loop dependant upon all the batch states and max number of allowable decode iterations
      while (batch_state == 1 && iter_batch < NbIter) {

          // Loop for each stream to set off sequential kernel execution by stream
          for (int k=0; k<stream_count; k++) {
            // Check if stream is still in active state
            if (stream_state[k] == 1){

              // Update VN to CN message array
              DataPassGB <<< ceil(N*numWords/(float)Block_size), Block_size, 0, stream[k] >>> (Dev_VtoC[k], Dev_CtoV[k], Dev_Receivedword[k], Dev_Interleaver, ColumnDegreeConst, N, NbBranch, iter_batch, numWords);
              
              // Update the CN to VN message array
              CheckPassGB<<< ceil(M*numWords/(float)Block_size), Block_size, 0, stream[k] >>> (Dev_CtoV[k], Dev_VtoC[k], M, NbBranch, RowDegreeConst, numWords); 

              //  Update the VN's (VN's are stored in the Decide array)
              APP_GB  <<< ceil(N*numWords/(float)Block_size), Block_size, 0, stream[k] >>> (Dev_Decide[k], Dev_CtoV[k], Dev_Receivedword[k], Dev_Interleaver, ColumnDegreeConst, N, M, NbBranch, numWords); 
              
              // Check to see if updated codeword has been recovered
              //ComputeSyndrome <<< ceil(M*numWords/(float)Block_size), Block_size, 0, stream[k] >>> (Dev_Decide[k], Dev_Mat, Dev_RowDegree, M, Dev_Syndrome[k], numWords); 
              ComputeSyndrome <<< 1, numWords, 0, stream[k] >>> (Dev_Decide[k], Dev_Mat, RowDegreeConst, M, Dev_Syndrome[k], numWords); 
              
              // Update host side memory for host controller error correction convergence check
              cudaMemcpyAsync(IsCodeword[k], Dev_Syndrome[k],  numWords * sizeof(int), cudaMemcpyDeviceToHost, stream[k]);

            }
          }

          // Sync Host and Device to ensure host side convergence checks are valid
          cudaDeviceSynchronize();

          // Check for codeword recovery and update stream states as neccissary
          for (int m=0; m<stream_count; m++) {
            if (stream_state[m] == 1) {
              
              
              // TMP CODE:  Verification !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
              /*
              printf("Checking iteration %d \n",iter_batch);
              for (int b=0; b<numWords; b++){
                printf("Codework check value for CW # %d in concatenated array =  %d  for Stream %d \n",b,IsCodeword[m][b],m);
              }
              */


              // Determine if batch of CWs in concatenated array are all valid
              // Perform Multi CW concatenated array check
              // Update iteration tracker for each word
              int concat_check = 1;
              for (int k=0; k<numWords; k++) {
                if (IsCodeword[m][k] == 1 && iter[m][k] == 0) {
                    iter[m][k] = iter_batch;
                }
                concat_check = concat_check && IsCodeword[m][k];
              }
              // Update stream state array and pull corrected concatenated array back to host
              if (concat_check == 1) {
                stream_state[m] = 0;
                cudaMemcpyAsync(Decide[m], Dev_Decide[m], numWords * N * sizeof(int), cudaMemcpyDeviceToHost, stream[m]);
              }
            }
          }

          // Update stream batch tracking variables for next round of streams
          iter_batch++;
          int state_check = 0;
          batch_state = 0;
          // Check to see if all streams converged to decoded codeword
          for (int m=0; m<stream_count; m++) {
            state_check ^= stream_state[m];
            if (state_check == 1) {
              batch_state = 1;
              break; 
            }
          }    
      }
      
      // Reconcile residual streams associated with codewords that were not recovered
      // and set the number of decoder iterations to the max allowable iterations
      for (int m=0; m<stream_count; m++) {
        if (stream_state[m] == 1) {
          cudaMemcpyAsync(Decide[m], Dev_Decide[m], numWords * N * sizeof(int), cudaMemcpyDeviceToHost, stream[m]);
          for (int k=0; k<numWords; k++) {
              if (IsCodeword[m][k] == 0) {
                  iter[m][k] = iter_batch;
              }
          }
        }
      }

      // Sync Host and Device to ensure error corrected words are available for stats processing
      cudaDeviceSynchronize();

      // ##############################################################
      //                         TIMER ENDS 
      end_time = clock();
      total_time = total_time + (float)(end_time - start_time) / CLOCKS_PER_SEC;
      // ##############################################################



      //============================================================================
  	  // Verification:  Uncomment for short runs
	    //============================================================================
      // Run H*CW syndrome check
      // Outer loop is for checking each stream
      // Inner loop is for checking the packed CWs within each stream
      /*
      int *parsedDecideCW;
      parsedDecideCW=(int *)calloc(1296,sizeof(int));
      for (int countBatch=0; countBatch<stream_count; countBatch++) {
          for (int countCW=0; countCW<numWords; countCW++){
            
            // Parse out CW and Error Corrected CW
            for (int countBit=0; countBit<1296; countBit++) {
              parsedDecideCW[countBit] = Decide[countBatch][countBit + countCW*N];
            }
            VerificationComputeSyndrome(parsedDecideCW,Mat,RowDegree,M);
            
            // Output uncorrupted codewords to file for verification purposes
            
            FILE *fptr1;
            fptr1 = (fopen("codewords_test_verification_pre_corrupt_02.txt", "a+"));
            // send codeword to output file
            for (k=0;k<N;k++) 
              fprintf(fptr1, "%u %s", Codeword[countBatch][k + countCW*N], "");
            fprintf(fptr1, "%s", "\n");
            fprintf(fptr1, "%s", "\n");
            fprintf(fptr1, "%s", "\n");
            fclose(fptr1);
            

            // Output corrected codewords to file for verification purposes
            FILE *fptr2;
            fptr2 = (fopen("codewords_test_verification_corrected_02.txt", "a+"));
            // send codeword to output file
            for (k=0;k<N;k++) 
              fprintf(fptr2, "%u %s", parsedDecideCW[k], "");
            fprintf(fptr2, "%s", "\n");
            fprintf(fptr2, "%s", "\n");
            fprintf(fptr2, "%s", "\n");
            fclose(fptr2);

            // Output message arrays to file for verification purposes
            FILE *fptr3;
            fptr3 = (fopen("CtoV_test_verification_02.txt", "a+"));
            // send codeword to output file
            for (k=0;k<NbBranch;k++) 
              fprintf(fptr3, "%u %s", CtoV[countBatch][k + countCW*NbBranch], "");
            fprintf(fptr3, "%s", "\n");
            fprintf(fptr3, "%s", "\n");
            fprintf(fptr3, "%s", "\n");
            fclose(fptr3);
          }
      }
      */
      
	    //============================================================================
  	  // Batch Compute Statistics
	    //============================================================================
      
      // update number of tested messages (aka frames)
      nbtestedframes += Batch_size*numWords;

      // Update total number of bit errors
      // Outer loop sweeps across the batch, inner loop calculates bit errors per message
	    // Initialize number of bit errors per codeword to zero
      for (int i=0; i<Batch_size; i++) {
        for (int j=0; j<numWords; j++) {
          NbError[i][j]=0;
        }
      }
      // bit error calcs
      for (int i=0; i<Batch_size; i++) {
        for (int j=0; j<numWords; j++) {
          for (int k=0; k<N; k++) {
            if (Decide[i][j+k] != Codeword[i][j+k]) 
              NbError[i][j]++;
          }
          NbBitError = NbBitError + NbError[i][j];
        }
      }
	    
      
      // Case Divergence
      for (int i=0; i<Batch_size; i++) {
        for (int j=0; j<numWords; j++) {
	        if (!(IsCodeword[i][j])) {
		        NiterMoy=NiterMoy+NbIter;
		        NbTotalErrors++;
          }
        }
	    }

	    // Update case divergence stats
      for (int i=0; i<Batch_size; i++) {
        for (int j=0; j<numWords; j++) {
          // Case Convergence to Right Codeword
          if ((IsCodeword[i][j]) && (NbError[i][j]==0)) { 
            NiterMax = max(NiterMax,iter[i][j]+1); 
            NiterMoy = NiterMoy+(iter[i][j]+1); 
          }
          // Case Convergence to Wrong Codeword
          if ((IsCodeword[i][j]) && (NbError[i][j]!=0)) {
            NiterMax = max(NiterMax,iter[i][j]+1); 
            NiterMoy = NiterMoy + (iter[i][j]+1);
            NbTotalErrors++; 
            NbUnDetectedErrors++;
            Dmin = min(Dmin,NbError[i][j]);
          }
        }
      }
      
      

      // Stopping Criterion 
	    if (NbTotalErrors >= NBframes) 
        break;

    }
    
    // Print final statistics for each alpha setting
    
    printf("%1.5f\t\t",alpha);
    printf("%10d (%1.16f)\t\t",NbBitError,(float)NbBitError/N/nbtestedframes);
    printf("%4d (%1.16f)\t\t",NbTotalErrors,(float)NbTotalErrors/nbtestedframes);
    printf("%10d\t\t",nbtestedframes);
    printf("%1.2f(%d)\t\t",(float)NiterMoy/nbtestedframes,NiterMax);
    printf("%d(%d)\n",NbUnDetectedErrors,Dmin);

    fprintf(f,"%1.5f\t\t",alpha);
    fprintf(f,"%10d (%1.8f)\t\t",NbBitError,(float)NbBitError/N/nbtestedframes);
    fprintf(f,"%4d (%1.8f)\t\t",NbTotalErrors,(float)NbTotalErrors/nbtestedframes);
    fprintf(f,"%10d\t\t",nbtestedframes);
    fprintf(f,"%1.2f(%d)\t\t",(float)NiterMoy/nbtestedframes,NiterMax);
    fprintf(f,"%d(%d)\n",NbUnDetectedErrors,Dmin);
    
    
  }

  //  Clean Memory
  for (int m=0; m<stream_count; m++) 
  {
    // Free pinned memory
    cudaFreeHost(Receivedword[m]);
    cudaFreeHost(Decide[m]);
    cudaFreeHost(IsCodeword[m]);
    cudaFreeHost(VtoC[m]);
    cudaFreeHost(CtoV[m]);

    // Free device memory
    cudaFree(Dev_Receivedword[m]);
    cudaFree(Dev_Decide[m]);
    cudaFree(Dev_Syndrome[m]);
    cudaFree(Dev_VtoC[m]);
    cudaFree(Dev_CtoV[m]);
  }

  cudaFree(Dev_Interleaver);
  cudaFree(Dev_ColumnDegree);
  cudaFree(Dev_RowDegree);
  cudaFree(Dev_Mat);

  free(ColumnDegree);
  free(RowDegree);
  free(FileName);
  free(FileMatrix);
  free(FileResult);
  free(name);
  free(Interleaver);
  free(NtoB);
  free(Mat);
  free(CtoV);
  free(VtoC);
  free(Codeword);
  free(Receivedword);
  free(Decide);
  free(IsCodeword);
  free(U);
  free(MatG);
  free(PermG);
  free(MatFull);
  free(Mat_flattened);
  
  fclose(f);
  
  // Update stats file for optimization data collection
  float minutesPerMillionWords = (total_time / 60) / ((float)total_words / 1000000);
  FILE *fptr5;
  fptr5 = (fopen("sweeping_studies_data_collection.txt", "a+"));
  fprintf(fptr5, "%d %d %d %f %d %f \n", Batch_size, Block_size, numWords, total_time, total_words, minutesPerMillionWords);
  fclose(fptr5);

  return(0);

}