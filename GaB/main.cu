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
		for (l=0;l<RowDegree[k];l++) Synd=Synd^Decide[Mat[k][l]];
		if (Synd==1) {
      printf("GPU Decoder failed decode process \n");
      break;
    }
	}
	return(1-Synd);
}

// Codeword generator Function (Batch Processing) #####################################################
int CodewordBatchGenerator(int **Codeword, int **Receivedword, int **MatG, int *PermG, float alpha, int rank, int N, int *U, int Batch_size) 
{
  int k, l, n;
  int bidx;
  
  for (bidx=0; bidx<Batch_size; bidx++){
    // Random generation and Encoding process
    for (k=0;k<rank;k++) 
      U[k]=0;
	  for (k=rank;k<N;k++) 
      U[k]=floor(drand48()*2);
	  for (k=rank-1;k>=0;k--) { 
      for (l=k+1;l<N;l++) 
        U[k]=U[k]^(MatG[k][l]*U[l]); 
    }
	  for (k=0;k<N;k++) 
      Codeword[bidx][PermG[k]]=U[k];

    // Add Noise and assign possibly corrupted Codeword to Receivedword
    for (n=0;n<N;n++)  
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
  int reserved_4 = std::stoi(argv[4]);
  int reserved_5 = std::stoi(argv[5]);
  int reserved_6 = std::stoi(argv[6]);
  /*
  printf("----------------------------------------\n");
  printf("Block size test: ");
  printf("%6d\n", block_size);
  */
  
  
  
  // Variables Declaration
  FILE *f;
  int Graine,NbIter,nbtestedframes,NBframes;
  float alpha_max, alpha_min,alpha_step,alpha,NbMonteCarlo;
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
  NbIter=200; 	            // Maximum nb of iterations
  alpha= 0.01;              // Channel probability of error
  NBframes=100;	            // Simulation stops when NBframes in error
  Graine=1;		            // Seed Initialization for Multiple Simulations

  // brkunl
  alpha_max= 0.0600;		    //Channel Crossover Probability Max and Min
  alpha_min= 0.0400;
  alpha_step=0.0100;
  
  // ----------------------------------------------------
  // Overrides for verification and testing runs
  // ----------------------------------------------------
  alpha= 0.04; 
  alpha_max = 0.04;
  alpha_min = 0.04;
  alpha_step = 0.04;
  NbMonteCarlo=100;
  
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
  int iter,numB;
  
  CtoV=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
      CtoV[k]=(int *)calloc(NbBranch,sizeof(int));
  
  VtoC=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
      VtoC[k]=(int *)calloc(NbBranch,sizeof(int));

  // debug code
  /*
  for (int dd=0; dd<NbBranch; dd++) {
    VtoC[0][dd] = 2;
  }
  */

  Codeword=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
        Codeword[k]=(int *)calloc(N,sizeof(int));

  Receivedword=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
        Receivedword[k]=(int *)calloc(N,sizeof(int));

  Decide=(int **)calloc(Batch_size,sizeof(int *));
  for (k=0;k<Batch_size;k++) 
        Decide[k]=(int *)calloc(N,sizeof(int));
 
  
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
  int NbUnDetectedErrors,NbError;
  int *energy;
  energy=(int *)calloc(N,sizeof(int));

  strcpy(FileName,FileResult);
  f=fopen(FileName,"w");
  /*
  fprintf(f,"-------------------------Gallager B--------------------------------------------------\n");
  fprintf(f,"alpha\t\tNbEr(BER)\t\tNbFer(FER)\t\tNbtested\t\tIterAver(Itermax)\t\tNbUndec(Dmin)\n");

  printf("-------------------------Gallager B  Parallel code--------------------------------------------------\n");
  printf("alpha\t\t\tNbEr(BER)\t\tNbFer(FER)\t\tNbtested\t\tIterAver(Itermax)\t\tNbUndec(Dmin)\n");
  */

  // ----------------------------------------------------
  // Constant Device Memory allocations
  // ----------------------------------------------------
  int *Dev_ColumnDegree, *Dev_RowDegree, *Dev_Interleaver, *Dev_Mat;

  if (cudaMalloc((void **) &Dev_Interleaver, NbBranch * sizeof(int)) != cudaSuccess) {
    printf("malloc error for *Dev_Interleaver \n");
    return 0;
  }
  if (cudaMalloc((void **) &Dev_ColumnDegree, N * sizeof(int)) != cudaSuccess) {
    printf("malloc error for *Dev_Interleaver \n");
    return 0;
  }
  if (cudaMalloc((void **) &Dev_RowDegree, M * sizeof(int)) != cudaSuccess) {
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
  //Copy column degree and row degree to Device 
  if (cudaMemcpy(Dev_ColumnDegree, ColumnDegree, N * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
    printf("data transfer error from host to device on Dev_Interleaver\n");
    return 0;
  }
  if (cudaMemcpy(Dev_RowDegree, RowDegree, M * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
    printf("data transfer error from host to device on Dev_Interleaver\n");
    return 0;
  }
  //Copy H-Matrix to Global memory 
  for (int i = 0; i < M ; i++) {
    if (cudaMemcpy((Dev_Mat+RowDegMax*i), *(Mat+i), RowDegMax * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
    printf("data transfer error from host to device on Dev_Mat\n");
    return 0;
    }
  }

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
  int *IsCodeword[stream_count];
  
  // Allocate memory on the GPU for each stream
  for (int m=0; m<stream_count; m++) {
    cudaMalloc((void **) &Dev_Receivedword[m], N * sizeof(int));
    cudaMalloc((void **) &Dev_Decide[m], N * sizeof(int));
    cudaMalloc((void **) &Dev_Syndrome[m], sizeof(int));
    cudaMalloc((void **) &Dev_VtoC[m], NbBranch * sizeof(int));
    cudaMalloc((void **) &Dev_CtoV[m], NbBranch * sizeof(int));
  }

  // Assign host side pinned memory where each word and initialized message arrays will be assigned to a stream
  for (int m=0; m<stream_count; m++) {

    cudaHostAlloc((void**) &Receivedword[m], N * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &Decide[m], N * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &IsCodeword[m], sizeof(int), cudaHostAllocDefault);
    
    // MIGHT BE ABLE TO INITIALIZE THESE ONCE SINCE THEY ARE THE SAME FOR EACH BATCH LAUNCH!!!!!!!!!!!!!!!!!!!!!
    cudaHostAlloc((void**) &VtoC[m], NbBranch * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc((void**) &CtoV[m], NbBranch * sizeof(int), cudaHostAllocDefault);
  }


  // Loop for different channel bit error rates
  for(alpha=alpha_max;alpha>=alpha_min;alpha-=alpha_step) {
    NiterMoy=0;NiterMax=0;
    Dmin=1e5;
    NbTotalErrors=0;NbBitError=0;
    NbUnDetectedErrors=0;NbError=0;

    // Main loop for max number of codeword simulations per bit error rate
    for (nb=0, nbtestedframes=0; nb<NbMonteCarlo; nb += Batch_size) {
      
      // Fill codeword pseudo buffer
      CodewordBatchGenerator(Codeword, Receivedword, MatG, PermG, alpha, rank, N, U, Batch_size);

      // Initialize stream state array (1 = not complete, 0 = complete)
      // and batch state detector
      for (int m=0; m<Batch_size; m++){
        stream_state[m] = 1;
      }
      int batch_state = 1;
      int iter=0;

      // invoke streams
      for (int m=0; m<stream_count; m++) {
        cudaStreamCreate(&stream[m]);
      }

      // Stream corrupted codeword and message arrays from pinned host memory to device
      for (int m=0; m<stream_count; m++) {
        
        

        if(cudaMemcpyAsync(Dev_Receivedword[m], Receivedword[m], N * sizeof(int), cudaMemcpyHostToDevice, stream[m]) != cudaSuccess){
          printf("Mem tx issue \n");
        }
        
        
        // I'm commenting the following out since syndrome does not need to be initialized with any particular value for kernel
        //cudaMemcpyAsync(Dev_Decide[m], Decide[m], N * sizeof(int), cudaMemcpyHostToDevice, stream[m]);
        //cudaMemcpyAsync(Dev_Syndrome[m], Syndrome[m], sizeof(int), cudaMemcpyHostToDevice, stream[m]);

        //cudaMemcpyAsync(Dev_VtoC[m], VtoC[m], NbBranch * sizeof(int), cudaMemcpyHostToDevice, stream[m]);
        cudaMemcpyAsync(Dev_CtoV[m], CtoV[m], NbBranch * sizeof(int), cudaMemcpyHostToDevice, stream[m]);
      }

            // Sync Host and Device
            cudaDeviceSynchronize();

      // Outer loop dependant upon all the batch states and max number of allowable decode iterations
      while (batch_state == 1 && iter < NbIter) {

          // Loop for each stream to set off sequential kernel execution by stream
          for (int k=0; k<stream_count; k++) {
            // Check if stream is still in active state
            if (stream_state[k] == 1){
              
              //debug code
              /*
              printf("Receivedword[0] is: %d \n",Receivedword[k][0]);
              printf("Receivedword[1] is: %d \n",Receivedword[k][1]);
              printf("Receivedword[2] is: %d \n",Receivedword[k][2]);
              */


              
              // Update VN to CN message array
              DataPassGB <<< ceil(N/32.0), 32, 0, stream[k] >>> (Dev_VtoC[k], Dev_CtoV[k], Dev_Receivedword[k], Dev_Interleaver, Dev_ColumnDegree, N, NbBranch, iter);  

              // debug code
              cudaMemcpyAsync(VtoC[k],Dev_VtoC[k], NbBranch * sizeof(int), cudaMemcpyDeviceToHost, stream[k]);


              // Update the CN to VN message array
              CheckPassGB<<< ceil(M/32.0), 32, 0, stream[k] >>> (Dev_CtoV[k], Dev_VtoC[k], M, NbBranch, Dev_RowDegree); 
              //  Update the VN's (VN's are stored in the Decide array)
              APP_GB  <<< ceil(N/32.0), 32, 0, stream[k] >>> (Dev_Decide[k], Dev_CtoV[k], Dev_Receivedword[k], Dev_Interleaver, Dev_ColumnDegree, N, M, NbBranch); 
              // Check to see if updated codeword has been recovered
              ComputeSyndrome <<< ceil(M/32.0), 32, 0, stream[k] >>> (Dev_Decide[k], Dev_Mat, Dev_RowDegree, M, Dev_Syndrome[k]); 
              // Update host side memory for host controller decoder convergence check
              cudaMemcpyAsync(IsCodeword[k], Dev_Syndrome[k],  sizeof(int), cudaMemcpyDeviceToHost, stream[k]);
              // Update most recent decoded codeword copy to host memory (Most messages will recover with no iterations)
              cudaMemcpyAsync(Decide[k], Dev_Decide[k], N * sizeof(int), cudaMemcpyDeviceToHost, stream[k]);
            }
          }

            // Sync Host and Device
            cudaDeviceSynchronize();

          // Sync active streams for host side checks and updates
          for (int m=0; m<stream_count; m++) {
            if (stream_state[m] == 1)
              cudaStreamSynchronize(stream[m]);
              cudaDeviceSynchronize();
              //printf("In CHECK loop \n");
          }

          // Check for codeword recovery and update stream states as neccissary
          for (int m=0; m<stream_count; m++) {
            if (stream_state[m] == 1) {
              if (*IsCodeword[m]) {
              //if (IsCodeword[m]) {
                // Update stream state array and kill stream
                stream_state[m] = 0;
                cudaStreamDestroy(stream[m]);
                // Debug code
                // Print out number of iterations decoder took
                printf("Codework check value =  %d  for Stream %d recovered in %d runs \n",*IsCodeword[m],m,iter);
                
                /*
                for (int j=0; j<N; j++) {
                  printf("%d ", Decide[m][j]);
                }
                printf("\n");
                printf("\n");
                */
                 /*
                printf("\n");
                printf("\n");
                for (int j=0; j<NbBranch; j++) {
                  printf("%d ", VtoC[m][j]);
                }         
                */       
              }
            }
          }

          // Update variables for next round of streams
          iter++;
          int state_check = 0;
          batch_state = 0;
          for (int m=0; m<stream_count; m++) {
            state_check ^= stream_state[m];
            if (state_check == 1) {
              batch_state = 1;
              break; 
            }
          }
              
      }
      
      // Kill any residual streams associated with codeword that was not recovered
      for (int m=0; m<stream_count; m++) {
        if (stream_state[m] == 1)
          cudaStreamDestroy(stream[m]);
      }
      
      // Sync Host and Device
      cudaDeviceSynchronize();



      //============================================================================
  	  // Verification:  Uncomment for short runs
	    //============================================================================
      /*
      // Output uncorrupted codewords to file for verification purposes
      FILE *fptr1;
      fptr1 = (fopen("codewords_test_verification_pre_corrupt_02.txt", "a+"));
      // send codeword to output file
      for (k=0;k<N;k++) 
        fprintf(fptr1, "%u %s", Codeword[k], "");
      fprintf(fptr1, "%s", "\n");
      fclose(fptr1);
            
      // Output decoded codewords to file for verification purposes
      FILE *fptr2;
      fptr2 = (fopen("codewords_test_verification_decoded_02.txt", "a+"));
      // send codeword to output file
      for (k=0;k<N;k++) 
        fprintf(fptr2, "%u %s", Decide[k], "");
      fprintf(fptr2, "%s", "\n");
      fclose(fptr2);
      */



      // Run H*CW syndrome check
      //VerificationComputeSyndrome(Decide[bidx],Mat,RowDegree,M);
      

      
	    //============================================================================
  	  // Compute Statistics
	    //============================================================================
      /*
      nbtestedframes++;
	    NbError=0;for (k=0;k<N;k++)  if (Decide[k]!=Codeword[k]) NbError++;
	    NbBitError=NbBitError+NbError;
	    // Case Divergence
	    if (!(*IsCodeword)) {
		    NiterMoy=NiterMoy+NbIter;
		    NbTotalErrors++;
	    }

	    // Case Convergence to Right Codeword
	    if ((*IsCodeword)&&(NbError==0)) { NiterMax=max(NiterMax,iter+1); NiterMoy=NiterMoy+(iter+1); }
	      // Case Convergence to Wrong Codeword
	      if ((*IsCodeword)&&(NbError!=0)) {
		      NiterMax=max(NiterMax,iter+1); NiterMoy=NiterMoy+(iter+1);
		      NbTotalErrors++; NbUnDetectedErrors++;
		      Dmin=min(Dmin,NbError);
	      }
	
      // Stopping Criterion (IMPORTANT TO ADD THIS BACK IN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
	    if (NbTotalErrors==NBframes) break;
      */
      }
    

    // Print final statistics
    /*
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
    */




  
  }

  //  Clean Memory
  for (int m=0; m<stream_count; m++) {
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

  fclose(f);
  return(0);

}