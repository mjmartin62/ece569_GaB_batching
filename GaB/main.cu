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

#define arrondi(x) ((ceil(x)-x)<(x-floor(x))?(int)ceil(x):(int)floor(x))
#define min(x,y) ((x)<(y)?(x):(y))
#define signf(x) ((x)>=0?0:1)
#define	max(x,y) ((x)<(y)?(y):(x))
#define SQR(A) ((A)*(A))
#define BPSK(x) (1-2*(x))
#define PI 3.1415926536


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
int main(int argc, char * argv[])
{
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



  strcpy(FileMatrix,argv[1]); 	// Matrix file
  strcpy(FileResult,argv[2]); 	// Results file
  //--------------Simulation input for GaB BF-------------------------
  NbMonteCarlo=100000;	    // Maximum nb of codewords sent
  NbIter=200; 	            // Maximum nb of iterations
  alpha= 0.01;              // Channel probability of error
  NBframes=100;	            // Simulation stops when NBframes in error
  Graine=1;		            // Seed Initialization for Multiple Simulations

    // brkunl
  alpha_max= 0.060;		    //Channel Crossover Probability Max and Min
  alpha_min= 0.020;
  alpha_step=0.010;


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

  printf("Matrix Loaded on Host side\n");

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

  printf("Graph Build \n");

  // ----------------------------------------------------
  // Decoder variables and memory allocation
  // ----------------------------------------------------
  int *CtoV,*VtoC,*Codeword,*Receivedword,*Decide,*U,l,kk;
  int iter,numB;
  CtoV=(int *)calloc(NbBranch,sizeof(int));
  VtoC=(int *)calloc(NbBranch,sizeof(int));
  Codeword=(int *)calloc(N,sizeof(int));
  Receivedword=(int *)calloc(N,sizeof(int));
  Decide=(int *)calloc(N,sizeof(int));
  U=(int *)calloc(N,sizeof(int));
  srand48(time(0)+Graine*31+113);


  //Allocate VtoC, CtoV, Receivedword, Interleaver location in the device
  int *Dev_VtoC, *Dev_CtoV, *Dev_Receivedword, *Dev_Interleaver, *Dev_ColumnDegree;
  int *Dev_Decide, *Dev_RowDegree;
  int *Dev_Mat, *Dev_Syndrome;
  if (cudaMalloc((void **) &Dev_VtoC, NbBranch * sizeof(int)) != cudaSuccess) {
  printf("malloc error for *Dev_VtoC \n");
  return 0;
  }
  if (cudaMalloc((void **) &Dev_CtoV, NbBranch * sizeof(int)) != cudaSuccess) {
  printf("malloc error for *Dev_CtoV \n");
  return 0;
  }
  if (cudaMalloc((void **) &Dev_Receivedword, N * sizeof(int)) != cudaSuccess) {
  printf("malloc error for *Dev_Receivedword \n");
  return 0;
  }
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
  if (cudaMalloc((void **) &Dev_Decide, N * sizeof(int)) != cudaSuccess) {
  printf("malloc error for *Dev_Decide \n");
  return 0;
  }

  if(cudaMalloc((void **) &Dev_Mat, M * RowDegMax * sizeof(int)) != cudaSuccess) {
  printf("malloc error for *Dev_Mat \n");
  return 0;
  }
  
  if(cudaMalloc((void **) &Dev_Syndrome, sizeof(int)) != cudaSuccess) {
  printf("malloc error for *Dev_Syndrome \n");
  return 0;
  }

  //Copy interleaver to Device 
   if (cudaMemcpy(Dev_Interleaver, Interleaver, NbBranch * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
   printf("data transfer error from host to device on Dev_Interleaver\n");
   return 0;
   }

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

    //Copy matrixB to device memory
     if (cudaMemcpy(Dev_VtoC, VtoC, NbBranch * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
     printf("data transfer error from host to device on deviceB\n");
     return 0;
      }
     if (cudaMemcpy(Dev_CtoV, CtoV, NbBranch * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
      printf("data transfer error from host to device on deviceB\n");
      return 0;
     }
      


   
 //  printf("\n Matrix copied ");
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
  int *IsCodeword = (int *)calloc(1, sizeof(int *)) ,nb;
  int NiterMoy,NiterMax;
  int Dmin;
  int NbTotalErrors,NbBitError;
  int NbUnDetectedErrors,NbError;
  int *energy;
  energy=(int *)calloc(N,sizeof(int));
 
  strcpy(FileName,FileResult);
  f=fopen(FileName,"w");
  fprintf(f,"-------------------------Gallager B--------------------------------------------------\n");
  fprintf(f,"alpha\t\tNbEr(BER)\t\tNbFer(FER)\t\tNbtested\t\tIterAver(Itermax)\t\tNbUndec(Dmin)\n");

  printf("-------------------------Gallager B  Parallel code--------------------------------------------------\n");
  printf("alpha\t\t\tNbEr(BER)\t\tNbFer(FER)\t\tNbtested\t\tIterAver(Itermax)\t\tNbUndec(Dmin)\n");

  for(alpha=alpha_max;alpha>=alpha_min;alpha-=alpha_step) {

  NiterMoy=0;NiterMax=0;
  Dmin=1e5;
  NbTotalErrors=0;NbBitError=0;
  NbUnDetectedErrors=0;NbError=0;

  //--------------------------------------------------------------
  // Main loop for max number of codeword simulations
  for (nb=0,nbtestedframes=0;nb<NbMonteCarlo;nb++)
  {
    //encoding
    for (k=0;k<rank;k++) U[k]=0;
	for (k=rank;k<N;k++) U[k]=floor(drand48()*2);
	for (k=rank-1;k>=0;k--) { for (l=k+1;l<N;l++) U[k]=U[k]^(MatG[k][l]*U[l]); }
	for (k=0;k<N;k++) Codeword[PermG[k]]=U[k];
	// All zero codeword
	//for (n=0;n<N;n++) { Codeword[n]=0; }

    // Add Noise and send possibly corrupted Codeword to Receivedword
    for (n=0;n<N;n++)  
      if (drand48()<alpha) 
        Receivedword[n]=1-Codeword[n]; 
      else 
        Receivedword[n]=Codeword[n];

	//============================================================================
 	// Decoder
	//============================================================================
	// Initialize the CN to VN message array to 0
  for (k=0; k<NbBranch; k++) {
    CtoV[k]=0;
  }

 


       if (cudaMemcpy(Dev_Receivedword, Receivedword, N * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess){
      printf("data transfer error from host to device on deviceB\n");
      return 0;
     }
  // Outer loop to limit (max of 100) the number of VN node updates thru parity checks
	for (iter=0;iter<NbIter;iter++)
	  {

          // Update VN to CN message array
        DataPassGB <<< ceil(N/32.0), 32 >>> (Dev_VtoC, Dev_CtoV, Dev_Receivedword, Dev_Interleaver, Dev_ColumnDegree,N,NbBranch, iter);  cudaDeviceSynchronize();

        // Update the CN to VN message array
        CheckPassGB<<< ceil(M/32.0), 32 >>> (Dev_CtoV, Dev_VtoC, M,NbBranch,Dev_RowDegree); cudaDeviceSynchronize();

        //  Update the VN's (VN's are stored in the Decide array)
        APP_GB  <<<ceil(N/32.0), 32>>> (Dev_Decide, Dev_CtoV, Dev_Receivedword, Dev_Interleaver, Dev_ColumnDegree, N,M, NbBranch); cudaDeviceSynchronize();
  
        // Check to see if updated codeword has been recovered
        ComputeSyndrome <<<ceil(M/32.0), 32>>>(Dev_Decide, Dev_Mat, Dev_RowDegree, M, Dev_Syndrome); cudaDeviceSynchronize();

         if (cudaMemcpy(IsCodeword, Dev_Syndrome,  sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess){
         printf("data transfer error from device to Dev Syndrome\n");
         return 0;
          }
          cudaDeviceSynchronize();
  
         //printf(" \n ISCodeword is : %d \n", *IsCodeword);
     
        if (*IsCodeword) 
          break;
	  }

        if (cudaMemcpy(Decide, Dev_Decide, N * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess){
         printf("data transfer error from device to host on Dev Decide\n");
         return 0;
         }
          cudaDeviceSynchronize();

	//============================================================================
  	// Compute Statistics
	//============================================================================
      nbtestedframes++;
	  NbError=0;for (k=0;k<N;k++)  if (Decide[k]!=Codeword[k]) NbError++;
	  NbBitError=NbBitError+NbError;
	// Case Divergence
	  if (!(*IsCodeword))
	  {
		  NiterMoy=NiterMoy+NbIter;
		  NbTotalErrors++;
	  }
	// Case Convergence to Right Codeword
	  if ((*IsCodeword)&&(NbError==0)) { NiterMax=max(NiterMax,iter+1); NiterMoy=NiterMoy+(iter+1); }
	// Case Convergence to Wrong Codeword
	  if ((*IsCodeword)&&(NbError!=0))
	  {
		  NiterMax=max(NiterMax,iter+1); NiterMoy=NiterMoy+(iter+1);
		  NbTotalErrors++; NbUnDetectedErrors++;
		  Dmin=min(Dmin,NbError);
	  }
	// Stopping Criterion
	 if (NbTotalErrors==NBframes) break;
  }
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
  fclose(f);
  return(0);
}


