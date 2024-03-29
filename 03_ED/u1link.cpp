#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

/* according to the rules of cpp, the variables are declared here 
 * and also in the header file as extern such that they are avl to
 * the other functions.
 */
int *next[2*DIM+1];
int *nextCHK[2*DIM+1];
int *chk2lin,*lin2chk;
int LX,LY,VOL,VOL2;
// Labels the winding number sectors
// lookup[LX+1+wx] refers to the wx winding number sector
int nWindSector;
int *lookup;
std::vector<WindNo> Wind;
unsigned int nWind;
double lam,Ti,Tf,dT;
int NTOT,NH;
std::vector<std::vector<bool>> basis;
std::vector<std::vector<bool>> basis_nonflip;
std::vector<std::vector<bool>> basis_flip;
int STORE_SVD;
int CHKDIAG;
int FBC;

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,q;
  int x,y;
  int wx,wy;
  int sector;
  extern void initneighbor(void);
  extern void conststatesEVEN(void);
  extern void conststatesODD(void);
  extern void printbasis(void);
  extern int** allocateint2d(int, int);
  extern void deallocateint2d(int**,int,int);
  extern int calc_WindNo(int);

  fptr = fopen("QUEUE","r");
  if(fptr == NULL){
      printf("could not open QUEUE FILE to open\n");
      exit(1);
  }
  fscanf(fptr,"%s %d\n",string,&LX);
  fscanf(fptr,"%s %d\n",string,&LY);
  fscanf(fptr,"%s %lf\n",string,&lam);
  fscanf(fptr,"%s %lf\n",string,&Ti);
  fscanf(fptr,"%s %lf\n",string,&Tf);
  fscanf(fptr,"%s %lf\n",string,&dT);
  fscanf(fptr,"%s %d\n",string,&LEN_A);
  fclose(fptr);
  if( LY !=2 ) { printf("Code does not work with LY>2. OBC is not implemented \n"); exit(0); }
  if(LX<LY) printf("Please make sure LX >= LY. Unforseen errors can occur otherwise. \n");
  VOL = LX*LY;
  VOL2 = VOL/2;

  // decide whether to check the results of the diagonalization 
  CHKDIAG=1;

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i] = (int *)malloc(VOL*sizeof(int)); 
    nextCHK[i] = (int *)malloc(VOL*sizeof(int));
  }

  /* Initialize checkerboard co-ordinates */
  lin2chk = (int *)malloc(VOL*sizeof(int));
  chk2lin = (int *)malloc(VOL*sizeof(int));
  initneighbor();
  
  /* Winding number sectors */
  lookup = (int *)malloc((2*(LX-1)+1)*sizeof(int));

  /* set up the flag for fixed boundary condition */
  FBC = 1;
  /* build basis states satisfying Gauss' Law */
  if((LX%2)==0) conststatesEVEN();
  else if((LX%2)==1) conststatesODD();

  totalH_diag();
  
  /* get number of winding number sectors */
  //nWind = calc_WindNo(LX);  
  //Wind.reserve(nWind); 
  //winding_no_decompose();
  // get the winding number sector (wx,wy)
  //if((LX%2)==0)      wx = -1;
  //else if((LX%2)==1) wx = 0; 
  //sector = lookup[LX-1+wx];
  //constH(sector);
  //std::cout<<"Back from the ED routine"<<std::endl;

  // calculate the expectation value of Oflip for every eigenstate 
  //calc_Oflip(sector);
  // real-time evolution of cartoon states and Locshmidt Echo 
  //evolve_cartoons(sector);
  // real-time evolution of Entanglement Entropy
  // for this, one needs to store the SVD coefficients
  //evolve_Eent(sector);
  // real-time correlation function
  //evolve_corrf1(sector);
  // calculate the Entanglement Entropy for the states
  // to only calculate the 
  //entanglementEntropy(sector);

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); free(nextCHK[i]); }
  free(chk2lin); free(lin2chk); 
  //free(lookup);
  //Wind.clear();

  return 0;
}
