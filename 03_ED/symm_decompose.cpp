#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

 /* Decompose the full Hilbert space using symmetry transformations */
 /* Decompose according to Winding Numbers */
 void winding_no_decompose(){
  extern void init_WindNo(std::vector<WindNo>&, int*);
  extern void initcheck_WindNo(std::vector<WindNo>&,int,int*);
  extern void display_basisNo(std::vector<WindNo>&,int);

  unsigned int i,j,count,p,p1,p2,p3,p4;
  int ix,iy,wx,wy;
  bool pxy,pyz,pzw,pwx;
  int count_flip;
  
  init_WindNo(Wind, lookup);
  // check the Winding number assignment
  //initcheck_WindNo(Wind,nWind,lookup);
  // scan basis states and pick out the Winding number
  for(i=0;i<basis.size();i++){

    // X-Winding (iy=0)
    // starting from ix=1, since ix=0 is fixed
    wx=0;
    for(ix=1;ix<LX;ix++){
      p=2*ix+1; // linear co-ordinate of the basis state; y-link
      if(basis[i][p]) wx++;
      else wx--;
    }
    /* copy to the relevant basis vector */
    count = lookup[LX-1+wx];
    Wind[count].nBasis++;
    Wind[count].basisVec.push_back(basis[i]);
  }
  display_basisNo(Wind,nWind);

  // sorting the basis states in each winding number sector
  //std::cout<<"Sorting the basis in each wind-no sector."<<std::endl;
  for(i=0; i<nWind; i++){
     Wind[i].sortbasis();
  }

  // compute the no of flippable plaquettes in the basis (ice) states
  // Note that this can only be calculated after sorting the basis!
  for(i=0; i<nWind; i++){
     Wind[i].flip_plaq();
  }

 }

/* Compute total number of sectors */
/* 2*Lx+1 sectors from -Lx, ..., 0, ..., Lx
 */
int calc_WindNo(int LX){
  int nWindSector=1 + 2*(LX-1); 
  printf("actual (Lx,Ly) extent of the lattice = (%d,%d)\n",LX,LY);
  printf("Effective winding number sectors = %d\n",nWindSector);
  return nWindSector;
}

/* initialize Winding numbers */
void init_WindNo(std::vector<WindNo> &Wind,int *lookup){
  int count,ix,iy;
  WindNo iWind;
  /* effective Lx for the lattice with FBC */
  int mLx;

  mLx   = LX-1; 
  count = 0;
  for(ix=-mLx; ix<=mLx; ix++){
      lookup[mLx+ix] = count; iWind.Wx = ix; iWind.Wy = 0;
      iWind.nBasis = 0; Wind.push_back(iWind); count++;
  }

}

void initcheck_WindNo(std::vector<WindNo> &Wind,int size,int *lookup){
  int i;
  std::cout<<"Checking the Winding Number initialization "<<std::endl;
  for(i=0;i<size;i++){
    printf(" count=%d;  Wx=% d; lookup[wx]=%d\n",i, Wind[i].Wx, lookup[LX-1+Wind[i].Wx]);
  }
}

void display_basisNo(std::vector<WindNo> &Wind,int size){
  int i,sum;
  sum=0;
  for(i=0;i<size;i++){
    //std::cout<<"count = "<<i<<std::endl;
    Wind[i].display();
    sum = sum + Wind[i].nBasis;
  }
  printf("Total states = %d\n",sum);
  NTOT=sum;
}

// sorts the basis states in a given Winding number sector
void WindNo::sortbasis(){
  std::sort(basisVec.begin(),basisVec.end());
}


// calculates the no of flippable plaquettes in each basis state
void WindNo::flip_plaq(){
    int i,j,p,p1,p2,p3,p4;
    int count_flip;
    bool pxy,pyz,pzw,pwx;

    //printf("In sector (% d); basis states = %ld \n",Wx,nBasis);
    // error message for zero basis
    //if(!nBasis){ std::cout<<" Skipping the sector with no basis! "<<std::endl; return; }
    for(i=0; i<nBasis; i++){
      // compute the diagonal term in the Hamiltonian
      /* a single plaquette is arranged as 
                pzw
             o-------o
             |       |
        pwx  |   p   |  pyz
             |       |
             o-------o
                pxy
      */
      count_flip=0;
      for(p=0;p<VOL;p++){
          p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
          pxy=basisVec[i][p1]; pyz=basisVec[i][p2]; pzw=basisVec[i][p3]; pwx=basisVec[i][p4];
          if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)) count_flip++;
      }
      nflip.push_back(count_flip);
      //std::cout<<"Basis state "<<i<<" = "<<count_flip<<std::endl;
    }
}

