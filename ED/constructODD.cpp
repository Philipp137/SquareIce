/* construct gauge inequivalent states of the U(1) model  on square lattice */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include <algorithm>
#include <iterator>
#include "define.h"
//#include <algorithm> //std::count
// Modified version of construct.cpp for gauge invariant state counting
// when the lattice is odd in extent. 

void conststatesODD(){
  extern void add(char*,int);
  extern void addATpos(char*,int);
  extern int checkGL(bool*,int);
  extern void storeconf(bool*,int*);
  extern void Pbasis(int);
  int count,fluxcount,count1;
  int flagFLUX;
  int base=6;
  int cSIZE = VOL2+1;
  long long int n,NST;
  char str[cSIZE];
  int linkset[VOL];

  bool conf[2*VOL];
  int flagGI,flagJUMP;
  int p,pp,l,cc;
  int q,xf,yf,xb,yb;
  int q0x,q0y,qx,qy;
  int Q,e1,e2,e3,e4;
  char cG,cG1,cG2;

  count = 0; fluxcount=0; 
  NST = pow(6,cSIZE);
  for(p=0;p<cSIZE;p++) str[p]='0';
  str[cSIZE]='\0';
 
  n=0;
  while(n<NST){
     
    // Eliminate states based on Gauss Law realizations on pair of neighboring sites; see notes below
    // i) (0,(LX-1)/2) 
    
    // i think this is wrong: p = (LX-1)/2; cG1 = str[p]; cG2 = str[0];
    // need to test the two methods; and determine which one is faster; do not delete
    //if( (cG1=='0') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){      add(str,cSIZE); n++; continue; }
    //else if( (cG1=='1') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='4') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='2') && ((cG2=='0')||(cG2!='1')||(cG2!='3')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='3') && ((cG2!='0')||(cG2!='1')||(cG2!='3')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='5') && ((cG2!='0')||(cG2!='1')||(cG2!='3')) ){ add(str,cSIZE); n++; continue; } 
    if( (cG1=='0') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){      
	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='1') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ 
	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='4') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ 
	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='2') && ((cG2=='0')||(cG2!='1')||(cG2!='3')) ){ 
	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='3') && ((cG2!='0')||(cG2!='1')||(cG2!='3')) ){ 
	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='5') && ((cG2!='0')||(cG2!='1')||(cG2!='3')) ){ 
	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; } 
    // and ii) (cSIZE-1, cSIZE)
    cG1 = str[cSIZE-1]; cG2 = str[cSIZE-2]; 
    // in this case, using the formula n=n+pow(6,cSIZE-1-p) is same as n++
    if( (cG1=='0') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){      add(str,cSIZE); n++; continue; }
    else if( (cG1=='1') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='4') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='2') && ((cG2=='0')||(cG2!='1')||(cG2!='3')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='3') && ((cG2!='0')||(cG2!='1')||(cG2!='3')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='5') && ((cG2!='0')||(cG2!='1')||(cG2!='3')) ){ add(str,cSIZE); n++; continue; } 
    // remove states that do not satisfy the OBC


    // print the state to check
    printf("state=%s\n",str);

    flagGI=1; flagJUMP=0;
    //for(l=0;l<VOL;l++) linkset[l]=0;
    

    if(flagJUMP==0){ add(str,cSIZE); n++; }
  }
  //printf("No of gauge invariant states: %d %ld\n",count,basis.size());
  //printf("Going to delete any duplicates\n");
  // delete duplicates
  //std::sort( basis.begin(),basis.end() );
  //basis.erase( unique(basis.begin(),basis.end() ), basis.end() );
  //count1 = basis.size();
  //printf("No of gauge invariant states after removal of duplicates: %d\n",count1);
  
  //Pbasis(count1);
  //NTOT = count1;
}


/*
 * ========================================================
 *  Notes about the state counting on lattices with Lx=odd
 * ========================================================
 */
