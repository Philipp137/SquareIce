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

    //printf("state being examined =%s; \n",str);

    // Eliminate states based on Gauss Law realizations on pair of neighboring sites; see notes below
    // i) locations (LX-1)/2 and 0 in str  (+x neighbor for y=1, x=1)
    p = (LX-1)/2; cG1 = str[p]; cG2 = str[0];
    // need to test the two methods; and determine which one is faster; do not delete
    //if(      (cG1=='0') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='1') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='2') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='3') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='4') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    //else if( (cG1=='5') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ add(str,cSIZE); n++; continue; } 
    if(      (cG1=='0') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){      
    	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='1') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ 
    	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='2') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ 
    	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='3') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ 
    	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='4') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ 
    	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    else if( (cG1=='5') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ 
    	    addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; } 
    // ii) locations cSIZE-1 and cSIZE-2 in str (+x neighbor for y=2, x=Lx-1)
    cG1 = str[cSIZE-1]; cG2 = str[cSIZE-2]; 
    // in this case, using the formula n=n+pow(6,cSIZE-1-p) is same as n++
    if(      (cG1=='0') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='1') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='2') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='3') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='4') && ((cG2=='2')||(cG2=='4')||(cG2=='5')) ){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='5') && ((cG2=='0')||(cG2=='1')||(cG2=='3')) ){ add(str,cSIZE); n++; continue; }
    // iii) locations (LX-1)/2 and cSIZE-1 in str (+y neighbor for x=Lx) 
    p = (LX-1)/2; q = cSIZE-1; cG1 = str[p]; cG2 = str[q];
    if(      (cG1=='2') &&  (cG2!='0')) {             add(str,cSIZE); n++; continue; }
    else if( (cG1=='3') && ((cG2!='1')||(cG2!='4'))){ add(str,cSIZE); n++; continue; }
    else if( (cG1=='5') && ((cG2!='1')||(cG2!='4'))){ add(str,cSIZE); n++; continue; }
    // the rest of the cases are removed by the conditions below
    // remove states that do not satisfy the OBC
    if(FBC){
      p = 0; cG = str[p]; // str[0] can only have realization 5 of the Gauss' Law
      if(cG!='5'){ addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
      p = (LX-1)/2; cG = str[p]; // str[p] has the +x link to be -1/2
      if((cG=='0')||(cG=='1')||(cG=='4')){ addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
      p = (LX+1)/2; cG = str[p]; // str[p] has the -x link to be +1/2
      if((cG=='2')||(cG=='4')||(cG=='5')){ addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
      p = cSIZE-1; cG = str[p];  // str[p] has the +x link to be +1/2
      if((cG=='2')||(cG=='3')||(cG=='5')){ addATpos(str,p); n=n+pow(6,cSIZE-1-p); continue; }
    }

    for(l=0;l<VOL;l++) linkset[l]=0;
    flagGI=1; flagJUMP=0;
    // assign the links touching at the vertices
    l=0; // start with the site 0; links fixed by FBC 
    xf=next[DIM+1][l]; xb=next[DIM-1][l]; 
    yf=next[DIM+2][l]; yb=next[DIM-2][l]; 
    q0x=2*l; q0y=2*l+1; qx=2*xb; qy=2*yf+1; 
    conf[q0x]=false; conf[q0y]=true; conf[qy]=true;
    linkset[l]=linkset[l]+3; linkset[xf]++; linkset[yf]++; linkset[yb]++;
    // the other sites
    for(l=1; l<VOL2; l++){
      cG=str[l]; q=chk2lin[l];
      //std::cout<<"chkbrd site ="<<l<<"; linear site ="<<q<<std::endl;
      xf=next[DIM+1][q]; xb=next[DIM-1][q];
      yf=next[DIM+2][q]; yb=next[DIM-2][q];
      // convert to 1d co-ordinates on the lattice configuration
      q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
      // from Fig 1 (left to right) in https://arxiv.org/pdf/1311.2459.pdf
      switch(cG){
        case '0' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=true;  conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '1' : conf[q0x]=true;  conf[q0y]=true;  conf[qx]=true;  conf[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '2' : conf[q0x]=false; conf[q0y]=false; conf[qx]=false; conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '3' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=true;  conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '4' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=false; conf[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '5' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=false; conf[qy]=true;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
      } //close switch
      // check Gauss' Law at this point on the odd sublattices
      //if(linkset[yb]==4){ flagGI = checkGL(conf,yb);
      //   if(flagGI==0){ addATpos(str,l); n=n+pow(6,cSIZE-1-l); flagJUMP=1; break; }
      //}
      //if(linkset[xb]==4){ flagGI = checkGL(conf,xb);    
      //   if(flagGI==0){ addATpos(str,l);  n=n+pow(6,cSIZE-1-l); flagJUMP=1; break; }
      //}
      //if(linkset[xf]==4){ flagGI = checkGL(conf,xf);
      //   if(flagGI==0){  addATpos(str,l); n=n+pow(6,cSIZE-1-l); flagJUMP=1; break; }
      //}
      //if(linkset[yf]==4){ flagGI = checkGL(conf,yf);
      //   if(flagGI==0){ addATpos(str,l);  n=n+pow(6,cSIZE-1-l); flagJUMP=1; break; } 
      //}
    }// close l
    // the last link remaining to be specified at linear site (x=Lx,y=2)
    p=VOL-1; q0x=2*p; xf=next[DIM+1][p]; conf[q0x]=true; 
    linkset[p]++; linkset[xf]++;
    if(flagGI){
      for(p=VOL2;p<VOL;p++){
         pp=chk2lin[p]; 
         flagGI = checkGL(conf,pp);
         if(flagGI==0) break;
      }
    }// if flagGI closure
    if(flagGI){
	    storeconf(conf,&count);
            // check to see if all links are set
            for(l=0;l<VOL;l++){
              std::cout<<"Site = "<<l<<"; #-links set  = "<<linkset[l]<<std::endl;
            }
            printf("state=%s; str[0]=%c \n",str,str[0]);
    }
    if(flagJUMP==0){ add(str,cSIZE); n++; }
  }
  printf("No of gauge invariant states: %d %ld\n",count,basis.size());
  printf("Going to delete any duplicates\n");
  // delete duplicates
  std::sort( basis.begin(),basis.end() );
  basis.erase( unique(basis.begin(),basis.end() ), basis.end() );
  count1 = basis.size();
  //printf("No of gauge invariant states after removal of duplicates: %d\n",count1);
  
  Pbasis(count1);
  NTOT = count1;
}


/*
 * ========================================================
 *  Notes about the state counting on lattices with Lx=odd
 * ========================================================
 * 1. Follow the same base 6 system to assign Gauss' Law on sites
 *    on the checkerboard even sites. However, since Lx=odd, this
 *    needs more specification to completely specify the links of
 *    a basis state.
 * 2. Consider the example of Lx=3 and Lx=5 with Ly=2 always. 
 *
 *      x===o===x===x      x===o===x===o===x===x
 *      |  2|  3|   |      |  3|   |  4|  5|   |  
 *      o===x===x===o      o===x===o===x===x===o
 *     0|   |  1|  0|     0|   |  1|   |  2|  0|
 *      x===o===x===x      x===o===x===o===x===x
 * site 0   1   2   0      0   1   2   3   4   0
 *
 *    The sites marked with x denote where the Gauss' Law is
 *    specified. The numbers on the top right corner of the 
 *    x-sites denote the checkerboard co-ordinates. 
 * 3. A state is completely specified by specifying the Gauss' Law
 *    realizations on each of these sites. We follow the convention:
 *       string-array index:  0 1 2 3 
 *       state-string      :  2 3 0 4 
 *    For the 5x2 lattice we have
 *       string-array index:  0 1 2 3 4 5
 *       state-string      :  2 3 4 0 5 0   
 * 4. Using this information, we can already rule out certain states.
 *    For example, the pair of location in the string ((Lx-1)/2,0) 
 *    and (Lx-1,Lx) only has following realization of Gauss' Law 
 *    shown in the table below. The second column lists the GL cases
 *    while the third one has the possible GL cases in the +x neighbor.  
 *      ____________________________________________________
 *      site  |  either (Lx-1)/2 or Lx-1   |  either Lx or 0
 *      ----------------------------------------------------
 *      GL    |   0                        |   0,1,3
 *            |   1                        |   0,1,3
 *            |   2                        |   2,4,5
 *            |   3                        |   2,4,5
 *            |   4                        |   0,1,3
 *            |   5                        |   2,4,5
 *      ------------------------------------------------------      
 * 5. There is the only vertical neighboring sites whose GL cases
 *    also matter, they at the string positions ((LX-1)/2,cSIZE-1);
 *    corresponding to y=1,2 and x=Lx. 
 *      ________________________________
 *      site  |  (Lx-1)/2   |  cSIZE-1
 *      --------------------------------
 *      GL    |  *0         |   0
 *            |  *1         |   1,4
 *            |   2         |   0
 *            |   3         |   1,4
 *            |  *4         |   0
 *            |   5         |   1,4
 *      --------------------------------   
 *      The cSIZE-1 is the +y neighbor for the site (Lx-1)/2. The 
 *      latter is decided also by the OBC case for the outgoing link
 *      for the site cSIZE-1. The GL cases marked by stars show the
 *      ones inconsistent with the OBC on the outgoing link for the
 *      site (Lx-1)/2.  
 *  6. As usual, the generated configurations are then checked by 
 *     checking the GL on the other vertices.
 */     
