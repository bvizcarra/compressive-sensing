/*
Begin tsimplexv2.cpp
*/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <sys/sysctl.h>
#include <sys/times.h>
#include <sys/vtimes.h>
#include <sys/resource.h>
#include "tsimplexv2-N3-M3-K1.h"
//#include "tsimplexv2-N10-M6-K1.h"
//#include "tsimplexv2-N20-M11-K2.h"
//#include "tsimplexv2-N30-M11-K2.h"
//#include "tsimplexv2-N40-M11-K3.h"
//#include "tsimplexv2-N80-M16-K3.h"
//#include "tsimplexv2-N256-M80-K8.h"
using namespace std;

int  IPOSV[MMAX], IZROV[NMAX];
int  i,j,ICASE;
REAL R;

void simp1(MAT,int,int *,int,int,int *,REAL *);
void simp2(MAT,int,int,int *,int,int *,int,REAL *);
void simp3(MAT,int,int,int,int);


void simplx(MAT a,int m,int n,int m1,int m2,int m3,int *icase,int *izrov, int *iposv)  { 
	/*---------------------------------------------------------------------------------------- 
	 USES simp1,simp2,simp3 
	 Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3, 
	 and output parameters a, icase, izrov, and iposv are described above (see reference). 
	 Parameters: MMAX is the maximum number of constraints expected; NMAX is the maximum number 
	 of variables expected; EPS is the absolute precision, which should be adjusted to the 
	 scale of your variables.                                    
	 -----------------------------------------------------------------------------------------*/
	int i,ip,ir,is,k,kh,kp,m12,nl1,nl2,l1[NMAX],l2[MMAX],l3[MMAX]; 
	REAL bmax,q1,EPS=1e-6; 
	if(m != m1+m2+m3) {
		printf(" Bad input constraint counts in simplx.\n");
		return;
	}	
	nl1=n; 
	for (k=1; k<=n; k++) { 
		l1[k]=k;     //Initialize index list of columns admissible for exchange.
		izrov[k]=k;  //Initially make all variables right-hand. 
	} 
	nl2=m;
	for (i=1; i<=m; i++) { 
		if (a[i+1][1] < 0.0) {
			printf(" Bad input tableau in simplx, Constants bi must be nonnegative.\n");
			return;
		}
		l2[i]=i;
		iposv[i]=n+i; 
		/*------------------------------------------------------------------------------------------------
		 Initial left-hand variables. m1 type constraints are represented by having their slackv ariable 
		 initially left-hand, with no artificial variable. m2 type constraints have their slack 
		 variable initially left-hand, with a minus sign, and their artificial variable handled implicitly 
		 during their first exchange. m3 type constraints have their artificial variable initially 
		 left-hand.     
		 ------------------------------------------------------------------------------------------------*/ 
	}
	for (i=1; i<=m2; i++) l3[i]=1;
	ir=0;
	if(m2+m3 == 0) goto e30; //The origin is a feasible starting solution. Go to phase two. 
	ir=1;
	for (k=1; k<=n+1; k++) { //Compute the auxiliary objective function. 
		q1=0.0; 
		for (i=m1+1; i<=m; i++) q1 += a[i+1][k]; 
		a[m+2][k]=-q1; 
	} 
e10: simp1(a,m+1,l1,nl1,0,&kp,&bmax);    //Find max. coeff. of auxiliary objective fn 
	if(bmax <= EPS && a[m+2][1] < -EPS) { 
		*icase=-1;    //Auxiliary objective function is still negative and can't be improved, 
		return;       //hence no feasible solution exists.
	}	
	else if (bmax <= EPS && a[m+2][1] <= EPS) { 
		//Auxiliary objective function is zero and can't be improved; we have a feasible starting vector. 
		//Clean out the artificial variables corresponding to any remaining equality constraints by 
		//goto 1's and then move on to phase two by goto 30.
		m12=m1+m2+1;
		if (m12 <= m)
			for (ip=m12; ip<=m; ip++) 
				if(iposv[ip] == ip+n) {   //Found an artificial variable for an equalityconstraint. 
					simp1(a,ip,l1,nl1,1,&kp,&bmax); 
					if(bmax > EPS) goto e1; //Exchange with column corresponding to maximum 
				}                         //pivot element in row. 
		ir=0;
		m12=m12-1;
		if (m1+1 > m12) goto e30; 
		for (i=m1+1; i<=m1+m2; i++)     //Change sign of row for any m2 constraints 
			if(l3[i-m1] == 1)             //still present from the initial basis. 
				for (k=1; k<=n+1; k++) 
					a[i+1][k] *= -1.0; 
		goto e30;                       //Go to phase two. 
	} 
	
	simp2(a,m,n,l2,nl2,&ip,kp,&q1); //Locate a pivot element (phase one). 
	
	if(ip == 0) {                         //Maximum of auxiliary objective function is 
		*icase=-1;                          //unbounded, so no feasible solution exists.
		return; 
	} 
e1: simp3(a,m+1,n,ip,kp); 
	//Exchange a left- and a right-hand variable (phase one), then update lists. 
	if(iposv[ip] >= n+m1+m2+1) { //Exchanged out an artificial variable for an 
		//equality constraint. Make sure it stays 
		//out by removing it from the l1 list. 
		for (k=1; k<=nl1; k++) 
			if(l1[k] == kp) goto e2; 
	e2: nl1=nl1-1; 
		for (is=k; is<=nl1; is++)  l1[is]=l1[is+1]; 
	}	
	else {
		if(iposv[ip] < n+m1+1) goto e20;
		kh=iposv[ip]-m1-n; 
		if(l3[kh] == 0) goto e20;  //Exchanged out an m2 type constraint. 
		l3[kh]=0;                  //If it's the first time, correct the pivot column 
		//or the minus sign and the implicit 
		//artificial variable. 
	}   
	a[m+2][kp+1] += 1.0; 
	for (i=1; i<=m+2; i++)  a[i][kp+1] *= -1.0; 
e20: is=izrov[kp];             //Update lists of left- and right-hand variables. 
	izrov[kp]=iposv[ip]; 
	iposv[ip]=is; 
	if (ir != 0) goto e10;       //if still in phase one, go back to 10. 
	//End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it. 
e30: simp1(a,0,l1,nl1,0,&kp,&bmax); //Test the z-row for doneness. 
	if(bmax <= EPS) {          //Done. Solution found. Return with the good news. 
		*icase=0; 
		return; 
	}
	simp2(a,m,n,l2,nl2,&ip,kp,&q1);   //Locate a pivot element (phase two). 
	if(ip == 0) {             //Objective function is unbounded. Report and return. 
		*icase=1; 
		return; 
	} 
	simp3(a,m,n,ip,kp);       //Exchange a left- and a right-hand variable (phase two), 
	goto e20;                 //update lists of left- and right-hand variables and 
}                           //return for another iteration.

// The preceding routine makes use of the following utility subroutines: 

void simp1(MAT a,int mm,int *ll,int nll,int iabf,int *kp,REAL *bmax) {
	//Determines the maximum of those elements whose index is contained in the supplied list 
	//ll, either with or without taking the absolute value, as flagged by iabf. 
	int k; 
	REAL test; 
	*kp=ll[1]; 
	*bmax=a[mm+1][*kp+1];
	if (nll < 2) return; 
	for (k=2; k<=nll; k++) { 
		if(iabf == 0) 
			test=a[mm+1][ll[k]+1]-(*bmax); 
		else
			test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax); 
		if(test > 0.0) { 
			*bmax=a[mm+1][ll[k]+1]; 
			*kp=ll[k]; 
		} 
	} 
	return; 
} 

void simp2(MAT a, int m,int n,int *l2,int nl2,int *ip,int kp,REAL *q1) { 
	REAL EPS=1e-6; 
	//Locate a pivot element, taking degeneracy into account. 
	int i,ii,k; 
	REAL q,q0,qp; 
	*ip=0; 
	if(nl2 < 1) return;
	for (i=1; i<=nl2; i++) 
		if (a[i+1][kp+1] < -EPS) goto e2; 
	return;  //No possible pivots. Return with message. 
e2: *q1=-a[l2[i]+1][1]/a[l2[i]+1][kp+1]; 
	*ip=l2[i];
	if (i+1 > nl2) return; 
	for (i=i+1; i<=nl2; i++) { 
		ii=l2[i];
		if(a[ii+1][kp+1] < -EPS) { 
			q=-a[ii+1][1]/a[ii+1][kp+1]; 
			if (q <  *q1) { 
				*ip=ii; 
				*q1=q; 
			}
			else if (q == *q1) {  //We have a degeneracy.
				for (k=1; k<=n; k++) { 
					qp=-a[*ip+1][k+1]/a[*ip+1][kp+1]; 
					q0=-a[ii+1][k+1]/a[ii+1][kp+1]; 
					if (q0 != qp) goto e6; 
				} 
			e6:     if (q0 < qp) *ip=ii; 
			} 
		} 
	} 
	return; 
} 

void simp3(MAT a,int i1,int k1,int ip,int kp) { 
	//Matrix operations to exchange a left-hand and right-hand variable (see text). 
	int ii,kk; 
	REAL piv; 
	piv=1.0/a[ip+1][kp+1];
	if (i1 >= 0) 
		for (ii=1; ii<=i1+1; ii++) 
			if (ii-1 != ip) { 
				a[ii][kp+1] *= piv; 
				for (kk=1; kk<=k1+1; kk++) 
					if (kk-1 != kp) 
						a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1]; 
			} 
	for (kk=1; kk<=k1+1; kk++) 
		if(kk-1 !=  kp) a[ip+1][kk] =-a[ip+1][kk]*piv; 
	a[ip+1][kp+1]=piv; 
	return; 
} 

void getSystemValues() {
        FILE* file = fopen("/proc/self/status", "r");
        //int result = -1;
        char line[50];
        while (fgets(line, 50, file) != NULL) {
		//printf("       Status file line - %s\n", line);
		//if (strncmp(line, "VmSize:", 7) == 0) printf("\n       Virtual memory size - %s", line);
		//if (strncmp(line, "VmHWM:", 6) == 0) printf("       Peak resident set size - %s", line);
		//if (strncmp(line, "VmRSS:", 6) == 0) printf("       Resident set size - %s", line);
		if (strncmp(line, "VmData:", 7) == 0) printf("\n       Size of data segment - %s", line);
		if (strncmp(line, "VmStk:", 6) == 0) printf("       Size of stack segment - %s", line);
		if (strncmp(line, "VmExe:", 6) == 0) printf("       Size of text segment - %s", line);
        //if (strncmp(line, "VmSize:", 7) == 0) result = parseLine(line);
                //break;
        }
        fclose(file);
}

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

int main() {
	
	printf("\n       Program simplex version 2 %d-%d-%d results:\n", N,M,k);

	double MSE0, Ps, SNR; // MSE
	double X_est[N];

	long int my_diff;
	timespec init_time, final_time;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &init_time);
	
	simplx(Phi,M,N,M1,M2,M3,&ICASE,IZROV,IPOSV);
	
	if (ICASE==0) {  //result ok.
		//printf("\n ECONOMIC FUNCTION OPTIMIZED: = %f\n\n", Phi[1][1]);
		
		for (i=1; i<=N; i++) {
			for (j=1; j<=M; j++)
				if (IPOSV[j] == i)  {
					//printf("  X%d = %f\n", i, A[j+1][1]);
					goto e3;
				}
			//printf("  X%d = %f\n", i, 0.0);
			e3:;}
	}
	else
		printf(" No solution (error code = %d).\n", ICASE);
	
	printf("       M: %2d\n       N: %2d\n       K: %2d\n\n", M,N,k);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &final_time);
	my_diff = diff(init_time,final_time).tv_nsec;
	//Display execution time
	printf("       Execution time: %ld nanoseconds\n\n",my_diff);
	
	/* ********************************************************
	   Computing the Signal-to-Noise Ratio, SNR = 10*log(Ps/MSE)
	********************************************************* */
	
	//First compute the Mean Square Error, MSE

	for (i=1; i<=N; i++) {
		X_est[i-1] = 0.0;
	}
	
	for (i=1; i<=N; i++) {
		for (j=1; j<=M+1; j++)
			if (IPOSV[j] == i)  {
				//printf("       Approx X[%d]: %f\n", i, A[j+1][1]);
				X_est[i-1] = Phi[j+1][1];
				goto e3a;
			}
	e3a: ;
	}


	for (i=0; i<N; i++) {
		printf("       [%d] Actual X:[%.16f] Approx X:[%.16f]\n", i, X_act[i],X_est[i]);
	}

	
	MSE0=0;
	for (i=1; i<=N; i++) {
		MSE0+=pow(X_act[i]-X_est[i],2);
		//printf("\n       MSE grows %.16f\n", MSE);
	}
	//MSE=MSE/NV;
	printf("\n       Mean Square Error, MSE:              %.16f\n", MSE0/N);
	
	//Now compute the Signal Power, Ps
	
	Ps = 0;
	for (i=1; i<=N; i++) {
		Ps+=pow(X_act[i],2);
		//printf("\n       Ps grows %.16f\n", Ps);
	}
	//Ps=Ps/NV;
	printf("       Signal Power, Ps:                    %.16f\n", Ps/N);
	
	//Compute SNR = 10*log(Ps/MSE)
	
	SNR = 10*log10(Ps/MSE0);
	printf("       Signal-to-Noise Ratio, SNR:          %.16f dB\n", SNR);
	
	//For memory help see: http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process

	//Get system values
	getSystemValues();
	
}

// end of file tsimplexv2.cpp
