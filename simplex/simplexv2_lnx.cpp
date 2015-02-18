/*
Begin simplexv2.cpp
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
#include "simplexv2-N3-M3-K1.h"
//#include "simplexv2-N10-M6-K1.h"
//#include "simplexv2-N20-M11-K2.h"
//#include "simplexv2-N30-M11-K2.h"
//#include "simplexv2-N40-M11-K3.h"
//#include "simplexv2-N80-M16-K3.h"
//#include "simplexv2-N256-M80-K8.h"
using namespace std;

int NOPTIMAL,P1,P2,XERR;

void Data() {
	//printf("\n RESULTS:\n");
}

void Pivot();
void Formula();
void Optimize();

void Simplex() {
e10: Pivot();
	Formula();
	Optimize();
	if (NOPTIMAL == 1) goto e10;
}

void Pivot() {
	
	double RAP,V,XMAX;
	int I,J;
	
	XMAX = 0.0;
	for(J=2; J<=NV+1; J++) {
		if (TS[1][J] > 0.0 && TS[1][J] > XMAX) {
			XMAX = TS[1][J];
			P2 = J;
		}
	}
	RAP = 999999.0;
	for (I=2; I<=NC+1; I++) {
		if (TS[I][P2] >= 0.0) goto e10;
		V = fabs(TS[I][1] / TS[I][P2]);
		if (V < RAP) {
			RAP = V;
			P1 = I;
		}
		e10:;}
	V = TS[0][P2]; TS[0][P2] = TS[P1][0]; TS[P1][0] = V;
}

void Formula() {;
	int I,J;
	
	for (I=1; I<=NC+1; I++) {
		if (I == P1) goto e70;
		for (J=1; J<=NV+1; J++) {
			if (J == P2) goto e60;
			TS[I][J] -= TS[P1][J] * TS[I][P2] / TS[P1][P2];
			e60:;}
		e70:;}
	TS[P1][P2] = 1.0 / TS[P1][P2];
	for (J=1; J<=NV+1; J++) {
		if (J == P2) goto e100;
		TS[P1][J] *= fabs(TS[P1][P2]);
		e100:;}
	for (I=1; I<=NC+1; I++) {
		if (I == P1) goto e110;
		TS[I][P2] *= TS[P1][P2];
		e110:;}
}   

void Optimize() {
	int I,J;
	for (I=2; I<=NC+1; I++)
		if (TS[I][1] < 0.0)  XERR = 1;
	NOPTIMAL = 0;
	if (XERR == 1)  return;
	for (J=2; J<=NV+1; J++)
		if (TS[1][J] > 0.0)  NOPTIMAL = 1;
}

void Results() {
	int I,J;
	
	if (XERR == 0) goto e30;
	printf(" NO SOLUTION.\n"); goto e100;
e30:for (I=1; I<=NV; I++)
	for (J=2; J<=NC+1; J++) {
		if (TS[J][0] != 1.0*I) goto e70;
		//printf("       VARIABLE #%d: %.8f\n", I, TS[J][1]);
		e70:  ;}
	//printf("\n       ECONOMIC FUNCTION OPTIMIZED: %f\n", TS[1][1]);
e100:printf("\n");
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

int main()  {

	printf("\n       Program simplex version 1 %d-%d-%d results:\n", NV,NC,k);
	int I,J;
	double MSE, Ps, SNR;
	double X_est[NV];
	
	timespec time1, time2;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);

	Data();
	Simplex();
	Results();
	
	printf("       M: %2d\n       N: %2d\n       K: %2d\n\n", NC,NV,k);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
	//Display execution time
	printf("       Execution time: %ld nanoseconds\n\n",diff(time1,time2).tv_nsec);
	
	/* ********************************************************
	   Computing the Signal-to-Noise Ratio, SNR = 10*log(Ps/MSE)
	********************************************************* */
	//First compute the Mean Square Error, MSE
	
	for (I=1; I<=NV; I++) {
		X_est[I-1] = 0.0;
	}

	for (I=1; I<=NV; I++) {
		for (J=2; J<=NC+1; J++) {
			if (TS[J][0] != 1.0*I) goto e70;
			//printf("       Approx X[%d]: %.16f\n", I-1, TS[J][1]);
			X_est[I-1] = TS[J][1];
		e70: ;
		}
	}
	
	for (I=0; I<NV; I++) {
		printf("       [%d] Actual X:[%.16f] Approx X:[%.16f]\n", I, X_act[I],X_est[I]);
	}

	MSE=0;
	for (I=1; I<=NV; I++) {
		MSE+=pow(X_act[I]-X_est[I],2);
		//printf("\n       MSE grows %.16f\n", MSE);
	}
	//MSE=MSE/NV;
	printf("\n       Mean Square Error, MSE:              %.16f\n", MSE/NV);
	
	//Now compute the Signal Power, Ps
	
	Ps = 0;
	for (I=1; I<=NV; I++) {
		Ps+=pow(X_act[I],2);
	}
	//Ps=Ps/NV;
	printf("       Signal Power, Ps:                    %.16f\n", Ps/NV);
	
	//Compute SNR = 10*log(Ps/MSE)
	
	SNR = 10*log10(Ps/MSE);
	printf("       Signal-to-Noise Ratio, SNR:          %.16f dB\n", SNR);

	//For memory help see: http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
	
	//Get system values
	getSystemValues();

	return 0;
	
}

//end of file simplexv2.cpp
