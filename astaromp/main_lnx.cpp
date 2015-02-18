//Begin AStarOMP algorithm

#include "AStarOMPBuilder.h"

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

using namespace std;

void getSystemValues() {
        FILE* file = fopen("/proc/self/status", "r");
        //int result = -1;
        char line[50];
        while (fgets(line, 50, file) != NULL) {
		//printf("       Status file line - %s\n", line);
		//if (strncmp(line, "VmSize:", 7) == 0) printf("\n       Virtual memory size - %s", line);
		//if (strncmp(line, "VmHWM:", 6) == 0) printf("       Peak resident set size - %s", line);
		//if (strncmp(line, "VmRSS:", 6) == 0) printf("       Resident set size - %s", line);
		if (strncmp(line, "VmData:", 7) == 0) printf("       Size of data segment - %s", line);
		if (strncmp(line, "VmStk:", 6) == 0) printf("       Size of stack segment - %s", line);
		if (strncmp(line, "VmExe:", 6) == 0) printf("       Size of text segment - %s", line);
        	//if (strncmp(line, "VmSize:", 7) == 0) result = parseLine(line);
                //break;
        }
        fclose(file);
}

int main(int argc, char** argv)
{
	string confFileName;
	AStarOMPBuilder myBuilder;
	confFileName = "config.txt";
	bool waitForUserAfterTerm;
	for(int i = 1;i<argc;i++)
	{ 
		if (!strcmp(argv[i],"-nW"))
		{
			waitForUserAfterTerm = false;
		}
		else
		{
			if (!strcmp(argv[i],"-c"))
			{
				if(argc>i+1)
				{
					confFileName = argv[i+1];
					if (confFileName[0] == '-')
					{
						cout<<"Invalid config file name: "<< confFileName <<endl;
						return 0;
					}
					i++;
				}
				else
				{
					cout<<"argument -c should be followed by the config file name!"<<endl;
					return 0;
				}
			}
			else
			{
				//				cout<<"Unknown Argument "<<argv[i]<<", Argument is ignored!"<<endl;
				cout<<"Invalid Argument "<<argv[i]<<endl<<endl;
				cout<<"Valid Arguments"<<endl;
				cout<<"      -c filename : use config file 'filename'"<<endl;
				cout<<"      -nW : do not wait for key press after termination of the search"<<endl<<endl;
				return 0;
			}
		}
	}
	
	
	if(!myBuilder.init(confFileName))
	{
		
		//system("pause");
		
		cout<<"Cannot initialize builder from config.txt file"<<endl;
		return 0;
		
	}
	else
	{
		myBuilder.run();
		if(myBuilder.getTargetVectorsProvided())
			myBuilder.printEvaluation();
		
		//system("pause");
		
		//Get system values
		getSystemValues();		

		cout<<"Finished..."<<endl;
		return 1;
	}
};
