# include "AStarOMPBuilder.h"

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

timespec time1, time2;

// constructor
AStarOMPBuilder::AStarOMPBuilder()
{
	mX = NULL;
	mY = NULL;
	mDict = NULL;

	mExRec = NULL;
	mNoExRecVec = NULL;	
	mErr = NULL;
	mNMSE = NULL;

	mBaseOMP = NULL;
	mBaseAStar = NULL;

	mBinOutput = true;
	mCompExactRec = true;

}

// function for initialization with config file
/// This function reads input parameters from pConfigFile and initializes A*OMP.
/// Target vectors, dictionary and measurements are read from the files given in pConfigFile.
int AStarOMPBuilder::init( std::string const pConfigFile )
{
	ConfigFile cf(pConfigFile);
	cout<<"Reading parameters from "<<pConfigFile<<"..."<<endl;
	mAlpha = (float) cf.Value("A*OMP_Parameters","alpha");  
	mBeta = (float) cf.Value("A*OMP_Parameters","beta");;  
	mI = (int) cf.Value("A*OMP_Parameters","I");  
	mB = (int) cf.Value("A*OMP_Parameters","B");  
	mP = (int) cf.Value("A*OMP_Parameters","P");  
	mK = (int) cf.Value("A*OMP_Parameters","K");  
	mEps = (float) cf.Value("A*OMP_Parameters","Eps");
	mInitPL = (int) cf.Value("A*OMP_Parameters","InitPL");
	mM = (int) cf.Value("Data_Parameters","M");  
	mN = (int) cf.Value("Data_Parameters","N");  
	mNoVectors = (int) cf.Value("Data_Parameters","NoVectors");  
	string AuxFuncMode = cf.Value("A*OMP_Parameters","myAuxiliaryFunctionMode");
	string DictFileName = cf.Value("Data_Parameters","DictFileName");
	string MeasurementsFileName = cf.Value("Data_Parameters","MeasurementsFileName");
	string TargetFileName = cf.Value("Data_Parameters","TargetFileName");
	mRecVectorsFileName = (string) cf.Value("OutputFiles","RecVectorsFileName");
	string ResultsFileName = cf.Value("OutputFiles","ResultFile");
	
	string dictMode = cf.Value("Data_Parameters","DictMode");

	if((int) cf.Value("Data_Parameters","ReadTargetVectors"))
		mTargetVectorsProvided = true;
	else 
		mTargetVectorsProvided = false;

	char* ResultsFile = new char[ResultsFileName.size()+1];
	strcpy(ResultsFile,ResultsFileName.c_str());
	mResultOfstream.open(ResultsFile, std::ios::trunc|ios::binary);
	if (mResultOfstream.fail()) 
	{       // files opened OK?
		cout << ResultsFileName<<"Cannot open result file..." << endl;
	}

	if(dictMode =="single")
	{
		mMultiDict = false;
	}
	else
	{
		if(dictMode =="multi")
		{
			mMultiDict = true;
		}
		else 
		{
			cout<<"Invalid DictMode in config file.should be single or multi."<<endl;
			cout<<"Terminating...";
			mResultOfstream<<"Invalid DictMode in config file... should be single or multi.";
			return 0;
		}
	}

	// set auxiliary function mode
	if(AuxFuncMode == "ADAP")
		mAuxiliaryFunctionMode = ADAP;
	if(AuxFuncMode== "ADAPMUL")
		mAuxiliaryFunctionMode = ADAPMUL;
	if(AuxFuncMode == "MUL")
		mAuxiliaryFunctionMode = MUL;

	// initialize data
	mDict = allocateFloatMatrix(mN,mM);
	mY = allocateFloatMatrix(mNoVectors,mM);

	// read the dictionary (phi)
	// read N vectors of size M (reads consecutive vectors from the binary file)
	if(dictMode =="single")
	{
		cout<<"Reading a single dictionary : ";
		if(!InitMatrixFromFile(mDict, DictFileName, mM, mN))
		{	
			cerr<<"Initialization Failed!..."<<endl;
			cerr<<"Terminating..."<<endl<<endl;
			if(mResultOfstream.is_open())
				mResultOfstream<<"Search terminated...";
			return 0;
		}
	}
	//else 
	//{
	//	cout<<"Individual dictionary for each measurement, dictionaries will be read online..."<<endl;
	//	char *filename = new char[DictFileName.size()+1];
	//	strcpy(filename,DictFileName.c_str());
	//	mDictIfstream.open(filename,ios::in | ios::binary);
	//	if(!mDictIfstream.is_open())
	//	{
	//		cout<<"cannot open dictionary file, terminating...";
	//		if(mResultOfstream.is_open())
	//		{
	//			mResultOfstream<<"Cannot open dictionary file. Search terminated.";
	//		}
	//		return 0;
	//	}	
	//}

	else 
	{
		cout<<"Individual dictionary for each measurement, dictionaries will be read online..."<<endl;
		char *filename = new char[DictFileName.size()+1];
		strcpy(filename,DictFileName.c_str());

		if(!strcmp(filename+DictFileName.size()-3,"bin")) 
		{
			mDictIfstream.open(filename,ios::in | ios::binary);
			mBinDictIfstream = true;
		}
		else
			if(!strcmp( filename+DictFileName.size()-3, "txt"))
			{
				mDictIfstream.open(filename);
				mBinDictIfstream = false;
			}
		if(!mDictIfstream.is_open())
		{
			cout<<"cannot open dictionary file, terminating...";
			if(mResultOfstream.is_open())
			{
				mResultOfstream<<"Cannot open dictionary file. Search terminated.";
			}
			return 0;
		}	
	}

	// read measurement vectors (y)
	// read noVectors vectors of size M (reads consecutive vectors from the binary file)
	cout<<"Reading measurement vectors : ";
	if(!InitMatrixFromFile(mY, MeasurementsFileName, mM, mNoVectors))
	{
		cerr<<"Initialization Failed!..."<<endl;
		cerr<<"Terminating..."<<endl<<endl;
		if(mResultOfstream.is_open())
			mResultOfstream<<"Search terminated...";
		return 0;
	}

	// read sparse target vectors (x)
	// read noVectors vectors of size N (reads consecutive vectors from the binary file)
	if(mTargetVectorsProvided)
	{
		mX =  allocateFloatMatrix(mNoVectors,mN);
		cout<<"Reading sparse target vectors : ";
		if(!InitMatrixFromFile(mX, TargetFileName, mN, mNoVectors))
		{
			deleteFloatMatrix(mX,mNoVectors);
			mX = NULL;
			mTargetVectorsProvided = false;
			cout<<"Continuing without sparse target vectors..."<<endl;
		}
		else
		{
			int inc = max(mNoVectors/20,1);
			for(int i = 0;  i<mNoVectors; i=i+inc)
			{
				if(countNonzeroElements(mX[i],mN) > mK)
				{
					mCompExactRec = false;
					cout<<"Target vectors have higher sparsity than K. Exact reconstruction rate will not be computed."<<endl;
					break;
				}
			}
		}
	}
	else
	{
		cout<<"No target vectors provided..."<<endl;
	}

	// for evaluation
	mExRec = new bool[mNoVectors];	
	mErr = new float[mN];
	mNMSE = new float[mNoVectors];
	mNoExRecVec = 0;
	mTime = 0;

	cout<<endl<<"Initializing A*OMP..."<<endl;
	mBaseOMP = new BaseOMP(mK,mM,mN,mEps,mInitPL);
	mBaseOMP->setDict(mDict);
	mBaseAStar = new BaseAStar(mB,mP,mI,mK,mN,mM,mAlpha, mBeta, mAuxiliaryFunctionMode);
	mBaseAStar->getAlgorithmInterface()->setProblem(mBaseOMP);

	if(mResultOfstream.is_open())
	{
		mResultOfstream<<"ALGORITHM PARAMETERS: \r"<<endl;
		mResultOfstream<<myIntend<<"Max. Non-zero components (K): "<<mK<<"\r"<<endl;
		mResultOfstream<<myIntend<<"Error Tolerance for termination (Eps): "<<mEps<<"\r"<<endl;
		mResultOfstream<<myIntend<<"Signal length (N):"<<mN<<"\r"<<endl;
		mResultOfstream<<myIntend<<"No observations (M):"<<mM<<"\r"<<endl;
		mResultOfstream<<myIntend<<"No Initial Branches (I):"<<mI<<"\r"<<endl;
		mResultOfstream<<myIntend<<"No Branches per Extension (B): "<<mB<<"\r"<<endl;
		mResultOfstream<<myIntend<<"No Maximum Paths in Stack (P): "<<mP<<"\r"<<endl;
	}
	cout<<myIntend<<"Max. Non-zero components (K): "<<mK<<endl;
	cout<<myIntend<<"Error Tolerance for termination (Eps): "<<mEps<<endl;
	cout<<myIntend<<"Signal length (N):"<<mN<<endl;
	cout<<myIntend<<"No observations (M):"<<mM<<endl;
	cout<<myIntend<<"No Initial Branches (I):"<<mI<<endl;
	cout<<myIntend<<"No Branches per Extension (B): "<<mB<<endl;
	cout<<myIntend<<"No Maximum Paths in Stack (P): "<<mP<<endl;
	switch(mAuxiliaryFunctionMode)
	{
	case MUL :
		{
			cout<<myIntend<<"Auxiliary Function Type: Multiplicative"<<endl;
			cout<<myIntend<<"Alpha: "<<mAlpha<<endl<<endl;
			if(mResultOfstream.is_open())
			{
				mResultOfstream<<myIntend<<"Auxiliary Function Type: Multiplicative\r"<<endl;
				mResultOfstream<<myIntend<<"Alpha: "<<mAlpha<<"\r"<<endl<<"\r"<<endl;
			}
			break;
		}

	case ADAP :
		{
			cout<<myIntend<<"Auxiliary Function Type: Adaptive-Additive"<<endl;
			cout<<myIntend<<"Beta: "<<mBeta<<endl<<endl;
			if(mResultOfstream.is_open())
			{
				mResultOfstream<<myIntend<<"Auxiliary Function Type: Adaptive-Additive\r"<<endl;
				mResultOfstream<<myIntend<<"Beta: "<<mBeta<<"\r"<<endl<<"\r"<<endl;
			}
			break;
		}
	case ADAPMUL :
		{
			cout<<myIntend<<"Auxiliary Function Type: Adaptive-Multiplicative"<<endl;
			cout<<myIntend<<"Alpha: "<<mAlpha<<endl<<endl;
			if(mResultOfstream.is_open())
			{
				mResultOfstream<<myIntend<<"Auxiliary Function Type: Adaptive-Multiplicative\r"<<endl;
				mResultOfstream<<myIntend<<"Alpha: "<<mAlpha<<"\r"<<endl<<"\r"<<endl;
			}
			break;
		}
	}

	char* RecVectorsFile = new char[mRecVectorsFileName.size()+1];
	strcpy(RecVectorsFile,mRecVectorsFileName.c_str());
	mRecVectOfstream.open(RecVectorsFile, std::ios::trunc|ios::binary);
	if (mRecVectOfstream.fail()) 
	{       // files opened OK?
		cerr<<"Failed to open "<< RecVectorsFile <<", terminating..." << endl;
		mResultOfstream<<"Failed to open "<< RecVectorsFile <<"\r"<<endl;
		cerr<<"Initialization Failed!..."<<endl;
		cerr<<"Terminating..."<<endl<<endl;
		return 0;
	}

	if(!strcmp(RecVectorsFile+mRecVectorsFileName.size()-3,"txt")) 
	{
		mBinOutput = false;
	}

	if(!mTargetVectorsProvided)
		mResultOfstream<<"No target vectors are provided. Evaluation cannot be performed...\r"<<endl<<endl;
	mNoIterations = 0;
	mNoEqBranch = 0;
	mNoBranchAdded = 0;
	mNoBranchIgnored = 0;
	mNoBranchReplaced = 0;

	return 1;

}

// function run AStarOMP algorithm
/// This function runs the A*OMP algorithm. Reconstruction is performed for all the test vectors loaded.
int AStarOMPBuilder::run()
{
	cout<<"Running A*OMP";
	for(int j = 0; j<mNoVectors;j++)
	{
		if(!(j%20))
		cout<<".";
		if(mMultiDict && mBinDictIfstream)
		{
		for(int i=0;i<mN;i++)
			mDictIfstream.read ((char*)mDict[i] , sizeof(float)*mM);
			mBaseOMP->setDict(mDict);
		}
		if(mMultiDict && !mBinDictIfstream)
		{
			freadtxt(&mDictIfstream, mDict, mN, mM);
			mBaseOMP->setDict(mDict);
		}

		mBaseOMP->sety(mY[j]);
		if(mBaseAStar->initialize())
		{
			//cout<<"Running AStar Search..."<<endl<<endl;

			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
			mBaseAStar->run();
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
			mTime = mTime + diff(time1,time2).tv_nsec;
			mTime = mTime/1000000000;
			
			// write reconstructed vector to ofstream
			float* mySol = (float*)mBaseAStar->getSolution();
			if (mRecVectOfstream.is_open())
			{
				if(mBinOutput)
					fwritebin(&mRecVectOfstream,mySol,mN);
				else
				{
					fwritetxt(&mRecVectOfstream,mySol,mN);
					mRecVectOfstream<<"\r"<<endl;
				}
			}
			
			if(mTargetVectorsProvided)
				evaluateSingleVector(j);
		}
		else
		{
			cerr<<endl<<"Initialization of BaseAStar failed in AStarOMPBuilder for vector "<<j<<endl;
			if(mResultOfstream.is_open())
			{
				mResultOfstream<<"Initialization of BaseAStar failed in AStarOMPBuilder for vector "<<j<<"\r"<<endl;
				mResultOfstream<<"Search terminated...";
			}
			return 0;
		}
	}
	cout<<endl<<endl<<"Reconstruction finished, reconstructed vectors written in "<<mRecVectorsFileName<<endl<<endl;
	if(mResultOfstream.is_open())
		mResultOfstream<<"Reconstruction successful, reconstructed vectors written in "<<mRecVectorsFileName<<"\r"<<endl;
	return 1;
}

// function to print evaluation results
/// This function prints the reconstruction results after A*OMP is run.
void AStarOMPBuilder::printEvaluation()
{
	//Evaluation
	// command line output
	cout<<endl<<"RECONSTRUCTION RESULTS:"<<endl;
	cout<<myIntend<<"No. Total Vectors: "<<mNoVectors<<endl;	
	cout<<myIntend<<"Average Normalized Mean Squared Error: "<<computeMean(mNMSE,mNoVectors)<<endl<<endl;
	if(mCompExactRec)
	{
		cout<<myIntend<<"No. Exactly Reconstructed Vectors: "<<mNoExRecVec<<endl;
		cout<<myIntend<<"Exact Reconstruction Rate: %"<<(float)100*mNoExRecVec/mNoVectors<<endl;
	}
	cout<<myIntend<<"A* Search Analysis Results: "<<endl;
	cout<<myIntend<<myIntend<<"Total Time: "<<mTime<<" sec."<<endl;
	cout<<myIntend<<myIntend<<"Statistics per vector:"<<endl;
	cout<<myIntend<<myIntend<<myIntend<<"Average Time: "<<mTime/mNoVectors<<" sec."<<endl;
	cout<<myIntend<<myIntend<<myIntend<<"No. Iterations: "<<(float)mNoIterations/mNoVectors<<endl;
	cout<<myIntend<<myIntend<<myIntend<<"No. Added Branches: "<<(float)mNoBranchAdded/mNoVectors<<endl;
	cout<<myIntend<<myIntend<<myIntend<<"No. Replaced Branches: "<<(float)mNoBranchReplaced/mNoVectors<<endl;
	cout<<myIntend<<myIntend<<myIntend<<"No. Equivalent Branches: "<<(float)mNoEqBranch/mNoVectors<<endl;
	cout<<myIntend<<myIntend<<myIntend<<"No. Ignored Branches: "<<(float)mNoBranchIgnored/mNoVectors<<endl<<endl;

	// file output
	if(mResultOfstream.is_open())
	{
		mResultOfstream<<"\r"<<endl<<"RECONSTRUCTION RESULTS:\r"<<endl;
		mResultOfstream<<myIntend<<"No. Total Vectors: "<<mNoVectors<<"\r"<<endl;
		mResultOfstream<<myIntend<<"Average Normalized Mean Squared Error: "<<computeMean(mNMSE,mNoVectors)<<"\r"<<endl;
		if(mCompExactRec)
		{
			mResultOfstream<<myIntend<<"No. Exactly Reconstructed Vectors: "<<mNoExRecVec<<"\r"<<endl;
			mResultOfstream<<myIntend<<"Exact Reconstruction Rate: %"<<(float)100*mNoExRecVec/mNoVectors<<"\r"<<endl;
		}
		mResultOfstream<<myIntend<<"A* Search Analysis Results: "<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<"Total Time: "<<mTime<<" sec."<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<"Statistics per vector:"<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<myIntend<<"Average Time: "<<mTime/mNoVectors<<" sec."<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<myIntend<<"No. Iterations: "<<(float)mNoIterations/mNoVectors<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<myIntend<<"No. Added Branches: "<<(float)mNoBranchAdded/mNoVectors<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<myIntend<<"No. Replaced Branches: "<<(float)mNoBranchReplaced/mNoVectors<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<myIntend<<"No. Equivalent Branches: "<<(float)mNoEqBranch/mNoVectors<<"\r"<<endl;
		mResultOfstream<<myIntend<<myIntend<<myIntend<<"No. Ignored Branches: "<<(float)mNoBranchIgnored/mNoVectors<<"\r"<<endl;
	}
}

// function to evaluate reconstruction of a single vector
/// This function evaluates A*OMP performance for each vector after its reconstruction.
/// @param VectorInd index of the test vector in the test vector matrix (mX).
void AStarOMPBuilder::evaluateSingleVector( int VectorInd )
{
	float* mySol = (float*)mBaseAStar->getSolution();
	SideInfo *mySideInfo = (SideInfo*)(mBaseAStar->getBestPath()->mSideInfo);
	//float *res = mySideInfo->mRes;

	// compute NMSE
	subtractVectorfromVector(mX[VectorInd],mySol,mErr,mN);
	mNMSE[VectorInd] = l2Norm(mErr,mN)/l2Norm(mX[VectorInd],mN);

	mNMSE[VectorInd] = mNMSE[VectorInd]*mNMSE[VectorInd];

	//if(mNMSE[VectorInd] > 0.00001)
	//	cout<<mNMSE[VectorInd]<<endl;

	if(mCompExactRec)
	{
		mExRec[VectorInd] = true;
		vector<unsigned int> IndList = mySideInfo->mIndList;

		for(int i = 0; i<mN; i++)
		{
			if(//(mX[VectorInd][i] == 0.0f && abs(mySol[i]) > 0.000001f)	||				// x is 0 but mySol is over tolerance!
				(mX[VectorInd][i] > 0.000002f && mySol[i] < 0.000002f) ||		// x over (+)tolerance, but not mySol
				(mX[VectorInd][i] < -0.000002f && mySol[i] > -0.000002f)	||		// x under (-)tolerance, but not mySol
				(abs(mX[VectorInd][i]) < 0.000001f && abs(mySol[i]) > 0.000002f))		// x within tolerance, but not mySol
			{		
				//cout<<(mX[VectorInd][i] > 0.000002f && mySol[i] < 0.000002f)<<endl;
				//cout<<(mX[VectorInd][i] < -0.000002f && mySol[i] > -0.000002f)<<endl;
				//cout<<(abs(mX[VectorInd][i]) == 0.0f && abs(mySol[i]) > 0.000002f)<<endl;
				//cout<<mX[VectorInd][i]<<" "<<mySol[i];
				//cout<<"Vector "<<VectorInd<<", NMSE = "<<mNMSE[VectorInd]<<endl;
				mExRec[VectorInd] = false;
				break;
			}
		}
		if(mExRec[VectorInd])
			mNoExRecVec++;
	}

	mNoIterations += mBaseAStar->getNoIterations();
	mNoEqBranch += mBaseAStar->getNoEqBranch(); 
	mNoBranchAdded += mBaseAStar->getNoBranchAdded();
	mNoBranchIgnored += mBaseAStar->getNoBranchIgnored();
	mNoBranchReplaced += mBaseAStar->getNoBranchReplaced();	
}

AStarOMPBuilder::~AStarOMPBuilder()
{
	if(mDict)
		deleteFloatMatrix(mDict,mN);
	if(mY)
		deleteFloatMatrix(mY,mNoVectors);
	if(mX)
		deleteFloatMatrix(mX,mNoVectors);
	if(mErr)
		delete mErr;
	if(mNMSE)
		delete mNMSE;
	if(mExRec)
		delete mExRec;
	if(mBaseAStar)
		delete mBaseAStar;
	if(mBaseOMP)
		delete mBaseOMP;
	if (mRecVectOfstream.is_open())
		mRecVectOfstream.close();

	if(mResultOfstream.is_open())
		mResultOfstream.close();

}

// function to read matrix of size pRows x pCols from the file pFileName)
/// This function reads a matrix of size pRows x pCols from the file pFileName and stores it into pMatrix. 
/// pFileName should either be a binary (.bin) file, in which the columns of the matrix are concatenated in floating point format, 
/// or a text (.txt) file, whose lines contain rows of the matrix.
/// @param pMatrix pointer to the matrix in which the read data will be stored
/// @param pFileName name of the file from which the vectors are read. 
/// @param pRows number of rows of the matrix
/// @param pCols number of columns of the matrix
/// @return 1 if file is read successfully, 0 if it fails.
int AStarOMPBuilder::InitMatrixFromFile(float** pMatrix, string &pFileName, int pRows, int pCols )
{
	char *filename = new char[pFileName.size()+1];
	strcpy(filename,pFileName.c_str());
	if(!strcmp(filename+pFileName.size()-3,"bin")) 
	{
		if(!freadbin(filename, pMatrix,pRows,pCols))
			return 0;
	}
	else
		if(!strcmp( filename+pFileName.size()-3, "txt"))
		{
			if(!freadtxt(filename,pRows,pCols,pMatrix))
				return 0;
		}
		else
		{
			cerr<<"Error reading "<<pFileName<<" in AStarOMPBuilder::InitMatrixFromFile : Unknown file type"<<endl
				<<"File should be either binary (.bin) or text (.txt)"<<endl;
			if(mResultOfstream.is_open())
			{
				mResultOfstream<<"Error reading "<<pFileName<<" in AStarOMPBuilder::InitMatrixFromFile : Unknown file type"<<endl
					<<"File should be either binary (.bin) or text (.txt)"<<endl;
			}
			return 0;
		}
		delete filename;
		return 1;
}

// Function to get mTargetVectorsProvided
/// @return mTargetVectorsProvided
bool AStarOMPBuilder::getTargetVectorsProvided()
{
	return mTargetVectorsProvided;
}

int AStarOMPBuilder::countNonzeroElements(float *pVector, int pSize)
{
	int K=0;
	for(int i = 0; i<pSize; i++)
	{
		if( abs(pVector[i]) > 0.0f )
			K++;
	}
	return K;
}

bool AStarOMPBuilder::isKSparse( float *pVector, int pSize, int pK )
{
	int K = 0;
	for(int i = 0; i<pSize; i++)
	{
		if( abs(pVector[i]) > 0.00000001f )
			K++;
	}
	if(K == pK)
		return true;
	else return false;
}
