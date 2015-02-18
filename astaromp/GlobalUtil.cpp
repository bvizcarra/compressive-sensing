#include "GlobalUtil.h"

// Function to read a float array from a binary file. 
/// This function reads the float array pFloatReadArray of size pN from the binary file "pFilename".
/// @param pFilename name of the binary file
/// @param pFloatReadArray	pointer to the float array which will be read
/// @param pN number of float values to be read
/// @return 1 if file is read successfully, 0 otherwise
int freadbin(const char * pFilename, float* pFloatReadArray, int pN)
{
	if ((pFloatReadArray != NULL) && (pFilename != NULL))
	{		
		ifstream file(pFilename, ios::in | ios::binary);
		if (file.fail())
		{       // files opened OK?
			std::cerr <<  pFilename<<" FILE OPEN FAILED" << endl;
			return 0;
		}
		file.read ((char*)pFloatReadArray , sizeof(float)*pN);
		cout <<pFilename<<" is successfully read"<<endl;
		file.close();
		return 1;
	}
	else
	{
		cerr<<"NULL CHECK IS FAILED IN  readFromFile"<<endl;
		return 0;
	}
}

// Function to read a float matrix from a binary file.
/// This function reads the float matrix pFloatReadMatrix of size pRows x pCols from the binary file "pFilename".
/// The matrix is read column by column from the file, i.e. it should be stored in a column-concatenated manner.
/// @param pFilename name of the binary file
/// @param pFloatReadMatrix	pointer to the matrix which will be read
/// @param pRows number of rows of the matrix
/// @param pCols number of columns of the matrix
/// @return 1 if file is read successfully, 0 otherwise
int freadbin( const char * pFilename, float** pFloatReadMatrix, int pRows, int pCols )
{
	if ((pFloatReadMatrix != NULL) && (pFilename != NULL))
	{		
		ifstream file(pFilename, ios::in | ios::binary);
		if (file.fail())
		{       // files opened OK?
			std::cerr <<  pFilename<<" FILE OPEN FAILED" << endl;
			return 0;
		}
		for(int i=0;i<pCols;i++)
			file.read ((char*)pFloatReadMatrix[i] , sizeof(float)*pRows);
		cout <<pFilename<<" is successfully read"<<endl;
		file.close();
		return 1;
	}
	else
	{
		cerr<<"NULL CHECK IS FAILED IN  readFromFile"<<endl;
		return 0;
	}
}

// Function to read a float matrix from a text file.
/// This function reads the float matrix pDstMatrix of size pRows x pCols from the binary file "pFilename".
/// The matrix is read column by column from the file, i.e. it should be stored in a column-concatenated manner.
/// @param pFilename name of the binary file
/// @param pDstMatrix pointer to the matrix which will be read
/// @param pNoRows number of rows of the matrix
/// @param pNoCols number of columns of the matrix
/// @return 1 if file is read successfully, 0 otherwise
int freadtxt(char* pFilename, int pNoRows, int pNoCols, float** pDstMatrix)
{
	ifstream myfile(pFilename);
	int rows = 0;
	int cols = 0;
	if(!myfile) //Always test the file open.
	{
		cout<<"Error opening output file :"<<pFilename<<endl;
		return 0;
	}
	while(!myfile.eof() && cols<pNoCols)
	{
		myfile>>pDstMatrix[cols][rows];
		rows++;
		if(rows == pNoRows)
		{
			rows=0;
			cols++;
		}
	}
	myfile.close();
	return 1;
}

// Function to display contents of a float matrix on command line. 
/// This function outputs the float matrix pMatrix of size pRows x pCols row by row on the command line.
/// @param pMatrix pointer to the matrix which will be output
/// @param pNoRows number of rows of the matrix
/// @param pNoCols number of columns of the matrix
void displayMatrix( float** pMatrix, int pNoRows, int pNoCols )
{
	for(int i = 0; i<pNoRows;i++)
	{
		for(int j = 0; j<pNoCols;j++)
		{
			cout<<pMatrix[i][j]<<" ";
		}
		cout<<endl;
	}
}

// Function to display contents of an integer matrix on command line. 
/// This function outputs the integer matrix pMatrix of size pRows x pCols row by row on the command line.
/// @param pMatrix pointer to the matrix which will be output
/// @param pNoRows number of rows of the matrix
/// @param pNoCols number of columns of the matrix
void displayMatrix( int** pMatrix, int pNoRows, int pNoCols )
{
	for(int i = 0; i<pNoRows;i++)
	{
		for(int j = 0; j<pNoCols;j++)
		{
			cout<<pMatrix[i][j]<<" ";
		}
		cout<<endl;
	}
}

// Function to delete a vector of pointers to float arrays.
/// @param pVector vector to be deleted
void deletevector( vector<float*> &pVector )
{
	for(int i=0; i<(int)pVector.size();i++)
	{
		if(pVector.at(i))
			delete pVector.at(i);
	}
}

// Function to delete a vector of pointers to integer arrays.
/// @param pVector vector to be deleted
void deletevector( vector<int*> &pVector )
{
	for(int i=0; i<(int)pVector.size();i++)
		delete pVector.at(i);
}

// Function to delete a vector of pointers of type void*.
/// @param pVector vector to be deleted
void deletevector( vector<void*> &pVector )
{
	for(int i=0; i<(int)pVector.size();i++)
		delete pVector.at(i);
	pVector.clear();
}

void deletevector( vector<unsigned int*> &pVector )
{
	for(int i=0; i<(int)pVector.size();i++)
		delete pVector.at(i);
	pVector.clear();
}

// Function to write a float array to a text file. 
/// This function writes the float array pWriteFloatArray of size pSize to the text file "pFilename".
/// @param pFilename name of the binary file
/// @param pWriteFloatArray	pointer to the float array which will be written
/// @param pSize size of the array
/// @return 1 if file is write successfully, 0 otherwise
int fwritetxt(const char * pFilename, float* pWriteFloatArray,int pSize)
{
	if ((pWriteFloatArray !=0 ) && (pFilename != NULL))
	{
		ofstream file;
		file.open(pFilename, std::ios::trunc|ios::binary);
		if (file.fail()) 
		{       // files opened OK?
			cerr << pFilename<<" FILE OPEN FAILED" << endl;
			return 1;
		}

		for(int i=0;i<pSize;i++)
			file<< *(pWriteFloatArray+i)<<"  ";
		file.close() ;
		return 1;
	}
	else
	{
		cerr<<"NULL CHECK FAILED IN  fwrite"<<endl;
		return 0;
	}
};

// Function to write a float array to a binary file. 
/// This function writes the float array pWriteFloatArray of size pSize to the binary file "pFilename".
/// @param pFilename name of the binary file
/// @param pWriteFloatArray	pointer to the float array which will be written
/// @param pSize size of the array
/// @return 1 if file is write successfully, 0 otherwise
int fwritebin(const char * pFilename, float* pWriteFloatArray,int pSize)
{
	if ((pWriteFloatArray !=0 ) && (pFilename != NULL))
	{
		ofstream file;
		file.open(pFilename, std::ios::trunc|ios::binary);
		if (file.fail()) 
		{       // files opened OK?
			cerr << pFilename<<" FILE OPEN FAILED" << endl;
			return 1;
		}

		file.write((char *)pWriteFloatArray, sizeof(float)*pSize);
		file.close() ;
		return 1;
	}
	else
	{
		cerr<<"NULL CHECK FAILED IN  fwrite"<<endl;
		return 0;
	}
};

// Function to write a float array to an ofstream in binary format. 
/// This function writes the float array pWriteFloatArray of size pSize to the binary ofstream pFile.
/// @param pFile binary ofstream to write the array
/// @param pWriteFloatArray	pointer to the float array which will be written
/// @param pSize size of the array
/// @return 1 if file is write successfully, 0 otherwise
int fwritebin(ofstream* pFile, float* pWriteFloatArray,int pSize)
{
	if ((pWriteFloatArray !=0 ) && (pFile != NULL))
	{
		pFile->write((char *)pWriteFloatArray, sizeof(float)*pSize);
		return 1;
	}
	else
	{
		cerr<<"NULL CHECK FAILED IN  fwrite"<<endl;
		return 0;
	}
};

// Function to write a float array to an ofstream in text format. 
/// This function writes the float array pWriteFloatArray of size pSize to the text ofstream pFile.
/// @param pFile text ofstream to write the array
/// @param pWriteFloatArray	pointer to the float array which will be written
/// @param pSize size of the array
/// @return 1 if file is write successfully, 0 otherwise
int fwritetxt(ofstream* pFile, float* pWriteFloatArray,int pSize)
{
	if ((pWriteFloatArray !=0 ) && (pFile != NULL))
	{
		for(int i=0;i<pSize;i++)
			*pFile<<*(pWriteFloatArray+i)<<"  ";
		return 1;
	}
	else
	{
		cerr<<"NULL CHECK FAILED IN  fwrite"<<endl;
		return 0;
	}
};

// Function to read a float matrix from a text file.
/// This function reads the float matrix pDstMatrix of size pRows x pCols from the binary file "pFilename".
/// The matrix is read column by column from the file, i.e. it should be stored in a column-concatenated manner.
/// @param pFile pointer to the text ifstream
/// @param pDstMatrix pointer to the matrix which will be read
/// @param pNoRows number of rows of the matrix
/// @param pNoCols number of columns of the matrix
/// @return 1 if file is read successfully, 0 otherwise
int freadtxt(ifstream *pFile, float **pDstMatrix ,int pNoCols, int pNoRows)
{
	if ((pDstMatrix !=0 ) && (pFile != NULL))
	{
		int rows = 0;
		int cols = 0;
		while(!pFile->eof() && cols<pNoCols)
		{
			*pFile>>pDstMatrix[cols][rows];
			rows++;
			if(rows == pNoRows)
			{
				rows=0;
				cols++;
			}
		}
	}
	else
	{
		cerr<<"NULL CHECK FAILED IN  fwrite"<<endl;
		return 0;
	}
	return 1;
}