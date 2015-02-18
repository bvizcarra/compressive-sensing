#include "VectorMath.h"

// Function to allocate a float matrix
/// This function allocates a float matrix of size pRows x pCols. It returns a pointer to 
/// a pointer array whose elements hold columns of the matrix.
/// @param pCols number of columns of the matrix
/// @param pRows number of rows of the matrix
/// @return pointer to the allocated matrix
float** allocateFloatMatrix( int pCols, int pRows )
{
	float** dstArray = new float*[pCols];
	for(int i = 0; i<pCols; i++)
	{
		dstArray[i] = new float[pRows];
		memset(dstArray[i],0,pRows*sizeof(float));
	};
	return dstArray;
}

// Function to allocate an integer matrix
/// This function allocates an integer matrix of size pRows x pCols. It returns a pointer to 
/// a pointer array whose elements hold columns of the matrix.
/// @param pCols number of columns of the matrix
/// @param pRows number of rows of the matrix
/// @return pointer to the allocated matrix
int** allocateIntMatrix(int pCols, int pRows)
{
	int** dstArray = new int*[pCols];
	for(int i = 0; i<pCols; i++)
	{
		dstArray[i] = new int[pRows];
		memset(dstArray[i],0,pRows*sizeof(int));
	};
	return dstArray;
}

// Function to delete a float matrix
/// This function deletes a float matrix with pCols columns.
/// @param pMatrix pointer to the matrix
/// @param pCols number of columns of the matrix
void deleteFloatMatrix( float** pMatrix, int pCols )
{
	if(pMatrix)
	{
		for(int i = 0; i<pCols; i++)
		{
			if(pMatrix[i])
				delete pMatrix[i];
		};
		delete pMatrix;
	}
}

// Function to delete an integer matrix
/// This function deletes an integer matrix with pCols columns.
/// @param pMatrix pointer to the matrix
/// @param pCols number of columns of the matrix
void deleteIntMatrix( int** pMatrix, int pCols )
{
	if(pMatrix)
	{
		for(int i = 0; i<pCols; i++)
		{
			if(pMatrix[i])
				delete pMatrix[i];
		};
		delete pMatrix;
	}
}

// Function to compute \f$l_2\f$ norm of a vector
/// This function returns \f$l_2\f$ norm of the vector pVector of size pSize.
/// @param pVector pointer to the array holding the vector
/// @param pSize size of pVector
/// @return \f$l_2\f$ norm of pVector
float l2Norm( float* pVector, int pSize )
{
	float res = 0;
	for(int i=0;i<pSize;i++)
		res += pow(pVector[i],2);

	return sqrt(res);
}

// Function to subtract one vector from another
/// This function subtracts pSrc2 from pSrc1 and returns the difference in pDst.
/// @param pSrc1 pointer to the minuend
/// @param pSrc2 pointer to the subtrahend
/// @param pDst pointer to the difference
/// @param pSize size of vectors
void subtractVectorfromVector( float* pSrc1, float* pSrc2, float* pDst,int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pDst[i] = pSrc1[i] - pSrc2[i];
	}
}

// Function to subtract one vector from another (in-place operation)
/// This function subtracts pSrc2 from pSrc1 and returns the difference in pSrc1.
/// @param pSrcDst pointer to the minuend, which is overwritten to contain the difference
/// @param pSrc2 pointer to the subtrahend
/// @param pSize size of vectors
void subtractVectorfromVector_I( float* pSrcDst, float* pSrc2, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pSrcDst[i] -= pSrc2[i];
	}

}
// Function to multiply a vector with a scalar
/// This function multiplies vector pSrc with the scalar pScalar and returns the result in pDst.
/// @param pSrc pointer to the source vector
/// @param pDst pointer to the destination vector
/// @param pScalar scalar value
/// @param pSize size of pSrc and pDst
void multVectorwithScalar( float *pSrc, float* pDst, float pScalar, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pDst[i] = pSrc[i]*pScalar;
	}
}

// Function to multiply a vector with a scalar
/// This function multiplies vector pSrc with the scalar pScalar and stores the result in pSrc.
/// @param pSrcDst pointer to the source-destination vector
/// @param pScalar scalar value
/// @param pSize size of pSrc and pDst
void multVectorwithScalar_I( float *pSrcDst, float pScalar, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pSrcDst[i] *= pScalar;
	}
}

// Function to sum two vectors (in-place operation)
/// This function adds pSrc2 to pSrcDst, overwriting the result on pSrcDst.
/// @param pSrcDst pointer to the source-destination vector
/// @param pSrc2 pointer to the second source vector
/// @param pSize size of vectors
void sumVectors_I( float *pSrcDst, float *pSrc2, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pSrcDst[i] += pSrc2[i];
	}
}

// Function to sum two vectors (in-place operation)
/// This function adds pSrc2 to pSrcDst, overwriting the result on pSrcDst.
/// @param pSrc1 pointer to the first source vector
/// @param pSrc2 pointer to the second source vector
/// @param pDst pointer to the destination vector
/// @param pSize size of vectors
void sumVectors( float *pSrc1, float *pSrc2, float* pDst, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pDst[i] = pSrc1[i] + pSrc2[i];
	}
}

// Function to divide a vector by a scalar (in-place operation)
/// This function divides the vector pSrcDst by the scalar value pScalar and overwrites the result to pSrcDst.
/// @param pSrcDst pointer to dividend and quotient
/// @param pScalar scalar divisor value
/// @param pSize size of vectors
void divideVectorByScalar_I( float *pSrcDst, float pScalar, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pSrcDst[i] = pSrcDst[i]/pScalar;
	}
}

// Function to divide a vector by a scalar (in-place operation)
/// This function divides the vector pSrcDst by the scalar value pScalar and overwrites the result to pSrcDst.
/// @param pSrc pointer to dividend
/// @param pDst pointer to quotient
/// @param pScalar scalar divisor value
/// @param pSize size of vectors
void divideVectorByScalar( float *pSrc, float* pDst, float pScalar, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pDst[i] = pSrc[i]/pScalar;
	}
}

// Function to copy a vector into another
/// This function copies the contents of pSrc into pDst.
/// @param pSrc pointer to the source vector
/// @param pDst pointer to the destination vector
/// @param pSize size of vectors
void copyVector( float* pSrc, float* pDst, int pSize )
{
	memcpy(pDst,pSrc,sizeof(float)*pSize);
}

// Function to copy a Matrix into another
/// This function copies the contents of pSrc into pDst.
/// @param pSrc pointer to the source matrix
/// @param pDst pointer to the destination matrix
/// @param pCols number of columns of the matrices
/// @param pRows number of rows of the matrices
void copyMatrix( float** pSrc, float** pDst, int pCols, int pRows )
{
	for(int i = 0;i<pCols; i++)
		memcpy(pDst[i],pSrc[i],sizeof(float)*pRows);
}


//map<float,int,greater<float>>* findMaxInnerProdIndList( float** pVectorArray, float* pVector, int pNoVectors, int pReturnSize, int pSize )
//{
//	map<float,int,greater<float>>* tempCost = new map<float,int,greater<float>>;
//	map<float,int,greater<float>>::iterator tempCostIter;
//	float tempCorr;
//	for(int i=0; i<pNoVectors;i++)
//	{
//		tempCorr = computeInnerProd(pVectorArray[i],pVector,pSize);
//		tempCost->insert(pair<float,int>(abs(tempCorr),i));
//		if((int)tempCost->size() > pReturnSize)
//		{
//			tempCostIter = tempCost->end();
//			tempCostIter--;
//			tempCost->erase(tempCostIter);
//		}
//	}
//	return tempCost;
//}

//// Function to find a sorted list of vectors among an array of vectors those have maximum inner-product to another vector
///// This function computes the inner products of pVector with the vectors stored in pVectorArray, and stores 
///// the indices of pReturnSize vectors which have maximum inner-product into pReturnList. Indices in pReturnList
///// are sorted wrt. decreasing inner-product.
///// @param pVectorArray pointer to the vector array
///// @param pVector pointer to the vector
///// @param pReturnList pointer to the array that stores returned index list
///// @param pNoVectors number of vectors in pVectorArray
///// @param pReturnSize number of vectors to be returned
///// @param pSize length of vectors
//void findMaxInnerProdIndList( float** pVectorArray, float* pVector,unsigned int* pReturnList ,int pNoVectors, int pReturnSize, int pSize )
//{
//	map<float,int,greater<float>> tempCost;
//	map<float,int,greater<float>>::iterator tempCostIter;
//	float tempCorr;
//	for(int i=0; i<pNoVectors;i++)
//	{
//		tempCorr = computeInnerProd(pVectorArray[i],pVector,pSize);
//		tempCost.insert(pair<float,int>(abs(tempCorr),i));
//		if((int)tempCost.size() > pReturnSize)
//		{
//			tempCostIter = tempCost.end();
//			tempCostIter--;
//			tempCost.erase(tempCostIter);
//		}
//	}
//	tempCostIter = tempCost.begin();
//	for(int i = 0; i<pReturnSize; i++, tempCostIter++)
//		pReturnList[i] = tempCostIter->second;
//}

// Function to compute inner-product of two vectors
/// This function returns the inner-product of pSrc1 and pScr2.
/// @param pSrc1 pointer to the first source vector
/// @param pSrc2 pointer to the second source vector
/// @param pSize length of vectors
/// @return inner-product of pSrc1 and pSrc2
float computeInnerProd( float* pSrc1, float* pSrc2, int pSize )
{
	float val = 0.0f;
	for(int i=0; i<pSize;i++)
		val += pSrc1[i]*pSrc2[i];
	return val;
}

/// Function to normalize a vector (in place)
/// This function normalizes pSrcDst (such that its \f$l_2\f$ norm is 1)
/// and stores the normalized vector in pSrcDst.
/// It returns \f$l_2\f$ norm of pSrcDst before normalization.
/// @param pSrcDst pointer to the source-destination vector
/// @param pSize length of pSrcDst
/// @return \f$l_2\f$ norm of pSrcDst before normalization
float normalizeVector_I( float* pSrcDst, int pSize )
{
	float normx = l2Norm(pSrcDst,pSize);
	divideVectorByScalar_I(pSrcDst,normx,pSize);
	return normx;
}

// Function to normalize a vector (in place)
/// This function normalizes pSrc (such that its \f$l_2\f$ norm is 1) 
/// and stores the normalized vector in pDst.
/// It returns \f$l_2\f$ norm of pSrcVec.
/// @param pSrc pointer to the source vector
/// @param pDst pointer to the (normalized) destination vector
/// @param pSize length of vectors
/// @return \f$l_2\f$ norm of pSrc
float normalizeVector( float* pSrc, float* pDst, int pSize )
{
	float normx = l2Norm(pSrc,pSize);
	divideVectorByScalar(pSrc,pDst,normx,pSize);
	return normx;
}

// Function to subtract product of a scalar and a vector from a vector (in-place operation)
/// This function multiplies pSrc2 with the scalar pScalar and subtracts the result from
/// pSrcDst, storing the result over pSrcDst.
/// @param pSrcDst pointer to the source-destination vector
/// @param pSrc2 pointer to the second source vector
/// @param pScalar scalar multiplier value
/// @param pSize length of vectors
void subtractProductScalarfromVector_I( float* pSrcDst, float* pSrc2, float pScalar, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pSrcDst[i] -= (pSrc2[i]*pScalar);
	}
}

// Function to subtract product of a scalar and a vector from a vector
/// This function multiplies pSrc2 with the scalar pScalar and subtracts the result from
/// pSrc1. Result is stored in pDst.
/// @param pSrc1 pointer to the first source  vector
/// @param pSrc2 pointer to the second source vector
/// @param pDst pointer to the destination vector
/// @param pScalar scalar multiplier value
/// @param pSize length of vectors
void subtractProductScalarfromVector( float* pSrc1, float* pSrc2, float *pDst, float pScalar, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pDst[i] = pSrc1[i] - (pSrc2[i]*pScalar);
	}
}

// Function to solve a linear system with upper triangular (UT) system matrix
/// This function solves x from a linear system Ax=b, where A is an upper-triangular(UT)
/// square matrix. 
/// @param pA pointer to the UT square matrix A
/// @param pb pointer to the vector b
/// @param px pointer to the vector x
/// @param pNoCols number of columns (and rows) of A, length of b and x
void linsolve_UT( float **pA, float *pb, float *px, int pNoCols )
{
	for(int i = pNoCols-1;i<0;i--)
	{
		//x(i) = b(i)/A(i,i);
		px[i] = pb[i]/pA[i][i];
		//b(1:i-1) = b(1:i-1) - A(1:i-1,i)*x(i);
		subtractProductScalarfromVector_I(pb, pA[i],px[i],i);
	}
	px[0] = pb[0]/pA[0][0];
}

// Function to subtract element-wise product of two vectors from a vector (in-place operation)
/// This function computed elementwise multiplication of pSrc1 and pSrc2 and subtracts the result from pSrcDst.
/// Result is written over pSrcDst.
/// @param pSrcDst pointer to the source-destination vector
/// @param pSrc1 pointer to the first source vector
/// @param pSrc2 pointer to the second source vector
/// @param pSize length of vectors
void subtractProductfromVector_I( float* pSrcDst, float* pSrc1, float* pSrc2, int pSize )
{
	for(int i=0;i<pSize;i++)
	{
		pSrcDst[i] -= (pSrc1[i]*pSrc2[i]);
	}
}

// Function to compute mean of a vector
/// This function returns mean value of pSrc.
/// @param pSrc pointer to the source vector
/// @param pSize length of pSrc
/// @return mean of pSrc
float computeMean( float *pSrc, int pSize )
{
	return sumVector(pSrc,pSize)/(float)pSize;
}

// Function to sum elements of a vector
/// This function returns the sum of elements of pSrc.
/// @param pSrc pointer to the source vector
/// @param pSize length of pSrc
/// @return sum of elements of pSrc
float sumVector( float *pSrc, int pSize )
{
	float sum = pSrc[0];
	for(int i = 1; i<pSize;i++)
		sum+=pSrc[i];
	return sum;
}
