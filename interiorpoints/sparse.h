#include "sort.h"

struct NRsparseCol							//sparse vector data structure
{
	Int nrows;							//number of rows
	Int nvals;							//maximum number of nonzeros
	VecInt row_ind;							//row indices of nonzeros
	VecDoub val;							//array of nonzero values

	NRsparseCol(Int m,Int nnvals) : nrows(m), nvals(nnvals),
	row_ind(nnvals,0),val(nnvals,0.0) {}				//constructor. initializes vector to zero.

	NRsparseCol() : nrows(0),nvals(0),row_ind(),val() {}		//default constructor

	void resize(Int m, Int nnvals) {
		nrows = m;
		nvals = nnvals;
		row_ind.assign(nnvals,0);
		val.assign(nnvals,0.0);
	}

};
struct NRsparseMat							//sparse matrix data structure for compressed column storage
{
	Int nrows;							//number of rows
	Int ncols;							//number of columns
	Int nvals;							//maximum number of nonzeros
	VecInt col_ptr;							//pointers to start of columns. length is ncols+1
	VecInt row_ind;							//row indices of nonzeros
	VecDoub val;							//array of nonzero values

	NRsparseMat();							//default constructor
	NRsparseMat(Int m,Int n,Int nnvals);				//constructor. initializes vector to zero
	VecDoub ax(const VecDoub &x) const;				//multiply A by a vector x[0…ncols-1]
	VecDoub atx(const VecDoub &x) const;				//multiply A transpose by a vector x[0…nrows-1]
	NRsparseMat transpose() const;					//form A transpose
};
NRsparseMat::NRsparseMat() : nrows(0),ncols(0),nvals(0),col_ptr(),
	row_ind(),val() {}
NRsparseMat::NRsparseMat(Int m,Int n,Int nnvals) : nrows(m),ncols(n),
	nvals(nnvals),col_ptr(n+1,0),row_ind(nnvals,0),val(nnvals,0.0) {}
VecDoub NRsparseMat::ax(const VecDoub &x) const {
	VecDoub y(nrows,0.0);
	for (Int j=0;j<ncols;j++) {
		for (Int i=col_ptr[j];i<col_ptr[j+1];i++)
			y[row_ind[i]] += val[i]*x[j];
	}
	return y;
}
VecDoub NRsparseMat::atx(const VecDoub &x) const {
	VecDoub y(ncols);
	for (Int i=0;i<ncols;i++) {
		y[i]=0.0;
		for (Int j=col_ptr[i];j<col_ptr[i+1];j++)
			y[i] += val[j]*x[row_ind[j]];
	}
	return y;
}
NRsparseMat NRsparseMat::transpose() const {
	Int i,j,k,index,m=nrows,n=ncols;
	NRsparseMat at(n,m,nvals);
	VecInt count(m,0);
	for (i=0;i<n;i++)
		for (j=col_ptr[i];j<col_ptr[i+1];j++) {
			k=row_ind[j];
			count[k]++;
		}
	for (j=0;j<m;j++)
		at.col_ptr[j+1]=at.col_ptr[j]+count[j];
	for(j=0;j<m;j++)
		count[j]=0;
	for (i=0;i<n;i++)
		for (j=col_ptr[i];j<col_ptr[i+1];j++) {
			k=row_ind[j];
			index=at.col_ptr[k]+count[k];
			at.row_ind[index]=i;
			at.val[index]=val[j];
			count[k]++;
		}
	return at;
}
struct ADAT {
	const NRsparseMat &a,&at;
	NRsparseMat *adat;

	ADAT(const NRsparseMat &A,const NRsparseMat &AT);
	void updateD(const VecDoub &D);
	NRsparseMat &ref();
	~ADAT();
};
ADAT::ADAT(const NRsparseMat &A,const NRsparseMat &AT) : a(A), at(AT) {
	Int h,i,j,k,l,nvals,m=AT.ncols;
	VecInt done(m);
	for (i=0;i<m;i++)
		done[i]=-1;
	nvals=0;
	for (j=0;j<m;j++) {
		for (i=AT.col_ptr[j];i<AT.col_ptr[j+1];i++) {
			k=AT.row_ind[i];
			for (l=A.col_ptr[k];l<A.col_ptr[k+1];l++) {
				h=A.row_ind[l];
				if (done[h] != j) {
					done[h]=j;
					nvals++;
				}
			}
		}
	}
	adat = new NRsparseMat(m,m,nvals);
	for (i=0;i<m;i++)
		done[i]=-1;
	nvals=0;
	for (j=0;j<m;j++) {
		adat->col_ptr[j]=nvals;
		for (i=AT.col_ptr[j];i<AT.col_ptr[j+1];i++) {
			k=AT.row_ind[i];
			for (l=A.col_ptr[k];l<A.col_ptr[k+1];l++) {
				h=A.row_ind[l];
				if (done[h] != j) {
					done[h]=j;
					adat->row_ind[nvals]=h;
					nvals++;
				}
			}
		}
	}
	adat->col_ptr[m]=nvals;
	for (j=0;j<m;j++) {
		i=adat->col_ptr[j];
		Int size=adat->col_ptr[j+1]-i;
		if (size > 1) {
			VecInt col(size,&adat->row_ind[i]);
			sort(col);
			for (k=0;k<size;k++)
				adat->row_ind[i+k]=col[k];
		}
	}
}
void ADAT::updateD(const VecDoub &D) {
	Int h,i,j,k,l,m=a.nrows,n=a.ncols;
	VecDoub temp(n),temp2(m,0.0);
	for (i=0;i<m;i++) {
		for (j=at.col_ptr[i];j< at.col_ptr[i+1];j++) {
			k=at.row_ind[j];
			temp[k]=at.val[j]*D[k];
		}
		for (j=at.col_ptr[i];j<at.col_ptr[i+1];j++) {
			k=at.row_ind[j];
			for (l=a.col_ptr[k];l<a.col_ptr[k+1];l++) {
				h=a.row_ind[l];
				temp2[h] += temp[k]*a.val[l];
			}
		}
		for (j=adat->col_ptr[i];j<adat->col_ptr[i+1];j++) {
			k=adat->row_ind[j];
			adat->val[j]=temp2[k];
			temp2[k]=0.0;
		}
	}
}
NRsparseMat & ADAT::ref() {
	return *adat;
}

ADAT::~ADAT() {
	delete adat;
}
