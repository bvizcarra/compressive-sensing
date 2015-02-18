extern "C" {
	#include "ldl.h"
	#include "amd.h"
}

struct NRldl {
// Interface between Numerical Recipes routine intpt and the required packages LDL and AMD
	Doub Info [AMD_INFO];
	Int lnz,n,nz;
	VecInt PP,PPinv,PPattern,LLnz,LLp,PParent,FFlag,*LLi;
	VecDoub YY,DD,*LLx;
	Doub *Ax, *Lx, *B, *D, *X, *Y;
	Int *Ai, *Ap, *Li, *Lp, *P, *Pinv, *Flag,*Pattern, *Lnz, *Parent;
	NRldl(NRsparseMat &adat);
	// Constructor only needs adat to have been declared with appropriate dimensions
	void order();
	// AMD ordering and LDL symbolic factorization. Only neds nonzero pattern of adat, not actual values
	void factorize();
	// Numerical factorization of matrix 
	void solve(VecDoub_O &y,VecDoub &rhs);
	// Solves for y given rhs. Can be invoked multiple times after a single call to factorize
	~NRldl();
};