/*
 Begin NRldlmain.cpp
 */
/*
 #include <stdio.h>
 #include <math.h>
 #include <time.h>
 #include <stdint.h>
 #include <sys/types.h>
 #include <sys/sysctl.h>
 #include <sys/resource.h>
 #include <mach/mach.h>
 #include <mach/vm_statistics.h>
 #include <mach/mach_types.h> 
 #include <mach/mach_init.h>
 #include <mach/mach_host.h>
 */
#include "nr3.h"
#include "sparse.h"
#include "NRldl.h"
#include "interiorv1-N3-M3-K1.h"
//#include "interiorv1-N10-M6-K1.h"
//#include "interiorv1-N20-M11-K2.h"
//#include "interiorv1-N30-M11-K2.h"
//#include "interiorv1-N40-M11-K3.h"
//#include "interiorv1-N80-M16-K3.h"
//#include "interiorv1-N256-M80-K8.h"

void init_Abc(Int M, Int N, NRsparseMat &a, const Int ptr_arr[], const Int row_arr[], const Doub A_arr[])
{
	a.nrows = M;
	a.ncols = N;
	a.nvals = M*N;
	VecInt v(11,ptr_arr), w(9,row_arr);
	VecDoub y(9,A_arr);
	a.col_ptr=v;
	a.row_ind=w;
	a.val=y;
}

Doub dotprod(VecDoub_I &x, VecDoub_I &y)
// Compute the dot product of two vectors, x dot y
{
	Doub sum=0.0;
	for (Int i=0;i<x.size();i++)
		sum += x[i]*y[i];
	return sum;
}

Int intpt(const NRsparseMat &a, VecDoub_I &b, VecDoub_I &c, VecDoub_O &x)
/*
 Interior point method for linear programming. On input a contains the coefficient matrix for the 
 constraints in the form of A dot x = b. The right-handed side of the contraints is input in b[0...m-1].
 The coefficients of the objective function to be minimized, c dot x, are input in c[0...n-1]. Note
 that c should generally be padded with zeros corresponding to the slack variables that extend
 the number of columns to be n. The function returns 0 if an optimal solution is found; 1 if 
 the problem is infeasible; 2 if the dual problem is infeasible, i.e., if the problem is unbounded
 or perhaps infeasiblele; and 3 if the number of iterations is exceeded.  The solution is returned in
 x[0...n-1].
 */
{
	const Int MAXITS=200;
	const Doub EPS=1.0e-6;
	const Doub SIGMA=0.9;
	const Doub DELTA=0.02;
	const Doub BIG=numeric_limits<Doub>::max();
	Int i,j,iter,status;
	Int m=a.nrows;
	Int n=a.ncols;
	VecDoub y(m),z(n),ax(m),aty(n),rp(m),rd(n),d(n),dx(n),dy(m),dz(n),
	rhs(m),tempm(m),tempn(n);
	NRsparseMat at=a.transpose();
	ADAT adat(a,at);
	NRldl solver(adat.ref());
	solver.order();
	Doub rpfact=1.0+sqrt(dotprod(b,b));
	Doub rdfact=1.0+sqrt(dotprod(c,c));
	for (j=0;j<n;j++) {
		x[j]=1000.0;
		z[j]=1000.0;
	}
	for (i=0;i<m;i++) {
		y[i]=1000.0;
	}
	Doub normrp_old=BIG;
	Doub normrd_old=BIG;
	cout << setw(4) << "iter" << setw(12) << "Primal obj." << setw(9) <<
	"||r_p||" << setw(13) << "Dual obj." << setw(11) << "||r_d||" <<
	setw(13) << "duality gap" << setw(16) << "normalized gap" << endl;
	cout << scientific << setprecision(4);
	for (iter=0;iter<MAXITS;iter++) {
		ax=a.ax(x);
		for (i=0;i<m;i++)
			rp[i]=ax[i]-b[i];
		Doub normrp=sqrt(dotprod(rp,rp))/rpfact;
		aty=at.ax(y);
		for (j=0;j<n;j++)
			rd[j]=aty[j]+z[j]-c[j];
		Doub normrd=sqrt(dotprod(rd,rd))/rdfact;
		Doub gamma=dotprod(x,z);
		Doub mu=DELTA*gamma/n;
		Doub primal_obj=dotprod(c,x);
		Doub dual_obj=dotprod(b,y);
		Doub gamma_norm=gamma/(1.0+abs(primal_obj));
	 	cout << setw(3) << iter << setw(12) << primal_obj << setw(12) <<
		normrp << setw(12) << dual_obj << setw(12) << normrd << setw(12)
		<< gamma << setw(12) << gamma_norm<<endl;
		if (normrp < EPS && normrd < EPS && gamma_norm < EPS)
			return status=0;
		if (normrp > 1000*normrp_old && normrp > EPS)
			return status=1;
		if (normrd > 1000*normrd_old && normrd > EPS)
			return status=2;
		for (j=0;j<n;j++)
			d[j]=x[j]/z[j];
		adat.updateD(d);
		solver.factorize();
		for (j=0;j<n;j++)
			tempn[j]=x[j]-mu/z[j]-d[j]*rd[j];
		tempm=a.ax(tempn);
		for (i=0;i<m;i++)
			rhs[i]=-rp[i]+tempm[i];
		solver.solve(dy,rhs);
		tempn=at.ax(dy);
		for (j=0;j<n;j++)
			dz[j]=-tempn[j]-rd[j];
		for (j=0;j<n;j++)
			dx[j]=-d[j]*dz[j]+mu/z[j]-x[j];
		Doub alpha_p=1.0;
		for (j=0;j<n;j++)
			if (x[j]+alpha_p*dx[j] < 0.0)
				alpha_p=-x[j]/dx[j];
		Doub alpha_d=1.0;
		for (j=0;j<n;j++)
			if (z[j]+alpha_d*dz[j] < 0.0)
				alpha_d=-z[j]/dz[j];
		alpha_p = MIN(alpha_p*SIGMA,1.0);
		alpha_d = MIN(alpha_d*SIGMA,1.0);
		for (j=0;j<n;j++) {
			x[j]+=alpha_p*dx[j];
			z[j]+=alpha_d*dz[j];
		}
		for (i=0;i<m;i++)
			y[i]+=alpha_d*dy[i];
		normrp_old=normrp;
		normrd_old=normrd;
	}
	return status=3;
}

Int main() {
	
	// Begin memory calculation
	/*
	 There are five types of memory pages in Mac OS X. They are as follows:
	 Wired pages that are locked in place and cannot be swapped out
	 Active pages that are loading into physical memory and would be relatively difficult to swap out
	 Inactive pages that are loaded into memory, but haven't been used recently and may not even be needed at all. These are potential candidates for swapping. This memory would probably need to be flushed.
	 Cached pages that have been some how cached that are likely to be easily reused. Cached memory probably would not require flushing. It is still possible for cached pages to be reactivated
	 Free pages that are completely free and ready to be used.
	 */
	
	// RAM Currently Used
	/*
	 vm_size_t page_size;
	 mach_port_t mach_port;
	 mach_msg_type_number_t count;
	 vm_statistics_data_t vm_stats;
	 long myFreeMemory, used_memory, myFreeMemoryEnd, used_memoryEnd;
	 
	 mach_port = mach_host_self();
	 count = sizeof(vm_stats) / sizeof(natural_t);
	 if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
	 KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO, 
	 (host_info_t)&vm_stats, &count))
	 {
	 myFreeMemory = (int64_t)vm_stats.free_count * (int64_t)page_size;
	 
	 used_memory = ((int64_t)vm_stats.active_count + 
	 (int64_t)vm_stats.inactive_count + 
	 (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
	 }
	 
	 //printf("       Free memory size: %ld\n", myFreeMemory);
	 //printf("       Used memory size: %ld\n", used_memory);
	 */
	
	Int i, my_Answer;
	NRsparseMat A;
	Doub time_end, MSE, Ps, SNR;
	
	/* #include
	 Doub X_act[3]= {0.000000,0.957167,0.000000};
	 const Int N=3, M=3, k=1;
	 const Doub my_A_arr[N*M] = {0.814724,0.905792,0.126987,0.913376,0.632359,
	 0.097540,0.278498,0.546882,0.957507}; // our A matrix in column-wise fashion
	 const Int my_A_ptr_arr[N+1] = {0,3,6,9}; // [0...N+1]
	 const Int my_A_row_arr[N*M] = {0,1,2,0,1,2,0,1,2};
	 const Doub my_b_arr[M] = {0.874253,0.605273,0.093362}; // our contraint vector
	 const Doub my_c_arr[N] = {1.000000,1.000000,1.000000}; // our objective function coefficients
	 */
	
	printf("Running program InteriorPoints %d-%d-%d\n", N,M,k);
	clock_t tStart = clock();
	VecDoub_I b(M,my_y_arr), c(N,my_c_arr);
	VecDoub_O my_x(N);
	init_Abc(M, N, A, my_Phi_ptr_arr, my_Phi_row_arr, my_Phi_arr);
	my_Answer = intpt(A, b, c, my_x);
	printf ("\nReturn status: %d (if 0, an optimal solution has been found)\n", my_Answer);
	
	//Display execution time
	time_end = (double)(clock() - tStart)/CLOCKS_PER_SEC;
	printf("\n       M: %2d\n       N: %2d\n       K: %2d\n\n", M,N,k);
	printf("       Execution time: %.9fs\n", time_end);
	
	/*
	 for (i=0; i<N; i++) {
	 printf("       [%d] Actual X:[%.16f] Approx X[%.16f]\n", i, X_act[i], my_x[i]);
	 }
	 */
	
	/******* Compute the Signal-to-Noise Ratio, SNR = 10*log(Ps/MSE) *******/
	
	//First compute the Mean Square Error, MSE
	MSE=0;
	for (i=0; i<N; i++) {
		MSE+=pow(X_act[i]-my_x[i],2);
		//printf("\n       MSE grows %.16f\n", MSE);
	}
	//MSE=MSE/NV;
	printf("\n       The MSE is   %.16f\n", MSE/N);
	
	//Now compute the Signal Power, Ps
	Ps = 0;
	for (i=0; i<N; i++) {
		Ps+=pow(X_act[i],2);
	}
	//Ps=Ps/NV;
	printf("       The Ps is    %.16f\n", Ps/N);
	
	//Compute SNR = 10*log(Ps/MSE)
	SNR = 10*log10(Ps/MSE);
	printf("       The SNR is %.16f dB\n", SNR);
	printf("\n");
	
	// End memory calculation
	/*
	 mach_port = mach_host_self();
	 count = sizeof(vm_stats) / sizeof(natural_t);
	 if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
	 KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO, 
	 (host_info_t)&vm_stats, &count))
	 {
	 myFreeMemoryEnd = (int64_t)vm_stats.free_count * (int64_t)page_size;
	 
	 used_memoryEnd = ((int64_t)vm_stats.active_count + 
	 (int64_t)vm_stats.inactive_count + 
	 (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
	 }
	 //printf("       Free memory used by my process: %ld\n", myFreeMemory-myFreeMemoryEnd);
	 printf("       Memory used by my process: %ld bytes\n", used_memoryEnd-used_memory);
	 
	 //For memory help see: http://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
	 
	 //Get total available RAM memory for Mac OSX
	 
	 int mib[2] = { CTL_HW, HW_MEMSIZE };
	 u_int namelen = sizeof(mib) / sizeof(mib[0]);
	 uint64_t size;
	 size_t len = sizeof(size);
	 
	 if (sysctl(mib, namelen, &size, &len, NULL, 0) < 0)
	 {
	 printf("System Error. ");
	 perror("sysctl");
	 }
	 else
	 {
	 printf("       Maximum available RAM = %llu bytes\n", size);
	 }
	 
	 // Virtual Memory Currently Used by my Process
	 
	 struct task_basic_info t_info;
	 mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
	 
	 if (KERN_SUCCESS != task_info(mach_task_self(),
								  TASK_BASIC_INFO, (task_info_t)&t_info, 
								  &t_info_count))
	{
		return -1;
	}
	
	printf("       Resident memory size: %ld\n", t_info.resident_size);
	printf("       Virtual memory size: %ld\n", t_info.virtual_size);
	
	*/
	
	return 0;
}

//end of file NRldlmain.cpp