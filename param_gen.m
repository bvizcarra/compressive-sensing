% Given vector x with sparsity K, we compute Phi*x = y
% where Phi is an MxN matrix, x is an Nx1 vector and y is 
% an Mx1 vector. Then we generate C++ embeddable code.

% *********** Begin program ***********
% Creating parameters to be used for all algorithms

clear all;
N = 3;
M = 3;
K = 1;
Phi = rand(M,N)        % Create a random Phi matrix holding values from 0-1
x = zeros(N, 1);       % Initialize x vector
v = randperm(N)        % Create random permutation of integers from 1 to N
                       % acting as indices of x
x(v(1:K)) = rand(1,K)  % Create sparse x vector
y = Phi*x              % Create y vector
TS = zeros(M+2,N+2);
TS(3:M+2,3:N+2)=-Phi;  % Place -Phi matrix into TS
TS(3:M+2,2)=y;         % Place y vector into TS
TS(2,3:N+2)=1.0000;    % Place C transpose array into TS
TS
Phi_col_wise = reshape(Phi, 1, N*M)

% *********** Create header file for TSIMPLEX algorithm ***********

fname = strcat('tsimplexv2-N', num2str(N),'-M',num2str(M),'-K',num2str(K),'.h');
fid = fopen(fname, 'wt');
fprintf(fid, '#define  MMAX  %d',M+2);
fprintf(fid, '\n#define  NMAX  %d',N+2);
fprintf(fid, '\n#define  REAL  double');
fprintf(fid, '\nint k=%d;', K);

fprintf(fid, '\ndouble X_act[%d]= {',N);
for i=1:N
    fprintf(fid, '%2.6f,', x(i));    
end   
fprintf(fid, '};');

fprintf(fid, '\nint N=%d, M=%d, M1=%d, M2=0, M3=0;',N,M,M);
fprintf(fid, '\ntypedef REAL MAT[MMAX][NMAX];');

fprintf(fid, '\ndouble Phi[%d][%d] = {\n',M+2,N+2); % Build array declaration

for i=1:M+1
    fprintf(fid, '  {');  
    for j=1:N+2
        fprintf(fid, '%2.6f,',TS(i,j));  % Build remaining rows of array
    end
    fprintf(fid, '},\n');
end

for i=M+2:M+2
    fprintf(fid, '  {');    
    for j=1:N+2
        fprintf(fid, '%2.6f,',TS(i,j));   % Build final row of array
    end
    fprintf(fid, '}\n');
end
 
fprintf(fid, '};\n');
fclose(fid);

% *********** Create header file for SIMPLEX algorithm ***********

for i=1:N
    TS(1,i+2)=i;  %Label top row
end
for j=i+1:M+N
    TS(j+2-N,1)=j;  %Label left column
end

fname = strcat('simplexv2-N', num2str(N),'-M',num2str(M),'-K',num2str(K),'.h');
fid = fopen(fname, 'wt');
fprintf(fid, 'int k=%d;', K);

fprintf(fid, '\ndouble X_act[%d]= {',N);
for i=1:N
    fprintf(fid, '%2.6f,', x(i));    
end   
fprintf(fid, '};');

fprintf(fid, '\nint NV=%d, N=NV;',N);
fprintf(fid, '\nint NC=%d, M=NC;',M);

fprintf(fid, '\ndouble TS[%d][%d] = {\n',M+2,N+2); % Build array declaration

for i=1:M+1
    fprintf(fid, '  {');  
    for j=1:N+2
        fprintf(fid, '%2.6f,',TS(i,j));  % Build remaining rows of array
    end
    fprintf(fid, '},\n');
end

for i=M+2:M+2
    fprintf(fid, '  {');    
    for j=1:N+2
        fprintf(fid, '%2.6f,',TS(i,j));   % Build final row of array
    end
    fprintf(fid, '}\n');
end
 
fprintf(fid, '};\n');
fclose(fid);

% *********** Create header file for INTERIOR POINTS algorithm ***********

fname = strcat('interiorv1-N', num2str(N),'-M',num2str(M),'-K',num2str(K),'.h');
fid = fopen(fname, 'wt');

fprintf(fid, 'Doub X_act[%d]= {',N);
for i=1:N
    fprintf(fid, '%2.6f,', x(i));    
end   
fprintf(fid, '};');

fprintf(fid, '\nconst Int N=%d, M=%d, k=%d, Vvar=%d, Wvar=%d, Yvar=%d;', N,M,K,N+1,N*M,N*M);

fprintf(fid, '\nconst Doub my_Phi_arr[N*M] = {'); % our Phi matrix in column-wise fashion
for i=1:N*M
    fprintf(fid, '%2.6f,', Phi_col_wise(i));    
end   
fprintf(fid, '};');

fprintf(fid, '\nconst Int my_Phi_ptr_arr[N+1] = {'); % [0...N+1]
for i=0:N
    fprintf(fid, '%d,', i*M);    
end  
fprintf(fid, '};');

fprintf(fid, '\nconst Int my_Phi_row_arr[N*M] = {'); 
for i=1:N
    for j=0:M-1
        fprintf(fid, '%d,', j);
    end
end  
fprintf(fid, '};');

fprintf(fid, '\nconst Doub my_y_arr[M] = {'); % our constraint vector
for i=1:M
    fprintf(fid, '%2.6f,', y(i));    
end   
fprintf(fid, '};');

fprintf(fid, '\nconst Doub my_c_arr[N] = {'); % our objective function coefficients
for i=1:N
    fprintf(fid, '1.000000,');    
end   
 
fprintf(fid, '};\n');
fclose(fid);

% *********** Create text files for A STAR OMP algorithm ***********

fname = strcat('astaromp-x-N', num2str(N),'-M',num2str(M),'-K',num2str(K),'.txt');
fid = fopen(fname, 'wt');
for i=1:N
    fprintf(fid, '  %2.6f', x(i));    
end 
fclose(fid);

fname = strcat('astaromp-y-N', num2str(N),'-M',num2str(M),'-K',num2str(K),'.txt');
fid = fopen(fname, 'wt');
for i=1:M
    fprintf(fid, '  %2.6f', y(i));    
end 
fclose(fid);

fname = strcat('astaromp-Phi-N', num2str(N),'-M',num2str(M),'-K',num2str(K),'.txt');
fid = fopen(fname, 'wt');
for i=1:N*M
    fprintf(fid, '  %2.6f', Phi_col_wise(i));    
end 
fclose(fid);

% *********** End program ***********

fprintf('Program %d-%d-%d executed successfully.\n',N,M,K); 
% plot(x);