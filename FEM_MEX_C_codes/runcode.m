function [Efield,x,te2p2,p2,te2te]=runcode(te2p,p,conductivity,rs,js,ro,FEMord)
%te2p is 3 by nte
%p is 3 by np
%conductivity nte by 1
%rs is 3 by nsource
%js is 3 by nsource
%ro is 3 by nobs
%FEMord is the order of the FEM
%step 1 geometric preprocessing
tic
[te2p2,p2]=femgenmesh_c(te2p,p,FEMord);
[te2te]=gente2te(te2p);
toc
%% Step 2 assemble FEM matrix
tic
A=femassemble(te2p2,p2,conductivity,FEMord);
toc
%% Step 3 generate right hand side of equation
tic
[rhs]=femgenrhs(te2p2,p2,conductivity,rs,js,FEMord);
toc
%% Step 4 delete one unknown and equation and define preconditioner
A=A(1:end-1,1:end-1);
PRECON=sparse(1:numel(A(:,1)),1:numel(A(:,1)),sqrt(full(diag(A))));
rhs=rhs(1:end-1);
%% Step 5 solve system of equations
tic
x=minres(A,rhs,10^-10,10000,PRECON,PRECON);
x(end+1)=0;
toc
%% Step 6 evaluate field at desired locations
Efield=FEMinterpolator(te2te,te2p2,p2,rs,js,x,ro,FEMord);




