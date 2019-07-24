function [Efield,x,te2p2,p2,rs,js]=runcode(te2p,p,rs,js,conductivity,ro,FEMord)
%four layer example
% The required data structures loaded from file are the following:
%       FEMord is the order of FEM desired
%       tetrahedron to node id (te2p): dimensions number of tetrahedrons by
%       4 nodes
%       node locations (p): dimensions number of nodes by 3 cartesian coordinate (x,y,z)
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): dimensions 3 vector weight (Jx,Jy,Jz) by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values
%       ro is 3 by number of observation points
%       Efield is the Efield at points ro
te2p=te2p';p=p';
[te2p2,p2]=femgenmesh(te2p,p,FEMord);%update tetrahedron mesh to have the correct order
%% Step 1 assemble FEM matrix
tic
A=femassemble(te2p2,p2,conductivity,FEMord);
toc
%% Step 2 generate right hand side of equation
tic
[rhs]=femgenrhs(te2p2,p2,conductivity,rs,js,FEMord);
toc
%% Step 3 delete one unknown and equation and define preconditioner
A=A(1:end-1,1:end-1);
PRECON=sparse(1:numel(A(:,1)),1:numel(A(:,1)),sqrt(full(diag(A))));
rhs=rhs(1:end-1);
%% Step 4 solve system of equations
x=minres(A,rhs,10^-10,10000,PRECON,PRECON);
x(end+1)=0;
%% Step 5 evaluate field at desired locations
Efield=FEMinterpolator(ro,x,te2p2,p2,rs,js,FEMord);
