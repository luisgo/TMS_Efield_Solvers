clear all
FEMord=1; %choose FEM order

%%%%%%%%%%%%%%%%%%%%THIS CODE CURRENTLY ONLY WORKS IN LINUX (IT WILL WORK
%%%%%%%%%%%%%%%%%%%%IF NOT IN LINUX USE SLOW CODES
%%%%%%%%%%%%%%%%%%%%FAST CODES WILL WORK IN MAC AND WINDOWS BY  JUNE 1 2019.
%%%%%%%%%%%%choose points where you want to evalute the E-field (ro): ro is
%%%%%%%%%%%%3 by number of evaluation points.
[X,Y,Z]=sphere(200);
ro(1,:)=0.07*X(:).';
ro(2,:)=0.07*Y(:).';
ro(3,:)=0.07*Z(:).';

%four layer example
% The required data structures loaded from file are the following:
%       tetrahedron to node id (te2p): dimensions 4 nodes by number of tetrahedrons
%       node locations (p): dimensions 3 cartesian coordinate (x,y,z) by
%       number of nodes
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): dimensions 3 vector weight (Jx,Jy,Jz) by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values

load fourlayersphereexample.mat te2p p rs js conductivity;
rs(3,:)=rs(3,:)-min(rs(3,:))+.102; %move coil to be above head



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
E=FEMinterpolator(ro,x,te2p2,p2,rs,js,FEMord);

%% Step 5 plot results
Emag=sqrt(sum(E.^2,1));%compute Emagnitude
surf(0.07*X,0.07*Y,0.07*Z,reshape(Emag,size(X)),'edgealpha',0,'facecolor','interp');
