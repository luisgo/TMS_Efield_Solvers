clear all
FEMord=1; %choose FEM order
omega=1;
reflev=0;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
te2p=te2p';p=p';
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


%%%%%%%%%%%%choose points where you want to evalute the E-field (ro): ro is
%%%%%%%%%%%%3 by number of evaluation points.
[X,Y,Z]=sphere(200);
ro1(1,:)=0.07*X(:).';
ro1(2,:)=0.07*Y(:).';
ro1(3,:)=0.07*Z(:).';
Z2=.07:-.001:-.01;
ro2(3,:)=Z2;
ro2(1:2,:)=0;



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
E1=FEMinterpolator(ro1,x,te2p2,p2,rs,js,FEMord);
E2=FEMinterpolator(ro2,x,te2p2,p2,rs,js,FEMord);
E3=FEMinterpolator(ro3,x,te2p2,p2,rs,js,FEMord);
E4=FEMinterpolator(ro4,x,te2p2,p2,rs,js,FEMord);

%% Step 6 plot results
subplot(2,1,1),
Emag1=sqrt(sum(E1.^2,1));%compute Emagnitude
surf(0.07*X,0.07*Y,0.07*Z,reshape(Emag1,size(X)),'edgealpha',0,'facecolor','interp');
axis equal
axis off
colorbar
subplot(2,1,2),
Emag2=sqrt(sum(E2.^2,1));%compute Emagnitude
plot(Z2,Emag2/Emag2(1));
%% Compare errors
Esoln2(isnan(Esoln2))=0;
norm(Esoln1(:)-E1(:))/norm(Esoln1(:))
norm(Esoln2(:)-E2(:))/norm(Esoln2(:))
norm(Esoln3(:)-E3(:))/norm(Esoln3(:))
norm(Esoln4(:)-E4(:))/norm(Esoln4(:))