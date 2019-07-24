function [Efield,t2p,p,epseff,xval]=runcode(te2p,p,rs,js,conductivity,ro)
%zeroth order BEM code
%       FEMord is the order of FEM desired
%       tetrahedron to node id (te2p): dimensions number of tetrahedrons by
%       4 nodes
%       node locations (p): dimensions number of nodes by 3 cartesian coordinate (x,y,z)
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): dimensions by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values
%       ro is 3 by number of observation points
%       Efield is the Efield at points ro


%% Step 1 Extract BEM surface data structures
tic
[p,~,nf,~,f2te,f2p,~,~,~,~,~,~,~,~]=createmeshbem_cl(te2p,p,conductivity);
[t2p,p,epseff]=extractsurfaces(p,f2te,conductivity,nf,f2p);
  nc=numel(rs(1,:));
  nt=numel(t2p(:,1));
  np=numel(p(:,1));1
%% Step 2 Solve BEM for charges on interfaces 
tic
xval=bemsolver(t2p',nt,p',np,epseff(:),rs,js,nc);
toc

%% Step 3 generate obseravtions
Efield=computeEfields(t2p',nt,p',np,epseff(:),rs,js,nc,xval,ro,numel(ro(:))/3);

