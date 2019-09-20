function [Efield,t2p,p,epseff,xval]=runcodeold(te2p,p,rs,js,conductivity,ro)
%zeroth order BEM code
%       FEMord is the order of FEM desired
%       tetrahedron to node id (te2p): dimensions 4 by number of tetrahedrons
%       node locations (p): dimensions 3 cartesian coordinate (x,y,z) by number of nodes
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): 3 dimensions by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values
%       ro is 3 by number of observation points
%       Efield is the Efield at points ro

%% Step 1 Extract BEM surface data structures
tic
[t2p,p,epseff]=genbemmesh(te2p,p,conductivity);
  nt=numel(t2p(1,:));
  np=numel(p(1,:));
  nc=numel(rs(1,:));
Generate_mesh_time=toc
%% Step 2 Generate BEM matrix data structures
tic
xval=bemsolver(t2p,nt,p,np,reg,rs,js,nc);
solve_time=toc
%% Step 3 generate obseravtions
tic
Efield=computeEfields(t2p,nt,p,np,epseff(:),rs,js,nc,xval,ro,numel(ro(:))/3);
Evaluate_field_time=toc


