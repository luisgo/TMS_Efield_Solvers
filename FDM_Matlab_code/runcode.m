function [Efield,rce,xval,ncell]=runcode(te2p,p,rs,js,conductivity,ro,dcell)
% The required data structures loaded from file are the following:
%       FEMord is the order of FEM desired
%       tetrahedron to node id (te2p): dimensions 4 nodes by number of tetrahedrons
%       node locations (p): dimensions 3 cartesian coordinate (x,y,z) by
%       number of nodes
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): dimensions 3 vector weight (Jx,Jy,Jz) by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values
%       ro is 3 by number of observation points
%       Efield is the Efield at points ro
%% Step 1 run FDM
tic
[Efield,RELRES,ITER,RESVEC,condu,Einc,rce,xval,ncell,dcell]=runfdlag(te2p,p,conductivity,rs,js,dcell);
toc
%% Step 2 evaluate field at desired locations
 Efield=fdinterpolatorlag(rce,xval,ro',ncell,dcell,rs,js);
