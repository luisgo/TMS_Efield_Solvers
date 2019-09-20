function [t2p,p2,epseff]=genbemmesh(te2p,p,condu)
%extracts tissue boundary faces from the tetrahedron mesh
%te2p is the 4 by number of tetrahedrons element array
%p is the 3 by number of points node location array
%condu is the number of tetrahedron by one conductivity array
%t2p is 3 by number of boundary triangles
%p2 is the boundary triangle point array
%epseff is the (conductivity_inner-conductivity_outer)/(conductivity_inner+conductivity_outer) array
te2p=te2p-1;
nte=numel(te2p)/4;
np=numel(p)/3;np2=0;
mex_id_ = 'o int = countfaces(i int, i int[xx], i double[x], io int[x])';
[nf, np2] = BEM(mex_id_, nte, te2p, condu, np2, 4, nte, nte, 1);
t2p=zeros([3,nf]);
p2=zeros([3,np2]);
epseff=zeros([nf,1]);
mex_id_ = 'makeBEMmesh(i int, i int[xx], i double[xx], i double[x], io int[xx], io double[xx], io double[x])';
[t2p, p2, epseff] = BEM(mex_id_, nte, te2p, p, condu, t2p, p2, epseff, 4, nte, 3, np, nte, 3, nf, 3, np2, nf);
t2p=t2p+1;

