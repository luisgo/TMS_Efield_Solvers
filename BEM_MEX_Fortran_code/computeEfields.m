function Etotal=computeEfields(t2p,nt,p,np,reg,rs,js,nc,xval,robs,nobs)
%nt number of triangles
%np number of points
%nc number of coil dipoles
%nobs number of observation points
%t2p triangle to points (3 by nt)
%p cartesian point coordinates (3 by np)
%epseff contrast function (eps2-eps1)/(eps2+eps1) (nt by 1)
%rs cartesian locations of electrical dipoles (3 by nc)
%js electrical dipole weight and direction (3 by nc)
%rs cartesian locations of electrical dipoles (3 by nc)
%js electrical dipole weight and direction (3 by nc)
%xval charge weights from BEM solution  (i.e. pulse basis expansion coefficients)
%robs locations for E-field evaluation
Etotal=zeros([3,nobs]);

%%%%%%%%%%%%modify the following at your own risk%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%See supplemental material for explanation of these parameters
%nquad(1) is the number of quadrature points on the far triangle interactions nquad(2) number of quadrature points on near interactions valid inputs are (1,3,16). Suggested is 1 point and 16 for near
nquad=[3,16];
%iprec(1) FMM precision flag for computing E-primary iprec(2) and FMM precision flag for matrix vector multiply. Values of iprec 0,1,2,3,4,and 5 typically provide errors of FMM accuracy of 2,3,6,9,12,and 14 digits, respectivelly. Recommended value iprec=[1,0] and [2,0] for error lower than 0.1%.
iprec=[1,0];
%dnear is number of average edge lengths of the near-field interactions. Our tests indicate that this is accurate enough for FEM simulations.
dnear=1;

% Note: treated as contiguous 1D array as far as C/FORTRAN is concerned:
mex_id_ = 'evaluateEfield(i int[xx], i int[x], i double[xx], i int[x], i double[x], i double[xx], i double[xx], i int[x], i double[x], i double[xx], i int[x], io double[xx], i int[x], i int[x], i double[x])';
[Etotal] = BEM(mex_id_, t2p, nt, p, np, reg, rs, js, nc, xval, robs, nobs, Etotal, nquad, iprec, dnear, 3, nt, 1, 3, np, 1, nt, 3, nc, 3, nc, 1, nt, 3, nobs, 1, 3, nobs, 2, 2, 1);
