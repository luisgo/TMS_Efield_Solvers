function [Eout]=computeEprimary(rs,js,nsource,robs,ntarget)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
Eout=zeros([3 ntarget]);
mex_id_ = 'Eprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx])';
[Eout] = FEM(mex_id_, nsource, ntarget, rs, js, robs, Eout, 3, nsource, 3, nsource, 3, ntarget, 3, ntarget);





