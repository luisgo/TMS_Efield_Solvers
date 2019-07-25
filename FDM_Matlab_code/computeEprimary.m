function [Eout]=computeEprimary(rs,js,nsource,robs,ntarget)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
Eout=zeros([3 ntarget]);
iprecEp=2;%FMM accuracy flag lower is faster but less accurate
mex_id_ = 'Eprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], i int)';
[Eout] = FEM(mex_id_, nsource, ntarget, rs, js, robs, Eout, iprecEp, 3, nsource, 3, nsource, 3, ntarget, 3, ntarget);





