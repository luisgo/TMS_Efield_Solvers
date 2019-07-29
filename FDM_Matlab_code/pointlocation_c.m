function [teid,bari]=pointlocation_c(robs,te2p,p,te2te,te_stride)
%compute Eobs at observation points (mostly for testing C++ to Fortran interface)
nobs=numel(robs)/3;
np=numel(p)/3;
nte=numel(te2p)/te_stride;
bari=zeros([4,nobs]);
teid=zeros([nobs,1]);
mex_id_ = 'pointLocation_f(i double[xx], i int, i int[xx], i int, i double[xx], i int, i int[xx], io int[x], io double[xx])';
[teid, bari] = FEM(mex_id_, robs, nobs, te2p, nte, p, te_stride, te2te, teid, bari, 3, nobs, te_stride, nte, 3, np, 4, nte, nobs, 4, nobs);
teid=teid+1;%change indexing from C style to matlab style

