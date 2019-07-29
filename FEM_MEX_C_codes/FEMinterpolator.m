function [Eout]=FEMinterpolator(te2te,te2p,p,rs,js,xval,robs,FEMord)
%compute Etotal at observation points robs 
nobs=numel(robs)/3;
np=numel(p)/3;
nsource=numel(rs)/3;
Eout=zeros([3 nobs]);
iprecEp=2;%FMM accuracy flag lower is faster but less accurate
if FEMord==1
nte=numel(te2p)/4;
mex_id_ = 'totalEfield(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[xx], i double[xx], i int, i double[x], io double[xx], i int)';
[Eout] = FEM(mex_id_, robs, nobs, te2te, te2p, nte, p, rs, js, nsource, xval, Eout, iprecEp, 3, nobs, 4, nte, 4, nte, 3, np, 3, nsource, 3, nsource, np, 3, nobs);
elseif FEMord==2
nte=numel(te2p)/10;
mex_id_ = 'totalEfield2nd(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[xx], i double[xx], i int, i double[x], io double[xx], i int)';
[Eout] = FEM(mex_id_, robs, nobs, te2te, te2p, nte, p, rs, js, nsource, xval, Eout, iprecEp, 3, nobs, 4, nte, 10, nte, 3, np, 3, nsource, 3, nsource, np, 3, nobs);
elseif FEMord==3
nte=numel(te2p)/20;
mex_id_ = 'totalEfield3rd(i double[xx], i int, i int[xx], i int[xx], i int, i double[xx], i double[xx], i double[xx], i int, i double[x], io double[xx], i int)';
[Eout] = FEM(mex_id_, robs, nobs, te2te, te2p, nte, p, rs, js, nsource, xval, Eout, iprecEp, 3, nobs, 4, nte, 20, nte, 3, np, 3, nsource, 3, nsource, np, 3, nobs);
end

