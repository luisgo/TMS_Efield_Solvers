function [te2p2,p2]=femgenmesh_c(te2p,p,FEMord)
% generates correct order mesh from 1st order mesh
te2p=te2p-1;
nte=numel(te2p(1,:));
np=numel(p(1,:));
if FEMord==1
te2p2=te2p;p2=p;
elseif FEMord==2
mex_id_ = 'o int = generateedge_ct(i int, i int[xx])';
[ne] = FEM(mex_id_, nte, te2p, 4, nte);
npne=np+ne;
te2p2=zeros([10,nte]);p2=zeros([3,npne]);
mex_id_ = 'generatemesh2nd(i int, i int[xx], i double[xx], i int, io int[xx], io double[xx])';
[te2p2, p2] = FEM(mex_id_, nte, te2p, p, np, te2p2, p2, 4, nte, 3, np, 10, nte, 3, npne);
elseif FEMord==3
mex_id_ = 'o int = generateedge_ct(i int, i int[xx])';
[ne] = FEM(mex_id_, nte, te2p, 4, nte);
mex_id_ = 'o int = generateface_ct(i int, i int[xx])';
[nf] = FEM(mex_id_, nte, te2p, 4, nte);
npne=np+2*ne+nf;
te2p2=zeros([20,nte]);p2=zeros([3,npne]);
mex_id_ = 'generatemesh3rd(i int, i int[xx], i double[xx], i int, i int, io int[xx], io double[xx])';
[te2p2, p2] = FEM(mex_id_, nte, te2p, p, np, ne, te2p2, p2, 4, nte, 3, np, 20, nte, 3, npne);
end

