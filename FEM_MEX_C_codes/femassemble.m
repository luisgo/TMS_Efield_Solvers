function A=femassemble(te2p,p,reg,FEMord)
%generates 1nd order FEM matrix
te2p=te2p-1;
nte=numel(te2p(1,:));
np=numel(p(1,:));
if FEMord==1
nte16=16*nte;
crvindex=zeros([3 nte16]);
val=zeros([1 nte16]);
mex_id_ = 'assemblemat(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindex, val] = FEM(mex_id_, nte, np, te2p, p, reg, crvindex, val, 4, nte, 3, np, nte, 3, nte16, nte16);
elseif FEMord==2
nte100=100*nte;
crvindex=zeros([3 nte100]);
val=zeros([1 nte100]);
mex_id_ = 'assemblemat2nd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindex, val] = FEM(mex_id_, nte, np, te2p, p, reg, crvindex, val, 10, nte, 3, np, nte, 3, nte100, nte100);
elseif FEMord==3
nte400=400*nte;
crvindex=zeros([3 nte400]);
val=zeros([1 nte400]);
mex_id_ = 'assemblemat3rd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindex, val] = FEM(mex_id_, nte, np, te2p, p, reg, crvindex, val, 20, nte, 3, np, nte, 3, nte400, nte400);
end
A=sparse(crvindex(1,:)',crvindex(2,:)',val(:));

