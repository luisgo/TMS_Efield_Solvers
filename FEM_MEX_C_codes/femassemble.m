function A=femassemble(te2p,p,reg,FEMord)
%generates FEM matrix
nte=numel(te2p(1,:));
np=numel(p(1,:));
if FEMord==1
nte16=16*nte;
crvindex=zeros([3 nte16]);
val=zeros([1 nte16]);
mex_id_ = 'assemblemat(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindex, val] = FEM(mex_id_, nte, np, te2p, p, reg, crvindex, val, 4, nte, 3, np, nte, 3, nte16, nte16);
elseif FEMord==2

batsiz=5*10^5;
nte100=100*nte;
crvindex=zeros([3 nte100]);
val=zeros([1 nte100]);
if nte<batsiz
mex_id_ = 'assemblemat2nd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindex, val] = FEM(mex_id_, nte, np, te2p, p, reg, crvindex, val, 10, nte, 3, np, nte, 3, nte100, nte100);
else
nbat=floor(nte/batsiz);
ct=0;
st=1;
en=batsiz;
for ibatch=1:nbat
nbat100=100*batsiz;
crvindexbat=zeros([3 nbat100]);
valbat=zeros([1 nbat100]);
mex_id_ = 'assemblemat2nd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindexbat, valbat] = FEM(mex_id_, batsiz, np, te2p(:,st:en), p, reg(st:en), crvindexbat, valbat,...
 10, batsiz, 3, np, batsiz, 3, nbat100, nbat100);
crvindex(:,ct+1:ct+nbat100)=crvindexbat;
val(ct+1:ct+nbat100)=valbat(:);
ct=ct+nbat100;
st=st+batsiz;
en=en+batsiz;
end
en=nte;
batsiz=en-st+1;
nbat100=100*batsiz;
crvindexbat=zeros([3 nbat100]);
valbat=zeros([1 nbat100]);
mex_id_ = 'assemblemat2nd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindexbat, valbat] = FEM(mex_id_, batsiz, np, te2p(:,st:en), p, reg(st:en), crvindexbat, valbat,...
 10, batsiz, 3, np, batsiz, 3, nbat100, nbat100);
crvindex(:,ct+1:ct+nbat100)=crvindexbat;
val(ct+1:ct+nbat100)=valbat(:);
end
elseif FEMord==3
batsiz=5*10^5;
nte400=400*nte;
crvindex=zeros([3 nte400]);
val=zeros([1 nte400]);
if nte<batsiz
mex_id_ = 'assemblemat3rd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindex, val] = FEM(mex_id_, nte, np, te2p, p, reg, crvindex, val, 20, nte, 3, np, nte, 3, nte400, nte400);

else
nbat=floor(nte/batsiz);
ct=0;
st=1;
en=batsiz;
for ibatch=1:nbat
nbat400=400*batsiz;
crvindexbat=zeros([3 nbat400]);
valbat=zeros([1 nbat400]);
mex_id_ = 'assemblemat3rd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindexbat, valbat] = FEM(mex_id_, batsiz, np, te2p(:,st:en), p, reg(st:en), crvindexbat, valbat,...
 20, batsiz, 3, np, batsiz, 3, nbat400, nbat400);
crvindex(:,ct+1:ct+nbat400)=crvindexbat;
val(ct+1:ct+nbat400)=valbat(:);
ct=ct+nbat400;
st=st+batsiz;
en=en+batsiz;
end
en=nte;
batsiz=en-st+1;
nbat400=400*batsiz;
crvindexbat=zeros([3 nbat400]);
valbat=zeros([1 nbat400]);
mex_id_ = 'assemblemat3rd(i int, i int, i int[xx], i double[xx], i double[x], io int[xx], io double[x])';
[crvindexbat, valbat] = FEM(mex_id_, batsiz, np, te2p(:,st:en), p, reg(st:en), crvindexbat, valbat,...
 20, batsiz, 3, np, batsiz, 3, nbat400, nbat400);
crvindex(:,ct+1:ct+nbat400)=crvindexbat;
val(ct+1:ct+nbat400)=valbat(:);
end
end
A=sparse(crvindex(1,:)',crvindex(2,:)',val(:));

