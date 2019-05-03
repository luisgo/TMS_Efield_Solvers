function A=femassembleslow(te2p2,p2,conductivity,FEMorder,nproc)
np=size(p2,2);
nte=size(te2p2,2);
fast_det2=@(pm,a,b) pm(1)*(a(2)*b(3)-a(3)*b(2))-pm(2)*(a(1)*b(3)-a(3)*b(1))+pm(3)*(a(1)*b(2)-a(2)*b(1));
fast_det=@(a,b) (a(2)*b(3)-a(3)*b(2))-(a(1)*b(3)-a(3)*b(1))+(a(1)*b(2)-a(2)*b(1));
aval=@(x,y,z) [ fast_det2([x(2),x(3),x(4)],[y(2),y(3),y(4)],[z(2),z(3),z(4)]),...
               -fast_det2([x(1),x(3),x(4)],[y(1),y(3),y(4)],[z(1),z(3),z(4)]),...
                fast_det2([x(1),x(2),x(4)],[y(1),y(2),y(4)],[z(1),z(2),z(4)]),...
               -fast_det2([x(1),x(2),x(3)],[y(1),y(2),y(3)],[z(1),z(2),z(3)])];
bval=@(x,y,z) [-fast_det([y(2),y(3),y(4)],[z(2),z(3),z(4)]),...
                fast_det([y(1),y(3),y(4)],[z(1),z(3),z(4)]),...
               -fast_det([y(1),y(2),y(4)],[z(1),z(2),z(4)]),...
                fast_det([y(1),y(2),y(3)],[z(1),z(2),z(3)])];
cval=@(x,y,z) [ fast_det([x(2),x(3),x(4)],[z(2),z(3),z(4)]),...
               -fast_det([x(1),x(3),x(4)],[z(1),z(3),z(4)]),...
                fast_det([x(1),x(2),x(4)],[z(1),z(2),z(4)]),...
               -fast_det([x(1),x(2),x(3)],[z(1),z(2),z(3)])];
dval=@(x,y,z) [-fast_det([x(2),x(3),x(4)],[y(2),y(3),y(4)]),...
                fast_det([x(1),x(3),x(4)],[y(1),y(3),y(4)]),...
               -fast_det([x(1),x(2),x(4)],[y(1),y(2),y(4)]),...
                fast_det([x(1),x(2),x(3)],[y(1),y(2),y(3)])];
%create element matrices;
a=zeros([4 nte]);
b=zeros([4 nte]);
c=zeros([4 nte]);
d=zeros([4 nte]);
vol=zeros([nte 1]);
te2p=te2p2(1:4,:);
for i=1:nte
a(:,i)=aval(p2(1,te2p(:,i)),p2(2,te2p(:,i)),p2(3,te2p(:,i)));
vol(i)=(a(1,i)+a(2,i)+a(3,i)+a(4,i))/6;
b(:,i)=bval(p2(1,te2p(:,i)),p2(2,te2p(:,i)),p2(3,te2p(:,i)));
c(:,i)=cval(p2(1,te2p(:,i)),p2(2,te2p(:,i)),p2(3,te2p(:,i)));
d(:,i)=dval(p2(1,te2p(:,i)),p2(2,te2p(:,i)),p2(3,te2p(:,i)));
end
clear te2p;
if FEMorder==1
%assemble matrix
col=ones([16*nte 1]);
row=ones([16*nte 1]);
val=zeros([16*nte 1]);
ct=1;
for i=1:nte
for j=1:4
    for k=1:4
        col(ct)=te2p2(j,i);
        row(ct)=te2p2(k,i);
        val(ct)=conductivity(i)*(b(j,i)*b(k,i)+c(j,i)*c(k,i)+d(j,i)*d(k,i))/(36*vol(i));
        ct=ct+1;
    end
end
end
A=sparse(col(:),row(:),val(:));
elseif FEMorder==2
load secondordint.mat matrixfunc; %precomputed integrals

col=ones([10,10,nte]);
row=ones([10,10,nte]);
val=zeros([10,10,nte]);
[XX,YY]=ndgrid(1:4,1:4);
matrixfunc2=matrixfunc;
if nproc~=1
mypool=parpool(nproc);
end
for i=1:nte
[col(:,:,i),row(:,:,i),val(:,:,i)]=...
    loopit(XX,YY,te2p2(:,i),b(:,i),c(:,i),d(:,i),...
    vol(i),conductivity(i),matrixfunc2,10);
end
if nproc~=1
delete(mypool)
end

A=sparse(col(:),row(:),val(:));
createintegral=toc
elseif FEMorder==3
    nbasis=20;
%% write matrix
np2=numel(p2(1,:));

load thirdordint.mat matrixfunc; %precomputed integrals
batchsiz=(1.5*10^6);
nbatch=floor(nte/batchsiz);
[XX,YY]=ndgrid(1:4,1:4);
matrixfunc2=matrixfunc;
for ibatch=1:nbatch+1
   if ibatch==nbatch+1
       ran=((ibatch-1)*batchsiz+1):nte;
   else
       ran=((ibatch-1)*batchsiz+1):(ibatch*batchsiz);
   end
ntemp=numel(ran);
col=ones([nbasis,nbasis,ntemp]);
row=ones([nbasis,nbasis,ntemp]);
val=zeros([nbasis,nbasis,ntemp]);
teran=te2p2(:,ran(:));
bran=b(:,ran(:));
cran=c(:,ran(:));
dran=d(:,ran(:));
volran=vol(ran(:));
condran=conductivity(ran(:));
if nproc~=1
mypool=parpool(nproc);
end
for i=1:ntemp
[col(:,:,i),row(:,:,i),val(:,:,i)]=...
    loopit(XX,YY,teran(:,i),bran(:,i),cran(:,i),dran(:,i),...
    volran(i),condran(i),matrixfunc2,nbasis);
end
if nproc~=1
delete(mypool)
end
if ibatch==1
A=sparse(col(:),row(:),val(:),np2,np2);
else
A=A+sparse(col(:),row(:),val(:),np2,np2);
end

end
end
end
function [col,row,val]=loopit(XX,YY,te2p2,b,c,d,vol,conductivity,matrixfunc,nbasis)
col=ones([nbasis nbasis]);
row=ones([nbasis nbasis]);
val=zeros([nbasis nbasis]);
dot=conductivity*(b(XX(:)).*b(YY(:))+...
     c(XX(:)).*c(YY(:))+...
     d(XX(:)).*d(YY(:)))/(6*vol);
for j=1:nbasis
    for k=1:nbasis
        col(k,j)=te2p2(j);
        row(k,j)=te2p2(k);
        val(k,j)=sum(matrixfunc(:,j,k).*dot);
    end
end

end

