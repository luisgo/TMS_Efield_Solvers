function [rhs]=femrhsgenslow(te2p,p,reg,rs2,js2,FEMord)
mu0=1.25663706*10^-6;
rs=rs2';js=-mu0*js2'/(4*pi);
p=p';te2p=te2p';
% generate surface triangles
if FEMord==1
nbasis=3;
Lmat=@(Lvals) cat(2,Lvals(1),Lvals(2),Lvals(3));
faces=[te2p(:,[2,3,4])
       te2p(:,[3,4,1])
       te2p(:,[4,1,2])
       te2p(:,[1,2,3])];
elseif FEMord==2
nbasis=6;
Lmat=@(Lvals) cat(2,...
Lvals(1)*(2*Lvals(1)-1),...
Lvals(2)*(2*Lvals(2)-1),...
Lvals(3)*(2*Lvals(3)-1),...
4*Lvals(1)*Lvals(2),...
4*Lvals(1)*Lvals(3),...
4*Lvals(2)*Lvals(3));
faces=[te2p(:,[2,3,4,8,9,10])
       te2p(:,[3,4,1,10,6,7])
       te2p(:,[4,1,2,7,9,5])
       te2p(:,[1,2,3,5,6,8])];

elseif FEMord==3
nbasis=10;
Lmat=@(Lvals) cat(2,...
Lvals(1)*(3*Lvals(1)-1)*(3*Lvals(1)-2)/2,...
Lvals(2)*(3*Lvals(2)-1)*(3*Lvals(2)-2)/2,...
Lvals(3)*(3*Lvals(3)-1)*(3*Lvals(3)-2)/2,...
9*Lvals(1)*Lvals(2)*(3*Lvals(1)-1)/2,...
9*Lvals(1)*Lvals(2)*(3*Lvals(2)-1)/2,...
9*Lvals(1)*Lvals(3)*(3*Lvals(1)-1)/2,...
9*Lvals(1)*Lvals(3)*(3*Lvals(3)-1)/2,...
9*Lvals(2)*Lvals(3)*(3*Lvals(2)-1)/2,...
9*Lvals(2)*Lvals(3)*(3*Lvals(3)-1)/2,...
27*Lvals(1)*Lvals(2)*Lvals(3)...
);

faces=[te2p(:,[2,3,4,11,12,13,14,15,16,20])
       te2p(:,[3,4,1,15,16,8,7,10,9,19])
       te2p(:,[4,1,2,10,9,14,13,5,6,18])
       te2p(:,[1,2,3,5,6,7,8,11,12,17])];
end
% for j=1:numel(faces(:,1))
% scatter3(p(faces(j,:),1),p(faces(j,:),2),p(faces(j,:),3));
% for i=1:nbasis
%     text(p(faces(j,i),1),p(faces(j,i),2),p(faces(j,i),3),num2str(i));
% end
% pause
% close all
% end
faces2=sort(faces,2);
[foo,ff,faces2]=unique(faces2,'rows');
f2p=faces(ff,:);
nf=max(faces2);
nte=length(te2p(:,1));
f2te=zeros([nf 2]);
te2f=zeros([nte 4]);
for i1=1:nte
    for i2=1:4
    te2f(i1,i2)=faces2(i1+(i2-1)*nte);
    if f2te(faces2(i1+(i2-1)*nte),1)==0
        f2te(faces2(i1+(i2-1)*nte),1)=i1;
    else
        f2te(faces2(i1+(i2-1)*nte),2)=i1;
    end
    end
end

nhat=zeros([nf 3]);
for i=1:nf
    v1=p(f2p(i,1),:)-p(f2p(i,3),:);
    v2=p(f2p(i,2),:)-p(f2p(i,3),:);
nhatt=[v1(2)*v2(3)-v1(3)*v2(2) ...
      v1(3)*v2(1)-v1(1)*v2(3) ...
      v1(1)*v2(2)-v1(2)*v2(1)];
  if te2f(f2te(i,1),1)==i
      freenode=te2p(f2te(i,1),1);
  elseif te2f(f2te(i,1),2)==i
      freenode=te2p(f2te(i,1),2);
  elseif te2f(f2te(i,1),3)==i
      freenode=te2p(f2te(i,1),3);
  elseif te2f(f2te(i,1),4)==i
      freenode=te2p(f2te(i,1),4);
  end
  outn=p(f2p(i,1),:)-p(freenode,:); %on surface of triangle so pointing in rhat direction

  nhat(i,:)=nhatt*sign(sum(nhatt(:).*outn(:)));%nhat points toward tetra
end
interface=1:nf;
del=zeros([nf 1]);
del(f2te(:,2)==0)=reg(f2te(f2te(:,2)==0,1));
del(f2te(:,2)~=0)=reg(f2te(f2te(:,2)~=0,1))...
    -reg(f2te(f2te(:,2)~=0,2));
interface=interface(del~=0);
ninterface=numel(interface);
mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
%%the six comes from basis function
rhs=zeros([max(te2p(:)) 1]);


nquad=1;
qwt=1;
qpt=[1/3 1/3 1/3];
[qwt,qpt,nquad]=integrationrules2d(FEMord+1);


ro=zeros([ninterface,nquad,3]);
for i=1:nquad
   for j=1:ninterface
       ro(j,i,:)=reshape(qpt(i,1)*p(f2p(interface(j),1),:)+...
                         qpt(i,2)*p(f2p(interface(j),2),:)+...
                         qpt(i,3)*p(f2p(interface(j),3),:),[1 1 3]);
   end
end

ro=reshape(ro,[ninterface*nquad 3]);
Ef2=zeros([ninterface nquad 3]);
nobs=numel(ro(:,1));
nsource=size(rs);
acc=0;
[Ux]=lfmm3dpart(acc,nsource(1),rs',1,js(:,1),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,0);
[Uy]=lfmm3dpart(acc,nsource(1),rs',1,js(:,2),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,0);
[Uz]=lfmm3dpart(acc,nsource(1),rs',1,js(:,3),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,0);
Ef2(:,:,1)=reshape(real(Ux.pottarg(:)),[ninterface nquad]);
Ef2(:,:,2)=reshape(real(Uy.pottarg(:)),[ninterface nquad]);
Ef2(:,:,3)=reshape(real(Uz.pottarg(:)),[ninterface nquad]);
ninterface
rhs=zeros([numel(p(:,1)) 1]);
for j=1:nquad
          ba=Lmat([qpt(j,1),qpt(j,2),qpt(j,3)]);
for i=1:ninterface
dadtnhat=del(interface(i))*...
            (Ef2(i,j,1)*nhat(interface(i),1)+...
             Ef2(i,j,2)*nhat(interface(i),2)+...
             Ef2(i,j,3)*nhat(interface(i),3));
    for k=1:nbasis
    rhs(f2p(interface(i),k))=rhs(f2p(interface(i),k))...
        + ba(k)*dadtnhat*qwt(j)/2;
    end
end
end
end
