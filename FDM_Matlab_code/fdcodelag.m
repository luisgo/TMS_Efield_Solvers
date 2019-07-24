function [A,rhs,tempind,precon]=fdcodelag(conductivity,n,dcell,vectpot)

basisfn=@(x,y,z)[(1-x)*(1-y)*(1-z);...
                 (x  )*(1-y)*(1-z);...
                 (1-x)*(  y)*(1-z);...
                 (x  )*(  y)*(1-z);...
                 (1-x)*(1-y)*(  z);...
                 (x  )*(1-y)*(  z);...
                 (1-x)*(  y)*(  z);...
                 (x  )*(  y)*(  z)];
gradfn=@(x,y,z)[ -(1-y)*(1-z),-(1-x)*(1-z),-(1-x)*(1-y);...
                  (1-y)*(1-z),-(x  )*(1-z),-(x  )*(1-y);...
                 -(  y)*(1-z), (1-x)*(1-z),-(1-x)*(  y);...
                  (  y)*(1-z), (x  )*(1-z),-(x  )*(  y);...
                 -(1-y)*(  z),-(1-x)*(  z), (1-x)*(1-y);...
                  (1-y)*(  z),-(x  )*(  z), (x  )*(1-y);...
                 -(  y)*(  z), (1-x)*(  z), (1-x)*(  y);...
                  (  y)*(  z), (x  )*(  z), (x  )*(  y)];
nquad=8;
precon=zeros([prod(n+1),1]);
qpt([1,3,5,7],1)=1/2-1/sqrt(3)/2;
qpt([2,4,6,8],1)=1/2+1/sqrt(3)/2;
qpt([1,2,5,6],2)=1/2-1/sqrt(3)/2;
qpt([3,4,7,8],2)=1/2+1/sqrt(3)/2;
qpt([1,2,3,4],3)=1/2-1/sqrt(3)/2;
qpt([5,6,7,8],3)=1/2+1/sqrt(3)/2;
qwt(1:8)=prod(dcell)/8;
lookup=zeros([8,8]);
inter=zeros([nquad,8]);
gr=zeros([8,3,nquad]);
for iquad=1:nquad
    gr(:,:,iquad)=gradfn(qpt(iquad,1),qpt(iquad,2),qpt(iquad,3))*qwt(iquad);
    inter(:,iquad)=basisfn(qpt(iquad,1),qpt(iquad,2),qpt(iquad,3));
for i=1:8
    for j=1:8
    lookup(i,j)=lookup(i,j)+sum(gr(i,:,iquad).*gr(j,:,iquad)./(dcell(:)').^2)/qwt(iquad); 
    end
end
end
ival=1:8;
[cc,rr]=ndgrid(ival,ival);
ti=mod(ival-1,2);
tj=mod(floor(ival/2.2),2);
tk=mod(floor(ival/5),2);
ishift=ti+tj*(n(1)+1)+tk*(n(1)+1)*(n(2)+1);
globind=@(i,j,k) i+(j-1)*(n(1)+1)+(k-1)*(n(1)+1)*(n(2)+1);
ishift2=ival(:)+(ival(:)-1)*8;
col=zeros([nnz(conductivity)*64,1]);
row=zeros([nnz(conductivity)*64,1]);
val=zeros([nnz(conductivity)*64,1]);
ct=1;
for k=1:n(3)
for j=1:n(2)
for i=1:n(1)
    if conductivity(i,j,k)~=0
    tempind=globind(i,j,k);
col(ct:ct+63)=tempind+ishift(cc(:));
row(ct:ct+63)=tempind+ishift(rr(:));
val(ct:ct+63)=conductivity(i,j,k)*lookup(:);
   ct=ct+64;
precon(tempind+ishift(:))=precon(tempind+ishift(:))+conductivity(i,j,k)*lookup(ishift2);
    end
end
end
end
A=sparse(col,row,val,prod(n+1),prod(n+1));
clear col row val;

rhs=zeros([prod(n+1),1]);
cc=1:8;
Einc(1:3)=0;
for k=1:n(3)
for j=1:n(2)
for i=1:n(1)
    if conductivity(i,j,k)~=0
for iquad=1:nquad
    for e=1:3
Einc(e)=inter(1,iquad)*vectpot(e,i,j,k)+inter(2,iquad)*vectpot(e,i+1,j,k)+...
        inter(3,iquad)*vectpot(e,i,j+1,k)+inter(4,iquad)*vectpot(e,i+1,j+1,k)+...
        inter(5,iquad)*vectpot(e,i,j,k+1)+inter(6,iquad)*vectpot(e,i+1,j,k+1)+...
        inter(7,iquad)*vectpot(e,i,j+1,k+1)+inter(8,iquad)*vectpot(e,i+1,j+1,k+1);
    end
    tempind=globind(i,j,k);
rhs(tempind+ishift(cc(:)))=rhs(tempind+ishift(cc(:)))+conductivity(i,j,k)*...
                          (Einc(1)*gr(cc(:),1,iquad)/dcell(1)+...
                           Einc(2)*gr(cc(:),2,iquad)/dcell(2)+...
                           Einc(3)*gr(cc(:),3,iquad)/dcell(3));
end
     end
end
end
end