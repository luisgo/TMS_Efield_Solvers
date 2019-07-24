function [t2p,p,epseff]=extractsurfaces(p,f2te,sigma,nf,f2p);
 epseff=sigma(f2te(:,1));
  epseff(f2te(:,2)~=0)=...
     (epseff(f2te(:,2)~=0)-sigma(f2te(f2te(:,2)~=0,2)))...
   ./(epseff(f2te(:,2)~=0)+sigma(f2te(f2te(:,2)~=0,2)));
 
epseff(f2te(:,2)==0)=1;
faces=1:nf;
faces=faces(epseff~=0);
epseff=epseff(epseff~=0);
nt=numel(faces);
t2p=f2p(faces,:);
[p2,pindnew,t2p]=unique(t2p(:));
t2p=reshape(t2p,[nt 3]);
p=p(p2,:);
end