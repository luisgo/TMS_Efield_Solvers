function [Eout,RELRES,ITER,RESVEC,conductivity,Einc,rce,xval,ncell,dcell]=runfdlag(te2p,p,cond,rs,js,dcell)
%%%%%%%%%%%%%%%% generate mesh
tic
p=p';te2p=te2p';
[XX,YY,ZZ]=ndgrid(min(p(:,1)):dcell(1):max(p(:,1)),...
                  min(p(:,2)):dcell(2):max(p(:,2)),...
                  min(p(:,3)):dcell(3):max(p(:,3)));% cell node locations
TR=triangulation(double(te2p),p);
nproc=1;%number of processors can be increased for faster speed here I default to one processor
if nproc>1
nst=floor(numel(XX)/(nproc-1));
for i=1:nproc-1
XXi{i}=XX((i-1)*nst+1:i*nst)+dcell(1)/2;
YYi{i}=YY((i-1)*nst+1:i*nst)+dcell(2)/2;
ZZi{i}=ZZ((i-1)*nst+1:i*nst)+dcell(3)/2;
end
i=nproc;
XXi{i}=XX((i-1)*nst+1:end)+dcell(1)/2;
YYi{i}=YY((i-1)*nst+1:end)+dcell(2)/2;
ZZi{i}=ZZ((i-1)*nst+1:end)+dcell(3)/2;
mypool=parpool(nproc);
parfor i=1:nproc
teidi{i}=pointLocation(TR,XXi{i}(:),YYi{i}(:),ZZi{i}(:));
end
delete(mypool);
teid=teidi{1}(:);
for i=2:nproc
teid=cat(1,teid(:),teidi{i}(:));
end
clear teidi;
else
teid=pointLocation(TR,XX(:)+dcell(1)/2,YY(:)+dcell(2)/2,ZZ(:)+dcell(3)/2);%sample conductivity at center of each cell
end
conductivity=zeros(size(teid));
conductivity(isnan(teid)==0)=cond(teid(isnan(teid)==0));
conductivity=reshape(conductivity,size(XX));
conductivity=conductivity(1:end-1,1:end-1,1:end-1);

ncell=size(conductivity);
rce(:,1)=XX(:);rce(:,2)=YY(:);rce(:,3)=ZZ(:);
clear TR XX YY ZZ;
genmesh=toc
tic
Einc=computeEprimary(rs,js,numel(rs)/3,rce',numel(rce)/3).';
size(Einc)
ncell
Einc=permute(reshape(Einc,[ncell(1)+1,ncell(2)+1,ncell(3)+1,3]),[4,1,2,3]);
genEinc=toc
tic
[A,rhs,sink,precon]=fdcodelag(conductivity,ncell,dcell,Einc);

%%%%%%%%%%%%%solve system
A(sink,:)=0;% delete one unknown, i.e. defined a ground V=0
A(:,sink)=0;
rhs(sink)=0; 
precon(sink)=0;
precon(precon~=0)=(1./precon(precon~=0));
%A(precon~=0,precon~=0)=diag(sparse(precon(precon~=0)))*A(precon~=0,precon~=0)*diag(sparse(precon(precon~=0))); %diagonally precondition to make sure that the system is properly scaled
%rhs=rhs.*precon;
genA=toc
tic

 [xval,FLAG,RELRES,ITER,RESVEC] =minres(A,rhs,10^-7,50000,@(x)matvec(x,precon)); %a tolerance of up to 10^-4 would probably be okay if it doesnt converge
whos
%xval=xval.*precon;
solvesys=toc
%%%%%%%%%%%%%%%%%%
tic
xval=reshape(xval,ncell+1);
% rv(:,3)=rce(:,3)+dcell(3)/2;
% rv(:,1)=rce(:,1)+dcell(1)/2;
% rv(:,2)=rce(:,2)+dcell(2)/2;
% Eout=fdinterpolatorlag(Einc,rce,xval,rv,ncell,dcell);
% Eout=reshape(Eout,[ncell(1)+1,ncell(2)+1,ncell(3)+1,3]);
% for i=1:3
% Eout(1:end-1,1:end-1,1:end-1,i)=1i*Eout(1:end-1,1:end-1,1:end-1,i).*(conductivity~=0);
% end
% Eout=Eout(1:end-1,1:end-1,1:end-1,:);
Eout=0;
rv=0;
postproc=toc
end
function y=matvec(x,precon)
y=x.*precon;
end
