clear all
omega=1;
reflev=0;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
[E1,t2p,p,epseff,xval]=runcode(te2p,p,rs,omega*js,conductivity,ro1);%generates Efield at ro1

np=numel(p)/3;nt=numel(t2p)/3;nc=numel(rs)/3;
E2=computeEfields(t2p',nt,p',np,epseff(:),rs,omega*js,nc,xval,ro2,numel(ro2(:))/3);
E3=computeEfields(t2p',nt,p',np,epseff(:),rs,omega*js,nc,xval,ro3,numel(ro3(:))/3);
E4=computeEfields(t2p',nt,p',np,epseff(:),rs,omega*js,nc,xval,ro4,numel(ro4(:))/3);

%% Compare errors
Esoln2(isnan(Esoln2))=0;
norm(Esoln1(:)-E1(:))/norm(Esoln1(:))
norm(Esoln2(:)-E2(:))/norm(Esoln2(:))
norm(Esoln3(:)-E3(:))/norm(Esoln3(:))
norm(Esoln4(:)-E4(:))/norm(Esoln4(:))