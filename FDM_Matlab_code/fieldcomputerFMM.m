function Ef=fieldcomputerFMM(ro,rv,jv,omega)
mu0=1.25663706*10^-6;
Ef=zeros(size(ro));
nobs=numel(ro(:,1));
nsource=size(rv);
acc=3;
[U]=lfmm3dpart(acc,nsource(1),rv',1,jv(:,3),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,0);
Ef(:,3)=-omega*mu0*U.pottarg(:)/(4*pi);
[U]=lfmm3dpart(acc,nsource(1),rv',1,jv(:,2),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,0);
Ef(:,2)=-omega*mu0*U.pottarg(:)/(4*pi);
[U]=lfmm3dpart(acc,nsource(1),rv',1,jv(:,1),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,0);
Ef(:,1)=-omega*mu0*U.pottarg(:)/(4*pi);
Ef=1i*real(Ef);
end