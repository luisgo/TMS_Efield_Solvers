function Efield=computefields(rs,js,ks,rhos,omega,lmax,rsphere,pt)
%rs dipole locations
%js electric dipole vector weights
%ks magnetic dipole vector weights
%rhos charge weights
%omega dA/dt scalar
%lmax is the number of spherical harmonics (20-30 is enough) i put 40 to
%get 12 digits
%rsphere radius of expansion currents rsphere has to separate pt from rs
%(i.e. rs must be outside a sphere with radius rsphere and all pt must be
%inside sphere with radius rsphere).
%pt points within rpshere where E-field has to be evaluated
%%%%%%%%%%this code is slow because the parameter orde is 200 (i.e. 400 order quadrature rule) 
%%%%%%%%%%, which is overkill. It could be set to 100.
mu0=1.2566370614*10^-6;
coeff=computemultipoles4(rs,js,ks,rhos,omega,lmax,rsphere);
[Emesh]=Eevaluator2(pt(:,1),pt(:,2),pt(:,3),rsphere,lmax);
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
stack=globind(lmax,lmax);
Efield=zeros(size(Emesh{1}));
for l=1:lmax
    for m=1:1:l;
        Efield=Efield+Emesh{globind(m,l)+stack}*coeff(globind(m,l)+stack);
    end
    for m=0:1:l;
        Efield=Efield+Emesh{globind(m,l)}*coeff(globind(m,l)); 
    end
end
Efield=-reshape(Efield,[numel(Efield)/3 3]);
end


function [Efield,Srpos,Srneg]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
stack=globind(lmax,lmax);
%change to spherical
nx=size(Xeval);
if numel(nx)<3;
nx(end+1:3)=1;
end
Reval=sqrt(Xeval(:).^2+Yeval(:).^2+Zeval(:).^2);
rad=Reval;
phi = atan2(Yeval(:),Xeval(:));
theta = acos(Zeval(:)./sqrt(Xeval(:).^2+Yeval(:).^2+Zeval(:).^2));

[Ypos,Yneg,Srpos,Srneg]=Yfunc2(lmax,theta,phi,rad);
Npt=numel(theta);

for l=0:lmax
    cons=-1/(2*l+1)/bigR^l;
    radrr(1,:,1)=cons*reshape(rad(:).^(l),[1,numel(rad),1]);
    Br=radrr.*Ypos{2}(globind(0,l),:,1);Bth=radrr.*Ypos{2}(globind(0,l),:,2);Bph=radrr.*Ypos{2}(globind(0,l),:,3);
        Efield{globind(0,l)}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(0,l)}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}=zeros([nx(1) nx(2) nx(3) 3]);
    
        Efield{globind(0,l)+stack}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);
    for m=1:l
    Br=radrr.*Ypos{2}(globind(m,l),:,1);Bth=radrr.*Ypos{2}(globind(m,l),:,2);Bph=radrr.*Ypos{2}(globind(m,l),:,3);
        Efield{globind(m,l)}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(m,l)}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);        
    Br=radrr.*Yneg{2}(globind(m,l),:,1);Bth=radrr.*Yneg{2}(globind(m,l),:,2);Bph=radrr.*Yneg{2}(globind(m,l),:,3);
        Efield{globind(m,l)+stack}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(m,l)+stack}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)+stack}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)+stack}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);        
    end
end
end




function coeff=computemultipoles4(rs,js,ks,rhos,omega,lmax,rsphere)
mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
%load quadrature
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
orde=200;
[qth,qwth]=lgwt(orde,0,pi);
[qph,qwph]=lgwt(orde,0,2*pi);
[theta,phi]=ndgrid(qth,qph);
[qwth,qwph]=ndgrid(qwth,qwph);
w=rsphere^2.*sin(theta(:)).*qwth(:).*qwph(:);
x=rsphere.*sin(theta(:)).*cos(phi(:));
y=rsphere.*sin(theta(:)).*sin(phi(:));
z=rsphere.*cos(theta(:));

ro(:,1)=x(:);
ro(:,2)=y(:);
ro(:,3)=z(:);
Ef=EfieldcomputerFMM(ro,rs,js,rhos,ks,omega,5);
Ef(:,1)=w(:).*Ef(:,1);
Ef(:,2)=w(:).*Ef(:,2);
Ef(:,3)=w(:).*Ef(:,3);

[Emat,rhat]=Eevaluator3(ro(:,1),ro(:,2),ro(:,3),lmax);
stack=globind(lmax,lmax);
coeff=zeros([2*stack 1]);
for l=1:lmax
    for m=1:1:l;
        ind=globind(m,l)+stack;
coeff(ind)=sum(sum(Emat{ind}(:).*Ef(:)));
    end
    for m=0:1:l;
        ind=globind(m,l);
coeff(ind)=sum(sum(Emat{ind}(:).*Ef(:)));
    end
end
coeff=coeff/rsphere^2;
end
function Ef=EfieldcomputerFMM(ro,rs2,js2,rhos2,ks2,omega,acc)
mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;

Ef=zeros(size(ro));
nobs=numel(ro(:,1));
js=-omega*mu0*js2/(4*pi);
ks=ks2/(4*pi);
rhos=rhos2/(4*pi*eps0);
rs=rs2;
nsource=size(rs);
[Ux]=lfmm3dpart(acc,nsource(1),rs',1,js(:,1),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
[Uy]=lfmm3dpart(acc,nsource(1),rs',1,js(:,2),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
[Uz]=lfmm3dpart(acc,nsource(1),rs',1,js(:,3),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
[Vx]=lfmm3dpart(acc,nsource(1),rs',1,ks(:,1),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
[Vy]=lfmm3dpart(acc,nsource(1),rs',1,ks(:,2),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
[Vz]=lfmm3dpart(acc,nsource(1),rs',1,ks(:,3),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
[Z]=lfmm3dpart(acc,nsource(1),rs',1,rhos(:),...
    0,zeros([1 nsource(1)]),zeros([3 nsource(1)]),0,0,nobs,ro',1,1);
Ef(:,3)=Uz.pottarg(1:nobs).'-(Vx.fldtarg(2,1:nobs)-Vy.fldtarg(1,1:nobs)).'+Z.fldtarg(3,1:nobs).';
Ef(:,2)=Uy.pottarg(1:nobs).'-(Vz.fldtarg(1,1:nobs)-Vx.fldtarg(3,1:nobs)).'+Z.fldtarg(2,1:nobs).';
Ef(:,1)=Ux.pottarg(1:nobs).'-(Vy.fldtarg(3,1:nobs)-Vz.fldtarg(2,1:nobs)).'+Z.fldtarg(1,1:nobs).';
end
function Ef=fieldcomputer(ro,rs,js,rhos,ks,omega,mu0,eps0)

nsour=size(rs);
Ef=zeros(size(ro));
for is=1:nsour(1)
Green(:,1)=ro(:,1)-rs(is,1);
Green(:,2)=ro(:,2)-rs(is,2);
Green(:,3)=ro(:,3)-rs(is,3);
bigR=sqrt(sum(Green.^2,2));
L_s=-1i*omega*mu0*1./bigR;
Green(:,1)=-Green(:,1)./bigR.^3;
Green(:,2)=-Green(:,2)./bigR.^3;
Green(:,3)=-Green(:,3)./bigR.^3;
Ef(:,1)=Ef(:,1)+js(is,1)*L_s-Green(:,1)*rhos(is)/eps0;
Ef(:,2)=Ef(:,2)+js(is,2)*L_s-Green(:,2)*rhos(is)/eps0;
Ef(:,3)=Ef(:,3)+js(is,3)*L_s-Green(:,3)*rhos(is)/eps0;
Ef(:,1)=Ef(:,1)-(Green(:,2)*ks(is,3)-Green(:,3)*ks(is,2));
Ef(:,2)=Ef(:,2)-(Green(:,3)*ks(is,1)-Green(:,1)*ks(is,3));
Ef(:,3)=Ef(:,3)-(Green(:,1)*ks(is,2)-Green(:,2)*ks(is,1));
end
Ef=Ef/(4*pi);
end
function [Efield,rhat,Srpos,Srneg]=Eevaluator3(Xeval,Yeval,Zeval,lmax)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
stack=globind(lmax,lmax);
%change to spherical
nx=size(Xeval);
if numel(nx)<3;
nx(end+1:3)=1;
end
Reval=sqrt(Xeval(:).^2+Yeval(:).^2+Zeval(:).^2);
rad=Reval;
phi = atan2(Yeval(:),Xeval(:));
theta = acos(Zeval(:)./sqrt(Xeval(:).^2+Yeval(:).^2+Zeval(:).^2));
rhat(:,1)=cos(phi).*sin(theta);
rhat(:,2)=sin(phi).*sin(theta);
rhat(:,3)=cos(theta);
[Ypos,Yneg,Srpos,Srneg]=Yfunc2(lmax,theta,phi,rad);
Npt=numel(theta);

Efield{1}=zeros([globind(lmax,lmax),Npt 3]);
for l=0:lmax
    cons=1/(2*l+1);
    radrr(1,:,1)=1/cons;
    %forces the fiels to match optimization
    Br=radrr.*Ypos{2}(globind(0,l),:,1);Bth=radrr.*Ypos{2}(globind(0,l),:,2);Bph=radrr.*Ypos{2}(globind(0,l),:,3);
        Efield{globind(0,l)}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(0,l)}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}=zeros([nx(1) nx(2) nx(3) 3]);
    
        Efield{globind(0,l)+stack}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);
    for m=1:l
    Br=radrr.*Ypos{2}(globind(m,l),:,1);Bth=radrr.*Ypos{2}(globind(m,l),:,2);Bph=radrr.*Ypos{2}(globind(m,l),:,3);
        Efield{globind(m,l)}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(m,l)}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);        
    Br=radrr.*Yneg{2}(globind(m,l),:,1);Bth=radrr.*Yneg{2}(globind(m,l),:,2);Bph=radrr.*Yneg{2}(globind(m,l),:,3);
        Efield{globind(m,l)+stack}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(m,l)+stack}(:,:,:,1)=reshape(Br(:).*sin(theta).*cos(phi)+...
            Bth(:).*cos(theta).*cos(phi)-Bph(:).*sin(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)+stack}(:,:,:,2)=reshape(Br(:).*sin(theta).*sin(phi)+...
            Bth(:).*cos(theta).*sin(phi)+Bph(:).*cos(phi),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)+stack}(:,:,:,3)=reshape(Br(:).*cos(theta)-...
            Bth(:).*sin(theta),[nx(1) nx(2) nx(3) 1]);        
    end
end
end
function [Ypos,Yneg,Srpos,Srneg]=Yfunc2(lmax,theta,phi,rad)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
Npt=length(theta(:));
theta=reshape(theta(:),[1 Npt]);
phi=reshape(phi(:),[1 Npt]);
rad=reshape(rad(:),[1 Npt]);
Ypos{2}=zeros([globind(lmax,lmax),Npt 3]);
Yneg{2}=zeros([globind(lmax,lmax),Npt 3]);
ct=1;
for m=1:lmax
Pvals=legendrep(m,lmax,cos(theta));
Pvals2=legendrep2(m,lmax,cos(theta)); %ensure no 0/0 at a small cost
Pprime=legendrepprime(m,lmax,cos(theta),Pvals);%this prime has the -sin(theta) integrated
for inloop=m:lmax
    l=inloop;
    harmind=globind(m,l);
    locind=inloop-m+1;
    Const=1/sqrt(l*(l+1))*sqrt(((2*l+1)/(2*pi)*factorial(l-m)/factorial(l+m)));    
Ypos{2}(harmind,:,2)=-Const*Pvals2(locind,:).*(-m*sin(m*phi)); %derivative w.r.t. phi
Ypos{2}(harmind,:,3)=(Const).*Pprime(locind,:).*cos(m*phi); %derivative w.r.t. theta
Yneg{2}(harmind,:,2)=Const*Pvals2(locind,:).*(m*cos(m*phi)); %derivative w.r.t. phi
Yneg{2}(harmind,:,3)=-(Const).*Pprime(locind,:).*sin(m*phi); %derivative w.r.t. theta

end
end
Pvals=legendrep(0,lmax,cos(theta));
Pprime=legendrepprime(0,lmax,cos(theta),Pvals);
for inloop=2:lmax+1
    l=inloop-1;
    harmind=globind(0,l);m=0;
    Const=1/sqrt(l*(l+1))*sqrt(((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m)));
Ypos{2}(harmind,:,3)=Const*Pprime(inloop,:);
Yneg{2}(harmind,:,3)=Ypos{2}(harmind,:,3);
end
[Srpos,Srneg]=Yval(lmax,theta,phi);
end
function [Ypos,Yneg]=Yval(lmax,theta,phi)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
Ypos=zeros([globind(lmax,lmax),length(theta(:))]);
Yneg=zeros([globind(lmax,lmax),length(theta(:))]);
for m=1:lmax
    legen=legendrep(m,lmax,cos(theta(:)));
    for innerloop=m:lmax
        l=innerloop;
        harmind=globind(m,l);
        Const=-1/sqrt(l*(l+1))*sqrt(((2*l+1)/(2*pi)*factorial(l-m)/factorial(l+m)));
        Ypos(harmind,:)=Const*legen(l-m+1,:).*cos(m*phi);
        Yneg(harmind,:)=-Const*legen(l-m+1,:).*sin(m*phi);
    end
end
legen=legendrep(0,lmax,cos(theta(:)));
    for innerloop=1:lmax
        l=innerloop;
        m=0;
        Const=-1/sqrt(l*(l+1))*sqrt(((2*l+1)/(2*pi)*factorial(l-m)/factorial(l+m)));
Ypos(l+1,:)=Const*legen(l+1,:)/sqrt(2);
Yneg(l+1,:)=Ypos(innerloop+1,:);
    end
end
function P=legendrep(m,l,x)
%computes all legendre polynomials in an array P^{m}_{m,m+1,...,l} 0=<m<=l
P=zeros([(l-m)+1 length(x(:))]);
x=reshape(x(:),[1 length(x(:))]);
fac2=prod(1:2:(2*m-1));
P(1,:)=(-1)^(m)*fac2*(1-x.^2).^(m/2); %P^m_l l=m
if m<l
P(2,:)=x.*(2*m+1).*P(1,:);%P^m_l l=m+1
end
for lval=m+2:l
P(lval-m+1,:)=((2*lval-1)*x.*P(lval-m,:)-(lval+m-1)*P(lval-m-1,:))/(lval-m); %P^m_l l=m+1+loopit
end

end
function P=legendrep2(m,l,x)
%computes all legendre polynomials in an array P^{m}_{m,m+1,...,l} 0=<m<=l
P=zeros([(l-m)+1 length(x(:))]);
x=reshape(x(:),[1 length(x(:))]);
fac2=prod(1:2:(2*m-1));
P(1,:)=(-1)^(m)*fac2*(1-x.^2).^(m/2-1/2);
if m<l
P(2,:)=x*(2*m+1).*P(1,:);
end
for lval=m+2:l
P(lval-m+1,:)=((2*lval-1)*x.*P(lval-m,:)-(lval+m-1)*P(lval-m-1,:))/(lval-m);
end
end
function Pprime=legendrepprime(m,l,x,P)
%computes all legendre polynomials in an array P'^{m}_{m,m+1,...,l} 0=<m<=l
Pprime=zeros([(l-m)+1 length(x(:))]);
x=reshape(x(:),[1 length(x(:))]);
fac2=prod(1:2:(2*m-1));
if m~=0
Pprime(1,:)=m/2*(-1)^(m)*fac2*(-2*x).*(1-x.^2).^(m/2-1/2); 
else
Pprime(1,:)=0;     
end
if m<l
Pprime(2,:)=(2*m+1).*P(1,:).*(1-x.^2).^(1/2)+(2*m+1).*x.*Pprime(1,:);
end
for lval=m+2:l
Pprime(lval-m+1,:)=((2*lval-1)*P(lval-m,:).*(1-x.^2).^(1/2)+...
    (2*lval-1)*x.*Pprime(lval-m,:)-(lval+m-1)*Pprime(lval-m-1,:))/(lval-m);
end
Pprime=-Pprime;%-(1-x.^2)^1/2 is negative sin

end
function [x,w]=lgwt(N,a,b)

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;


end

