clear all
%%%define observation points and dA/dt=omega
omega=1;
reflev=0;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
js=omega*js;
% The required data structures loaded from file are the following:
%       tetrahedron to node id (te2p): dimensions number of tetrahedrons by
%       4 nodes
%       node locations (p): dimensions number of nodes by 3 cartesian coordinate (x,y,z)
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): dimensions by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values
%Esoln1,Esoln2,Esoln3 are analytical Efields for testcase comparison
%ro3 are observation points for ro3
[X,Y,Z]=sphere(200);
ro1(1,:)=0.07*X(:).';
ro1(2,:)=0.07*Y(:).';
ro1(3,:)=0.07*Z(:).';
Z2=.07:-.001:-.01;
ro2(3,:)=Z2;
ro2(1:2,:)=0;

%% Step 1 Extract BEM surface data structures
tic
[p,~,nf,~,f2te,f2p,~,~,~,~,~,~,~,~]=createmeshbem_cl(te2p,p,conductivity);
[t2p,p,epseff]=extractsurfaces(p,f2te,conductivity,nf,f2p);
  nc=numel(rs(1,:));
  nt=numel(t2p(:,1));
  np=numel(p(:,1));1
%% Step 2 Solve BEM for charges on interfaces 
tic
xval=bemsolver(t2p',nt,p',np,epseff(:),rs,js,nc);
toc

%% Step 3 generate obseravtions
E1=computeEfields(t2p',nt,p',np,epseff(:),rs,js,nc,xval,ro1,numel(ro1(:))/3);
E2=computeEfields(t2p',nt,p',np,epseff(:),rs,js,nc,xval,ro2,numel(ro2(:))/3);
E3=computeEfields(t2p',nt,p',np,epseff(:),rs,js,nc,xval,ro3,numel(ro3(:))/3);
E4=computeEfields(t2p',nt,p',np,epseff(:),rs,js,nc,xval,ro4,numel(ro4(:))/3);

%% Step 4 plot results
subplot(2,1,1),
Emag1=sqrt(sum(E1.^2,1));%compute Emagnitude
surf(0.07*X,0.07*Y,0.07*Z,reshape(Emag1,size(X)),'edgealpha',0,'facecolor','interp');
axis equal
axis off
colorbar
subplot(2,1,2),
Emag2=sqrt(sum(E2.^2,1));%compute Emagnitude
plot(Z2,Emag2/Emag2(1));
%% Compare errors
Esoln2(isnan(Esoln2))=0;
norm(Esoln1(:)-E1(:))/norm(Esoln1(:))
norm(Esoln2(:)-E2(:))/norm(Esoln2(:))
norm(Esoln3(:)-E3(:))/norm(Esoln3(:))
norm(Esoln4(:)-E4(:))/norm(Esoln4(:))