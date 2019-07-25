clear all
omega=1;
reflev=0;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
dcell=[.003;.003;.003];%%choose grid spacing (will work with non-cubic voxels, however, has not been tested

te2p=te2p';p=p';
%four layer example
% The required data structures loaded from file are the following:
%       FEMord is the order of FEM desired
%       tetrahedron to node id (te2p): dimensions number of tetrahedrons by
%       4 nodes
%       node locations (p): dimensions number of nodes by 3 cartesian coordinate (x,y,z)
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): dimensions 3 vector weight (Jx,Jy,Jz) by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values


%%%%%%%%%%%%choose points where you want to evalute the E-field (ro): ro is
%%%%%%%%%%%%3 by number of evaluation points.
[X,Y,Z]=sphere(200);
ro1(1,:)=0.07*X(:).';
ro1(2,:)=0.07*Y(:).';
ro1(3,:)=0.07*Z(:).';
Z2=.07:-.001:-.01;
ro2(3,:)=Z2;
ro2(1:2,:)=0;



%% Step 1 run FDM
js=js*omega;
tic
[Efield,RELRES,ITER,RESVEC,condu,Einc,rce,xval,ncell,dcell]=runfdlag(te2p,p,conductivity,rs,js,dcell);
toc
%% Step 2 evaluate field at desired locations
 E1=fdinterpolatorlag(rce,xval,ro1',ncell,dcell,rs,js);
 E2=fdinterpolatorlag(rce,xval,ro2',ncell,dcell,rs,js);
 E3=fdinterpolatorlag(rce,xval,ro3',ncell,dcell,rs,js);
 E4=fdinterpolatorlag(rce,xval,ro4',ncell,dcell,rs,js);

%% Step 3 plot results
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