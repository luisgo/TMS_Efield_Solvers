clear all
omega=1;
reflev=0;
omega=1;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
dcell=[.003;.003;.003];%%choose grid spacing (will work with non-cubic voxels, however, has not been tested
[E1,rce,xval,ncell]=runcode(te2p,p,rs,omega*js,conductivity,ro1,dcell);%generates Efield at ro1
 E2=fdinterpolatorlag(rce,xval,ro2',ncell,dcell,rs',omega*js');
 E3=fdinterpolatorlag(rce,xval,ro3',ncell,dcell,rs',omega*js');
 E4=fdinterpolatorlag(rce,xval,ro4',ncell,dcell,rs',omega*js');
 
%% Compare errors
Esoln2(isnan(Esoln2))=0;
norm(Esoln1(:)-E1(:))/norm(Esoln1(:))
norm(Esoln2(:)-E2(:))/norm(Esoln2(:))
norm(Esoln3(:)-E3(:))/norm(Esoln3(:))
norm(Esoln4(:)-E4(:))/norm(Esoln4(:))