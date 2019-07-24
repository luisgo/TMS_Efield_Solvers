clear all
FEMord=1; %choose FEM order
omega=1;
reflev=0;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
js=js*omega;
[E1,x,te2p2,p2]=runcode(te2p,p,rs,omega*js,conductivity,ro1,FEMord);%generates Efield at ro1
E2=FEMinterpolator(ro2,x,te2p2,p2,rs,omega*js,FEMord);
E3=FEMinterpolator(ro3,x,te2p2,p2,rs,omega*js,FEMord);
E4=FEMinterpolator(ro4,x,te2p2,p2,rs,omega*js,FEMord);
%% Compare errors
Esoln2(isnan(Esoln2))=0;
norm(Esoln1(:)-E1(:))/norm(Esoln1(:))
norm(Esoln2(:)-E2(:))/norm(Esoln2(:))
norm(Esoln3(:)-E3(:))/norm(Esoln3(:))
norm(Esoln4(:)-E4(:))/norm(Esoln4(:))