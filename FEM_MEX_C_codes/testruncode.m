clear all
FEMord=3; %choose FEM order
omega=1;
reflev=0;
load(strcat('../samplescenario',num2str(reflev),'.mat'));
js=js*omega;
[E1,x,te2p2,p2,te2te]=runcode(te2p',p',conductivity,rs,js,ro1,FEMord);%generates Efield at ro1



%% Step 2 evaluate field at aditional desired locations
E2=FEMinterpolator(te2te,te2p2,p2,rs,js,x,ro2,FEMord);
E3=FEMinterpolator(te2te,te2p2,p2,rs,js,x,ro3,FEMord);
E4=FEMinterpolator(te2te,te2p2,p2,rs,js,x,ro4,FEMord);

%% plot results
Esoln2(isnan(Esoln2))=0;

subplot(1,3,1),
[X,Y,Z]=sphere(200);
Emag1=sqrt(sum(E1.^2,1));%compute Emagnitude
surf(0.07*X,0.07*Y,0.07*Z,reshape(Emag1,size(X)),'edgealpha',0,'facecolor','interp');
axis equal
axis off
colorbar
axis square
view([0 0 1])
title('E-field magnitude (V/m) on 70 mm Surface')
subplot(1,3,2),
Emag2=sqrt(sum(E2.^2,1));%compute Emagnitude
Emagsoln2=sqrt(sum(Esoln2.^2,1));%compute Emagnitude
plot(ro2(3,:),Emag2,'b');
hold on
plot(ro2(3,:),Emagsoln2,'--r');
axis square
xlabel('z (x=0;y=0)')
ylabel('E-field magnitude (V/m)')
legend('FEM','Analytical')
xlim([-0.01 0.07])
subplot(1,3,3),
semilogy(ro2(3,:),Emag2,'b');
hold on
semilogy(ro2(3,:),Emagsoln2,'--r');
semilogy(ro2(3,:),abs((Emag2(:)-Emagsoln2(:))),'g');
xlabel('z (x=0;y=0)')
ylabel('E-field magnitude (V/m)')
legend('FEM','Analytical','Difference')
xlim([-0.01 0.07])
axis square
%% Compare errors
norm(Esoln1(:)-E1(:))/norm(Esoln1(:))
norm(Esoln2(:)-E2(:))/norm(Esoln2(:))
norm(Esoln3(:)-E3(:))/norm(Esoln3(:))
norm(Esoln4(:)-E4(:))/norm(Esoln4(:))