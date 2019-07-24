function [te2p2,p2]=femgenmesh(te2p,p,FEMord)
nte=numel(te2p(1,:));
np=numel(p(1,:));
if FEMord==1
te2p2=te2p;p2=p;
elseif FEMord==2
%% define edge elements
te2p2=zeros([10 nte]);
te2p2(1:4,:)=te2p;
faces=[te2p([1,2],:)';%5
       te2p([1,3],:)';%6
       te2p([1,4],:)';%7
       te2p([2,3],:)';%8
       te2p([2,4],:)';%9
       te2p([3,4],:)'];%10
faces=sort(faces,2);
[~,ix,jx]=unique(faces,'rows','stable');
e2p=faces(ix,:);
te2p2(5:10,:)=reshape(jx+np,[nte 6])';
%% define interppoints
ne=numel(e2p(:,1));
p2=zeros([3 np+ne]);
p2(:,1:np)=p;
p2(:,np+1:np+ne)=(p(:,e2p(:,1))+p(:,e2p(:,2)))/2;
elseif FEMord==3
nbasis=20;
te2p2=zeros([nbasis nte]);
%% define node elements
te2p2(1:4,:)=te2p;
%% define edge elements
em=[1,2;1,3;1,4;2,3;2,4;3,4];
faces=[te2p(em(1,:),:)';
   te2p(em(2,:),:)';
   te2p(em(3,:),:)';
   te2p(em(4,:),:)';
   te2p(em(5,:),:)';
   te2p(em(6,:),:)'];
faces=sort(faces,2);
[~,ix,jx]=unique(faces,'rows');
e2p=faces(ix,:);
ne=numel(e2p(:,1));
jx=reshape(jx,[nte 6]);
for j=1:6% ensure that edge elements are arranged properly
for i=1:nte
if te2p(em(j,1),i)==e2p(jx(i,j),1)
te2p2([2*j-1+4;2*j+4],i)=[jx(i,j)+np;jx(i,j)+np+ne];
elseif te2p(em(j,2),i)==e2p(jx(i,j),1)
te2p2([2*j-1+4;2*j+4],i)=[jx(i,j)+np+ne;jx(i,j)+np];
end
end
end
%% define triangle elements
faces=[te2p([1,2,3],:)';
   te2p([1,2,4],:)';
   te2p([1,3,4],:)';
   te2p([2,3,4],:)'];
faces=sort(faces,2);
[~,ix,jx]=unique(faces,'rows');
t2p=faces(ix,:);
nt=numel(t2p(:,1));
te2p2([17 18 19 20],:)=reshape(jx+np+2*ne,[nte 4])';
%% define interpolation points
ne=numel(e2p(:,1));
p2=zeros([3 np+ne]);
p2(:,1:np)=p;
p2(:,np+1:np+ne)=(2*p(:,e2p(:,1))/3+p(:,e2p(:,2))/3);
p2(:,np+ne+1:np+2*ne)=(p(:,e2p(:,1))/3+2*p(:,e2p(:,2))/3);
p2(:,np+2*ne+1:np+2*ne+nt)=(p(:,t2p(:,1))+p(:,t2p(:,2))+...
                        p(:,t2p(:,3)))/3;
end

