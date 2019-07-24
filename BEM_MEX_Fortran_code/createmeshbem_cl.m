 function [p,nte,nf,te2p,f2te,f2p,te2f,tevol,tecen,fvol,fcen,nhat,triaflat,sigma]=createmeshbem(te2p,p,condu)
%INPUT:
%GMSH .msh file location for a tetrahedral mesh

%OUTPUT:
%p  points
%nte number of  tetra
%nf number of faces
%ne number of edges
% te2p tetra 2 point array
% f2te face to tetra array
% f2p face to point array
% te2f tetra to face array
% t2p triangle to point array
% e2t edge to triangle array
% e2p edge to point array
% t2e triangle to edge array
% tevol volume of each tetra
% fvol area of all faces (internal and boundary)
% nhat outward pointing normals for all facets
% tarea area of each triangle
% elength length of each edge
% vavg average volume of each tetra
% region of each tetra

faces=[te2p(:,[2,3,4])
       te2p(:,[3,4,1])
       te2p(:,[4,1,2])
       te2p(:,[1,2,3])];
    faces=sort(faces,2);
[foo,foo,faces]=unique(faces,'rows');
nf=max(faces);
nte=length(te2p(:,1))
f2te=zeros([nf 2]);
f2p=zeros([nf 3]);
te2f=zeros([nte 4]);
for i1=1:nte
    for i2=1:4
    te2f(i1,i2)=faces(i1+(i2-1)*nte);
    if f2te(faces(i1+(i2-1)*nte),1)==0
        f2te(faces(i1+(i2-1)*nte),1)=i1;
        if i2==1
        f2p(faces(i1+(i2-1)*nte),1)=te2p(i1,2);
        f2p(faces(i1+(i2-1)*nte),2)=te2p(i1,3);
        f2p(faces(i1+(i2-1)*nte),3)=te2p(i1,4);
        elseif i2==2
        f2p(faces(i1+(i2-1)*nte),1)=te2p(i1,3);
        f2p(faces(i1+(i2-1)*nte),2)=te2p(i1,4);
        f2p(faces(i1+(i2-1)*nte),3)=te2p(i1,1);
        elseif i2==3
        f2p(faces(i1+(i2-1)*nte),1)=te2p(i1,4);
        f2p(faces(i1+(i2-1)*nte),2)=te2p(i1,1);
        f2p(faces(i1+(i2-1)*nte),3)=te2p(i1,2);
        elseif i2==4
        f2p(faces(i1+(i2-1)*nte),1)=te2p(i1,1);
        f2p(faces(i1+(i2-1)*nte),2)=te2p(i1,2);
        f2p(faces(i1+(i2-1)*nte),3)=te2p(i1,3);
        end
    else
        f2te(faces(i1+(i2-1)*nte),2)=i1;
    end
    end
end
tevol=zeros([nte 1]);
tecen=zeros([nte 3]);
for tes=1:nte
    v1=p(te2p(tes,1),:)-p(te2p(tes,4),:);
    v2=p(te2p(tes,2),:)-p(te2p(tes,4),:);
    v3=p(te2p(tes,3),:)-p(te2p(tes,4),:);
    tecen(tes,:)=(p(te2p(tes,1),:)+p(te2p(tes,2),:)+p(te2p(tes,3),:)+p(te2p(tes,4),:))/4;
    cr=[v2(2)*v3(3)-v2(3)*v3(2); ...
        v2(3)*v3(1)-v2(1)*v3(3); ...
        v2(1)*v3(2)-v2(2)*v3(1)];
    tevol(tes)=abs(sum(v1(:).*cr(:)))/6;
end
fvol=zeros([nf 1]);
nhat=zeros([nf 3]);
fcen=zeros([nf 3]);
for i=1:nf
    v1=p(f2p(i,1),:)-p(f2p(i,3),:);
    v2=p(f2p(i,2),:)-p(f2p(i,3),:);
nhatt=[v1(2)*v2(3)-v1(3)*v2(2) ...
      v1(3)*v2(1)-v1(1)*v2(3) ...
      v1(1)*v2(2)-v1(2)*v2(1)];
  fvol(i)=norm(nhatt)/2;  
      if te2f(f2te(i,1),1)==i
      freenode=te2p(f2te(i,1),1);
  elseif te2f(f2te(i,1),2)==i
      freenode=te2p(f2te(i,1),2);
  elseif te2f(f2te(i,1),3)==i
      freenode=te2p(f2te(i,1),3);
  elseif te2f(f2te(i,1),4)==i
      freenode=te2p(f2te(i,1),4);
  end
  outn=p(f2p(i,1),:)-p(freenode,:); %on surface of triangle so pointing in rhat direction
if sign(sum(nhatt(:).*outn(:)))<0
   f2p(i,:)=f2p(i,[1,3,2]); %reorient triangle
end
  nhat(i,:)=nhatt/(2*fvol(i))*sign(sum(nhatt(:).*outn(:)));%nhat points toward tetra
  fcen(i,:)=(p(f2p(i,1),:)+p(f2p(i,2),:)+p(f2p(i,3),:))/3;
end

triaflat(:,3,:)=reshape(p(f2p(:,3),:)',[3 1 nf]);
triaflat(:,2,:)=reshape(p(f2p(:,2),:)',[3 1 nf]);
triaflat(:,1,:)=reshape(p(f2p(:,1),:)',[3 1 nf]);
sigma=condu(:);

end


function tri=surftri(p,t)
%SURFTRI Find surface triangles from tetrahedra mesh
%   TRI=SURFTRI(P,T)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% Form all faces, non-duplicates are surface triangles
faces=[t(:,[1,2,3]);
       t(:,[1,2,4]);
       t(:,[1,3,4]);
       t(:,[2,3,4])];
node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
tri=faces(ix(qx),:);
node4=node4(ix(qx));

% Orientation
v1=p(tri(:,2),:)-p(tri(:,1),:);
v2=p(tri(:,3),:)-p(tri(:,1),:);
v3=p(node4,:)-p(tri(:,1),:);
ix=find(dot(cross(v1,v2,2),v3,2)>0);
tri(ix,[2,3])=tri(ix,[3,2]);
end

