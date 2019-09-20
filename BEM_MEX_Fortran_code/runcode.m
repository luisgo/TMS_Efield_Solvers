function [Efield,t2p,p,epseff,xval]=runcode(te2p,p,rs,js,conductivity,ro)
%zeroth order BEM code
%       FEMord is the order of FEM desired
%       tetrahedron to node id (te2p): dimensions 4 by number of tetrahedrons
%       node locations (p): dimensions 3 cartesian coordinate (x,y,z) by number of nodes
%       electric dipole source locations (rs): dimensions 3 cartesian coordinate (x,y,z) by
%       number of dipole locations
%       electric dipole source weights (js): 3 dimensions by
%       number of dipole locations
%       conductivity (conductivity): nte dimensional column vector of conductivity values
%       ro is 3 by number of observation points
%       Efield is the Efield at points ro

%% Step 1 Extract BEM surface data structures
tic
[t2p,p,epseff]=genbemmesh(te2p,p,conductivity);
  nt=numel(t2p(1,:));
  np=numel(p(1,:));
  nc=numel(rs(1,:));
Generate_mesh_time=toc
%% Step 2 Generate BEM matrix data structures
tic

v1=p(:,t2p(2,:))-p(:,t2p(1,:));
v2=p(:,t2p(3,:))-p(:,t2p(1,:));
normal=cross(v1,v2,1);
area=sqrt(normal(1,:).^2+normal(2,:).^2+normal(3,:).^2);
normal(1,:)=normal(1,:)./area;
normal(2,:)=normal(2,:)./area;
normal(3,:)=normal(3,:)./area;
area=area/2;

%See supplemental material for explanation of these parameters
%nquad(1) is the number of quadrature points on the far triangle interactions nquad(2) number of quadrature points on near interactions valid inputs are (1,3,16). Suggested is 1 point and 16 for near
nquad=[3,16];
%iprec(1) FMM precision flag for computing E-primary iprec(2) and FMM precision flag for matrix vector multiply. Values of iprec 0,1,2,3,4,and 5 typically provide errors of FMM accuracy of 2,3,6,9,12,and 14 digits, respectivelly. Recommended value iprec=[1,0] and [2,0] for error lower than 0.1%.
iprec=[1,0];
%dnear is number of average edge lengths of the near-field interactions. Our tests indicate that this is accurate enough for FEM simulations.
dnear=0.5;
iter=0;
relres=1;
del=dnear*sqrt(mean(area(:))*2);
mex_id_ = 'generategroupct(i int[xx], i int[x], i double[xx], i int[x], i double[x], o int[x])';
[ncol] = BEM(mex_id_, t2p, nt, p, np, del, 3, nt, 1, 3, np, 1, 1, 1);

mex_id_ = 'generategroupmatlab(i int[xx], i int[x], i double[xx], i int[x], i double[x], o int[x], o int[x], i int[x])';
[col, row] = BEM(mex_id_, t2p, nt, p, np, del, ncol, 3, nt, 1, 3, np, 1, 1, ncol, ncol, 1);
row=row(col~=0);
col=col(col~=0);
ncols=numel(col);
nt3=nquad(1)*nt;
mex_id_ = 'createbemdatastructmatlab(i int[xx], i int[x], i double[xx], i int[x], i double[x], i double[xx], i int[x], i int[x], i int[x], o double[x], o double[xx], i int[x], o int[x], o int[x], o double[x])';
[potout, triafl, colch, rowch, volch] = BEM(mex_id_, t2p, nt, p, np, area, normal, row, col, ncols, nquad, 3, nt, 1, 3, np, 1, nt, 3, nt, ncols, ncols, 1, ncols, 3, nt3, 2, nt3, nt3, nt3);
col=col(abs(potout)>10^-12*max(abs(potout)));
row=row(abs(potout)>10^-12*max(abs(potout)));
potout=potout(abs(potout)>10^-12*max(abs(potout)));
Anear=sparse(col,row,potout,nt,nt);
clear col row potout;
Ach=sparse(colch,rowch,volch);
clear colch rowch volch;
Time_for_matrix=toc
%% step 3 compute rhs
tic
mex_id_ = 'computerhsmatlab(i int[xx], i int[x], i double[xx], i int[x], i double[x], i double[xx], i double[xx], i int[x], i int[x], o double[x])';
[rhs] = BEM(mex_id_, t2p, nt, p, np, epseff, rs, js, nc, iprec, 3, nt, 1, 3, np, 1, nt, 3, nc, 3, nc, 1, 2, nt);

Time_for_RHS=toc
%% step 4 solve system of equations
epseff=epseff(:)./area(:)/(4*pi);
xval=tfqmr(@(x)matvec(x,triafl,nquad(1)*nt,Ach,Anear,normal,epseff), ...
rhs,10^-7,100,[],[],rhs);

BEM_matrix_time=toc

%% Step 5 generate obseravtions
tic
Efield=computeEfields(t2p,nt,p,np,epseff(:),rs,js,nc,xval,ro,numel(ro(:))/3);
Evaluate_field_time=toc

end
function y=matvec(x,triafl,nfmm,Ach,Anear,nhat,epseff)
chspace=Ach'*x;
iprec=[0 0];
mex_id_ = 'fmmchargematlab(i double[x], i double[xx], i int[x], o double[xx], i int[x])';
[fld] = BEM(mex_id_, chspace, triafl, nfmm, iprec, nfmm, 3, nfmm, 1, 3, nfmm, 2);
y=Anear*x;
y=y+(nhat(1,:)').*(Ach*(fld(1,:)'));
y=y+(nhat(2,:)').*(Ach*(fld(2,:)'));
y=y+(nhat(3,:)').*(Ach*(fld(3,:)'));
y=epseff.*y;
y=x*0.5-y;
end



