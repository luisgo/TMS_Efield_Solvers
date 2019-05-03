function [U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg)
%LFMM3DPART Laplace particle target FMM in R^3.
%
% Laplace FMM in R^3: evaluate all pairwise particle
% interactions (ignoring self-interactions) and interactions with targets.
%
% [U]=LFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC);
%
% [U]=LFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD);
%
% [U]=LFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET);
%
% [U]=LFMM3DPART(IPREC,NSOURCE,SOURCE,...
%         IFCHARGE,CHARGE,IFDIPOLE,DIPSTR,DIPVEC,IFPOT,IFFLD,...
%         NTARGET,TARGET,IFPOTTARG,IFFLDTARG);
%
%
% This subroutine evaluates the Laplace potential and field due
% to a collection of charges and dipoles. We use (1/r) for the 
% Green's function, without the (1/4 pi) scaling. 
% Self-interactions are not-included.
%
% Input parameters:
% 
% iprec - FMM precision flag
%
%             -2 => tolerance =.5d0   =>  
%             -1 => tolerance =.5d-1  =>  1 digit 
%              0 => tolerance =.5d-2  =>  2 digits
%              1 => tolerance =.5d-3  =>  3 digits
%              2 => tolerance =.5d-6  =>  6 digits
%              3 => tolerance =.5d-9  =>  9 digits
%              4 => tolerance =.5d-12 => 12 digits
%              5 => tolerance =.5d-15 => 15 digits
%
% nsource - number of sources
% source - real (3,nsource): source locations
% ifcharge - charge computation flag
%
%         0 => do not compute
%         1 => include charge contribution
% 
% charge - complex (nsource): charge strengths 
% ifdipole - dipole computation flag
%
%         0 => do not compute
%         1 => include dipole contributions
% 
% dipole - complex (nsource): dipole strengths
% dipvec - real (3,source): dipole orientation vectors
%
% ifpot - potential computation flag, 1 => compute the potential, otherwise no
% iffld - field computation flag, 1 => compute the field, otherwise no
%
% ntarget - number of targets
% target - real (3,ntarget): target locations
%
% ifpottarg - target potential computation flag, 
%      1 => compute the target potential, otherwise no
% iffldtarg - target field computation flag, 
%      1 => compute the target field, otherwise no
%
% Output parameters: 
%
% U.pot - complex (nsource) - potential at source locations
% U.fld - complex (3,nsource) - field (i.e. -gradient) at source locations
% U.pottarg - complex (ntarget) - potential at target locations
% U.fldtarg - complex (3,ntarget) - field (i.e. -gradient) at target locations
%
% U.ier - error return code
%
%             ier=0     =>  normal execution
%             ier=4     =>  cannot allocate tree workspace
%             ier=8     =>  cannot allocate bulk FMM  workspace
%             ier=16    =>  cannot allocate mpole expansion workspace in FMM
%

if( nargin == 8 ) 
  ifpot = 1;
  iffld = 1;
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 10 ) 
  ntarget = 0;
  target = zeros(3,1);
  ifpottarg = 0;
  iffldtarg = 0;
end

if( nargin == 12 ) 
  ifpottarg = 1;
  iffldtarg = 1;
end

ifcharge = double(ifcharge); ifdipole = double(ifdipole);
ifpot = double(ifpot); iffld = double(iffld);
ifpottarg = double(ifpottarg); iffldtarg = double(iffldtarg);

pot=0;
fld=zeros(3,1);
pottarg=0;
fldtarg=zeros(3,1);

if( ifpot == 1 ), pot=zeros(1,nsource)+1i*zeros(1,nsource); end;
if( iffld == 1 ), fld=zeros(3,nsource)+1i*zeros(3,nsource); end;
if( ifpottarg == 1 ), pottarg=zeros(1,ntarget)+1i*zeros(1,ntarget); end;
if( iffldtarg == 1 ), fldtarg=zeros(3,ntarget)+1i*zeros(3,ntarget); end;

ier=0;

if( ntarget == 0 ) 
mex_id_ = 'lfmm3dpartself(io int[x], i int[x], i int[x], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[])';
[ier, pot, fld] = fmm3d_r2012a(mex_id_, ier, iprec, nsource, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, pot, iffld, fld, 1, 1, 1, 3, nsource, 1, 1, 3, nsource, 1, 1);
else
mex_id_ = 'lfmm3dparttarg(io int[x], i int[x], i int[x], i double[xx], i int[x], i dcomplex[], i int[x], i dcomplex[], i double[xx], i int[x], io dcomplex[], i int[x], io dcomplex[], i int[x], i double[], i int[x], io dcomplex[], i int[x], io dcomplex[])';


nsa=2^23-2;nsas=2^23-2;
nst=floor(ntarget/nsa);
nsts=floor(nsource/nsas);
for i=1:nst+1
if i==nst+1
nran=((i-1)*nsa+1):ntarget;
else
nran=((i-1)*nsa+1):(i*nsa);
end
nran=nran(:);
ntarget1=numel(nran);
if ntarget1==0
continue;
end
target1=target(:,nran);
pottarg1=0;fldtarg1=zeros(3,1);
pottarg2=0;fldtarg2=zeros(3,1);
if( ifpottarg == 1 ), pottarg1=zeros(1,ntarget1)+1i*zeros(1,ntarget1); end;
if( ifpottarg == 1 ), pottarg2=zeros(1,ntarget1)+1i*zeros(1,ntarget1); end;
if( iffldtarg == 1 ), fldtarg1=zeros(3,ntarget1)+1i*zeros(3,ntarget1); end;
if( iffldtarg == 1 ), fldtarg2=zeros(3,ntarget1)+1i*zeros(3,ntarget1); end;

for j=1:nsts+1%source loop
if j==nsts+1%batch ranges
nrans=((j-1)*nsas+1):nsource;
else
nrans=((j-1)*nsas+1):(j*nsas);
end
nsource1=numel(nrans);
if nsource1==0
continue;
end
if j==1
[ier, pot, fld, pottarg1, fldtarg1] = fmm3d_r2012a(mex_id_, ier, iprec, nsource1, source(:,nrans),...
 ifcharge, charge(nrans), ifdipole, dipstr(nrans), dipvec(:,nrans), ifpot, pot, iffld, fld,...
 ntarget1, target1, ifpottarg, pottarg1, iffldtarg, fldtarg1, 1, 1, 1, 3, nsource1, 1, 1, 3, nsource1, 1, 1, 1, 1, 1);
else
[ier, pot, fld, pottarg2, fldtarg2] =  fmm3d_r2012a(mex_id_, ier, iprec, nsource1, source(:,nrans),...
 ifcharge, charge(nrans), ifdipole, dipstr(nrans), dipvec(:,nrans), ifpot, pot, iffld, fld,...
 ntarget1, target1, ifpottarg, pottarg2, iffldtarg, fldtarg2, 1, 1, 1, 3, nsource1, 1, 1, 3, nsource1, 1, 1, 1, 1, 1);
pottarg1=pottarg1+pottarg2;
fldtarg1=fldtarg1+fldtarg2;
end
end
if( ifpottarg == 1 ),pottarg(:,nran)=pottarg1; end;
if( iffldtarg == 1 ),fldtarg(:,nran)=fldtarg1; end;
end

end

if( ifpot == 1 ), U.pot=pot; end
if( iffld == 1 ), U.fld=fld; end
if( ifpottarg == 1 ), U.pottarg=pottarg; end
if( iffldtarg == 1 ), U.fldtarg=fldtarg; end
U.ier=ier;



