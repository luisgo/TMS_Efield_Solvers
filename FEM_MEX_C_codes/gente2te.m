function [te2te]=gente2te(te2p)
%computes array of face connectivity for point location via walk. 
te2p=te2p-1;
nte=numel(te2p)/4;
te2te=-ones([4,nte]);
mex_id_ = 'generatetetrapaths(i int, i int[xx], io int[xx])';
[te2te] = FEM(mex_id_, nte, te2p, te2te, 4, nte, 4, nte);


