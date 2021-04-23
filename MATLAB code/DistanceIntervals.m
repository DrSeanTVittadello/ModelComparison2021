%%% Calculates the distance between two multisets of persistence intervals associated with two labelled simplicial complexes.

%%% Implementation:
    % SOURCE: Vittadello, S. T. and Stumpf, M. P. H., Model comparison via simplicial complexes and persistent homology, arXiv preprint: 2012.13039, 2020.

    % INPUT: 'L1' is the first multiset of persistence intervals.
    % INPUT: 'L2' is the second multiset of persistence intervals.
    
    % OUTPUT: the distance between `L1' and `L2'.

%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


function [dist] = DistanceIntervals(L1,L2)

L1 = L1(~cellfun('isempty',L1));
L2 = L2(~cellfun('isempty',L2));

indices1 = zeros(0,1);
for i=1:size(L1,1)
    indices1 = [indices1;L1{i,1}(:,1);L1{i,1}(~isinf(L1{i,1}(:,2)),2)];
end
indices1 = unique(indices1);

indices2 = zeros(0,1);
for i=1:size(L2,1)
    indices2 = [indices2;L2{i,1}(:,1);L2{i,1}(~isinf(L2{i,1}(:,2)),2)];
end
indices2 = unique(indices2);

dist = 0;
if ~isempty(setxor(indices1,indices2))
    dist = size(setxor(indices1,indices2),1);
end
