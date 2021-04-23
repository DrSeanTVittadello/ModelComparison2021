%%% Calculates the distance between two labelled simplicial complexes.

%%% Implementation:
    % SOURCE: Vittadello, S. T. and Stumpf, M. P. H., Model comparison via simplicial complexes and persistent homology, arXiv preprint: 2012.13039, 2020.

    % INPUT: 'simplices1' is the first simplicial complex.
    % INPUT: 'simplices2' is the second simplicial complex.
    
    % OUTPUT: the distance between `simplices1' and `simplices2'.

%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


function [dist] = DistanceComplexes(simplices1,simplices2)

size1 = size(simplices1,1);
size2 = size(simplices2,1);

dist = 0;
for i=1:max(size1,size2)
    if i <= min(size1,size2)
        sym_diff = setxor(simplices1{i,1},simplices2{i,1},'rows');
        if ~isempty(sym_diff)
            dist = dist + size(sym_diff,1);
        end
    elseif i > size1
        dist = dist + size(simplices2{i,1},1);
    else
        dist = dist + size(simplices1{i,1},1);
    end
end
