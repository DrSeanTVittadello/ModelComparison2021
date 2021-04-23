%%% Converts the coordinates of a simplicial complex from linear coordinates in a vector to coordinates in a cell array.

%%% Implementation:
    % INPUT: 'x' is the scalar linear coordinate of a simplicial complex in a vector.
    % INPUT: 'num_simplices' is a row vector containing the number of simplices in each dimension of the simplicial complex.
    
    % OUTPUT: the coordinates of the simplicial complex in a cell array.

%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


function [cell_coord] = ConvertLinear(x,num_simplices)
    intervals = zeros(1,size(num_simplices,2)+1);
    for i=1:size(num_simplices,2)
        intervals(1,i+1) = sum(num_simplices(1:i));
    end
    a = find((x>intervals(1,1:end-1))&(x<=intervals(1,2:end)));
    cell_coord = [a,x-intervals(1,a)];
end