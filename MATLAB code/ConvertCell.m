%%% Converts the coordinates of a simplicial complex from coordinates in a cell array to linear coordinates in a vector.

%%% Implementation:
    % INPUT: 'w' is the 1x2 vector of coordinates of a simplicial complex in a cell.
    % INPUT: 'num_simplices' is a row vector containing the number of simplices in each dimension of the simplicial complex.
    
    % OUTPUT: the linear coordinates of the simplicial complex in a vector.

%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


function [lin_coord] = ConvertCell(w,num_simplices)
    intervals = zeros(1,size(num_simplices,2)+1);
    for i=1:size(num_simplices,2)
        intervals(1,i+1) = sum(num_simplices(1:i));
    end
    
    if w(1,1)==1
        lin_coord = w(1,2);
    else
        lin_coord = sum(num_simplices(1,1:w(1,1)-1),2) + w(1,2);
    end
end