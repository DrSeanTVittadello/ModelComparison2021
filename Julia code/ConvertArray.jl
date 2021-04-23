## Converts the coordinates of a simplicial complex from coordinates in an array to linear coordinates in a vector.

## Implementation:
    # INPUT: 'w' is the 2x1 vector of coordinates of a simplicial complex in an array.
    # INPUT: 'num_simplices' is a vector containing the number of simplices in each dimension of the simplicial complex.
    
    # OUTPUT: the linear coordinates of the simplicial complex in a vector.

## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 23 April 2021.


function ConvertArray(w,num_simplices)
    intervals = zeros(Int64,1,size(num_simplices,1)+1)
    for i=1:size(num_simplices,1)
        intervals[i+1] = sum(num_simplices[1:i])
    end

    if w[1]==1
        lin_coord = w[2]
    else
        lin_coord = sum(num_simplices[1:w[1]-1]) + w[2]
    end
end