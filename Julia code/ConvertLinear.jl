## Converts the coordinates of a simplicial complex from linear coordinates in a vector to coordinates in an array.

## Implementation:
    # INPUT: 'x' is the scalar linear coordinate of a simplicial complex in a vector.
    # INPUT: 'num_simplices' is a vector containing the number of simplices in each dimension of the simplicial complex.
    
    # OUTPUT: the coordinates of the simplicial complex in an array.

## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 23 April 2021.


function ConvertLinear(x,num_simplices)
    intervals = zeros(Int64,1,size(num_simplices,1)+1)
    for i=1:size(num_simplices,1)
        intervals[i+1] = sum(num_simplices[1:i])
    end
    
    a = intersect(findall(<(x),intervals[1,1:length(intervals)-1]),findall(>=(x),intervals[1,2:length(intervals)]))
    a = a[1]

    return reshape([a,x-intervals[a]],1,2)
end