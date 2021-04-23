## Calculates the distance between two labelled simplicial complexes.

## Implementation:
    # SOURCE: Vittadello, S. T. and Stumpf, M. P. H., Model comparison via simplicial complexes and persistent homology, arXiv preprint: 2012.13039, 2020.

    # INPUT: 'simplices1' is the first simplicial complex.
    # INPUT: 'simplices2' is the second simplicial complex.
    
    # OUTPUT: the distance between `simplices1' and `simplices2'.

## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 23 April 2021.


function DistanceComplexes(simplices1,simplices2)
    size1 = size(simplices1,1)
    size2 = size(simplices2,1)
    
    dist = 0
    for i=1:max(size1,size2)
        simplices_rows1 = Vector()
        if i <= size1
            for j=1:size(simplices1[i],1)
                push!(simplices_rows1,simplices1[i][[j],:])
            end
        end

        simplices_rows2 = Vector()
        if i <= size2
            for j=1:size(simplices2[i],1)
                push!(simplices_rows2,simplices2[i][[j],:])
            end
        end

        if isempty(simplices_rows1) && !isempty(simplices_rows2)
            dist = dist + size(simplices_rows2,1)
        elseif !isempty(simplices_rows1) && isempty(simplices_rows2)
            dist = dist + size(simplices_rows1,1)
        elseif !isempty(simplices_rows1) && !isempty(simplices_rows2)
            dist = dist + size(setdiff(simplices_rows1,simplices_rows2),1) + size(setdiff(simplices_rows2,simplices_rows1),1)
        end
    end
    return dist
end
