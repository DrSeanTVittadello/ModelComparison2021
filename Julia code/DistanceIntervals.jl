## Calculates the distance between two multisets of persistence intervals associated with two labelled simplicial complexes.

## Implementation:
    # SOURCE: Vittadello, S. T. and Stumpf, M. P. H., Model comparison via simplicial complexes and persistent homology, arXiv preprint: 2012.13039, 2020.

    # INPUT: 'L1' is the first multiset of persistence intervals.
    # INPUT: 'L2' is the second multiset of persistence intervals.
    
    # OUTPUT: the distance between `L1' and `L2'.

## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 21 April 2021.


function DistanceIntervals(L1,L2)
    indices1 = Vector()
    indices2 = Vector()

    for i=1:size(L1,1)
        indices1 = vcat(indices1,filter!(x->x!==Inf,L1[i][:]))
    end
    indices1 = convert(Vector{Int64},unique(indices1))

    for i=1:size(L2,1)
        indices2 = vcat(indices2,filter!(x->x!==Inf, L2[i][:]))
    end
    indices2 = convert(Vector{Int64},unique(indices2))
    
    dist = size(setdiff(indices1,indices2),1) + size(setdiff(indices2,indices1),1)

    return dist
end
