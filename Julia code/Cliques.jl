## Calculates the higher-dimensional cliques of the given simplicial complex.
    # Given a simplicial complex with dimension d, Cliques.jl calculates the cliques in dimension d+1, which is 'min_dim',
    # and then iteratively in each higher dimension until a dimension is reached which is either equal to 'max_dim' or for
    # which no more cliques exist.

## Implementation:
    # INPUT: 'simplices' contains the simplices of the complex.
    # INPUT: 'min_dim' is the minimum dimension clique.
    # INPUT: 'max_dim' is the maximum dimension clique.

    # OUTPUT: 'simplices' with the higher-dimensional cliques appended.

## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 23 April 2021.


using Combinatorics


function Cliques(simplices,min_dim,max_dim)
    for k=min_dim:max_dim
        simplices_check = reshape([],0,k+1) # possible k-simplices.
        for i=1:size(simplices[k],1) # (k-1)-simplices.
            diff_vertices = setdiff(simplices[1],simplices[k][i,:]) # vertices not in the (k-1)-simplex.
            for j=1:size(diff_vertices,1)
                simplices_check = vcat(simplices_check,[reshape(simplices[k][i,:],1,k) diff_vertices[j]])
            end
        end
        simplices_check = sort(simplices_check,dims=2)
        
        # Remove duplicate simplices from 'simplices_check'.
        simplices_temp = reshape([],0,k+1)
        for i=1:size(simplices_check,1)
            simplex = simplices_check[[i],:]
            if isempty(simplices_temp)
                simplices_temp = vcat(simplices_temp,simplex)
            elseif findfirst(t->all(s->simplices_temp[[t],s] == simplex[[1],s],1:size(simplices_temp,2)),1:size(simplices_temp,1)) === nothing
                simplices_temp = vcat(simplices_temp,simplex)
            end
        end
        simplices_check = simplices_temp

        next_dim = 0 # check for cliques in the next dimension? (Yes = 1, No = 0).
        for m=1:size(simplices_check,1)
            faces_check = collect(combinations(simplices_check[m,:],k)) # check for the (k-1)-faces.

            # Find (k-1)-simplices in 'faces_check' that aren't simplices in the given simplicial complex.
            faces_check = transpose(hcat(faces_check...))
            faces_rows = Vector()
            for i=1:size(faces_check,1)
                push!(faces_rows,faces_check[i,:])
            end
            simplices_rows = Vector()
            for i=1:size(simplices[k],1)
                push!(simplices_rows,simplices[k][i,:])
            end

            if isempty(setdiff(faces_rows,simplices_rows)) # true when all of the proper faces of the simplex are in the given simplicial complex.
                next_dim = 1 # check for cliques in the next dimension.
                if size(simplices,1) < k+1 # true when there are currently no k-dimensional simplices in the simplicial complex.
                    push!(simplices,reshape(convert(Array{Int64,1},simplices_check[m,:]),1,:))
                else
                    simplices[k+1] = vcat(simplices[k+1],reshape(convert(Array{Int64,1},simplices_check[m,:]),1,:))      
                end
            end
        end
        if next_dim == 0 # don't check for cliques in the next dimension.
            break
        end
    end
    return filter!(!isempty,simplices)
end
