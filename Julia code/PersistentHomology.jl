## Calculates the persistent homology of the given simplicial complex, and outputs the corresponding multiset of persistence intervals.

## Implementation:
    # SOURCE (persistent-homology algorithm): Zomorodian, A. & Carlsson, G., Computing persistent homology, Discrete & Computational Geometry, 2005, 33, 249--274.

    # INPUT: 'dim_complex' is the dimension of the simplicial complex.
    # INPUT: 'num_simplices' is the number of simplices in each dimension of the simplicial complex.
    # INPUT: 'simplices' contains the simplices of the simplicial complex.
    # INPUT: 'filtration_weight' contains the filtration weight (or filtration index) of each simplex of the simplicial complex.

    # OUTPUT: array 'L' containing the persistence intervals of the homology classes in each dimension.

## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 23 April 2021.


using Combinatorics
include("ConvertLinear.jl")
include("ConvertArray.jl")


function PersistentHomology(dim_complex,num_simplices,simplices,filtration_weight)
    ## Initial calculations.
        # Sort each simplex by the component 0-simplices to facilitate identification of the subcomplexes.
        simplices_sort = simplices
        for i=2:dim_complex+1
            simplices_sort[i] = sort(simplices[i],dims=2)
        end
        
        # Order the simplices in each dimension by increasing index of filtration.
        filtration_weight_sort = filtration_weight
        for i=1:dim_complex+1
            faces_indices = sortperm(vec(filtration_weight[i]))
            filtration_weight_sort[i] = filtration_weight[i][faces_indices,1]
            simplices_sort[i] = simplices_sort[i][faces_indices,:]
        end

        # The dimension of each simplex.
        dim_simplices = [zeros(Int8,num_simplices[1],1)]
        for i=2:dim_complex+1
            dim_simplices = vcat(dim_simplices,[(i-1)*ones(Int8,num_simplices[i],1)]) 
        end
        
    ## Calculate the persistence intervals.
        # Note that the simplices and associated data are ordered in the arrays first by
        # increasing dimension and then by increasing filtration index within each dimension.

        # Persistence intervals for each dimension.
        L = Array[[] for i=1:dim_complex+1]

        # Nonempty positions in 'chain' indicate that the corresponding simplex is a pivot in the boundary matrix,
        # and the position contains the boundary chain to which it belongs. The boundary chain is the boundary of a 
        # simplex with linear index (w.r.t. 'simplices') in 'chain_index'.
        chain = Array[[] for i=1:sum(num_simplices)] # boundary chains.
        
        # Index of the simplex corresponding to the boundary chain.
        chain_index = Array[Int64[] for i=1:sum(num_simplices)]

        # Marked simplices correspond to non-pivot columns in the boundary matrix (1/0).
        marked = Vector(undef,dim_complex+1)
        for i=1:dim_complex+1
            marked[i] = vec(zeros(Int64,size(simplices_sort[i],1),1))
        end
        
        # Consider each simplex w.r.t. the total ordering of the simplices until the end of the filtration is reached.
        for j=1:sum(num_simplices)
            # Convert the linear index j to an index in the array 'simplices'.
            index_array = reshape(ConvertLinear(j,num_simplices),1,2)
            
            dim_simplex = dim_simplices[index_array[1]][index_array[2]] # dimension of the current simplex.
            
            # Find the faces of the current simplex.
            faces = []
            if dim_simplex >= 1
                simplex = simplices_sort[index_array[1]][[index_array[2]],:]
                faces = collect(combinations(simplex,dim_simplex)) # Each row is a face.
            end
            
            # Remove the unmarked terms in 'faces'.
            faces_indices2 = zeros(Int64,size(faces,1),2)
            for t=1:size(faces,1)
                faces_indices2[t,:] = [dim_simplex findfirst(k->all(s->faces[t][s] == simplices_sort[dim_simplex][k,s],1:size(simplices_sort[dim_simplex],2)),1:size(simplices_sort[dim_simplex],1))]
            end
            
            keep = ones(Int64,size(faces,1),1)
            for i=1:size(faces,1)
                if marked[faces_indices2[i,1]][faces_indices2[i,2]] != 1
                    keep[i] = 0
                end
            end

            faces_temp = []
            for i=1:size(faces,1)
                if keep[i]==1
                    push!(faces_temp,faces[i])
                end
            end
            faces = faces_temp

            # Find the face of the simplex with the maximum index (possible pivot).
            faces_indices3 = []
            while !isempty(faces)
                faces_indices3 = zeros(Int64,size(faces,1),1)

                for t=1:size(faces,1)
                    faces_indices3[t] = findfirst(k->all(s->faces[t][s] == simplices_sort[dim_simplex][k,s],1:size(simplices_sort[dim_simplex],2)),1:size(simplices_sort[dim_simplex],1))
                end
                max_index = argmax(vec(faces_indices3))
                
                # Determine the index of the face in the array 'simplices'.
                face_index = [dim_simplex faces_indices3[max_index]]
                
                # If a pivot doesn't exist in the corresponding row of the boundary matrix then break from the while loop.
                if isempty(chain[ConvertArray(face_index,num_simplices)])
                    break
                end
                
                # If a pivot exists in the corresponding row of the boundary matrix, then use the inverse of its field coefficient
                # to eliminate the row from the chain (easy for our field Z/2Z).
                pivot = chain[ConvertArray(face_index,num_simplices)]
                for i=1:size(pivot,1)
                    term = findfirst(k->all(s->pivot[i][s] == faces[k][s],1:dim_simplex),1:size(faces,1))
                    if term !== nothing
                        faces = faces[1:size(faces,1) .!= term,:]
                        faces_indices3 = faces_indices3[1:size(faces_indices3,1) .!= term,:]
                    else
                        faces_temp = Vector(undef,size(faces,1)+1)
                        for w=1:size(faces,1)
                            faces_temp[w] = faces[w]
                        end
                        faces_temp[size(faces,1)+1] = pivot[i]
                        faces = faces_temp
                    end
                end
            end
            
            # If the simplex has empty boundary.
            if isempty(faces)                
                marked[index_array[1]][index_array[2]] = 1
               
            # Otherwise the chain 'faces' corresponds to a pivot column in the boundary matrix
            # and the term with the maximum index is the pivot.
            else
                max_index = argmax(vec(faces_indices3))
                face_index = [dim_simplex faces_indices3[max_index]]
                chain[ConvertArray(face_index,num_simplices)] = faces
                chain_index[ConvertArray(face_index,num_simplices)] = [j]

                # Add finite persistence intervals corresponding to the homology group of order dim_simplex-1.
                if isempty(L[dim_simplex])
                    L[dim_simplex] = [filtration_weight_sort[face_index[1]][face_index[2]],filtration_weight_sort[index_array[1]][index_array[2]]]
                else
                    L[dim_simplex] = vcat(L[dim_simplex],[filtration_weight_sort[face_index[1]][face_index[2]],filtration_weight_sort[index_array[1]][index_array[2]]])
                end
            end
        end
        
    # Infinite persistence intervals.
    # Add infinite persistence intervals corresponding to the homology group of order dim_simplex.
    for j=1:sum(num_simplices)
        index_array = ConvertLinear(j,num_simplices) # Convert linear index j to an index in the array 'simplices'.
        dim_simplex = dim_simplices[index_array[1]][index_array[2]]

        if marked[index_array[1]][index_array[2]] == 1 && isempty(chain[j])
            if isempty(L[dim_simplex+1])
                L[dim_simplex+1] = [filtration_weight_sort[index_array[1]][index_array[2]],Inf]
            else
                L[dim_simplex+1] = vcat(L[dim_simplex+1],[filtration_weight_sort[index_array[1]][index_array[2]],Inf])
            end
        end
    end
    
    for i=1:size(L,1)
        L[i] = transpose(reshape(L[i],2,:))
    end
    
    return L
end
