## Main file for model comparison via simplicial complexes and persistent homology.

## Implementation:
    # SOURCE (model-comparison methodology): Vittadello, S. T. and Stumpf, M. P. H., Model comparison via simplicial complexes and persistent homology, arXiv preprint: 2012.13039, 2020.
    # SOURCE (persistent-homology algorithm): Zomorodian, A. & Carlsson, G., Computing persistent homology, Discrete & Computational Geometry, 2005, 33, 249--274.

    # INPUT: 'SC_dist' indicates whether to calculate the distances directly between the simplicial complexes ("Y" or "N");
    #        the distances between the simplicial complexes are contained in the upper triangular matrix 'dist_complexes'.
    # INPUT: 'PI_dist' indicates whether to calculate the distances between the persistence intervals ("Y" or "N");
    #        the distances between the persistence intervals are contained in the upper triangular matrix 'dist_intervals'.
    # INPUT: 'model_components' contains all of the model components required by the models for comparison.
    # INPUT: 'dim_ref' is the maximum possible dimension for the reference complex (set equal to `Inf' if no maximum is required).
    # INPUT: 'num_models' is the number of models to compare.
    # INPUT: for each model the following is required -
             # INPUT: 'model_num' is the model number, beginning at 1 and incrementing by 1 for subsequent models.
             # INPUT: 'model_components_i' contains all of the model components for model number i.
             # INPUT: 'simplices' contains the simplices of the simplicial complex associated to the model;
             #        the 0-simplices (or vertices) are determined by the specified components of the model in 'model_components_i';
             #        input the simplices of dimension 1 and higher, as required;
             #        if simplices up to dimension d are given, then the function Cliques.jl may be employed to incrementally obtain 
             #        the cliques of dimensions higher than d.

    # OUTPUT (when SC_dist = "Y"): 'dist_complexes' containing the distances between the simplicial complexes associated with the
    #                              models for comparison.
    # OUTPUT (when PI_dist = "Y"): array 'L_all' containing the persistence intervals of the homology classes in each dimension.
    # OUTPUT (when PI_dist = "Y"): 'dist_intervals' containing the distances between the persistence intervals associated with the
    #                              models for comparison.
    
## Author: Sean T. Vittadello.
## Affiliation: The University of Melbourne.
## Date: 23 April 2021.


include("Cliques.jl")
include("PersistentHomology.jl")
include("DistanceComplexes.jl")
include("DistanceIntervals.jl")


## Input.
SC_dist = "Y"
PI_dist = "Y"

model_components = ["Morphogen 1","Diffusion 1","Degradation 1","Production 1","Basal production 1","Influx 1","Outflux 1","Morphogen 1 bound","Morphogen 2","Diffusion 2","Degradation 2","Basal production 2","Influx 2","Morphogen 3","Diffusion 3","Degradation 3","Production 3","Substrate 1","Diffusion of Substrate 1","Degradation of Substrate 1","Basal production of Substrate 1","Modulator 1","Diffusion of Modulator 1","Degradation of Modulator 1","Production of Modulator 1","Annihilation between Morphogens 1 and 2","Self-activation of Morphogen 1","Activation of Morphogen 2 by Morphogen 1","Inhibition of Morphogen 1 by Morphogen 2","Inhibition of Morphogen 1 by Morphogen 3","Inhibition of Morphogen 3 by Morphogen 1","Inhibition of an inhibition by Morphogen 2","Production of Morphogen 2 by Morphogen 1","Depletion of Substrate 1 by Morphogen 1","Modulation of Diffusion 1 by Modulator 1","Modulation of Degradation 1 by Modulator 1","Inhibition of Modulator 1 by Morphogen 1","Adsorption of Morphogen 1","Desorption of Morphogen 1 bound","Monotonic gradient","Oscillatory gradient","Local scale-invariance","Global scale-invariance"]

dim_ref = Inf

num_models = 9


## Initialisation.
dim_max = Array[[] for i=1:num_models] # dimension of the simplex spanned by the vertices for each model.
simplices_all = [] # contains the simplicial complexes associated with the models.
L_all = [] # contains the persistence intervals for the simplicial complexes associated with the models.


## Input the simplicial complex of each model for comparison (here, positional-information (PI) and Turing-pattern (TP) models).
# Model PI 1 - linear gradient.
    model_num = 1

    model_components_1 = ["Morphogen 1","Diffusion 1","Influx 1","Outflux 1","Monotonic gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_1),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 6;1 7;1 40;1 43;2 6;2 7;2 40;2 43;6 40;6 43;7 40;7 43;40 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model PI 2 - synthesis-diffusion-degradation (SDD)
    model_num = 2

    model_components_2 = ["Morphogen 1","Diffusion 1","Degradation 1","Influx 1","Monotonic gradient"]

    vertices = findfirst.(isequal.(model_components_2),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 6;1 40;2 3;2 6;2 40;3 40;6 40]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model PI 3 - Opposing gradients
    model_num = 3

    model_components_3 = ["Morphogen 1","Diffusion 1","Degradation 1","Influx 1","Morphogen 2","Diffusion 2","Degradation 2","Influx 2","Monotonic gradient","Local scale-invariance"]

    vertices = findfirst.(isequal.(model_components_3),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 6;1 40;1 42;2 3;2 6;2 40;2 42;3 40;3 42;6 40;6 42;9 10;9 11;9 13;9 40;9 42;10 11;10 13;10 40;10 42;11 40;11 42;13 40;13 42;40 42]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model PI 4 - Annihilation
    model_num = 4

    model_components_4 = ["Morphogen 1","Diffusion 1","Degradation 1","Influx 1","Morphogen 2","Diffusion 2","Degradation 2","Influx 2","Annihilation between Morphogens 1 and 2","Monotonic gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_4),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 6;1 26;1 40;1 43;2 3;2 6;2 26;2 40;2 43;3 40;3 43;6 40;6 43;9 10;9 11;9 13;9 26;9 40;9 43;10 11;10 13;10 26;10 40;10 43;11 40;11 43;13 40;13 43;26 40;26 43;40 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model PI 5 - Induction-contraction (active modulation)
    model_num = 5

    model_components_5 = ["Morphogen 1","Diffusion 1","Degradation 1","Influx 1","Modulator 1","Diffusion of Modulator 1","Degradation of Modulator 1","Production of Modulator 1","Modulation of Diffusion 1 by Modulator 1","Modulation of Degradation 1 by Modulator 1","Inhibition of Modulator 1 by Morphogen 1","Monotonic gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_5),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 6;1 37;1 40;1 43;2 3;2 6;2 35;2 40;2 43;3 36;3 40;3 43;6 40;6 43;22 23;22 24;22 25;22 35;22 36;22 40;22 43;23 24;23 25;23 40;23 43;24 40;24 43;25 37;25 40;25 43;35 40;35 43;36 40;36 43;37 40;37 43;40 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# # Model TP 1 - Activator-inhibitor
    model_num = 6

    model_components_6 = ["Morphogen 1","Diffusion 1","Degradation 1","Basal production 1","Morphogen 2","Diffusion 2","Degradation 2","Basal production 2","Self-activation of Morphogen 1","Activation of Morphogen 2 by Morphogen 1","Inhibition of Morphogen 1 by Morphogen 2","Oscillatory gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_6),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 5;1 27;1 28;1 29;1 41;1 43;2 3;2 5;2 27;2 28;2 29;2 41;2 43;3 41;3 43;5 27;5 41;5 43;9 10;9 11;9 12;9 28;9 29;9 41;9 43;10 11;10 12;10 28;10 29;10 41;10 43;11 41;11 43;12 41;12 43;27 41;27 43;28 29;28 41;28 43;29 41;29 43;41 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model TP 2 - Substrate depletion
    model_num = 7

    model_components_7 = ["Morphogen 1","Diffusion 1","Degradation 1","Basal production 1","Substrate 1","Diffusion of Substrate 1","Degradation of Substrate 1","Basal production of Substrate 1","Self-activation of Morphogen 1","Depletion of Substrate 1 by Morphogen 1","Oscillatory gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_7),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 5;1 27;1 34;1 41;1 43;2 3;2 5;2 27;2 34;2 41;2 43;3 41;3 43;5 27;5 41;5 43;18 19;18 20;18 21;18 34;18 41;18 43;19 20;19 21;19 41;19 43;20 41;20 43;21 41;21 43;27 41;27 43;34 41;34 43;41 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model TP 3 - Inhibition of an inhibition
    model_num = 8

    model_components_8 = ["Morphogen 1","Diffusion 1","Degradation 1","Production 1","Morphogen 2","Diffusion 2","Degradation 2","Morphogen 3","Diffusion 3","Degradation 3","Production 3","Inhibition of Morphogen 1 by Morphogen 3","Inhibition of Morphogen 3 by Morphogen 1","Inhibition of an inhibition by Morphogen 2","Production of Morphogen 2 by Morphogen 1","Oscillatory gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_8),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 4;1 31;1 33;1 41;1 43;2 3;2 4;2 30;2 31;2 32;2 33;2 41;2 43;3 41;3 43;4 30;4 31;4 32;4 33;4 41;4 43;9 10;9 11;9 32;9 33;9 41;9 43;10 11;10 30;10 31;10 32;10 33;10 41;10 43;11 41;11 43;14 15;14 16;14 17;14 30;14 41;14 43;15 16;15 17;15 30;15 31;15 32;15 33;15 41;15 43;16 41;16 43;17 30;17 31;17 32;17 33;17 41;17 43;30 41;30 43;31 32;31 41;31 43;32 41;32 43;33 41;33 43;41 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


# Model TP 4 - Modulation
    model_num = 9

    model_components_9 = ["Morphogen 1","Diffusion 1","Degradation 1","Basal production 1","Morphogen 1 bound","Morphogen 2","Diffusion 2","Degradation 2","Basal production 2","Self-activation of Morphogen 1","Activation of Morphogen 2 by Morphogen 1","Inhibition of Morphogen 1 by Morphogen 2","Adsorption of Morphogen 1","Desorption of Morphogen 1 bound","Oscillatory gradient","Global scale-invariance"]

    vertices = findfirst.(isequal.(model_components_9),(model_components,))

    dim_max[model_num] = [size(vertices,1)-1] # dimension of the simplex spanned by the vertices.

    # 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = Array[reshape(vertices,:,1),[1 2;1 3;1 5;1 27;1 28;1 29;1 38;1 39;1 41;1 43;2 3;2 5;2 27;2 28;2 29;2 38;2 39;2 41;2 43;3 8;3 41;3 43;5 27;5 41;5 43;8 27;8 28;8 29;8 38;8 39;8 41;8 43;9 10;9 11;9 12;9 28;9 29;9 41;9 43;10 11;10 12;10 27;10 28;10 29;10 38;10 39;10 41;10 43;11 41;11 43;12 41;12 43;27 29;27 41;27 43;28 29;28 41;28 43;29 41;29 43;38 41;38 43;39 41;39 43;41 43]]
    simplices = Cliques(simplices,size(simplices,1),size(simplices[1],1)-1) # find cliques, as required.
    push!(simplices_all,simplices)


## Calculations.
# Distance between the simplicial complexes associated with the models.
if SC_dist == "Y"
    dist_complexes = zeros(Int64,num_models,num_models)
    for i=1:num_models
        for j=(i+1):num_models
            dist_complexes[i,j] = DistanceComplexes(simplices_all[i],simplices_all[j])
        end
    end
end


# Distance between the multisets of persistence intervals associated with the models.
if PI_dist == "Y"
    ## Reference simplicial complex.
        # Maximum dimension of all simplicial representations of models.
        maximum_dim = 0
        for i=1:num_models
            maximum_dim = max(maximum_dim,size(simplices_all[i],1)-1)
        end

        # All required simplices of the reference complex.
        simplices_ref = [] # here we only need the simplices associated with the simplicial representations, rather than the full
                           # reference complex, since we can use any flat filtration to calculate the distances between the models;
                           # the required indices of the flat filtration are based on dimension - extension to the full reference
                           # complex, if required, provides multiple options for the filtration.
        for i=1:num_models
            for j=1:size(simplices_all[i],1)
                if size(simplices_ref,1) >= j
                    simplices_ref[j] = vcat(simplices_ref[j],simplices_all[i][j])
                else
                    push!(simplices_ref,simplices_all[i][j])
                end
            end
        end

        # Find the unique simplices (optionally ordered lexicographically in each dimension).
        for i=1:size(simplices_ref,1)
            simplices_temp = reshape([],0,i)
            for j=1:size(simplices_ref[i],1)
                simplex = reshape(simplices_ref[i][j,:],1,:)
                if isempty(simplices_temp)
                    simplices_temp = vcat(simplices_temp,simplex)
                elseif findfirst(t->all(s->simplices_temp[[t],s] == simplex[[1],s],1:size(simplices_temp,2)),1:size(simplices_temp,1)) === nothing
                    simplices_temp = vcat(simplices_temp,simplex)
                end
            end
            simplices_ref[i] = sortslices(simplices_temp,dims=1) # optionally ordered lexicographically.
        end

        # The flat filtration on the reference simplicial complex (multiple options).
        num_ref = [size(simplices_ref[i],1) for i=1:maximum_dim+1]
        filtration_weight_ref = [Int64[] for i=1:size(num_ref,1)]
        filtration_weight_ref[1] = 1:num_ref[1]
        for i=2:size(num_ref,1)
            filtration_weight_ref[i] = sum(num_ref[1:i-1])+1:sum(num_ref[1:i])
        end

    ## Persistence intervals for each model, based on the induced filtrations
        for i=1:num_models
            simplices = simplices_all[i] # simplices for model i.
            filtration_weight = [] # contains the induced-filtration weights for the given model.
            for j=1:size(simplices,1) # simplices of dimension j-1.
                weights = reshape(Int64[],0,1) # temporarily store the filtration weights for dimension j-1.
                for k=1:size(simplices[j],1) # each simplex of dimension j-1.
                    simplex = simplices[j][[k],:]
                    index = findfirst(t->all(s->simplices_ref[j][t,:][s] == simplex[[1],s][1],1:size(simplex,2)),1:size(simplices_ref[j],1))
                    weights = vcat(weights,filtration_weight_ref[j][index])
                end
                push!(filtration_weight,weights)
            end
       
            num_simplices = [size(simplices[i],1) for i=1:size(simplices,1)]
            L = PersistentHomology(size(simplices_all[i],1)-1,num_simplices,simplices_all[i],filtration_weight) # calculate the persistent homology.
            L = filter!(!isempty, L) # remove the empty arrays corresponding to dimensions with no homology classes.
            push!(L_all,L) # add the multiset of persistence intervals for the current model.
        end
    
    ## Distance between the multisets of persistence intervals associated with the models.
        dist_intervals = zeros(Int64,num_models,num_models)
        for i=1:num_models
            for j=(i+1):num_models
                dist_intervals[i,j] = DistanceIntervals(L_all[i],L_all[j])
            end
        end
end
