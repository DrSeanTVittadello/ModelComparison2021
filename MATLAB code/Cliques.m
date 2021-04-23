%%% Calculates the higher-dimensional cliques of the given simplicial complex.
    % Given a simplicial complex with dimension d, Cliques.m calculates the
    % cliques in dimension d+1, which is 'min_dim', and then iteratively in
    % each higher dimension until a dimension is reached which is either equal
    % to 'max_dim' or for which no more cliques exist.

%%% Implementation:
    % INPUT: 'simplices' contains the simplices of the complex.
    % INPUT: 'min_dim' is the minimum dimension clique.
    % INPUT: 'max_dim' is the maximum dimension clique.

    % OUTPUT: 'simplices' with the higher-dimensional cliques appended.

%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


function [simplices] = Cliques(simplices,min_dim,max_dim)
for k=min_dim:max_dim
    simplices_check = zeros(0,k+1); % possible k-simplices.
    for i=1:size(simplices{k,1},1) % (k-1)-simplices.
        diff_vertices = setdiff(simplices{1,1}',simplices{k,1}(i,:)); % vertices not in the (k-1)-simplex.
        diff_vertices = unique(diff_vertices);
        for j=1:size(diff_vertices,2)
            simplices_check = [simplices_check;[simplices{k,1}(i,:) diff_vertices(1,j)]];
        end
    end
    simplices_check = sort(simplices_check,2);
    simplices_check = unique(simplices_check,'rows');
    
    next_dim = 0; % check for cliques in the next dimension? (Yes = 1, No = 0).
    for m=1:size(simplices_check,1)
        faces_check = nchoosek(simplices_check(m,:),k); % check for the (k-1)-faces.
        
        % Find (k-1)-simplices in 'faces_check' that aren't simplices in the given simplicial complex.
        if isempty(setdiff(faces_check,simplices{k,1},'rows')) % true when all of the proper faces of the simplex are in the given simplicial complex.
            next_dim = 1; % check for cliques in the next dimension.
            if nnz(~cellfun(@isempty,simplices(:,1))) < k+1 % true when there are currently no k-dimensional simplices in the simplicial complex.
                simplices{k+1,1} = simplices_check(m,:);
            else
                simplices{k+1,1} = [simplices{k+1,1};simplices_check(m,:)];
            end
        end
    end
    if next_dim==0 % don't check for cliques in the next dimension.
        break
    end
end
simplices = simplices(~cellfun('isempty',simplices));
