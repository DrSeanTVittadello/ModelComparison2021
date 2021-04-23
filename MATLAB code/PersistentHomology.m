%%% Calculates the persistent homology of the given simplicial complex, and outputs the corresponding multiset of persistence intervals.

%%% Implementation:
    % SOURCE (persistent-homology algorithm): Zomorodian, A. & Carlsson, G., Computing persistent homology, Discrete & Computational Geometry, 2005, 33, 249--274.

    % INPUT: 'dim_complex' is the dimension of the simplicial complex.
    % INPUT: 'num_simplices' is the number of simplices in each dimension of the simplicial complex.
    % INPUT: 'simplices' contains the simplices of the simplicial complex.
    % INPUT: 'filtration_weight' contains the filtration weight (or filtration index) of each simplex of the simplicial complex.

    % OUTPUT: array 'L' containing the persistence intervals of the homology classes in each dimension.

%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


function [L] = PersistentHomology(dim_complex,num_simplices,simplices,filtration_weight)
%% Initial calculations.
% Sort each simplex by the component 0-simplices to facilitate
% identification of the subcomplexes.
simplices_sort{1,1} = simplices{1,1};
for i=2:dim_complex+1
    simplices_sort{i,1} = sort(simplices{i,1},2);
end

% Order the simplices in each dimension by increasing index of filtration.
for i=1:dim_complex+1
    [~,faces_indices] = sort(filtration_weight{i,1},1);
    filtration_weight{i,:} = filtration_weight{i,:}(faces_indices,:);
    simplices_sort{i,:} = simplices_sort{i,1}(faces_indices,:);
end

% The dimension of each simplex.
dim_simplices = cell(dim_complex+1,1);
for i=1:dim_complex+1
    dim_simplices{i,1} = (i-1)*ones(num_simplices(1,i),1);
end

%% Calculate the persistence intervals.
% Note that the simplices and associated data are ordered in the cells by
% increasing dimension and then by increasing filtration index within each
% dimension.

% Persistence intervals for each dimension.
L = cell(dim_complex+1,1);

% Nonempty positions in 'chain' indicate that the corresponding simplex is 
% a pivot in the boundary matrix, and the position contains the boundary 
% chain to which it belongs. The boundary chain is the boundary of a 
% simplex with index (w.r.t. the cell 'simplices_sort') in 'chain_index'.
chain = cell(sum(num_simplices),1); % boundary chains.

% Index of the simplex corresponding to the boundary chain.
chain_index = cell(dim_complex+1,1);

% Marked simplices correspond to non-pivot columns in the boundary
% matrix (1/0).
marked = cell(dim_complex+1,1);
for i=1:dim_complex+1
    marked{i,1} = 0.*simplices_sort{i,1}(:,1);
end

% Consider each simplex w.r.t. the total ordering of the simplices until
% the end of the filtration is reached.
for j=1:sum(num_simplices)
    % Convert the linear index j to an index in the cell 'simplices_sort'.
    index_cell = ConvertLinear(j,num_simplices);
    
    % Dimension of the current simplex.
    dim_simplex = dim_simplices{index_cell(1,1),1}(index_cell(1,2),1);
    
    % Find the faces of the current simplex.
    if dim_simplex==0
        faces = [];
    else
        simplex = simplices_sort{index_cell(1,1),1}(index_cell(1,2),:);
        faces = sort(nchoosek(simplex,dim_simplex)); % each row is a face.
    end
    
    %% Remove the unmarked terms in 'faces'
    faces_indices2 = zeros(size(faces,1),2);
    for t=1:size(faces,1)
         faces_indices2(t,:) = [dim_simplex,find(ismember(simplices_sort{dim_simplex,1},faces(t,:),'rows'))];
    end
    
    keep = ones(size(faces,1),1);
    for i=1:size(faces,1)
        if marked{faces_indices2(i,1),1}(faces_indices2(i,2,1))~=1
            keep(i,1) = 0;
        end
    end
    faces = faces(find(keep),:);

    %% Find the face of the simplex with the maximum index (possible pivot).
    flag = 0;
    while ~isempty(faces) && flag==0
        faces_indices3 = zeros(size(faces,1),1);
        for t=1:size(faces,1)
             faces_indices3(t,1) = find(ismember(simplices_sort{dim_simplex,1},faces(t,:),'rows'));
        end
        [~,max_index] = max(faces_indices3);
        
        % Determine the index of the face in the cell 'simplices_sort'.
        if dim_simplex==1
            face_index = [1,faces_indices3(max_index,1)];
        else
            face_index = [dim_simplex,faces_indices3(max_index,1)];
        end
        
        % If a pivot doesn't exist in the corresponding row of the boundary
        % matrix then terminate the while loop.
        if isempty(chain{ConvertCell(face_index,num_simplices),1})
            flag = 1;
        end
        
        % If a pivot exists in the corresponding row of the boundary matrix,
        % then use the inverse of its field coefficient to eliminate the row
        % from the chain (easy for our field Z/2Z).
        if flag~=1
            pivot = chain{ConvertCell(face_index,num_simplices),1};
            for i=1:size(pivot,1)
                term = find(ismember(faces,pivot(i,:),'rows'),1);
                if ~isempty(term) 
                    faces(term,:) = [];
                    faces_indices3(term,:) = [];
                else   
                    faces = cat(1,faces,pivot(i,:));  
                end
            end
        end
    end
    
    % If the simplex has empty boundary.
    if isempty(faces)
        marked{index_cell(1,1),1}(index_cell(1,2),1) = 1;
        
    % Otherwise the chain 'faces' corresponds to a pivot column in the
    % boundary matrix and the term with the maximum index is the pivot.
    else
        [~,max_index] = max(faces_indices3);
        face_index = [dim_simplex,faces_indices3(max_index,1)]; % w.r.t. the cell 'simplices_sort'.
        
        chain{ConvertCell(face_index,num_simplices),1} = faces;
        chain_index{face_index(1,1),1}(face_index(1,2),:) = j;
        
        % Add finite persistence intervals corresponding to the homology 
        % group of order dim_simplex-1.
        L{dim_simplex,1} = cat(1,L{dim_simplex,1},[filtration_weight{face_index(1,1),1}(face_index(1,2),:),filtration_weight{index_cell(1,1),1}(index_cell(1,2),:)]);
    end
end

%% Infinite persistence intervals
% Add infinite persistence intervals corresponding to the homology group
% of order dim_simplex.
for j=1:sum(num_simplices)
    index_cell = ConvertLinear(j,num_simplices); % convert linear index j to an index in the cell 'simplices_sort'.
    dim_simplex = dim_simplices{index_cell(1,1),1}(index_cell(1,2),1);
    if marked{index_cell(1,1),1}(index_cell(1,2),1)==1 && isempty(chain{j,1})
        L{dim_simplex+1,1} = cat(1,L{dim_simplex+1,1},[filtration_weight{index_cell(1,1),1}(index_cell(1,2),:),Inf]);
    end
end
