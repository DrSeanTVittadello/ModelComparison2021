%% Main file for model comparison via simplicial complexes and persistent homology.

%%% Implementation:
    % SOURCE (model-comparison methodology): Vittadello, S. T. and Stumpf, M. P. H., Model comparison via simplicial complexes and persistent homology, arXiv preprint: 2012.13039, 2020.
    % SOURCE (persistent-homology algorithm): Zomorodian, A. & Carlsson, G., Computing persistent homology, Discrete & Computational Geometry, 2005, 33, 249--274.

    % INPUT: 'SC_dist' indicates whether to calculate the distances directly between the simplicial complexes ('Y' or 'N');
    %        the distances between the simplicial complexes are contained in the upper triangular matrix 'dist_complexes'.
    % INPUT: 'PI_dist' indicates whether to calculate the distances between the persistence intervals ('Y' or 'N');
    %        the distances between the persistence intervals are contained in the upper triangular matrix 'dist_intervals'.
    % INPUT: 'model_components' contains all of the model components required by the models for comparison.
    % INPUT: 'dim_ref' is the maximum possible dimension for the reference complex (set equal to `inf' if no maximum is required).
    % INPUT: 'num_models' is the number of models to compare.
    % INPUT: for each model the following is required -
             % INPUT: 'model_num' is the model number, beginning at 1 and incrementing by 1 for subsequent models.
             % INPUT: 'model_components_i' contains all of the model components for model number i.
             % INPUT: 'simplices' contains the simplices of the simplicial complex associated to the model;
             %        the 0-simplices (or vertices) are determined by the specified components of the model in 'model_components_i';
             %        input the simplices of dimension 1 and higher, as required;
             %        if simplices up to dimension d are given, then the function Cliques.m may be employed to incrementally obtain 
             %        the cliques of dimensions higher than d.

    % OUTPUT (when SC_dist = 'Y'): 'dist_complexes' containing the distances between the simplicial complexes associated with the
    %                              models for comparison.
    % OUTPUT (when PI_dist = 'Y'): array 'L_all' containing the persistence intervals of the homology classes in each dimension.
    % OUTPUT (when PI_dist = 'Y'): 'dist_intervals' containing the distances between the persistence intervals associated with the
    %                              models for comparison.
    
%%% Author: Sean T. Vittadello.
%%% Affiliation: The University of Melbourne.
%%% Date: 23 April 2021.


%% Input.
SC_dist = 'Y';
PI_dist = 'Y';

model_components = {'Morphogen 1','Diffusion 1','Degradation 1','Production 1','Basal production 1','Influx 1','Outflux 1','Morphogen 1 bound','Morphogen 2','Diffusion 2','Degradation 2','Basal production 2','Influx 2','Morphogen 3','Diffusion 3','Degradation 3','Production 3','Substrate 1','Diffusion of Substrate 1','Degradation of Substrate 1','Basal production of Substrate 1','Modulator 1','Diffusion of Modulator 1','Degradation of Modulator 1','Production of Modulator 1','Annihilation between Morphogens 1 and 2','Self-activation of Morphogen 1','Activation of Morphogen 2 by Morphogen 1','Inhibition of Morphogen 1 by Morphogen 2','Inhibition of Morphogen 1 by Morphogen 3','Inhibition of Morphogen 3 by Morphogen 1','Inhibition of an inhibition by Morphogen 2','Production of Morphogen 2 by Morphogen 1','Depletion of Substrate 1 by Morphogen 1','Modulation of Diffusion 1 by Modulator 1','Modulation of Degradation 1 by Modulator 1','Inhibition of Modulator 1 by Morphogen 1','Adsorption of Morphogen 1','Desorption of Morphogen 1 bound','Monotonic gradient','Oscillatory gradient','Local scale-invariance','Global scale-invariance'};

dim_ref = inf;

num_models = 9;


%% Initialisation.
dim_max = zeros(1,num_models); % dimension of the simplex spanned by the vertices for each model.
simplices_all = cell(0,num_models); % contains the simplicial complexes associated with the models.
L_all = cell(0,num_models); % contains the persistence intervals for the simplicial complexes associated with the models.


%% Input the simplicial complex of each model for comparison (here, positional-information (PI) and Turing-pattern (TP) models).
    %% Model PI 1 - linear gradient.
    model_num = 1;

    model_components_1 = {'Morphogen 1','Diffusion 1','Influx 1','Outflux 1','Monotonic gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_1));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices.

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,6;1,7;1,40;1,43;2,6;2,7;2,40;2,43;6,40;6,43;7,40;7,43;40,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max,dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model PI 2 - SDD.
    model_num = 2;

    model_components_2 = {'Morphogen 1','Diffusion 1','Degradation 1','Influx 1','Monotonic gradient'};

    vertices = find(ismember(model_components,model_components_2));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices.

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,6;1,40;2,3;2,6;2,40;3,40;6,40];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model PI 3 - Opposing gradients.
    model_num = 3;

    model_components_3 = {'Morphogen 1','Diffusion 1','Degradation 1','Influx 1','Morphogen 2','Diffusion 2','Degradation 2','Influx 2','Monotonic gradient','Local scale-invariance'};

    vertices = find(ismember(model_components,model_components_3));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,6;1,40;1,42;2,3;2,6;2,40;2,42;3,40;3,42;6,40;6,42;9,10;9,11;9,13;9,40;9,42;10,11;10,13;10,40;10,42;11,40;11,42;13,40;13,42;40,42];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model PI 4 - Annihilation.
    model_num = 4;

    model_components_4 = {'Morphogen 1','Diffusion 1','Degradation 1','Influx 1','Morphogen 2','Diffusion 2','Degradation 2','Influx 2','Annihilation between Morphogens 1 and 2','Monotonic gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_4));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,6;1,26;1,40;1,43;2,3;2,6;2,26;2,40;2,43;3,40;3,43;6,40;6,43;9,10;9,11;9,13;9,26;9,40;9,43;10,11;10,13;10,26;10,40;10,43;11,40;11,43;13,40;13,43;26,40;26,43;40,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model PI 5 - Induction-contraction (active modulation).
    model_num = 5;

    model_components_5 = {'Morphogen 1','Diffusion 1','Degradation 1','Influx 1','Modulator 1','Diffusion of Modulator 1','Degradation of Modulator 1','Production of Modulator 1','Modulation of Diffusion 1 by Modulator 1','Modulation of Degradation 1 by Modulator 1','Inhibition of Modulator 1 by Morphogen 1','Monotonic gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_5));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,6;1,37;1,40;1,43;2,3;2,6;2,35;2,40;2,43;3,36;3,40;3,43;6,40;6,43;22,23;22,24;22,25;22,35;22,36;22,40;22,43;23,24;23,25;23,40;23,43;24,40;24,43;25,37;25,40;25,43;35,40;35,43;36,40;36,43;37,40;37,43;40,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model T 1 - Activator-inhibitor.
    model_num = 6;

    model_components_6 = {'Morphogen 1','Diffusion 1','Degradation 1','Basal production 1','Morphogen 2','Diffusion 2','Degradation 2','Basal production 2','Self-activation of Morphogen 1','Activation of Morphogen 2 by Morphogen 1','Inhibition of Morphogen 1 by Morphogen 2','Oscillatory gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_6));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,5;1,27;1,28;1,29;1,41;1,43;2,3;2,5;2,27;2,28;2,29;2,41;2,43;3,41;3,43;5,27;5,41;5,43;9,10;9,11;9,12;9,28;9,29;9,41;9,43;10,11;10,12;10,28;10,29;10,41;10,43;11,41;11,43;12,41;12,43;27,41;27,43;28,29;28,41;28,43;29,41;29,43;41,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model T 2 - Substrate depletion.
    model_num = 7;

    model_components_7 = {'Morphogen 1','Diffusion 1','Degradation 1','Basal production 1','Substrate 1','Diffusion of Substrate 1','Degradation of Substrate 1','Basal production of Substrate 1','Self-activation of Morphogen 1','Depletion of Substrate 1 by Morphogen 1','Oscillatory gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_7));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,5;1,27;1,34;1,41;1,43;2,3;2,5;2,27;2,34;2,41;2,43;3,41;3,43;5,27;5,41;5,43;18,19;18,20;18,21;18,34;18,41;18,43;19,20;19,21;19,41;19,43;20,41;20,43;21,41;21,43;27,41;27,43;34,41;34,43;41,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model T 3 - Inhibition of an inhibition.
    model_num = 8;

    model_components_8 = {'Morphogen 1','Diffusion 1','Degradation 1','Production 1','Morphogen 2','Diffusion 2','Degradation 2','Morphogen 3','Diffusion 3','Degradation 3','Production 3','Inhibition of Morphogen 1 by Morphogen 3','Inhibition of Morphogen 3 by Morphogen 1','Inhibition of an inhibition by Morphogen 2','Production of Morphogen 2 by Morphogen 1','Oscillatory gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_8));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,4;1,31;1,33;1,41;1,43;2,3;2,4;2,30;2,31;2,32;2,33;2,41;2,43;3,41;3,43;4,30;4,31;4,32;4,33;4,41;4,43;9,10;9,11;9,32;9,33;9,41;9,43;10,11;10,30;10,31;10,32;10,33;10,41;10,43;11,41;11,43;14,15;14,16;14,17;14,30;14,41;14,43;15,16;15,17;15,30;15,31;15,32;15,33;15,41;15,43;16,41;16,43;17,30;17,31;17,32;17,33;17,41;17,43;30,41;30,43;31,32;31,41;31,43;32,41;32,43;33,41;33,43;41,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


    %% Model T 4 - Modulation.
    model_num = 9;

    model_components_9 = {'Morphogen 1','Diffusion 1','Degradation 1','Basal production 1','Morphogen 1 bound','Morphogen 2','Diffusion 2','Degradation 2','Basal production 2','Self-activation of Morphogen 1','Activation of Morphogen 2 by Morphogen 1','Inhibition of Morphogen 1 by Morphogen 2','Adsorption of Morphogen 1','Desorption of Morphogen 1 bound','Oscillatory gradient','Global scale-invariance'};

    vertices = find(ismember(model_components,model_components_9));

    dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)

    % 0-simplices represent model components; 1-simplices indicate interconnections between the model components.
    simplices = cell(0,1);
    simplices{1,1} = sortrows(vertices');
    simplices{2,1} = [1,2;1,3;1,5;1,27;1,28;1,29;1,38;1,39;1,41;1,43;2,3;2,5;2,27;2,28;2,29;2,38;2,39;2,41;2,43;3,8;3,41;3,43;5,27;5,41;5,43;8,27;8,28;8,29;8,38;8,39;8,41;8,43;9,10;9,11;9,12;9,28;9,29;9,41;9,43;10,11;10,12;10,27;10,28;10,29;10,38;10,39;10,41;10,43;11,41;11,43;12,41;12,43;27,29;27,41;27,43;28,29;28,41;28,43;29,41;29,43;38,41;38,43;39,41;39,43;41,43];
    simplices = Cliques(simplices,nnz(~cellfun(@isempty,simplices)),min(dim_max(1,model_num),dim_ref));

    simplices_all{1,model_num} = simplices;


%% Calculations
    %% Distance between the simplicial complexes associated with the models.
    if SC_dist == 'Y'
        dist_complexes = zeros(num_models,num_models);
        for i=1:num_models
            for j=(i+1):num_models
                dist_complexes(i,j) = DistanceComplexes(simplices_all{1,i},simplices_all{1,j});
            end
        end
    end


    %% Distance between the multisets of persistence intervals associated with the models.
    if PI_dist == 'Y'
        %% Reference simplicial complex.
        % Maximum dimension of all simplicial representations of models.
        maximum_dim = 0;
        for i=1:num_models
            maximum_dim = max(maximum_dim,nnz(~cellfun(@isempty,simplices_all{:,i}))-1);
        end

        % All required simplices of the reference complex.
            % Here we only need the simplices associated with the simplicial representations, rather than the full
            % reference complex, since we can use any flat filtration to calculate the distances between the models;
            % the required indices of the flat filtration are based on dimension - extension to the full reference
            % complex, if required, provides multiple options for the filtration.
        simplices_ref = cell(maximum_dim+1,1);
        for i=1:num_models
            for p=1:size(simplices_all{1,i},1)
                simplices_ref{p,1} = vertcat(simplices_ref{p,1},simplices_all{1,i}{p,:});
            end     
        end
        for r=1:size(simplices_ref,1)
            simplices_ref{r,1} = unique(simplices_ref{r,1},'rows'); % optionally ordered lexicographically in each dimension.
        end

        % The flat filtration on the reference simplicial complex (multiple options).
        num_ref = zeros(1,maximum_dim+1);
        for i=1:maximum_dim+1
            num_ref(1,i) = size(simplices_ref{i,1},1);
        end

        filtration_weight_ref = cell(maximum_dim+1,1);
        filtration_weight_ref{1,1} = [1:num_ref(1,1)]';
        for i=2:size(num_ref,2)
            filtration_weight_ref{i,1} = [sum(num_ref(1,1:i-1))+1:sum(num_ref(1:i))]';
        end


        %% Persistence intervals for each model, based on the induced filtrations.
        for i=1:num_models
            simplices = simplices_all{:,i}; % simplices for model i.
            filtration_weight = cell(size(simplices)); % contains the filtration weights for the given model.
            for j=1:size(simplices,1) % simplices of dimension j-1.
                filtration_weight{j,1} = zeros(size(simplices{j,1},1),1);
                for k=1:size(simplices{j,1},1) % each simplex of dimension j-1.
                    index = find(ismember(simplices_ref{j,1},simplices{j,1}(k,:),'rows'));
                    filtration_weight{j,1}(k,1) = filtration_weight_ref{j,1}(index,1);
                end
            end

            num_simplices = zeros(1,size(simplices,1));
            for m=1:size(simplices,1)
                num_simplices(1,m) = size(simplices{m,1},1);
            end

            [L] = PersistentHomology(nnz(~cellfun(@isempty,simplices(:,1)))-1,num_simplices,simplices,filtration_weight);

            L = L(~cellfun('isempty',L)); % remove any empty cells from the multiset of persistence intervals L.

            L_all{1,i} = L;
        end

        %% Distance between the multisets of persistence intervals associated with the models.
        dist_intervals = zeros(num_models,num_models);
        for i=1:num_models
            for j=(i+1):num_models
                dist_intervals(i,j) = DistanceIntervals(L_all{:,i},L_all{:,j});
            end
        end
    end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST

% model_num = 1;
% 
% simplices{1,model_num} = [1;2;3;4];% 0-simplices
% simplices{2,model_num} = [2,3;1,2;1,3;1,4;2,4;3,4]; % 1-simplices, indicating interconnections between the model components
% simplices{3,model_num} = [1,2,3;1,2,4;1,3,4;2,3,4];
% 
% dim_max(1,model_num) = size(simplices{1,1},1)-1; % dimension of the simplex spanned by the vertices (maximum dimension)
% 
% simplices_all = Cliques(simplices(:,model_num),nnz(~cellfun(@isempty,simplices(:,model_num))),min(dim_max,dim_ref));
% 
% n = nnz(~cellfun(@isempty,simplices(:,model_num)));
% for i=n+1:size(simplices_all,1)
%     simplices(i,model_num) = simplices_all(i,1);
% end

%%%%%%%%%%%%%%%%%%%



%% 5-cell
% dim_max = 4; % dimension of the simplex spanned by the vertices (maximum dimension)
% 
% simplices = cell(dim_max+1,1); % contains the simplices of the simplicial representation
% simplices{1,1} = [1;2;3;4;5];% 0-simplices
% % simplices{2,1} = [1,2;1,3;1,4;1,5;2,3;2,4;2,5;3,4;3,5;4,5]; % 1-simplices, indicating interconnections between the model components
% simplices{2,1} = [1,2;1,4;2,4]; % 1-simplices, indicating interconnections between the model components
% 
% simplices = Cliques(simplices,2,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_num = 1;
% 
% vertices = [1,2,3,4];
% 
% dim_max(1,model_num) = size(vertices,2)-1; % dimension of the simplex spanned by the vertices (maximum dimension)
% 
% % simplices(:,model_num) = cell(dim_max+1,model_num); % contains the simplices of the simplicial representation
% simplices{1,model_num} = sortrows(vertices');% 0-simplices
% simplices{2,model_num} = [1,2;1,3;1,4;2,3;2,4;3,4];
% simplices{3,model_num} = [1,2,3;1,2,4;1,3,4;2,3,4];
% simplices_all = Cliques(simplices(:,model_num),nnz(~cellfun(@isempty,simplices(:,model_num))),min(dim_max,dim_ref));





%% %%%%%%%%%%%%%%%%%%%%%%%% Plots
% num = 6; % specify the model number
% intervals = L(:,num); % persistence intervals for the specified model
% intervals = intervals(~cellfun('isempty',intervals));clc
% 
% for i=1:size(intervals,1)
%     intervals{i,1} = sortrows(intervals{i,1});
% end
% 
% figure
% tiledlayout(size(intervals,1),1);
% x = 1:sum(num_ref,2);
% plot_colours = 'brgcmbrg';
% group_labels = ["$H_0 \qquad$","$H_1 \qquad$","$H_2 \qquad$","$H_3 \qquad$","$H_4 \qquad$","$H_5 \qquad$","$H_6 \qquad$","$H_7 \qquad$","$H_8 \qquad$"];
% for i=1:size(intervals,1)
%     index = 1;
%     if ~isempty(intervals{i,1})
%         nexttile
%         for j=1:size(intervals{i,1},1)
%             if intervals{i,1}(j,1)==intervals{i,1}(j,2)
%                 plot(intervals{i,1}(j,1),index,'Marker','.','MarkerSize',16,'Color',plot_colours(1,i));
%                 index = index+1;
%             elseif intervals{i,1}(j,1) < min(intervals{i,1}(j,2),max(x))
%                 plot([intervals{i,1}(j,1):min(intervals{i,1}(j,2),max(x))],index*ones(1,min(intervals{i,1}(j,2),max(x))-intervals{i,1}(j,1)+1),'LineStyle','-','LineWidth',2,'Color',plot_colours(1,i));
%                 index = index+1;
%             end
%             hold on
%         end
%         xlabel('Filtration index');
% %         xlim([0,sum(num_ref,2)+20]);
%         xlim([0,900]);
% %         ylim([0,size(intervals{i,1},1)+1]);
% 
% %         ylim([0,54]);
% 
% 
% %         xlim([0,WeightCritical(end,1)]);
% %         xlim([0,ceil(sum(num_skeleton,2)/100)*100]);
%     %%%%%%%%%%
% %     if i==1
% %         xlim([0 10*max(L{1,1}(~isinf(L{1,1}(:,2)),2))]);
% %     end
%     %%%%%%%%%%%%
%         yticks([])
%         yl = ylabel(group_labels(1,i),'Interpreter','latex','Rotation',0);
%         set(gca,'FontSize',20)
%     end
% end
% set(gca,'FontSize',20)