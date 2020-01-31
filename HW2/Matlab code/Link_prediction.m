    close all
clear all
clc

G = importdata('sport_edges.csv', ',',1);
T = readtable('sport_ids.csv');
T{:,'id'}=T{:,'id'}+1; % all indexes must be >=1
T = sortrows(T,'new_id','ascend');
T{5258:end,'new_id'}=T{5258:end,'new_id'}-1; % fixing indexes 5258:end

% adjacency matrix
G.data(:,:) = G.data(:,:)+1;
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
N = length(A); % # of nodes

%% Global Link-prediction techniques

G = digraph(A, 'OmitSelfLoops');
N = numnodes(G);
E = numedges(G);
Edges = G.Edges{1:E,1};
din = indegree(G);

comparisons = 1e5; % # number of comparisons performed to compute AUC
n_iterations = 10;
AUC = zeros(n_iterations,2);
precision = zeros(n_iterations,2);

for t = 1:n_iterations
	test_fraction = 0.8;
	i_test = sort(randperm(E, round(test_fraction*E)), 'ascend'); % randomly choose edges in the test set
	test_set = Edges(i_test, :); %test set
    i_probe = setdiff(1:E, i_test);
    probe_set = Edges(i_probe, :); %probe set
	
	% graph corresponding to the test set
    G_T = digraph(zeros(N));
    G_T = addedge(G_T, test_set(:,1), test_set(:,2), ones(1, size(test_set, 1))); % add the links of the test set
    A_T = full(adjacency(G_T));
    dinT = indegree(G_T);
	
	%global techiques
    beta = 1e-2;
    S_Katz = inv(eye(N)-beta*A_T)-eye(N); % Katz
    S_LP = A_T^2 + beta*A_T^3; % Local Path
    
	% current AUC
    AUC_t = zeros(2,1);
    AUC_t(1) = computeAUC(S_Katz,A,comparisons,N,probe_set);
    AUC_t(2) = computeAUC(S_LP,A,comparisons,N,probe_set);
    AUC(t,:) = AUC_t; % store current AUC 
    
    % precision
    L = 100;
    precision_t = zeros(2,1);
    precision_t(1) = computePR(S_Katz,L,test_set,probe_set);
    precision_t(2) = computePR(S_LP,L,test_set,probe_set);
    precision(t,:) = precision_t; % store current precision
end

avgAUC = mean(AUC);
avgPR = mean(precision);

%% Display results ________________________________________________________

disp(['AUC for Katz: ' num2str(avgAUC(1))]);
disp(['Precision for Katz: ' num2str(avgPR(1))]);
disp(['AUC for Local Path ' num2str(avgAUC(2))]);
disp(['Precision for Local Path: ' num2str(avgPR(2))]);

%% [Functions definition]

% AUC _____________________________________________________________________
function auc = computeAUC(S,A,num_comparisons,n,probe_set)
    k = 0;
    auc = 0;
    
    while(k < num_comparisons)
        % randomly select element in the probe set
        i = randi(size(probe_set,1));
        a1 = probe_set(i,1);
        a2 = probe_set(i,2);
        
        b1 = randi(n);
        b2 = randi(n);
        if(A(b1,b2) == 0)
            k = k+1;
            if S(a1,a2) > S(b1,b2) 
                auc = auc+1;
            elseif S(a1,a2) == S(b1,b2)
                auc = auc+0.5;
            end
        end
    end
    auc = auc/num_comparisons; %normalization
end
% _________________________________________________________________________

% Prec ____________________________________________________________________
function prec = computePR(S,L,test_set,probe_set)
    % no self loops
    S1 = S-diag(diag(S));

    % only consider the complementary of the test set
    for k = 1:size(test_set,2)
        S1(test_set(k,1),test_set(k,2)) = 0;
    end

    prec = 0;
    for k = 1:L
        top = max(max(S1));
        [x,y] = find(S1==top,1);

        % remove current max from S1
        S1(x,y) = 0;

        % if the maximum value is in the probe set update prec
        if (ismember([x,y], probe_set, 'rows')) || (ismember([y,x], probe_set, 'rows'))
            prec = prec+1;
        end
    end
    prec = prec/L;
end
% _________________________________________________________________________