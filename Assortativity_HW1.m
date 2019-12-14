close all
clear 
clc

G = importdata('CA-GrQc.txt', '\t', 4);
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);
clear G;

%remove nodes which are NOT connected
pos = find(sum(A)~=0);
A = A(pos,pos);
N = size(A,1);

%% Assortativity __________________________________________________________
all_d = sum(A,2);
k_tmp = (A*all_d)./all_d; % averages number of neighbours
u = unique(all_d); 
k_nn = zeros(1, length(u)); 

% extract averages for each k
for k = 1:length(u)
    k_nn(k) = mean(k_tmp(all_d==u(k))); 
end

p = polyfit(log(u'),log(k_nn),1); % linear fitting
disp(['- Assortativity factor                     : ---> ' num2str(p(1)) ' <---'])

% Show results ____________________________________________________________
figure(4)
loglog(all_d, k_tmp, 'g.');
hold on
loglog(u,exp(p(2)+log(u)*p(1)),'r-');
loglog(u,k_nn,'k.');
hold off
grid
xlabel('k')
ylabel('k_{nn}')
title('Assortativity')

%% Cutoffs ________________________________________________________________
gamma = 2.592; %% estimated in the Degree_Distribution_HW1 file
S_cut = sqrt(mean(all_d)*N); % structural
N_cut = N^(1/(gamma-1)); % natural (k_min factor is negligible)

disp(['- Structural cutoff: ---> ' num2str(S_cut) ' <---'])
disp(['- Natural cutoff:    ---> ' num2str(N_cut) ' <---'])

if(S_cut < N_cut && gamma<3)
    disp('   ==> structural disassortativity')
end
%% Assortativity under Randomization ______________________________________

% Network rewiring
degree_preservation_A = sparse(A); % matrix a that will be rewired randomly

% Execution of this section may vary on the runs. Might take too long or
% give NaN as result. This is due to the randomness under which the
% rewiring is performed. It might be necessary to run the script multiple 
% times to get a correct value of the new assortativity factor. 
% Values I found: {0.524; 0.523; 0.526} ___________________________________
for i = 2:2:N
    Z1 = i-1;
    Z2 = i;
    links_1 = find(degree_preservation_A(:,Z1));
    links_2 = find(degree_preservation_A(:,Z2));
    while (true)
        C1 = datasample(links_1,1,'Replace',false);
        C2 = datasample(links_2,1,'Replace',false);
        if not((ismember(C1,links_2)) && (ismember(C2,links_1)))
            break;
        end
    end
    degree_preservation_A(Z1,C1) = 0;
    degree_preservation_A(C1,Z1) = 0;
    degree_preservation_A(Z2,C2) = 0;
    degree_preservation_A(C2,Z2) = 0;
    degree_preservation_A(Z1,C2) = 1;
    degree_preservation_A(C2,Z1) = 1;
    degree_preservation_A(Z2,C1) = 1;
    degree_preservation_A(C1,Z2) = 1;    
end %%_____________________________________________________________________

% New assortativity computation, after randomization
new_d = sum(degree_preservation_A,2);
k_tmp_ = (degree_preservation_A*new_d)./new_d; % averages # of neighbours
uns = unique(new_d); 
k_nn_new = zeros(1, length(uns)); 

% extract averages for each k
for k = 1:length(uns)
    k_nn_new(k) = mean(k_tmp_(all_d==uns(k))); 
end

p_new = polyfit(log(uns'),log(k_nn_new),1); % linear fitting
disp(['- Assortativity under randomization factor : ---> ' num2str(p_new(1)) ' <---'])

% Show results ____________________________________________________________
% figure(5)
% loglog(uns,k_nn_new,'b.');
% hold on
% loglog(u, k_nn, 'r.');
% legend('After randomization', 'Original Network', 'Location','northwest')
% hold off
% grid
% xlabel('k')
% ylabel('k_{nn}')
% title('Assortativity under randomization')