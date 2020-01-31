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

%% Preprocessing __________________________________________________________
A = 1*(A+A'>0); % undirected network

% Remove self loops _______________________________________________________
for i=1:N
    A(i,i)=0;
end

% Remove nodes which are NOT connected ____________________________________

pos = find(sum(A)~=0);
A = A(pos,pos);

% Remove dead ends  _______________________________________________________
exit = false;

 while (~exit)
    pos = find(sum(A)~=0);
    A = A(pos,pos);
    exit = isempty(find(sum(A)==0, 1));
end

N = size(A,1); % update N

% Find the largest connected component ____________________________________
e1 = [1;zeros(N-1,1)];
exit = false;
while(~exit)
    e1_old = e1;
    e1 = 1*(A*e1>0);
    exit = (sum(e1-e1_old)==0);
end
pos = find(e1);
A = A(pos,pos);
N = size(A,1); % update N

%% ________________________________________________________________________
% computing nodes rank %
c = 0.85; % damping factor
q = ones(N,1)/N; % normalized teleportation vector
M = A*sparse(diag(1./sum(A))); % normalized M
r = sparse((eye(N)-c*M)/(1-c))\q;
r = r/sum(r);
T = sortrows(T,'new_id','ascend');

%% SimRank ________________________________________________________________
% [max, argmax] = max(r); 
q_simrank = zeros(N,1);
q_simrank(1823) = 1; % Mirza Teletovic's page

c = 0.85; % damping factor
M = A*sparse(diag(1./sum(A))); % normalized M
sim_r = sparse((eye(N)-c*M)/(1-c))\q_simrank; % SimRank teleportation vector
sim_r = sim_r/sum(sim_r);

rankT = table(sim_r);
T = [T rankT];% add the simRank score as the fourth column
T = sortrows(T,'sim_r','descend'); % sort the dataset 
T(2:11,:)  % display the top-10 similar pages
%%
