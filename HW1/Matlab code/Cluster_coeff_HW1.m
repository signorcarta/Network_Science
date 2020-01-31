close all
clear all
clc

G = importdata('CA-GrQc.txt', '\t', 4);

% adjacency matrix
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);

clear G;

%% Clustering coefficient _________________________________________________
cluster_coeffs = zeros(1, N);

for i=find(diag(A)) % Diagonal cleaning
    A(i,i)=0;
end

d = full(sum(A,1)); 
for i=find(d>1)
    neigh_indexes = find(A(i, :)); % indexed of nodes connected to i
    L = A(neigh_indexes, neigh_indexes); % a factor 2 is already included 
                                         % here, since the graph is undirected
    neigh_connections_i = sum(L(:)); % # of links between neighs
    max_possible_neighs_connections_i = (d(i) * (d(i)-1));
    cluster_coeffs(i) = neigh_connections_i / max_possible_neighs_connections_i;
end

%average clustering coefficient of degrees
average_C_coeff = mean(cluster_coeffs);
disp([ 'average clustering coefficient: ---> ' num2str(average_C_coeff) ' <---']);

%Show result _____________________________________________________________
% figure(3);
% plot(log(d), cluster_coeffs, 'x');
% grid
% xlabel('log(d_{i})')
% ylabel('C_{i}');
% title('Clustering coefficients')

